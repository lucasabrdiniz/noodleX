#!/usr/bin/env python3
"""
noodlex.py: B-factor "putty" (worm) representation for UCSF ChimeraX
=====================================================================

Draws a per-residue backbone tube whose radius AND color track a scalar
value (by default, the B-factor column of the input PDB) -- the same idea
as PyMOL's `cartoon putty`.

Why this exists instead of ChimeraX's native `worm` / `cartoon byattribute`
command (available since ChimeraX 1.8):
  - it needs no ChimeraX version check -- it works on any ChimeraX (or even
    Chimera classic) that can open a .bild file;
  - no bundle/plugin install: like porcupineX, this just emits a plain
    .bild file of primitives (`.sphere` / `.cylinder`) that you `open`
    directly, so it drops into a scripted/batch pipeline with zero setup;
  - the radius/color mapping is a few lines of Python you can edit freely,
    instead of being locked to whatever curve the built-in command offers.

If you're working interactively in ChimeraX and don't need any of that,
the native command is less code:
    worm bfactor #1
    color byattribute bfactor #1 palette blue:white:red

Usage:
    python3 noodlex.py input.pdb output.bild [options]

Dependencies: numpy, MDAnalysis. matplotlib (+ optionally seaborn) only if
--plot is used.
"""

import argparse
import csv as _csv
import glob
import os
import platform
import shutil
import subprocess

import numpy as np
import MDAnalysis as mda


# --------------------------------------------------------------------------
# color helpers
# --------------------------------------------------------------------------

_NAMED_COLORS = {
    "blue": (0.0, 0.0, 1.0),
    "white": (1.0, 1.0, 1.0),
    "red": (1.0, 0.0, 0.0),
    "green": (0.0, 1.0, 0.0),
    "yellow": (1.0, 1.0, 0.0),
    "cyan": (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "black": (0.0, 0.0, 0.0),
    "orange": (1.0, 0.6, 0.0),
    "gray": (0.5, 0.5, 0.5),
    "grey": (0.5, 0.5, 0.5),
}


def _parse_color(s):
    """Parse a named color, '#rrggbb', or 'r,g,b' (0-1 floats) into an
    (r, g, b) tuple."""
    s = s.strip()
    if s.lower() in _NAMED_COLORS:
        return _NAMED_COLORS[s.lower()]
    if s.startswith("#"):
        s = s[1:]
        r = int(s[0:2], 16) / 255.0
        g = int(s[2:4], 16) / 255.0
        b = int(s[4:6], 16) / 255.0
        return (r, g, b)
    parts = [float(p) for p in s.split(",")]
    assert len(parts) == 3, f"bad color spec: {s!r}"
    return tuple(parts)


_PALETTE_PRESETS = {
    # AlphaFold DB / ColabFold pLDDT coloring (low -> high confidence),
    # as a continuous gradient across the same 4 reference colors the
    # official (discrete, binned) scheme uses. pLDDT is always 0-100, so
    # --palette alphafold defaults --value-min/--value-max to that range
    # (see main()) unless you override them.
    "alphafold": ["#FF7D45", "#FFDB13", "#65CBF3", "#0053D6"],
}


def _parse_palette(s):
    """Parse a ':'-separated list of colors (>=2 stops), or a preset name
    from _PALETTE_PRESETS, into RGB tuples."""
    if s in _PALETTE_PRESETS:
        s = ":".join(_PALETTE_PRESETS[s])
    stops = [_parse_color(p) for p in s.split(":")]
    assert len(stops) >= 2, "palette needs at least 2 color stops"
    return stops


def _rgb_to_hex(rgb):
    r, g, b = (max(0, min(255, round(c * 255))) for c in rgb)
    return f"#{r:02x}{g:02x}{b:02x}"


def _lerp(t, a, b):
    return a + (b - a) * t


def _lerp_color(t, stops):
    """Piecewise-linear interpolation of t (0-1) across N color stops."""
    if t <= 0.0:
        return stops[0]
    if t >= 1.0:
        return stops[-1]
    n = len(stops) - 1
    seg = t * n
    i = min(int(seg), n - 1)
    local_t = seg - i
    a, b = stops[i], stops[i + 1]
    return tuple(_lerp(local_t, a[k], b[k]) for k in range(3))


def _value_to_t(v, vmin, vmax):
    if vmax <= vmin:
        return 0.0
    return min(1.0, max(0.0, (v - vmin) / (vmax - vmin)))


def _resolve_value_range(values, value_min, value_max, clip_percentile):
    """Pick (vmin, vmax) for the value->radius/color mapping. An explicit
    --value-min/--value-max always wins; otherwise --clip-percentile
    derives the range from percentiles of the data (robust to a single
    outlier stretching the whole scale), falling back to plain min/max."""
    if clip_percentile is not None:
        lo, hi = clip_percentile
        vmin = value_min if value_min is not None else float(np.percentile(values, lo))
        vmax = value_max if value_max is not None else float(np.percentile(values, hi))
    else:
        vmin = value_min if value_min is not None else float(values.min())
        vmax = value_max if value_max is not None else float(values.max())
    return vmin, vmax


# --------------------------------------------------------------------------
# geometry: catmull-rom spline through the CA (or chosen atom) trace
# --------------------------------------------------------------------------

def _catmull_rom_centripetal(p0, p1, p2, p3, local_fracs, alpha=0.5):
    """Barry-Goldman formulation of Catmull-Rom, parameterized by chord
    length**alpha. alpha=0.5 ("centripetal") avoids the loops/overshoot
    that plain uniform Catmull-Rom (alpha=0) produces on sharply-curved or
    unevenly-spaced control points -- exactly the case for a protein
    backbone with tight turns. `local_fracs` are fractions in [0, 1)
    across the p1->p2 segment; the t1/t2 parameter bounds they map into
    are derived here, from these same four points, so there's a single
    source of truth for the parameterization."""
    def _seg(ti, pa, pb):
        d = np.linalg.norm(pb - pa)
        return ti + (d ** alpha if d > 1e-9 else 1e-6)

    t0 = 0.0
    t1 = _seg(t0, p0, p1)
    t2 = _seg(t1, p1, p2)
    t3 = _seg(t2, p2, p3)

    out = []
    for f in local_fracs:
        t = t1 + f * (t2 - t1)
        a1 = p0 * ((t1 - t) / (t1 - t0)) + p1 * ((t - t0) / (t1 - t0))
        a2 = p1 * ((t2 - t) / (t2 - t1)) + p2 * ((t - t1) / (t2 - t1))
        a3 = p2 * ((t3 - t) / (t3 - t2)) + p3 * ((t - t2) / (t3 - t2))
        b1 = a1 * ((t2 - t) / (t2 - t0)) + a2 * ((t - t0) / (t2 - t0))
        b2 = a2 * ((t3 - t) / (t3 - t1)) + a3 * ((t - t1) / (t3 - t1))
        c = b1 * ((t2 - t) / (t2 - t1)) + b2 * ((t - t1) / (t2 - t1))
        out.append(c)
    return out


def _smooth_chain(points, values, subdiv):
    """Subdivide a (N,3) point chain + (N,) values with a centripetal
    Catmull-Rom spline, linearly interpolating values along the way.
    subdiv<=1 returns the input unchanged."""
    n = len(points)
    if n < 2 or subdiv <= 1:
        return np.asarray(points, dtype=float), np.asarray(values, dtype=float)

    pts, vals = [], []
    for i in range(n - 1):
        p0 = points[i - 1] if i - 1 >= 0 else points[i] - (points[i + 1] - points[i])
        p1 = points[i]
        p2 = points[i + 1]
        p3 = points[i + 2] if i + 2 < n else points[i + 1] + (points[i + 1] - points[i])
        v1, v2 = values[i], values[i + 1]
        local_fracs = [s / subdiv for s in range(subdiv)]
        seg_pts = _catmull_rom_centripetal(p0, p1, p2, p3, local_fracs)
        pts.extend(seg_pts)
        for f in local_fracs:
            vals.append(_lerp(f, v1, v2))
    pts.append(points[-1])
    vals.append(values[-1])
    return np.array(pts), np.array(vals)


def _contiguous_runs(resids):
    """Split a sorted array of resids into index runs wherever there's a
    gap (missing residue) so we don't draw a tube across a chain break."""
    runs = []
    start = 0
    for i in range(1, len(resids)):
        if resids[i] - resids[i - 1] != 1:
            runs.append((start, i))
            start = i
    runs.append((start, len(resids)))
    return runs


# --------------------------------------------------------------------------
# .bild writing
# --------------------------------------------------------------------------

def _write_run_bild(fh, points, radii, colors):
    """Write one contiguous chain run as spheres (joints) + cylinders
    (segments) so bends don't show a gap between primitives."""
    n = len(points)
    if n == 0:
        return
    if n == 1:
        r, g, b = colors[0]
        fh.write(f".color {r:.3f} {g:.3f} {b:.3f}\n")
        fh.write(f".sphere {points[0][0]:.3f} {points[0][1]:.3f} {points[0][2]:.3f} {radii[0]:.3f}\n")
        return

    for i in range(n):
        r, g, b = colors[i]
        fh.write(f".color {r:.3f} {g:.3f} {b:.3f}\n")
        p = points[i]
        fh.write(f".sphere {p[0]:.3f} {p[1]:.3f} {p[2]:.3f} {radii[i]:.3f}\n")
        if i < n - 1:
            p2 = points[i + 1]
            r_avg = 0.5 * (radii[i] + radii[i + 1])
            fh.write(
                f".cylinder {p[0]:.3f} {p[1]:.3f} {p[2]:.3f} "
                f"{p2[0]:.3f} {p2[1]:.3f} {p2[2]:.3f} {r_avg:.3f}\n"
            )


def _write_legend_bar(fh, bbox_min, bbox_max, stops, vmin, vmax, scale_min, scale_max, n=20):
    """Append a straight reference bar next to the structure, built with
    the same primitives as the tube: radius sweeps scale_min->scale_max
    and color sweeps the palette, in lockstep, so it doubles as a key for
    both the thickness and the color mapping. Plain .bild has no text
    primitive that reliably renders in 3D, so the numeric endpoints
    (vmin/vmax) aren't labeled here -- when opened via --open, the
    generated .cxc adds a proper labeled `key` legend on top of this."""
    bbox_min = np.asarray(bbox_min, dtype=float)
    bbox_max = np.asarray(bbox_max, dtype=float)
    extents = bbox_max - bbox_min
    long_axis = int(np.argmax(extents))
    offset_axis = (long_axis + 1) % 3
    margin = 0.15 * extents[long_axis] + scale_max * 3

    start = bbox_min.copy()
    start[offset_axis] -= margin
    end = start.copy()
    end[long_axis] = bbox_min[long_axis] + extents[long_axis] * 0.6

    bar_pts = np.linspace(start, end, n)
    t = np.linspace(0.0, 1.0, n)
    radii = scale_min + t * (scale_max - scale_min)
    colors = [_lerp_color(ti, stops) for ti in t]
    _write_run_bild(fh, bar_pts, radii, colors)


# --------------------------------------------------------------------------
# main pipeline
# --------------------------------------------------------------------------

def noodle_bild(
    pdb_path,
    bild_path,
    atom="CA",
    chains=None,
    frame=0,
    scale_min=0.3,
    scale_max=2.0,
    value_min=None,
    value_max=None,
    clip_percentile=None,
    palette="blue:white:red",
    smooth=4,
    legend=False,
    csv_path=None,
    plot_path=None,
    label_top=0,
):
    """Returns (vmin, vmax, stops) actually used for the value->radius/color
    mapping, so callers (e.g. --open's `key` legend) can match it exactly."""
    u = mda.Universe(pdb_path)
    u.trajectory[frame]

    sel = u.select_atoms(f"protein and name {atom}")
    assert len(sel) > 0, f"no '{atom}' atoms found on the protein backbone"

    try:
        chain_ids = np.asarray(sel.chainIDs)
    except (AttributeError, mda.exceptions.NoDataError):
        chain_ids = np.asarray(sel.segids)

    if chains:
        wanted = set(c.strip() for c in chains.split(","))
        mask = np.isin(chain_ids, list(wanted))
        sel = sel[mask]
        chain_ids = chain_ids[mask]
        assert len(sel) > 0, f"no atoms left after filtering chains={chains!r}"

    values = sel.tempfactors.astype(float)
    vmin, vmax = _resolve_value_range(values, value_min, value_max, clip_percentile)

    stops = _parse_palette(palette) if isinstance(palette, str) else list(palette)

    bbox_min = sel.positions.min(axis=0)
    bbox_max = sel.positions.max(axis=0)

    all_resids, all_resnames, all_chains, all_values, all_radii = [], [], [], [], []

    with open(bild_path, "w") as fh:
        for cid in dict.fromkeys(chain_ids):  # preserve first-seen order
            idx = np.where(chain_ids == cid)[0]
            order = np.argsort(sel.resids[idx])
            idx = idx[order]

            resids = sel.resids[idx]
            pts = sel.positions[idx]
            vals = values[idx]

            for start, end in _contiguous_runs(resids):
                run_pts = pts[start:end]
                run_vals = vals[start:end]
                s_pts, s_vals = _smooth_chain(run_pts, run_vals, smooth)

                t = np.array([_value_to_t(v, vmin, vmax) for v in s_vals])
                radii = scale_min + t * (scale_max - scale_min)
                colors = [_lerp_color(ti, stops) for ti in t]

                _write_run_bild(fh, s_pts, radii, colors)

                # bookkeeping (per-residue, not per spline subdivision)
                run_resids = resids[start:end]
                run_resnames = sel.resnames[idx][start:end]
                t_res = np.array([_value_to_t(v, vmin, vmax) for v in run_vals])
                radii_res = scale_min + t_res * (scale_max - scale_min)
                all_resids.extend(run_resids.tolist())
                all_resnames.extend(run_resnames.tolist())
                all_chains.extend([cid] * len(run_resids))
                all_values.extend(run_vals.tolist())
                all_radii.extend(radii_res.tolist())

        if legend:
            _write_legend_bar(fh, bbox_min, bbox_max, stops, vmin, vmax, scale_min, scale_max)

    if csv_path:
        _write_csv(csv_path, all_chains, all_resids, all_resnames, all_values, all_radii)
    if plot_path:
        _write_plot(plot_path, all_resids, all_values, label_top=label_top)

    return vmin, vmax, stops


def _write_csv(csv_path, chains, resids, resnames, values, radii):
    with open(csv_path, "w", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(["chain", "resid", "resname", "value", "radius"])
        for c, i, n, v, r in zip(chains, resids, resnames, values, radii):
            w.writerow([c, i, n, f"{v:.3f}", f"{r:.3f}"])


def _write_plot(plot_path, resids, values, label_top=0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
        sns.set_style("whitegrid")
    except ImportError:
        pass

    resids = np.asarray(resids)
    values = np.asarray(values)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(resids, values, color="#d92626", lw=1.5)
    ax.set_xlabel("Residue")
    ax.set_ylabel("value (B-factor)")
    ax.set_title("noodleX profile")

    if label_top > 0:
        order = np.argsort(values)[::-1][:label_top]
        for i in order:
            ax.annotate(
                str(resids[i]), (resids[i], values[i]),
                textcoords="offset points", xytext=(0, 5), fontsize=7,
            )

    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)


# --------------------------------------------------------------------------
# --open: generate a viewer .cxc and launch ChimeraX with it
# --------------------------------------------------------------------------

def _find_chimerax():
    """Best-effort search for a ChimeraX executable/app. Returns a path
    suitable for launching, or None if nothing was found."""
    env = os.environ.get("CHIMERAX_APP")
    if env:
        return env
    if platform.system() == "Darwin":
        candidates = sorted(glob.glob("/Applications/ChimeraX*.app"), reverse=True)
        if candidates:
            return candidates[0]
    exe = shutil.which("chimerax") or shutil.which("ChimeraX")
    return exe


def _write_viewer_cxc(cxc_path, pdb_path, bild_path, stops, vmin, vmax):
    """Write a small ChimeraX command script that opens the structure +
    the tube and adds a proper labeled `key` legend (the .bild
    file's own legend bar has no text -- BILD has no reliable 3D text
    primitive -- so the readable min/max labels live here instead)."""
    key_parts = []
    n = len(stops)
    for i, c in enumerate(stops):
        v = vmin + (i / (n - 1)) * (vmax - vmin)
        key_parts.append(f"{_rgb_to_hex(c)}:{v:.2f}")
    key_cmd = "key " + " ".join(key_parts) + " pos 0.62,0.06 size 0.18,0.025 fontSize 16"

    lines = [
        "close session",  # start clean even if ChimeraX is already running with other models open
        f"open {os.path.abspath(pdb_path)}",
        "hide atoms",
        "color /A gray(200)",
        f"open {os.path.abspath(bild_path)}",
        "lighting soft",
        "set bgColor white",
        key_cmd,
        "view",
    ]
    with open(cxc_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _launch_chimerax(cxc_path):
    app = _find_chimerax()
    if app is None:
        print(
            "[noodlex] Could not find ChimeraX automatically. Open manually:\n"
            f"  open -a ChimeraX {cxc_path}   (macOS)\n"
            f"  chimerax {cxc_path}           (Linux/Windows)"
        )
        return
    try:
        if platform.system() == "Darwin" and app.endswith(".app"):
            subprocess.Popen(["open", "-a", app, cxc_path])
        else:
            subprocess.Popen([app, cxc_path])
        print(f"[noodlex] Launched ChimeraX with {cxc_path}")
    except OSError as e:
        print(f"[noodlex] Failed to launch ChimeraX automatically ({e}). Open manually:\n  {cxc_path}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdb", help="input structure (PDB) with B-factor values to visualize")
    p.add_argument("bild", help="output .bild path to open in ChimeraX")
    p.add_argument("--atom", default="CA", help="backbone atom to trace (default: CA)")
    p.add_argument("--chain", default=None, help="comma-separated chain IDs to include (default: all)")
    p.add_argument("--frame", type=int, default=0, help="model/frame index to read (default: 0)")
    p.add_argument("--scale-min", type=float, default=0.3, help="tube radius at the lowest value (default: 0.3)")
    p.add_argument("--scale-max", type=float, default=2.0, help="tube radius at the highest value (default: 2.0)")
    p.add_argument("--value-min", type=float, default=None, help="fix the low end of the value range instead of using the data min (or percentile)")
    p.add_argument("--value-max", type=float, default=None, help="fix the high end of the value range instead of using the data max (or percentile)")
    p.add_argument("--clip-percentile", default=None, metavar="LOW,HIGH", help="derive the value range from data percentiles instead of raw min/max, robust to a single outlier stretching the whole scale, e.g. '5,95' (0-100 each)")
    p.add_argument("--palette", default="blue:white:red", help="':'-separated color stops (each a name, '#rrggbb', or 'r,g,b'), or a preset name ('alphafold') (default: blue:white:red)")
    p.add_argument("--smooth", type=int, default=4, help="spline subdivisions per residue-residue segment; 1 disables smoothing (default: 4)")
    p.add_argument("--legend", action="store_true", help="append a straight color/radius reference bar to the .bild (unlabeled -- prefer --open for a labeled 2D key)")
    p.add_argument("--open", action="store_true", help="write a companion .cxc viewer script (structure + tube + labeled key legend) and launch ChimeraX with it")
    p.add_argument("--csv", default=None, help="write a per-residue value/radius table to this CSV path")
    p.add_argument("--plot", default=None, help="write a PNG profile plot to this path (needs matplotlib)")
    p.add_argument("--label-top", type=int, default=0, help="annotate the N highest-value residues on --plot")
    args = p.parse_args()

    clip_percentile = None
    if args.clip_percentile:
        lo, hi = (float(x) for x in args.clip_percentile.split(","))
        assert 0 <= lo < hi <= 100, "--clip-percentile needs 0 <= LOW < HIGH <= 100"
        clip_percentile = (lo, hi)

    value_min, value_max = args.value_min, args.value_max
    if args.palette == "alphafold" and value_min is None and value_max is None:
        value_min, value_max = 0.0, 100.0  # pLDDT is always on a fixed 0-100 scale

    vmin, vmax, stops = noodle_bild(
        args.pdb,
        args.bild,
        atom=args.atom,
        chains=args.chain,
        frame=args.frame,
        scale_min=args.scale_min,
        scale_max=args.scale_max,
        value_min=value_min,
        value_max=value_max,
        clip_percentile=clip_percentile,
        palette=args.palette,
        smooth=args.smooth,
        legend=args.legend,
        csv_path=args.csv,
        plot_path=args.plot,
        label_top=args.label_top,
    )

    if args.open:
        cxc_path = os.path.splitext(args.bild)[0] + "_view.cxc"
        _write_viewer_cxc(cxc_path, args.pdb, args.bild, stops, vmin, vmax)
        _launch_chimerax(cxc_path)


if __name__ == "__main__":
    main()
