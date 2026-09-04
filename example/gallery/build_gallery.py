"""One-off driver that generated the README gallery renders. Not part of
the public CLI -- just calls noodlex's functions directly so each example
gets its own clean .cxc (no gray reference chain, just tube + key) and its
own independent ChimeraX window.

Regenerate + reopen everything with:
    python3 example/gallery/build_gallery.py
    open -n -a ChimeraX example/gallery/showcase_1_default_view.cxc
    open -n -a ChimeraX example/gallery/showcase_2_clipped_view.cxc
    open -n -a ChimeraX example/gallery/showcase_3_alphafold_view.cxc
then File > Save Image (or just a screenshot) in each window.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO_ROOT)
import noodlex as px

BASE = os.path.dirname(os.path.abspath(__file__))
CRAMBIN = os.path.join(REPO_ROOT, "example", "1crn.pdb")
P53 = os.path.join(BASE, "af_p53.pdb")  # AF-P04637-F1 (human p53), fetched from alphafold.ebi.ac.uk


def write_clean_viewer(name, bild_path, vmin, vmax, stops):
    cxc_path = os.path.join(BASE, f"{name}_view.cxc")
    key_parts = []
    n = len(stops)
    for i, c in enumerate(stops):
        v = vmin + (i / (n - 1)) * (vmax - vmin)
        key_parts.append(f"{px._rgb_to_hex(c)}:{v:.2f}")
    key_cmd = "key " + " ".join(key_parts) + " pos 0.62,0.06 size 0.18,0.025 fontSize 16"
    lines = [
        "close session",
        f"open {bild_path}",
        "lighting soft",
        "set bgColor white",
        key_cmd,
        "view",
    ]
    with open(cxc_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return cxc_path


examples = []

# 1. hero shot / baseline: default blue:white:red, data-driven range
bild = os.path.join(BASE, "showcase_1_default.bild")
vmin, vmax, stops = px.noodle_bild(CRAMBIN, bild)
examples.append(write_clean_viewer("showcase_1_default", bild, vmin, vmax, stops))

# 2. outlier robustness: same data, percentile-clipped range
bild = os.path.join(BASE, "showcase_2_clipped.bild")
vmin, vmax, stops = px.noodle_bild(CRAMBIN, bild, clip_percentile=(5, 95))
examples.append(write_clean_viewer("showcase_2_clipped", bild, vmin, vmax, stops))

# 3. AlphaFold pLDDT preset on a real AlphaFold model (p53, mixed confidence)
bild = os.path.join(BASE, "showcase_3_alphafold.bild")
vmin, vmax, stops = px.noodle_bild(
    P53, bild, palette="alphafold", value_min=0.0, value_max=100.0, smooth=3
)
examples.append(write_clean_viewer("showcase_3_alphafold", bild, vmin, vmax, stops))

for e in examples:
    print(e)
