"""
I-1 cold-open — transitive-inference RELATION CHAIN (the shown premises), talk schematic.

The "here's the ranking" reveal shown BEFORE the race bracket. Four disc-golf-player
photo-placeholder boxes in the TRUE order FB · AG · DG · SM, with a large ">" between
each adjacent pair — the three results you are actually SHOWN:

    FB > AG ,  AG > DG ,  DG > SM

(gold ">" = a seen/given result). The bracket that follows re-seeds these four so its
two semifinals are the NON-adjacent inferences and only the final is one of these
shown premises. Boxes are PHOTO PLACEHOLDERS — drop the four headshots in Canva.

Dark ink on a LIGHT slide: black text/borders, gray initial hints, gold ">". Transparent PDF.
Output: data/outputs/figures/talk/ti_relation_chain.pdf
QA (white-matte PNG, opt-in): python fig_ti_relation_chain.py --qa /path/to/qa.png
"""
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _ti_faces import load_face, place_face  # noqa: E402

use_lrg_style()

# ---- palette (light slide, dark ink; matches the race bracket) ------------
DARK = "#1a1a1a"     # box borders, names
GRAY = "#6f757c"     # placeholder initial hints
GOLD = "#f2c14e"     # ">" = a shown/seen result

# true rank order (left = highest)
NAMES = ["Battiston", "Gabrielli", "Garlaschelli", "Meloni"]
INIT = ["FB", "AG", "DG", "SM"]
BW = BH = 1.6                     # square photo-placeholder box
GAP = 1.5                         # room for the ">" between boxes
XS = [i * (BW + GAP) for i in range(4)]


def main():
    fig, ax = plt.subplots(figsize=(10.5, 3.1))
    ax.set_aspect("equal")
    ax.axis("off")

    for x, name, ini in zip(XS, NAMES, INIT):
        img = load_face(ini)
        patch = FancyBboxPatch(
            (x, 0), BW, BH,
            boxstyle="round,pad=0.02,rounding_size=0.14",
            linewidth=1.8, edgecolor=DARK,
            facecolor=("none" if img is not None else "white"), zorder=2)
        ax.add_patch(patch)
        if img is not None:
            place_face(ax, img, x + BW / 2, BH / 2, BW, patch, zorder=1.7)
        else:
            ax.text(x + BW / 2, BH / 2, ini, ha="center", va="center", color=GRAY,
                    fontsize=15, alpha=0.75, style="italic")       # photo hint
        ax.text(x + BW / 2, -0.30, name, ha="center", va="top", color=DARK,
                fontsize=13)

    # ">" between adjacent boxes = the three SHOWN premises
    for i in range(3):
        xc = XS[i] + BW + GAP / 2
        ax.text(xc, BH / 2, ">", ha="center", va="center", color=GOLD,
                fontsize=34, fontweight="bold")

    ax.set_xlim(-0.3, XS[-1] + BW + 0.3)
    ax.set_ylim(-0.98, BH + 0.30)

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "ti_relation_chain.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    print("wrote", out)

    if "--qa" in sys.argv:
        qa = sys.argv[sys.argv.index("--qa") + 1]
        fig.savefig(qa, facecolor="white", bbox_inches="tight", dpi=130)
        print("wrote QA", qa)
    plt.close(fig)


if __name__ == "__main__":
    main()
