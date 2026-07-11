"""
Slide 03 (sEEG dataset) — the four-phase acquisition protocol as a horizontal timeline.

No paper describes OUR protocol, so it is a schematic. A left->right time arrow carries four
rounded phase blocks in sequence, coloured REST (cool) vs TASK (warm):

    rest_pre  ->  task_learn  ->  task_test  ->  rest_post
    (baseline)   (adjacent      (non-adjacent   (post-task
                  premises)      inferences)     rest)

The two middle blocks are bracketed as the transitive-inference task (learn adjacent premises,
then test the non-adjacent inferences — the disc-golf exercise scaled to many items). The outer
two are rest. Chevrons mark the ordering; a time axis runs underneath.

DURATIONS: left OFF. Our recordings are heterogeneous per patient (a slide talking point) and the
acquisition-protocol lengths are not yet pinned (Lane R: Ricci et al. 2023, or the collaborator).
The block below is wire-ready — fill DURATIONS with confirmed minutes and they render as "~X min"
inside each block. rest is ~10 min; do NOT invent task_learn / task_test until sourced.

Dark ink on a LIGHT slide (matches the slide-02 TI cold-open figures): black text, cool/warm
blocks, amber task bracket. Transparent vector PDF.
Output: data/outputs/figures/talk/ti_protocol_timeline.pdf
QA (white-matte PNG, opt-in): python fig_ti_protocol_timeline.py --qa /path/to/qa.png
"""
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

use_lrg_style()

# ---- palette (light slide, dark ink; shared with the TI cold-open figures) ----
DARK = "#1a1a1a"      # tokens, chevrons, time axis, bracket label
GRAY = "#6f757c"      # descriptors, time label
REST_FILL, REST_EDGE = "#dce4ec", "#5b6b7a"   # cool = rest
TASK_FILL, TASK_EDGE = "#f7e4bb", "#b9832a"   # warm = task
AMBER_TXT = "#9a6a12"                          # task-bracket label

# ---- protocol ------------------------------------------------------------------
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
KIND = {"rest_pre": "rest", "task_learn": "task", "task_test": "task", "rest_post": "rest"}
DESC = {
    "rest_pre": "baseline rest",
    "task_learn": "learn adjacent premises",
    "task_test": "test non-adjacent inferences",
    "rest_post": "post-task rest",
}
# wire-ready: fill with confirmed minutes (Ricci 2023 / collaborator) -> renders "~X min".
# rest is ~10 min; leave task_* as None until sourced (do not invent).
DURATIONS = {ph: None for ph in PHASES}

# ---- geometry ------------------------------------------------------------------
BW, BH, GAP = 2.6, 1.3, 1.15
XS = [i * (BW + GAP) for i in range(4)]          # left edges
CX = [x + BW / 2 for x in XS]                     # centres
RIGHT = XS[-1] + BW


def _block(ax, x, kind):
    fill, edge = (REST_FILL, REST_EDGE) if kind == "rest" else (TASK_FILL, TASK_EDGE)
    ax.add_patch(FancyBboxPatch(
        (x, 0), BW, BH, boxstyle="round,pad=0.02,rounding_size=0.12",
        linewidth=2.0, edgecolor=edge, facecolor=fill, zorder=2))


def _chevron(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch(
        (x0, y), (x1, y), arrowstyle="-|>", mutation_scale=18,
        color=DARK, lw=2.2, shrinkA=0, shrinkB=0, zorder=3))


def _bracket(ax, x0, x1, y, label):
    tick = 0.13
    ax.plot([x0, x1], [y, y], color=TASK_EDGE, lw=2.0, zorder=3)
    ax.plot([x0, x0], [y, y - tick], color=TASK_EDGE, lw=2.0, zorder=3)
    ax.plot([x1, x1], [y, y - tick], color=TASK_EDGE, lw=2.0, zorder=3)
    ax.text(0.5 * (x0 + x1), y + 0.13, label, ha="center", va="bottom",
            color=AMBER_TXT, fontsize=12, fontweight="bold")


def draw():
    fig, ax = plt.subplots(figsize=(11.5, 3.5))
    ax.set_aspect("equal")
    ax.axis("off")

    for x, ph in zip(XS, PHASES):
        _block(ax, x, KIND[ph])
        cx = x + BW / 2
        dur = DURATIONS[ph]
        ty = BH / 2 + (0.14 if dur else 0.0)                 # nudge up if a duration sits below
        ax.text(cx, ty, ph, ha="center", va="center", color=DARK,
                fontsize=13.5, fontweight="bold")
        if dur:
            ax.text(cx, BH / 2 - 0.30, f"~{dur} min", ha="center", va="center",
                    color=GRAY, fontsize=9.5)
        ax.text(cx, -0.30, DESC[ph], ha="center", va="top", color=GRAY,
                fontsize=10.5, style="italic")

    # chevrons in the gaps (ordering)
    for i in range(3):
        gl, gr = XS[i] + BW, XS[i + 1]
        _chevron(ax, gl + 0.20, gr - 0.20, BH / 2)

    # task bracket over the two middle blocks
    _bracket(ax, XS[1], XS[2] + BW, BH + 0.30, "transitive-inference task")

    # time axis underneath
    y_axis = -1.05
    ax.add_patch(FancyArrowPatch(
        (-0.3, y_axis), (RIGHT + 0.45, y_axis), arrowstyle="-|>", mutation_scale=20,
        color=DARK, lw=2.0, shrinkA=0, shrinkB=0, zorder=1))
    ax.text(0.5 * RIGHT, y_axis - 0.34, "time", ha="center", va="top",
            color=GRAY, fontsize=10.5, style="italic")

    ax.set_xlim(-0.6, RIGHT + 0.7)
    ax.set_ylim(-1.75, BH + 0.85)
    return fig


def main():
    fig = draw()
    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "ti_protocol_timeline.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    print("wrote", out)

    if "--qa" in sys.argv:
        qa = sys.argv[sys.argv.index("--qa") + 1]
        fig.savefig(qa, facecolor="white", bbox_inches="tight", dpi=130)
        print("wrote QA", qa)
    plt.close(fig)


if __name__ == "__main__":
    main()
