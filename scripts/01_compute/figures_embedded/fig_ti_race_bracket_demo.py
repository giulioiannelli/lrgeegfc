"""
I-1 cold-open — transitive-inference DISC-GOLF race bracket, round-by-round (talk schematic).

Sport = disc golf (1-v-1, so an individual ranking makes sense). A knockout bracket
over the four Lipari guests, re-seeded so BOTH first-round matchups are NON-adjacent
(inferred) pairs of the order FB > AG > DG > SM, and the FINAL is the one adjacent
result you were actually SHOWN:

    Semifinal 1:  FB vs DG -> FB   INFERRED (FB>AG>DG, never shown)      [amber]
    Semifinal 2:  AG vs SM -> AG   INFERRED (AG>DG>SM, never shown)      [amber]
    Final:        FB vs AG -> FB   SEEN     (FB>AG is a shown result)    [gold]

Champion = FB (Battiston). The whole tournament is decided by two inferences and one
seen fact — the point of the cold open.

THREE ROUND-STATES (one PDF each) so the deck builds it in the Canva reveal:
    round0  full bracket drawn in GRAY, unresolved, semis tagged "inferred"  ("who wins these matches?")
    round1  both semifinals resolve -> FB, AG advance (amber); final still gray  ("who wins the final?")
    round2  final resolves -> champion FB ★ (gold); tag "seen (encoded): FB > AG"
The whole skeleton is gray at every round so the shape is anticipated; resolved matches
are over-drawn in colour (STRAIGHT lines). Winners advance as photo-placeholder chips
(initials) — drop the headshots in Canva, same as the leaves.

Dark ink on a LIGHT slide: black text, dark-gray skeleton, amber/gold accents. Transparent PDF.
Output: data/outputs/figures/talk/ti_race_bracket_round{0,1,2}.pdf  (+ ti_race_bracket.pdf = round2)
QA (white-matte PNG, opt-in): python fig_ti_race_bracket_demo.py --qa /path/to/qa_prefix
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

# ---- palette (light slide, dark ink) --------------------------------------
DARK = "#1a1a1a"     # names, box borders, prompts, tags, champion label
GRAY = "#6f757c"     # unresolved bracket skeleton + placeholder initial hints
AMBER = "#e08d00"    # INFERRED matchup (both semifinals: never shown)
GOLD = "#f2c14e"     # SEEN matchup (the final) + champion star

# ---- geometry -------------------------------------------------------------
# leaves left->right = FB, DG, AG, SM  (matchups pair NON-adjacent ranks)
LEAF_X = [1.0, 3.0, 5.0, 7.0]
LEAF_NAMES = ["Battiston", "Garlaschelli", "Gabrielli", "Meloni"]
LEAF_INIT = ["FB", "DG", "AG", "SM"]
BW = BH = 1.4                     # photo-placeholder box (square headshot)
Y_SEMI = 3.0                      # semifinal bar height
Y_FINAL = 4.7                     # final bar height
Y_CHAMP = 5.55                    # champion chip baseline
MID1 = 0.5 * (LEAF_X[0] + LEAF_X[1])      # 2.0  (SF1 winner column)
MID2 = 0.5 * (LEAF_X[2] + LEAF_X[3])      # 6.0  (SF2 winner column)
CHAMP_X = 0.5 * (MID1 + MID2)             # 4.0
Y_SEMI_CHIP = Y_SEMI - 0.62               # advancing-winner chip below each semi bar


# ---- drawing primitives ---------------------------------------------------
def _line(ax, p0, p1, color, lw, z=1):
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw,
            solid_capstyle="round", zorder=z)


def _U(ax, xL, xR, y0, ybar, color, lw, z=1):
    """Inverted-U connector: two verticals rising to a shared horizontal bar."""
    _line(ax, (xL, y0), (xL, ybar), color, lw, z)
    _line(ax, (xR, y0), (xR, ybar), color, lw, z)
    _line(ax, (xL, ybar), (xR, ybar), color, lw, z)


def _semis(ax, color, lw, z=1):
    _U(ax, LEAF_X[0], LEAF_X[1], BH, Y_SEMI, color, lw, z)
    _U(ax, LEAF_X[2], LEAF_X[3], BH, Y_SEMI, color, lw, z)


def _final(ax, color, lw, z=1):
    _U(ax, MID1, MID2, Y_SEMI, Y_FINAL, color, lw, z)
    _line(ax, (CHAMP_X, Y_FINAL), (CHAMP_X, Y_CHAMP - 0.34), color, lw, z)


def _leaf_boxes(ax):
    for x, name, ini in zip(LEAF_X, LEAF_NAMES, LEAF_INIT):
        img = load_face(ini)
        patch = FancyBboxPatch(
            (x - BW / 2, 0), BW, BH,
            boxstyle="round,pad=0.02,rounding_size=0.12",
            linewidth=1.8, edgecolor=DARK,
            facecolor=("none" if img is not None else "white"), zorder=2)
        ax.add_patch(patch)
        if img is not None:
            place_face(ax, img, x, BH / 2, BW, patch, zorder=1.7)
        else:
            ax.text(x, BH / 2, ini, ha="center", va="center", color=GRAY,
                    fontsize=13, alpha=0.75, style="italic")        # photo hint
        ax.text(x, -0.30, name, ha="center", va="top", color=DARK, fontsize=12)


def _chip(ax, x, y, ini, edge, size=0.64, fs=12.5):
    """Advancing-winner chip — the guest's headshot behind a coloured border."""
    img = load_face(ini)
    patch = FancyBboxPatch(
        (x - size / 2, y - size / 2), size, size,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=2.2, edgecolor=edge,
        facecolor=("none" if img is not None else "white"), zorder=4)
    ax.add_patch(patch)
    if img is not None:
        place_face(ax, img, x, y, size, patch, zorder=3.7)
    else:
        ax.text(x, y, ini, ha="center", va="center", color=DARK,
                fontsize=fs, fontweight="bold", zorder=5)


def _prompt(ax, text):
    ax.text(CHAMP_X, Y_CHAMP + 0.10, text, ha="center", va="center", color=DARK,
            fontsize=12.5, style="italic", linespacing=1.1)


# ---- round-by-round figure ------------------------------------------------
def draw(rnd):
    fig, ax = plt.subplots(figsize=(9.5, 6.6))
    ax.set_aspect("equal")
    ax.axis("off")

    _leaf_boxes(ax)
    _semis(ax, GRAY, 1.8, z=1)                              # full skeleton, gray
    _final(ax, GRAY, 1.8, z=1)

    if rnd == 0:                                            # all unresolved
        for xx in (MID1, MID2):
            ax.text(xx, Y_SEMI_CHIP, "?", ha="center", va="center", color=GRAY,
                    fontsize=15)
        _prompt(ax, "who wins\nthese matches?")

    if rnd >= 1:                                            # semifinals resolved
        _semis(ax, AMBER, 3.0, z=3)
        _chip(ax, MID1, Y_SEMI_CHIP, "FB", AMBER)
        _chip(ax, MID2, Y_SEMI_CHIP, "AG", AMBER)
        for xx in (MID1, MID2):
            ax.text(xx, Y_SEMI_CHIP - 0.58, "★ inferred", ha="center", va="top",
                    color=DARK, fontsize=9.0)

    if rnd == 1:                                            # final still unresolved
        _prompt(ax, "who wins\nthe final?")

    if rnd >= 2:                                            # final resolved + champion
        _final(ax, GOLD, 3.0, z=3)
        ax.text(CHAMP_X, 0.5 * (Y_SEMI + Y_FINAL), "seen (encoded):\nFB > AG",
                ha="center", va="center", color=DARK, fontsize=9.5,
                fontweight="bold", linespacing=1.15)
        ax.plot(CHAMP_X, Y_CHAMP + 0.70, marker="*", markersize=20, color=GOLD,
                markeredgecolor=DARK, markeredgewidth=0.7, zorder=6)
        _chip(ax, CHAMP_X, Y_CHAMP, "FB", GOLD, size=0.78, fs=15)

    ax.set_xlim(-0.2, 8.2)
    ax.set_ylim(-0.95, Y_CHAMP + 1.25)
    return fig


def main():
    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    qa_prefix = sys.argv[sys.argv.index("--qa") + 1] if "--qa" in sys.argv else None

    for rnd in range(3):
        fig = draw(rnd)
        out = outdir / f"ti_race_bracket_round{rnd}.pdf"
        fig.savefig(out, transparent=True, bbox_inches="tight")
        if rnd == 2:                                        # alias: full bracket
            fig.savefig(outdir / "ti_race_bracket.pdf", transparent=True,
                        bbox_inches="tight")
        if qa_prefix:
            fig.savefig(f"{qa_prefix}_round{rnd}.png", facecolor="white",
                        bbox_inches="tight", dpi=130)
        plt.close(fig)
        print(f"wrote {out.name}")
    print("wrote ti_race_bracket.pdf (= round2)")


if __name__ == "__main__":
    main()
