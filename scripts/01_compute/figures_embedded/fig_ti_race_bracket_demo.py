"""
I-1 cold-open — transitive-inference RACE BRACKET (talk schematic).

A knockout bracket over the four Lipari guests whose first-round matchups are
the NON-adjacent (inferred) pairs of the order FB > AG > DG > SM:

    Semifinal 1 (HERO):  Battiston vs Garlaschelli -> Battiston  (inferred FB>AG>DG)
    Semifinal 2:         Gabrielli vs Meloni       -> Gabrielli  (inferred AG>DG>SM)
    Final:               Battiston vs Gabrielli    -> Battiston  (shown FB>AG)

Champion = Battiston; Meloni out in round 1. The accented left semifinal
(Battiston vs Garlaschelli) is the "you were never shown this" inference.

The leaf boxes are PHOTO PLACEHOLDERS — drop the four headshots in Canva.
Vector PDF, transparent background, light ink → sits on a DARK slide.

Output: data/outputs/figures/talk/ti_race_bracket.pdf
QA (dark-matte PNG, opt-in): python fig_ti_race_bracket_demo.py --qa /path/to/qa.png
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

use_lrg_style()

# ---- palette (dark-slide, light ink) --------------------------------------
INK = "#e6edf3"      # bracket lines, box borders, neutral names
MUTED = "#8b949e"    # placeholder initials, secondary text
ACCENT = "#e08d00"   # HERO node: the inferred, never-shown matchup (FB vs DG)
GOLD = "#f2c14e"     # champion star

# ---- geometry -------------------------------------------------------------
# leaves left->right = FB, DG, AG, SM  (matchups pair NON-adjacent ranks)
LEAF_X = [1.0, 3.0, 5.0, 7.0]
LEAF_NAMES = ["Battiston", "Garlaschelli", "Gabrielli", "Meloni"]
LEAF_INIT = ["FB", "DG", "AG", "SM"]
BW = BH = 1.4                     # photo-placeholder box (square headshot)
Y_SEMI = 3.0                      # semifinal bar height
Y_FINAL = 4.7                     # final bar height
Y_CHAMP = 5.15                    # champion label baseline

fig, ax = plt.subplots(figsize=(9.5, 6.2))
ax.set_aspect("equal")
ax.axis("off")


def bracket(xL, xR, y_from, y_bar, color, lw):
    """Inverted-U connector: two verticals rising to a shared horizontal bar."""
    ax.plot([xL, xL], [y_from, y_bar], color=color, lw=lw, solid_capstyle="round")
    ax.plot([xR, xR], [y_from, y_bar], color=color, lw=lw, solid_capstyle="round")
    ax.plot([xL, xR], [y_bar, y_bar], color=color, lw=lw, solid_capstyle="round")


# --- bracket lines ---------------------------------------------------------
bracket(LEAF_X[0], LEAF_X[1], BH, Y_SEMI, ACCENT, 2.8)   # semi 1 (HERO)
bracket(LEAF_X[2], LEAF_X[3], BH, Y_SEMI, INK, 2.0)      # semi 2
mid1 = 0.5 * (LEAF_X[0] + LEAF_X[1])                     # 2.0
mid2 = 0.5 * (LEAF_X[2] + LEAF_X[3])                     # 6.0
bracket(mid1, mid2, Y_SEMI, Y_FINAL, INK, 2.0)           # final
champ_x = 0.5 * (mid1 + mid2)                            # 4.0
ax.plot([champ_x, champ_x], [Y_FINAL, Y_CHAMP - 0.10], color=INK, lw=2.0)

# --- leaf photo-placeholder boxes ------------------------------------------
for x, name, ini in zip(LEAF_X, LEAF_NAMES, LEAF_INIT):
    ax.add_patch(FancyBboxPatch(
        (x - BW / 2, 0), BW, BH,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        linewidth=1.8, edgecolor=INK, facecolor="none"))
    ax.text(x, BH / 2, ini, ha="center", va="center", color=MUTED,
            fontsize=13, alpha=0.5, style="italic")            # placeholder hint
    ax.text(x, -0.30, name, ha="center", va="top", color=INK, fontsize=12)

# --- winner labels (inside each U, so the rising final lines stay clean) ----
ax.text(mid1, Y_SEMI - 0.34, "Battiston", ha="center", va="top",
        color=ACCENT, fontsize=12.5, fontweight="bold")       # semi-1 winner (hero)
ax.text(mid1, Y_SEMI - 0.74, "★ inferred", ha="center", va="top",
        color=ACCENT, fontsize=9.0)                            # never-shown tag
ax.text(mid2, Y_SEMI - 0.34, "Gabrielli", ha="center", va="top",
        color=INK, fontsize=12)                               # semi-2 winner

# --- champion --------------------------------------------------------------
ax.plot(champ_x, Y_CHAMP + 0.74, marker="*", markersize=18,
        color=GOLD, markeredgecolor=INK, markeredgewidth=0.6)
ax.text(champ_x, Y_CHAMP, "Battiston", ha="center", va="bottom",
        color=INK, fontsize=15, fontweight="bold")

ax.set_xlim(-0.2, 8.2)
ax.set_ylim(-0.95, Y_CHAMP + 1.25)

out = FIGURES_ROOT / "talk" / "ti_race_bracket.pdf"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, transparent=True, bbox_inches="tight")
print("wrote", out)

if "--qa" in sys.argv:                                        # dark-matte QA only
    qa = sys.argv[sys.argv.index("--qa") + 1]
    fig.savefig(qa, facecolor="#0b0e17", bbox_inches="tight", dpi=130)
    print("wrote QA", qa)

plt.close(fig)
