#!/usr/bin/env python3
r"""Talk slide 06 — the Ellenbogen 2007 bridge: ONSET (ours) vs MATURATION (theirs).

A two-lane "level-of-analysis bridge", NOT a shared axis. The earlier version pinned our
recording WINDOW onto Ellenbogen's 20-min TEST bar — welding a continuous neural read
(a span, no % coordinate) onto a between-subjects delayed behavioural result (a value on a
'% correct' axis), which hid task_test and manufactured a brain-behaviour comparison the
project forbids. This redraw fixes that by keeping the two studies in physically separate
rows joined only by a labelled relationship.

Layout — a broken "time after learning" x-axis (minutes | hours):
  TOP    = THEIR behaviour (Ellenbogen & Walker 2007, PNAS 104:7723, healthy adults, n=56).
           % correct axis. Premises are held at once (~90%); the INFERRED order is at chance
           at the 20-min test (52%) and only becomes significant after hours (75%); sleep
           lifts the hardest pair B>E (69% wake -> 93% sleep). Their offline delay is drawn
           EMPTY (no test given) — the binding is a black box, bracketed only by a later test.
  BOTTOM = OUR protocol (intracranial sEEG, epilepsy patients, n=10). A schematic with NO
           y-quantity: rest_pre -> task_learn -> task_test -> REST_POST. task_test is drawn
           explicitly BEFORE rest_post (the patient has already attempted the inference before
           the window we analyse); rest_post is the teal read-window (wake, minutes) where the
           network trace begins to form. Our lane ENDS at rest_post — the hours/sleep region
           under their matured effect is deliberately empty ("beyond our recorded window").

The one honest bridge: same offline period, two vantage points — they read its OUTPUT by
testing AFTER it; we record the PROCESS beginning INSIDE it. No shared subjects, no shared
axis, no brain-behaviour correlation. Ellenbogen is external motivation, not our data.

House rules: use_lrg_style(), transparent PDF (+ .png sibling for the Canva talk export),
no suptitle, no rasterisation, plt.close. Light deck -> dark ink. NOT the band palette
(schematic + external behavioural data, not band-indexed).

Usage:
    /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python \
        scripts/07_figures/talk_fig_ellenbogen_offline_inference.py
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.lines import Line2D

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "talk"

# ---- palette (hardcoded by design — this is schematic + external data, NOT band-indexed) --
INK = "#1f2a37"          # structure / dark ink for a light deck
CHANCE = "#9aa0a8"       # chance reference
C_SHOWN = "#3a6ea5"      # premises (what was SHOWN) — slate blue [THEIR data only]
C_INFER = "#e0713a"      # inference (what was WORKED OUT) — coral   [THEIR data only]
C_SLEEP = "#c8531f"      # the sleep-dependent B>E endpoint — deep coral [THEIR data only]
C_OURS = "#0f766e"       # OUR neural read window — deep teal [OUR lane only]
BAND = "#e7b24d"         # shared offline-period band — pale amber (very low alpha)
REST_F, REST_E = "#dce4ec", "#5b6b7a"   # house rest-block fill / edge
TASK_F, TASK_E = "#f7e4bb", "#b9832a"   # house task-block fill / edge
TEAL_F = "#d6ece9"       # light teal fill for the rest_post read-window

# ---- reported numbers (Ellenbogen & Walker 2007) ------------------------------------------
PREM_20, INF_20 = 90.0, 52.0            # 20-min group: premises held, inference at chance
PREM_H, INF_H = 89.0, 75.0              # 12/24 h: premises held, inference significant
BE_WAKE, BE_SLEEP = 69.0, 93.0          # hardest pair B>E, 12 h awake vs asleep
CHANCE_LVL = 50.0


def _break_marks(ax, side):
    """Draw a small '//' axis-break glyph at the inner (broken) edge of *ax*."""
    x = 1.0 if side == "right" else 0.0
    for dx in (-0.007, 0.007):
        ax.plot([x + dx - 0.009, x + dx + 0.009], [0.5 - 0.05, 0.5 + 0.05],
                transform=ax.transAxes, color=INK, lw=1.1, clip_on=False, zorder=30)


def _block(ax, x0, x1, yc, h, face, edge, label, tcolor, *, lw=1.6, fs=8.6):
    ax.add_patch(Rectangle((x0, yc - h / 2), x1 - x0, h, facecolor=face, edgecolor=edge,
                           lw=lw, zorder=4))
    ax.text((x0 + x1) / 2, yc, label, ha="center", va="center", fontsize=fs, color=tcolor,
            fontweight="bold", zorder=6, linespacing=1.0)


def _chevron(ax, x, yc):
    ax.annotate("", xy=(x + 0.7, yc), xytext=(x - 0.7, yc),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3), zorder=5)


def _sprout(ax, x0, base):
    """Tiny emerging dendrogram inside rest_post — the qualitative trace-onset glyph."""
    xs = [x0, x0 + 1.1, x0 + 2.2]
    top = base + 0.03
    for xx in xs:
        ax.plot([xx, xx], [base, top], color=C_OURS, lw=1.3, zorder=7)
    h1 = base + 0.10
    ax.plot([xs[0], xs[0]], [top, h1], color=C_OURS, lw=1.3, zorder=7)
    ax.plot([xs[1], xs[1]], [top, h1], color=C_OURS, lw=1.3, zorder=7)
    ax.plot([xs[0], xs[1]], [h1, h1], color=C_OURS, lw=1.3, zorder=7)
    m1 = (xs[0] + xs[1]) / 2
    h2 = base + 0.17
    ax.plot([m1, m1], [h1, h2], color=C_OURS, lw=1.3, zorder=7)
    ax.plot([xs[2], xs[2]], [top, h2], color=C_OURS, lw=1.3, zorder=7)
    ax.plot([m1, xs[2]], [h2, h2], color=C_OURS, lw=1.3, zorder=7)


def main() -> None:
    use_lrg_style()
    fig = plt.figure(figsize=(11.0, 6.3))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.0], height_ratios=[1.5, 1.0],
                          wspace=0.06, hspace=0.62, left=0.075, right=0.985,
                          top=0.80, bottom=0.11)
    axTL = fig.add_subplot(gs[0, 0])                      # their behaviour, minutes
    axTR = fig.add_subplot(gs[0, 1], sharey=axTL)         # their behaviour, hours
    axBL = fig.add_subplot(gs[1, 0], sharex=axTL)         # our protocol, minutes
    axBR = fig.add_subplot(gs[1, 1], sharex=axTR)         # our protocol, hours

    axTL.set_xlim(-6.5, 26)
    axTR.set_xlim(8, 27)
    axTL.set_ylim(40, 100)
    for ax in (axBL, axBR):
        ax.set_ylim(0, 1)

    # ---- shared offline-period band (both rows point at ONE period) -----------------------
    # kept VERY faint so it reads as a soft backdrop, never competing with the coral/teal data
    for ax, (x0, x1) in [(axTL, (2, 26)), (axTR, (8, 27)),
                         (axBL, (2, 26)), (axBR, (8, 27))]:
        ax.axvspan(x0, x1, color=BAND, alpha=0.045, lw=0, zorder=0)

    # ===================================================================================
    # TOP-LEFT — their behaviour, the 20-min group (minutes region)
    # ===================================================================================
    for ax in (axTL, axTR):
        ax.axhline(CHANCE_LVL, color=CHANCE, lw=1.1, ls=(0, (4, 3)), zorder=2)
    axTL.text(-6.2, 51.2, "chance", color=CHANCE, fontsize=8, va="bottom")

    # their offline delay = EMPTY (no test during it — the black box)
    axTL.add_patch(Rectangle((1, 41.5), 19, 4.0, facecolor="none", hatch="////",
                             edgecolor="0.62", lw=0.8, zorder=1))
    axTL.text(10.5, 43.5, "offline delay — no test given", ha="center", va="center",
              fontsize=8, color="0.5", style="italic")

    # anchor: the only shared event
    axTL.axvline(0, color=INK, lw=1.2, zorder=3)
    axTL.text(0, 1.03, "learning\nacquired", transform=axTL.get_xaxis_transform(),
              ha="center", va="bottom", fontsize=8, color=INK, fontweight="bold",
              linespacing=1.0)

    # 20-min markers: premises held, inference at chance (HOLLOW = not yet significant)
    axTL.plot([20], [PREM_20], "o", ms=9, mfc=C_SHOWN, mec="white", mew=1.2, zorder=6)
    axTL.plot([20], [INF_20], "o", ms=9, mfc="none", mec=C_INFER, mew=1.9, zorder=6)
    axTL.text(20, PREM_20 + 3.2, "premises", ha="center", va="bottom", color=C_SHOWN,
              fontsize=9, fontweight="bold")
    axTL.text(18.6, INF_20 + 2.6, "inference\nstill at chance", ha="right", va="bottom",
              color=C_INFER, fontsize=8.3, fontweight="bold", linespacing=1.0)
    axTL.text(20, -0.055, "20 min", transform=axTL.get_xaxis_transform(), ha="center",
              va="top", fontsize=8, color=INK)

    axTL.set_ylabel("% correct")
    axTL.set_yticks([50, 75, 100])
    axTL.set_xticks([])
    for sp in ("top", "right"):
        axTL.spines[sp].set_visible(False)

    # ===================================================================================
    # TOP-RIGHT — their behaviour matures (hours region)
    # ===================================================================================
    axTR.plot([12, 24], [PREM_H, PREM_H], "-o", color=C_SHOWN, lw=1.6, ms=8, mfc=C_SHOWN,
              mec="white", mew=1.2, zorder=6)
    axTR.plot([12, 24], [INF_H, INF_H], "-o", color=C_INFER, lw=1.6, ms=8, mfc=C_INFER,
              mec="white", mew=1.2, zorder=6)
    axTR.text(24, INF_H - 2.4, "inference\nnow significant", ha="right", va="top",
              color=C_INFER, fontsize=8.3, fontweight="bold", linespacing=1.0)

    # B>E sleep split, docked at its true temporal home (12 h), offset right of the mean marker
    xbe = 13.6
    axTR.plot([xbe], [BE_WAKE], "s", ms=8, mfc="none", mec=C_SLEEP, mew=1.8, zorder=6)
    axTR.plot([xbe], [BE_SLEEP], "s", ms=8, mfc=C_SLEEP, mec="white", mew=1.0, zorder=6)
    axTR.annotate("", xy=(xbe, BE_SLEEP - 1.4), xytext=(xbe, BE_WAKE + 1.4),
                  arrowprops=dict(arrowstyle="-|>", color=C_SLEEP, lw=1.6), zorder=5)
    axTR.text(xbe + 0.7, (BE_WAKE + BE_SLEEP) / 2, "sleep\nlifts B>E", ha="left",
              va="center", color=C_SLEEP, fontsize=8, fontweight="bold", linespacing=1.0)
    axTR.text(xbe, BE_WAKE - 3.0, "B>E (hardest)", ha="center", va="top", color=C_SLEEP,
              fontsize=7.4)

    for xh, lab in [(12, "12 h"), (24, "24 h")]:
        axTR.text(xh, -0.055, lab, transform=axTR.get_xaxis_transform(), ha="center",
                  va="top", fontsize=8, color=INK)
    axTR.set_xticks([])
    axTR.tick_params(labelleft=False)
    for sp in ("top", "right", "left"):
        axTR.spines[sp].set_visible(False)

    # cross-break guides (dashed) so the rise reads as one trajectory — NEVER a solid tie
    fig.add_artist(ConnectionPatch(xyA=(20, INF_20), coordsA=axTL.transData,
                                   xyB=(12, INF_H), coordsB=axTR.transData,
                                   color=C_INFER, lw=1.3, ls="--", alpha=0.6, zorder=3))
    fig.add_artist(ConnectionPatch(xyA=(20, PREM_20), coordsA=axTL.transData,
                                   xyB=(12, PREM_H), coordsB=axTR.transData,
                                   color=C_SHOWN, lw=1.3, ls="--", alpha=0.45, zorder=3))

    # ===================================================================================
    # BOTTOM-LEFT — our protocol (minutes region), NO y-quantity
    # ===================================================================================
    # lane identity, at the very top of the bottom row (clear of the blocks)
    axBL.text(0.0, 0.99, "OUR PARADIGM   ·   epilepsy patients (n = 10)   ·   intracranial network",
              transform=axBL.transAxes, ha="left", va="top", fontsize=8.5, color=C_OURS,
              fontweight="bold")

    yc, h = 0.52, 0.26
    _block(axBL, -5.8, -2.2, yc, h, REST_F, REST_E, "rest_pre", REST_E)
    _block(axBL, -1.4, 2.6, yc, h, TASK_F, TASK_E, "task_learn", TASK_E)
    _block(axBL, 3.4, 7.6, yc, h, TASK_F, TASK_E, "task_test", TASK_E)
    # rest_post: label placed on the LEFT half so the trace glyph has the right half
    axBL.add_patch(Rectangle((9.0, yc - h / 2), 13.0, h, facecolor=TEAL_F, edgecolor=C_OURS,
                             lw=2.4, zorder=4))
    axBL.text(12.4, yc, "rest_post", ha="center", va="center", fontsize=9.4, color=C_OURS,
              fontweight="bold", zorder=6)
    for xg in (-1.8, 3.0, 8.3):
        _chevron(axBL, xg, yc)

    # anchor continues through the schematic (task_learn straddles it)
    axBL.axvline(0, color=INK, lw=1.2, zorder=3)

    # task_test precedes the window we read — the ordering the old figure hid.
    # Annotated BELOW its block (mirrors the rest_post caption) to clear the identity line.
    axBL.annotate("", xy=(5.5, yc - h / 2 - 0.01), xytext=(5.5, 0.31),
                  arrowprops=dict(arrowstyle="-|>", color=C_INFER, lw=1.4), zorder=5)
    axBL.text(5.5, 0.285, "inference attempted here\n(before our window)", ha="center",
              va="top", fontsize=7.8, color=C_INFER, fontweight="bold", linespacing=1.0)

    # emerging-structure glyph inside rest_post (right half) + honest one-liners below
    _sprout(axBL, 16.4, 0.44)
    axBL.text(19.4, 0.53, "trace\nbegins", ha="left", va="center", fontsize=7.6,
              color=C_OURS, fontweight="bold", linespacing=1.0)
    axBL.text(15.5, 0.30, "our sEEG read window · wake · minutes", ha="center", va="top",
              fontsize=8.4, color=C_OURS, fontweight="bold")
    axBL.text(15.5, 0.185, "network state — not a behavioural score", ha="center", va="top",
              fontsize=7.8, color="0.5", style="italic")

    axBL.get_yaxis().set_visible(False)
    axBL.set_xticks([])
    for sp in ("top", "right", "left", "bottom"):
        axBL.spines[sp].set_visible(False)

    # ===================================================================================
    # BOTTOM-RIGHT — our lane is BLANK here (we claim onset, not the matured endpoint)
    # ===================================================================================
    axBR.text(0.5, 0.55, "beyond our\nrecorded window", transform=axBR.transAxes,
              ha="center", va="center", fontsize=8.6, color="0.55", style="italic",
              linespacing=1.1)
    axBR.text(0.5, 0.30, "(maturation not measured)", transform=axBR.transAxes,
              ha="center", va="center", fontsize=7.6, color="0.6", style="italic")
    axBR.get_yaxis().set_visible(False)
    axBR.set_xticks([])
    for sp in ("top", "right", "left", "bottom"):
        axBR.spines[sp].set_visible(False)

    # ---- axis-break glyphs (the sole guard against reading minutes ≈ hours) ----------------
    _break_marks(axTL, "right")
    _break_marks(axTR, "left")
    _break_marks(axBL, "right")
    _break_marks(axBR, "left")

    # ---- zone headers (carry the onset-vs-maturation thesis) -------------------------------
    axTL.text(0.5, 1.20, "minutes · wake", transform=axTL.transAxes, ha="center",
              va="bottom", fontsize=8.5, color="0.45")
    axTL.text(0.5, 1.11, "ONSET — offline binding begins", transform=axTL.transAxes,
              ha="center", va="bottom", fontsize=10, color=INK, fontweight="bold")
    axTR.text(0.5, 1.20, "hours · sleep", transform=axTR.transAxes, ha="center",
              va="bottom", fontsize=8.5, color="0.45")
    axTR.text(0.5, 1.11, "MATURATION — inference expressed", transform=axTR.transAxes,
              ha="center", va="bottom", fontsize=10, color=INK, fontweight="bold")

    # ---- the ONE bridge: a labelled relationship set between the rows (never a shared axis) -
    fig.text(0.52, 0.395,
             "same offline period — they read its OUTPUT by testing AFTER it  ·  "
             "we record the PROCESS beginning INSIDE it",
             ha="center", va="center", fontsize=9.6, color=INK, fontweight="bold")
    fig.text(0.52, 0.368,
             "no shared subjects   ·   no shared axis   ·   no brain–behaviour correlation",
             ha="center", va="center", fontsize=8, color="0.5")

    # ---- lane identity + attribution -------------------------------------------------------
    fig.text(0.078, 0.905, "THEIR PARADIGM   ·   behaviour · healthy adults (n = 56)",
             ha="left", va="bottom", fontsize=8.5, color=INK, fontweight="bold")
    fig.text(0.985, 0.905,
             "Ellenbogen & Walker 2007, PNAS 104:7723 — motivation, not our data",
             ha="right", va="bottom", fontsize=7.5, color="0.5", style="italic")

    # ---- legend (3 entries) ----------------------------------------------------------------
    handles = [
        Line2D([0], [0], marker="o", color="none", mfc=C_SHOWN, mec="white", ms=8,
               label="premises (shown)"),
        Line2D([0], [0], marker="o", color="none", mfc=C_INFER, mec="white", ms=8,
               label="inference (worked out)"),
        Line2D([0], [0], marker="s", color="none", mfc=TEAL_F, mec=C_OURS, mew=1.6, ms=8,
               label="our sEEG read window"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.015), ncol=3,
               frameon=False, fontsize=9)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    pdf = OUTDIR / "slide06_ellenbogen_offline_inference.pdf"
    fig.savefig(pdf, transparent=True, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTDIR / "slide06_ellenbogen_offline_inference.png", transparent=True,
                bbox_inches="tight", pad_inches=0.04, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf}  (+ .png sibling)", flush=True)


if __name__ == "__main__":
    main()
