#!/usr/bin/env python3
r"""fig:enc_anatomy — the §2 "anatomy" compound (Fig 4): three region COLUMNS.

COMPOUND, built entirely in matplotlib (no LaTeX minipage). Three region columns side by
side (landscape); within each column the lollipop (system enrichment) sits ABOVE its glass
brain, plus ONE shared bottom legend:

  a  encoding component anchors in orbitofrontal cortex, beta          -- fig_arc_c_...
  b  a focal low-gamma memory trace on cingulate the whole brain misses -- fig_arc_e_...
  c  the inference-specific component leans to the cingulate (lead)     -- fig_arc_f_...

Each column is drawn by the single-figure scripts' own helpers (lollipop + nilearn glass
brain), imported here and composed onto ONE flat matplotlib Figure; nothing is re-derived.
Flat (not SubFigures) because nilearn.plot_glass_brain resolves `figure=` to the root
Figure and lays its brain axes in ROOT coordinates -- a SubFigure's local frame is ignored,
so every brain would stack at one rect. Rects here are therefore in figure coordinates.
Columns share one legend (they share the enriched / depleted / n.s. / excluded vocabulary);
the column-specific meaning is on each column's caption. Lollipop fonts are shrunk here so
the long system labels and axis title fit the narrower column; the single-figure scripts
keep their own (wider) sizes.

Reads : (delegated) data/audit/inference_localization_rhosym/*
Writes: data/preprint/figures/results_section2/fig_arc_enc_anatomy.pdf
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import band_color, use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
import fig_arc_c_encoding_ofc_anchor as arc_c  # noqa: E402
import fig_arc_e_lowgamma_cingulate as arc_e  # noqa: E402
import fig_arc_f_inference_cingulate as arc_f  # noqa: E402

OUT = ROOT / "data/preprint/figures/results_section2/fig_arc_enc_anatomy.pdf"

# three columns (figure coords, 0=left): a | b | c. Within a column the lollipop is on top,
# the glass brain below, and a coloured caption under the brain; one legend spans the bottom.
LM, RM = 0.028, 0.995
COLW = (RM - LM) / 3.0
LOLLI_BOT, LOLLI_H = 0.555, 0.370            # lollipop band (top); top=0.925 clears the chip
LOLLI_LPAD, LOLLI_RPAD = 0.058, 0.010        # left pad = room for system labels
BRAIN_BOT, BRAIN_H = 0.220, 0.285            # glass-brain band (bottom), snug to content
CAP_Y, LET_Y, CHIP_Y = 0.205, 0.980, 0.962   # caption under the brain; band chip top-left
LEG_Y = 0.140                                # legend just under the captions (tight crop below)
LOLLI_XLAB_FS, LOLLI_YTICK_FS = 9.0, 9.5     # shrunk to fit the narrower column

# each column is a SINGLE band's localization (not a per-band sweep): the y-axis lists
# anatomical SYSTEMS, so the band is stated as a chip at the column top, coloured by the
# canonical band palette (beta = green, low-gamma = blue, same identity as Fig 3).
BAND_OF_COL = {"a": "beta", "b": "low_gamma", "c": "beta"}
BAND_HZ = {"beta": "13-30 Hz", "low_gamma": "30-80 Hz"}


def _col(i):
    """Figure-coord rects for column i (0=a, 1=b, 2=c)."""
    colL = LM + i * COLW
    lolli = [colL + LOLLI_LPAD, LOLLI_BOT, COLW - LOLLI_LPAD - LOLLI_RPAD, LOLLI_H]
    brain = (colL + 0.003, BRAIN_BOT, COLW - 0.010, BRAIN_H)
    return dict(lolli=lolli, brain=brain, cx=colL + 0.5 * COLW, letx=colL + 0.004)


def _fit_lollipop(ax):
    """Shrink the lollipop's axis-title + system labels for the narrow column."""
    ax.xaxis.label.set_size(LOLLI_XLAB_FS)
    ax.tick_params(axis="y", labelsize=LOLLI_YTICK_FS)


def _letter(fig, ch, letx):
    fig.text(letx, LET_Y, rf"$\mathbf{{{ch}}}$", fontsize=20, va="top", ha="left",
             fontweight="bold")
    band = BAND_OF_COL[ch]
    txt = f"{BRAIN_BAND_TEX_DICT[band]}  {BAND_HZ[band]}"
    fig.text(letx + 0.026, CHIP_Y, txt, fontsize=12.5, va="center", ha="left",
             fontweight="bold", color=band_color(band, 0.70),
             bbox=dict(boxstyle="round,pad=0.30",
                       facecolor=to_rgba(band_color(band, 1.0), 0.15),
                       edgecolor=band_color(band, 0.70), linewidth=1.1))


def draw_col_c(fig, col, coords, systems):
    """a — encoding -> OFC anchor (beta)."""
    inc, exc = arc_c.load_sys("include"), arc_c.load_sys("exclude")
    ax = fig.add_axes(col["lolli"])
    arc_c.draw_lollipop(ax, inc, exc)
    _fit_lollipop(ax)
    if arc_c._BRAIN_OK:
        arc_c.draw_brain(fig, col["brain"], coords, systems, inc)
        fig.text(col["cx"], CAP_Y, r"encoding $\rightarrow$ OFC", ha="center",
                 va="top", fontsize=12.5, color=arc_c.C_ENR, fontweight="bold")
    _letter(fig, "a", col["letx"])


def draw_col_e(fig, col, coords, systems):
    """b — focal low-gamma memory trace on cingulate."""
    inc, exc = arc_e.load("include"), arc_e.load("exclude")
    sig_win = bool(inc.loc[arc_e.WINNER, "bh_q"] < arc_e.Q_SIG)
    ax = fig.add_axes(col["lolli"])
    arc_e.draw_lollipop(ax, inc, exc)
    _fit_lollipop(ax)
    if arc_e._BRAIN_OK:
        arc_e.draw_brain(fig, col["brain"], coords, systems, sig_win)
        fig.text(col["cx"], CAP_Y, r"memory (encoding) $\rightarrow$ cingulate, low-$\gamma$",
                 ha="center", va="top", fontsize=11.5, color=arc_e.C_WIN, fontweight="bold")
    _letter(fig, "b", col["letx"])


def draw_col_f(fig, col, coords, systems):
    """c — inference-specific component concentrates in the cingulate (beta)."""
    inc, exc = arc_f.load_sys("include"), arc_f.load_sys("exclude")
    ax = fig.add_axes(col["lolli"])
    arc_f.draw_lollipop(ax, inc, exc)
    _fit_lollipop(ax)
    if arc_f._BRAIN_OK:
        arc_f.draw_brain(fig, col["brain"], coords, systems, inc)
        fig.text(col["cx"], CAP_Y, r"inference $\rightarrow$ cingulate",
                 ha="center", va="top", fontsize=12.5, color=arc_f.C_INF, fontweight="bold")
    _letter(fig, "c", col["letx"])


def shared_legend(fig):
    # one row, palette-correct: each column's positive marker gets its own colour
    # (a green encoding / b blue low-gamma / c red-orange inference), plus depleted /
    # n.s. / replicate. No provisional diamond (the length demotion was withdrawn).
    handles = [
        Line2D([], [], color=arc_c.C_ENR, marker="o", ls="", mfc=arc_c.C_ENR,
               mec="black", mew=1.4, ms=11,
               label=r"enriched (encoding), BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color=arc_e.C_WIN, marker="o", ls="", mfc=arc_e.C_WIN,
               mec="black", mew=1.2, ms=11,
               label=r"within-region trace (low-$\gamma$), $q<0.05$ ($\ast$)"),
        Line2D([], [], color=arc_f.C_INF, marker="o", ls="", mfc=arc_f.C_INF,
               mec="black", mew=1.4, ms=11,
               label=r"inference concentration, $q<0.05$ ($\ast$)"),
        Line2D([], [], color=arc_c.C_DEP, marker="o", ls="", mfc=arc_c.C_DEP,
               mec="black", mew=1.1, ms=10, label=r"depleted, BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color="0.4", marker="o", ls="", mfc="white", mec="0.4",
               mew=1.4, ms=9, label=r"n.s. ($q\geq0.05$)"),
        Line2D([], [], color="0.4", marker="o", ls="", mfc="none", mec="0.4",
               mew=0.9, ms=7, alpha=0.7, label="epilepsy / SOZ-excluded replicate"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, LEG_Y),
               ncol=6, frameon=False, fontsize=9.3, handletextpad=0.5, columnspacing=1.3)


def main():
    coords, systems = (arc_c.pooled_coords_systems() if arc_c._BRAIN_OK else (None, None))

    fig = plt.figure(figsize=(18.0, 7.4))
    draw_col_c(fig, _col(0), coords, systems)
    draw_col_e(fig, _col(1), coords, systems)
    draw_col_f(fig, _col(2), coords, systems)
    shared_legend(fig)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True, pad_inches=0.05)
    plt.close(fig)
    print("fig:enc_anatomy — compound (3 columns: encoding OFC | low-gamma cingulate | "
          "inference lead), landscape")
    print(f"  brain machinery: {'ON' if arc_c._BRAIN_OK else 'OFF (lollipops only)'}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
