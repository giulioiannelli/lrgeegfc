#!/usr/bin/env python3
r"""Per-patient own-null cophenetic-trace forest across all six bands (companion to fig:trace_a).

Reference style (data/preprint/figures/all_bands/patients_forest_plot) adapted to the
rho_sym estimator + house palette. One panel per band, frequency-ordered; each of the
10 patients is a row with its OWN matched-strength null (grey p5-p95 bar), a stem to
its rho_sym, and a dot: filled = clears its own null, open = within null, red = anti.
The per-band verdict glyph (corner) marks cohort-wide trace / present-but-split / absent.

This is the rigorous per-patient view behind fig:trace_a — here color and null are
always consistent because every patient carries its own bar.

Reads : data/audit/rho_sym_gate/per_patient_per_band.csv
        data/audit/rho_sym_gate/cohort_summary.csv
Writes: data/preprint/figures/results_section1/fig_perpatient_bands_forest.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import band_marker, use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

GATE = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
SUMMARY = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
OUT = ROOT / "data/preprint/figures/results_section1/fig_perpatient_bands_forest.pdf"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

C_TRACE, C_ANTI, C_NULL = "#3a9a4f", "#d1352b", "#c2c2c2"
C_SPLIT, C_ABSENT = "#d99a1f", "#9a9a9a"


def _verdict(p):
    if p < 0.05:
        return "clears"
    if p > 0.95:
        return "anti"
    return "within"


def _band_glyph(p):
    if p < 0.05:
        return "full", C_TRACE          # cohort-wide trace
    if p < 0.15:
        return "left", C_SPLIT          # present but split
    return "none", C_ABSENT             # absent


def draw(target):
    g = pd.read_csv(GATE)
    s = pd.read_csv(SUMMARY).set_index("band")
    ymap = {p: len(COHORT) - 1 - i for i, p in enumerate(COHORT)}

    # sizes tuned for the compound's ~3-inch-wide tile (the assembler does NOT rescale
    # this panel): small dots so each stem to zero stays visible, legible patient labels,
    # per-band marker so a stray dot is never confounded, and a compact legend that stays
    # inside the panel width (it used to spill left over panel c).
    gs = target.add_gridspec(2, 3, left=0.105, right=0.985, top=0.965, bottom=0.185,
                             hspace=0.10, wspace=0.10)
    axes = gs.subplots(sharex=True, sharey=True)

    for ax, band in zip(axes.ravel(), BAND_ORDER):
        gb = g[g.band == band].set_index("patient")
        mk = band_marker(band)
        ax.axvline(0.0, color="0.55", lw=0.7, ls="--", zorder=1)
        for pat in COHORT:
            if pat not in gb.index:
                continue
            r = gb.loc[pat]
            y = ymap[pat]
            ax.plot([r.surr_p5, r.surr_p95], [y, y], color=C_NULL, lw=2.8,
                    solid_capstyle="round", zorder=2)                       # own MS null p5-p95
            v = _verdict(r.obs_p_one_sided)
            # colour follows the tree direction: positive=green, negative=red (an open red
            # marker = leaned anti but n.s.)
            sign_col = C_TRACE if r.obs_rho >= 0 else C_ANTI
            ax.plot([0.0, r.obs_rho], [y, y], color=sign_col, lw=1.0, zorder=3)
            if v == "within":
                ax.scatter(r.obs_rho, y, s=15, marker=mk, facecolor="white",
                           edgecolor=sign_col, linewidth=0.9, zorder=4)
            else:
                fill = C_TRACE if v == "clears" else C_ANTI
                ax.scatter(r.obs_rho, y, s=15, marker=mk, color=fill,
                           edgecolor="white", linewidth=0.4, zorder=4)
        # band label INSIDE each axis, lower-right, colour-coded: green=cohort-wide trace,
        # amber=present-but-split, grey=absent
        _, gcol = _band_glyph(s.loc[band, "gate_p_sym"])
        ax.text(0.965, 0.05, BRAIN_BAND_TEX_DICT[band], transform=ax.transAxes,
                fontsize=11, color=gcol, fontweight="bold", ha="right", va="bottom")

    for ax in axes.ravel():
        ax.set_xlim(-0.55, 1.0)
        ax.set_ylim(-0.6, len(COHORT) - 0.4)
    for ax in axes[:, 0]:
        ax.set_yticks(list(ymap.values()))
        ax.set_yticklabels([p.replace("Pat_", "") for p in COHORT], fontsize=7.5)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\rho^{\mathrm{coph}}$", fontsize=8)
    for ax in axes.ravel():
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=5.8, pad=1)

    handles = [
        Line2D([0], [0], color=C_NULL, lw=2.8, solid_capstyle="round",
               label="matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=4.5,
               label="clears null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec=C_TRACE, mew=1.0, ms=4.5,
               label="within null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_ANTI, mec="white", ms=4.5,
               label="anti"),
    ]
    # single row along the very bottom, clear of the x-axis labels above it
    target.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.012),
                  ncol=4, frameon=False, fontsize=6.0, handletextpad=0.35,
                  columnspacing=0.9, labelspacing=0.3)


def main():
    # portrait, matched to the compound's right-column tile so the standalone previews
    # exactly what the assembler embeds (this panel is authored at tile size, not rescaled).
    fig = plt.figure(figsize=(3.4, 4.7))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
