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
Writes: data/reports/results_section1/fig_perpatient_bands_forest.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

GATE = ROOT / "data/audit/rho_sym_gate/per_patient_per_band.csv"
SUMMARY = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
OUT = ROOT / "data/reports/results_section1/fig_perpatient_bands_forest.pdf"

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


def main():
    g = pd.read_csv(GATE)
    s = pd.read_csv(SUMMARY).set_index("band")
    ymap = {p: len(COHORT) - 1 - i for i, p in enumerate(COHORT)}

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.2), sharex=True, sharey=True)

    for ax, band in zip(axes.ravel(), BAND_ORDER):
        gb = g[g.band == band].set_index("patient")
        ax.axvline(0.0, color="0.55", lw=0.9, ls="--", zorder=1)
        for pat in COHORT:
            if pat not in gb.index:
                continue
            r = gb.loc[pat]
            y = ymap[pat]
            # own matched-strength null (p5-p95)
            ax.plot([r.surr_p5, r.surr_p95], [y, y], color=C_NULL, lw=5.5,
                    solid_capstyle="round", zorder=2)
            v = _verdict(r.obs_p_one_sided)
            # colour follows the DIRECTION of the tree: positive=green, negative=red,
            # regardless of significance (an empty red dot = leaned anti but n.s.)
            sign_col = C_TRACE if r.obs_rho >= 0 else C_ANTI
            ax.plot([0.0, r.obs_rho], [y, y], color=sign_col, lw=1.8, zorder=3)
            if v == "within":
                ax.scatter(r.obs_rho, y, s=95, facecolor="white", edgecolor=sign_col,
                           linewidth=1.7, zorder=4)
            else:
                fill = C_TRACE if v == "clears" else C_ANTI
                ax.scatter(r.obs_rho, y, s=95, color=fill, edgecolor="white",
                           linewidth=0.7, zorder=4)

        # per-band verdict glyph + label, lower-right
        fill, gcol = _band_glyph(s.loc[band, "gate_p_sym"])
        ax.plot(0.80, 0.8, marker=MarkerStyle("o", fillstyle=fill), ms=17,
                mfc=gcol, mec=gcol, mew=1.8, transform=ax.get_yaxis_transform(),
                clip_on=False, zorder=5)
        ax.text(0.93, 0.8, BRAIN_BAND_TEX_DICT[band], transform=ax.get_yaxis_transform(),
                fontsize=20, color=gcol, ha="center", va="center")

    for ax in axes.ravel():
        ax.set_xlim(-0.55, 1.0)
        ax.set_ylim(-0.6, len(COHORT) - 0.4)
    for ax in axes[:, 0]:
        ax.set_yticks(list(ymap.values()))
        ax.set_yticklabels([p.replace("Pat_", "Pat ") for p in COHORT], fontsize=11.5)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\rho^{\mathrm{coph}}$", fontsize=15)
    for ax in axes.ravel():
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=12)

    handles = [
        Line2D([0], [0], color=C_NULL, lw=5.5, solid_capstyle="round",
               label="matched-strength null (p5–p95)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec="white", ms=10,
               label="patient clears null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec=C_TRACE, mew=1.7, ms=10,
               label="patient within null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_ANTI, mec="white", ms=10,
               label="anti"),
        Line2D([0], [0], marker="o", ls="", mfc=C_TRACE, mec=C_TRACE, ms=12,
               label="band: cohort-wide trace"),
        Line2D([0], [0], marker=MarkerStyle("o", fillstyle="left"), ls="",
               mfc=C_SPLIT, mec=C_SPLIT, ms=12, label="band: present but split"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec=C_ABSENT, mew=1.7, ms=12,
               label="band: absent"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=4, frameon=False, fontsize=10.5, handletextpad=0.5, columnspacing=1.6)

    fig.subplots_adjust(left=0.07, right=0.98, top=0.97, bottom=0.16, hspace=0.12, wspace=0.08)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
