#!/usr/bin/env python3
r"""fig:trace_c — the multiscale hierarchy is a band-specific filter (Results §1, para c).

CORE MESSAGE: read through raw edges (x-axis), every band shows a similar blunt drift
back toward the task state — six bands piled into a narrow strip, theta as high as beta,
and only the edge-level gate barely separating them. Read through the LRG cophenetic
hierarchy (y-axis), that generic overlap is stripped: every band falls BELOW the diagonal
(the hierarchy sees less tracking than the raw edges claim), and only beta stays high and
close to the diagonal — its edge overlap was genuinely hierarchical. theta collapses
through zero; delta, which cleared the raw edge gate, is filtered out. The single band
that clears the cophenetic gate but NOT the raw one (beta) is REVEALED by the hierarchy;
the band that clears raw but not cophenetic (delta) is FILTERED by it.

x = raw edge overlap rho (cohort median); y = cophenetic trace rho_sym (cohort median).
Diagonal = equal tracking at both levels. Drop below diagonal = attenuation by hierarchy.
Filled marker = clears the cophenetic gate; thick ring = clears the raw edge gate.

Reads : data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv   (edge-level)
        data/audit/rho_sym_gate/cohort_summary.csv                        (cophenetic)
Writes: data/preprint/figures/results_section1/fig_trace_c_raw_vs_multiscale.pdf
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
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

RAW = ROOT / "data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv"
COPH = ROOT / "data/audit/rho_sym_gate/cohort_summary.csv"
OUT = ROOT / "data/preprint/figures/results_section1/fig_trace_c_raw_vs_multiscale.pdf"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
COL = {"beta": "#2166ac", "alpha": "#e08214"}
GREY = "#8a8a8a"
# manual label anchor offsets (dx, dy, ha) for the crowded grey cluster
LOFF = {
    "delta": (0.004, -0.028, "left"),
    "theta": (0.006, -0.006, "left"),
    "alpha": (0.010, 0.004, "left"),
    "beta": (0.012, 0.006, "left"),
    "low_gamma": (-0.004, 0.028, "right"),
    "high_gamma": (0.010, 0.006, "left"),
}


def main():
    raw = pd.read_csv(RAW).set_index("band")
    coph = pd.read_csv(COPH).set_index("band")

    x = {b: float(raw.loc[b, "obs_median_rho"]) for b in BANDS}
    y = {b: float(coph.loc[b, "obs_median_rho_sym"]) for b in BANDS}
    xp = {b: float(raw.loc[b, "paired_wilcoxon_p"]) for b in BANDS}       # raw gate
    yp = {b: float(coph.loc[b, "gate_p_sym"]) for b in BANDS}            # coph gate

    fig, ax = plt.subplots(figsize=(7.4, 6.8))

    lo, hi = -0.09, 0.30
    # diagonal: equal tracking at both representational levels
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=1.3, ls="--", zorder=1)
    ax.text(0.283, 0.283, "equal\ntracking", color="0.5", fontsize=10.5,
            ha="left", va="center", rotation=0)
    ax.axhline(0.0, color="0.8", lw=0.9, zorder=0)

    for b in BANDS:
        col = COL.get(b, GREY)
        # drop connector: from the diagonal (x,x) down to the cophenetic reading (x,y)
        ax.plot([x[b], x[b]], [x[b], y[b]], color=col, lw=1.4, alpha=0.45, zorder=2)
        raw_sig = xp[b] < 0.05
        coph_sig = yp[b] < 0.05
        ax.scatter(x[b], y[b], s=260 if b in COL else 190,
                   facecolor=col if coph_sig else "white",
                   edgecolor=col, linewidth=2.8 if raw_sig else 1.1, zorder=5)
        dx, dy, ha = LOFF[b]
        ax.text(x[b] + dx, y[b] + dy, BRAIN_BAND_TEX_DICT[b], color=col,
                fontsize=17 if b in COL else 14, ha=ha, va="center",
                fontweight="bold" if b in COL else "normal", zorder=6)

    # story callouts
    ax.annotate("retained:\nedge overlap was\ngenuinely hierarchical",
                xy=(x["beta"], y["beta"]), xytext=(0.115, 0.245),
                fontsize=10.5, color=COL["beta"], ha="center", va="center",
                arrowprops=dict(arrowstyle="->", color=COL["beta"], lw=1.3))
    ax.annotate("generic drift\nfiltered out",
                xy=(x["theta"], y["theta"]), xytext=(0.205, -0.055),
                fontsize=10.5, color="0.35", ha="center", va="center",
                arrowprops=dict(arrowstyle="->", color="0.5", lw=1.2))

    ax.set_xlim(lo, hi + 0.02)
    ax.set_ylim(-0.10, 0.30)
    ax.set_xlabel(r"raw edge overlap  $\rho$   (cohort median)", fontsize=15)
    ax.set_ylabel(r"cophenetic trace  $\rho^{\mathrm{coph}}$   (cohort median)", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.text(-0.145, 1.02, r"$\mathbf{c}$", transform=ax.transAxes, fontsize=17,
            va="bottom", ha="left")

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="0.4", ms=12,
               label="clears cophenetic gate"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.1, ms=12,
               label="does not"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.3", mew=2.8, ms=12,
               label="clears raw edge gate (thick ring)"),
        Line2D([0], [0], color="0.6", lw=1.3, ls="--", label="equal tracking (diagonal)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=2, frameon=False, fontsize=10.5, handletextpad=0.4,
               columnspacing=1.6)

    fig.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.20)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:trace_c — raw edge vs cophenetic hierarchy (scatter)\n")
    print(f"{'band':11s} {'raw_rho':>8s} {'raw_p':>7s} {'coph_rho':>9s} {'coph_p':>7s} {'drop':>7s}")
    for b in BANDS:
        print(f"{b:11s} {x[b]:+8.3f} {xp[b]:7.3f} {y[b]:+9.3f} {yp[b]:7.3f} {x[b]-y[b]:+7.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
