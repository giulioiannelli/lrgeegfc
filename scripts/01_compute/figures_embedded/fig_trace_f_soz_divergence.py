#!/usr/bin/env python3
r"""fig:trace_f — over the seizure-onset zone, beta and alpha diverge (Results §1, para f).

CORE MESSAGE: the two bands that carry the trace treat the diseased core oppositely.
beta — the focal, OFC-anchored band — is carried by healthy coupling and LEAVES the SOZ
out: its SOZ-internal persistence is generic (p=0.93) and does not clear the strength null.
alpha — the diffuse band with no cortical address — does the REVERSE: its persistence
CONCENTRATES in SOZ-internal coupling (p<0.001) and STRENGTHENS under the strength-matched
surrogate (p=0.005), a genuine recruitment of the diseased core. The two lines cross: beta
steers around the epileptogenic tissue, alpha pulls it in.

Marker filled = clears matched-strength null; star = clears AND strengthens under it.

Reads : data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv
Writes: data/reports/results_section1/fig_trace_f_soz_divergence.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv"
OUT = ROOT / "data/reports/results_section1/fig_trace_f_soz_divergence.pdf"

C_B, C_A = "#2166ac", "#e08214"


def main():
    d = pd.read_csv(SRC)
    d = d[(d.stratify == "epi")]

    def get(band, cls):
        r = d[(d.band == band) & (d.config == cls)].iloc[0]
        return float(r.cohort_obs_rho_sym), float(r.ms_wilcoxon_p)

    b_out, b_out_p = get("beta", "nonepi_nonepi")
    b_in, b_in_p = get("beta", "epi_epi")
    a_out, a_out_p = get("alpha", "nonepi_nonepi")
    a_in, a_in_p = get("alpha", "epi_epi")

    fig, ax = plt.subplots(figsize=(7.6, 6.0))

    def _mk(x, y, p, col, star=False):
        if star:
            ax.scatter(x, y, s=340, color=col, edgecolor="white", lw=1.0,
                       marker="*", zorder=6)
        elif p < 0.05:
            ax.scatter(x, y, s=150, color=col, edgecolor="white", lw=1.0, zorder=6)
        else:
            ax.scatter(x, y, s=150, facecolor="white", edgecolor=col, lw=2.2, zorder=6)

    ax.plot([0, 1], [b_out, b_in], color=C_B, lw=3.4, zorder=4)
    ax.plot([0, 1], [a_out, a_in], color=C_A, lw=3.4, zorder=4)
    _mk(0, b_out, b_out_p, C_B)
    _mk(1, b_in, b_in_p, C_B)
    _mk(0, a_out, a_out_p, C_A)
    _mk(1, a_in, a_in_p, C_A, star=True)

    # endpoint band labels
    ax.text(-0.05, b_out, r"$\beta$", color=C_B, fontsize=19, ha="right", va="center", fontweight="bold")
    ax.text(-0.05, a_out, r"$\alpha$", color=C_A, fontsize=19, ha="right", va="center", fontweight="bold")
    ax.text(1.05, b_in, r"$\beta$", color=C_B, fontsize=19, ha="left", va="center", fontweight="bold")
    ax.text(1.05, a_in, r"$\alpha$", color=C_A, fontsize=19, ha="left", va="center", fontweight="bold")

    # callouts
    ax.annotate("recruits the SOZ\nstrengthens under null ($p=0.005$)",
                xy=(1, a_in), xytext=(0.42, 0.40), fontsize=11, color=C_A, ha="center",
                arrowprops=dict(arrowstyle="->", color=C_A, lw=1.4))
    ax.annotate("steers around the SOZ\ngeneric ($p=0.93$)", xy=(1, b_in),
                xytext=(0.60, 0.135), fontsize=11, color=C_B, ha="center",
                arrowprops=dict(arrowstyle="->", color=C_B, lw=1.4))

    ax.set_xlim(-0.32, 1.42)
    ax.set_ylim(0.0, 0.46)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["non-SOZ\ncoupling", "SOZ-internal\ncoupling"], fontsize=14)
    ax.set_ylabel(r"held trace  $\rho^{\mathrm{coph}}$  (cohort)", fontsize=15)
    ax.tick_params(axis="y", labelsize=13)
    ax.text(0.0, 1.02, r"$\mathbf{f}$", transform=ax.transAxes, fontsize=17,
            va="bottom", ha="left")

    handles = [
        Line2D([0], [0], marker="*", ls="", mfc="0.4", mec="white", ms=17,
               label="clears null & strengthens"),
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="white", ms=11,
               label="clears matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=2.2, ms=11,
               label="does not clear"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=3, frameon=False, fontsize=10.5, handletextpad=0.4, columnspacing=1.5)

    fig.subplots_adjust(left=0.13, right=0.92, top=0.94, bottom=0.19)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:trace_f — SOZ beta/alpha divergence\n")
    print(f"  beta : non-SOZ {b_out:+.3f} (p={b_out_p:.3f}) -> SOZ {b_in:+.3f} (p={b_in_p:.3f})")
    print(f"  alpha: non-SOZ {a_out:+.3f} (p={a_out_p:.3f}) -> SOZ {a_in:+.3f} (p={a_in_p:.3f})")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
