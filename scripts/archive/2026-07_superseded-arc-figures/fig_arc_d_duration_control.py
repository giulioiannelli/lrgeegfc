#!/usr/bin/env python3
r"""fig:arc_d — the inference component does not track recording length (Results §2, R2.4).

CORE MESSAGE: the test phase runs 1.4-2.5x longer than the learning phase, and the
inference contrast f = D_test - D_learn inherits that asymmetry -- so a component that
simply grew with recording length would be the obvious confound. It does not. Across
patients the beta inference-specific component is UNCORRELATED with the test/learn
duration ratio (Spearman rho = +0.10, p = 0.78), while the bands whose components DO
track length -- alpha (rho = +0.55) and high-gamma (rho = +0.70, p = 0.025) -- carry no
inference-specific trace at all (both inference-null, R2.1). The rhythms that follow
recording length and the rhythm that carries the inference are DISJOINT.

Reads : data/audit/consolidation_arc_rhosym/arc_per_patient.csv   (T_infspec_pe / patient)
        data/audit/inference_localization/length_control.csv      (test/learn ratio)
Writes: data/reports/results_section2/fig_arc_d_duration_control.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

ARC = ROOT / "data/audit/consolidation_arc_rhosym/arc_per_patient.csv"
LEN = ROOT / "data/audit/inference_localization/length_control.csv"
OUT = ROOT / "data/reports/results_section2/fig_arc_d_duration_control.pdf"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
INFER_BAND = "beta"                       # the sole inference-carrying band (R2.1)
TRACKERS = ("alpha", "high_gamma")        # the length-tracking bands (paragraph)
C_INFER, C_TRACK, C_NEU = "#d1491c", "#6a51a3", "#b9b9b9"


def band_color(b: str) -> str:
    if b == INFER_BAND:
        return C_INFER
    if b in TRACKERS:
        return C_TRACK
    return C_NEU


def main():
    arc = pd.read_csv(ARC)
    ratio = (pd.read_csv(LEN).drop_duplicates("pat").set_index("pat")["ratio"])

    # per-band Spearman(T_infspec_pe, test/learn ratio)
    rows = {}
    for b in BANDS:
        sub = arc[arc.band == b].set_index("patient")["T_infspec_pe"]
        common = sub.index.intersection(ratio.index)
        rho, p = spearmanr(sub.loc[common], ratio.loc[common])
        rows[b] = dict(rho=float(rho), p=float(p),
                       x=ratio.loc[common].to_numpy(),
                       y=sub.loc[common].to_numpy())

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.2, 5.6),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.15]))

    # -- panel a: per-band Spearman rho lollipop --------------------------------
    order = sorted(BANDS, key=lambda b: rows[b]["rho"])       # ascending -> tracker on top
    y0 = np.arange(len(order))
    axa.axvline(0, color="0.55", lw=1.0, zorder=1)
    for y, b in zip(y0, order):
        col = band_color(b)
        rho, p = rows[b]["rho"], rows[b]["p"]
        axa.plot([0, rho], [y, y], color=col, lw=2.4, zorder=2,
                 solid_capstyle="round")
        sig = p < 0.05
        axa.scatter([rho], [y], s=150 if b == INFER_BAND else 120, marker="o",
                    facecolors=col if sig else "white", edgecolors=col,
                    linewidths=1.8, zorder=4)
        if sig:
            axa.annotate(r"$\ast$", (rho, y), textcoords="offset points",
                         xytext=(9, 3), fontsize=13, color=col)
        axa.text(rho + (0.03 if rho >= 0 else -0.03), y + 0.34,
                 rf"$\rho={rho:+.2f}$", color=col, fontsize=10.5,
                 ha="left" if rho >= 0 else "right", va="center")
    axa.set_yticks(y0)
    axa.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in order], fontsize=17)
    for tick, b in zip(axa.get_yticklabels(), order):
        tick.set_color(band_color(b))
        if b in (INFER_BAND,) + TRACKERS:
            tick.set_fontweight("bold")
    axa.set_ylim(-0.7, len(order) - 0.3)
    axa.set_xlim(-0.35, 0.9)
    axa.set_xlabel(r"Spearman $\rho$   ($T_{\mathrm{inf}\cdot e}$  vs  "
                   r"test/learn length ratio)", fontsize=12)
    axa.tick_params(axis="y", length=0)
    axa.spines[["top", "right"]].set_visible(False)
    axa.text(0.0, 1.02, r"$\mathbf{a}$", transform=axa.transAxes, fontsize=17,
             va="bottom", fontweight="bold")

    # -- panel b: beta (flat) vs high_gamma (steep) per-patient scatter ---------
    axb.axhline(0, color="0.6", lw=0.9, ls=":", zorder=1)
    for b, mk in ((INFER_BAND, "o"), ("high_gamma", "^")):
        col = band_color(b)
        x, y = rows[b]["x"], rows[b]["y"]
        axb.scatter(x, y, s=95, marker=mk, facecolors=col, edgecolors="white",
                    linewidths=0.9, alpha=0.9, zorder=4,
                    label=rf"{BRAIN_BAND_TEX_DICT[b]}  $\rho={rows[b]['rho']:+.2f}$, "
                          rf"$p={rows[b]['p']:.2f}$")
        # visual OLS guide line (Spearman is the reported stat, not this fit)
        sl, ic = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        axb.plot(xs, sl * xs + ic, color=col, lw=2.0,
                 ls="-" if b != INFER_BAND else "--", alpha=0.85, zorder=3)
    axb.set_xlabel(r"test/learn recording-length ratio", fontsize=12)
    axb.set_ylabel(r"inference-specific persistence  $T_{\mathrm{inf}\cdot e}$",
                   fontsize=12)
    axb.spines[["top", "right"]].set_visible(False)
    axb.text(0.0, 1.02, r"$\mathbf{b}$", transform=axb.transAxes, fontsize=17,
             va="bottom", fontweight="bold")
    axb.legend(loc="upper left", frameon=False, fontsize=10.5, handletextpad=0.4)

    # figure-level key for panel a categories
    handles = [
        Line2D([0], [0], color=C_INFER, lw=2.4, marker="o", mfc="white",
               mec=C_INFER, ms=10, label="inference-carrying band (flat vs length)"),
        Line2D([0], [0], color=C_TRACK, lw=2.4, marker="o", mfc=C_TRACK,
               mec=C_TRACK, ms=10,
               label=r"length-tracking band ($\ast$ = $p<0.05$; inference-null)"),
        Line2D([0], [0], color=C_NEU, lw=2.4, marker="o", mfc="white", mec=C_NEU,
               ms=10, label="other band"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False, fontsize=10, handletextpad=0.5,
               columnspacing=1.6)

    fig.subplots_adjust(left=0.11, right=0.975, top=0.92, bottom=0.20, wspace=0.28)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_d — duration control (T_infspec_pe vs test/learn length ratio)\n")
    print(f"  ratio range {ratio.min():.2f}-{ratio.max():.2f}")
    for b in order[::-1]:
        cat = ("INFERENCE" if b == INFER_BAND else
               "tracker" if b in TRACKERS else "other")
        print(f"  {b:11s} rho={rows[b]['rho']:+.2f} p={rows[b]['p']:.3f}  [{cat}]")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
