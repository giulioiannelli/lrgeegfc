#!/usr/bin/env python3
"""preprint_40 — from the verified marker to a deployable probability: the propagator+strength
compound and the leave-one-patient-out P(SOZ). Companion to preprint_39 (the core marker).

Four panels, all from audit_103/104 result CSVs (read-only):

(a) Two populations (the motivation): per-patient enrichment (lift@5) of the propagator vs
    node strength. Above the diagonal = propagator wins (SOZ = diffusion community,
    responders); below = strength wins (SOZ = hubs, Pat_10/15). They are complementary.
(b) Compound comparison: cohort-median lift@5 for each propagator/strength combination. The
    seed-regime "adaptive" weight (hub-like seeds -> weight strength) wins.
(c) Coverage gain: per-patient off-shaft AUC, propagator alone -> adaptive compound; the
    compound rescues a hub-patient (8/10 -> 9/10) without hurting responders.
(d) Calibrated probability: leave-one-patient-out logistic P(SOZ) vs the observed SOZ
    fraction (reliability curve). On the diagonal = honest probability; the high-P end is
    well-calibrated, the values are honestly modest (SOZ are rare off-shaft).

PDF only, full vector, transparent. House style: stats to stdout, no in-axes stat text.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

CP = ROOT / "data/audit/epi_marker_compound"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

RESP_THR = 0.568
C_AFF = "#1b7837"      # propagator / affinity
C_STR = "#b2182b"      # strength
C_SWITCH = "#2166ac"   # hard regime switch (the honest winner)
C_BLEND = "#e08214"    # adaptive blend (cautionary — damages responders)
C_OTHER = "#9e9e9e"
COMP_ORDER = ["strength", "mean", "affinity", "union_max", "adaptive", "switch"]
COMP_LABEL = {"strength": "strength", "affinity": "propagator", "union_max": "union",
              "mean": "mean", "adaptive": "blend", "switch": "switch\n(seed-regime)"}
COMP_COL = {"strength": C_STR, "affinity": C_AFF, "union_max": C_OTHER,
            "mean": C_OTHER, "adaptive": C_BLEND, "switch": C_SWITCH}


def panel_a(ax, pp):
    aff = pp[pp.compound == "affinity"].set_index("patient")
    strg = pp[pp.compound == "strength"].set_index("patient")
    lim = 16.0
    ax.plot([0, lim], [0, lim], "-", color="0.7", lw=0.9, zorder=0)
    for pat in aff.index:
        x = float(strg.loc[pat, "lift_at_5"])
        y = float(aff.loc[pat, "lift_at_5"])
        resp = float(aff.loc[pat, "auc"]) >= RESP_THR
        ax.scatter(x, y, s=52, c=(C_AFF if resp else C_STR), edgecolor="black",
                   linewidth=0.4, zorder=3)
        ax.annotate(pat.replace("Pat_", ""), (x, y), textcoords="offset points",
                    xytext=(3, 3), fontsize=6)
    ax.set_xlim(-0.6, lim)
    ax.set_ylim(-0.6, lim)
    ax.set_xlabel("strength enrichment (lift@5)")
    ax.set_ylabel("propagator enrichment (lift@5)")
    ax.text(0.03, 0.95, "propagator wins\n(SOZ = community)", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.4, color=C_AFF)
    ax.text(0.97, 0.05, "strength wins\n(SOZ = hubs)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.4, color=C_STR)
    ax.text(0.02, 1.03, "a  two complementary populations", transform=ax.transAxes,
            va="bottom", fontsize=8.3)


def panel_b(ax, coh):
    d = coh.set_index("compound").loc[COMP_ORDER]
    x = np.arange(len(COMP_ORDER))
    ax.axhspan(0.487, 0.568, color="0.86", zorder=0)                 # label-shuffle null band
    ax.axhline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    ax.bar(x, d.med_auc.values, width=0.62, color=[COMP_COL[c] for c in COMP_ORDER],
           edgecolor="black", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([COMP_LABEL[c] for c in COMP_ORDER], fontsize=6.6)
    ax.set_ylabel("off-shaft AUC (cohort median)")
    ax.set_ylim(0.45, 0.80)
    ax.text(0.02, 1.03, "b  compound comparison (AUC)", transform=ax.transAxes,
            va="bottom", fontsize=8.3)


def panel_c(ax, pp):
    aff = pp[pp.compound == "affinity"].set_index("patient")["auc"]
    blend = pp[pp.compound == "adaptive"].set_index("patient")["auc"]
    sw = pp[pp.compound == "switch"].set_index("patient")["auc"]
    order = aff.sort_values().index
    y = np.arange(len(order))
    ax.axvspan(0.487, 0.568, color="0.86", zorder=0)                 # null band
    ax.axvline(0.5, color="black", ls="--", lw=0.8, zorder=1)
    for yi, pat in zip(y, order):
        ax.plot([float(aff[pat]), float(sw[pat])], [yi, yi], "-", color="0.7", lw=0.8, zorder=2)
        ax.scatter(float(aff[pat]), yi, s=40, c=C_AFF, edgecolor="black", linewidth=0.3, zorder=3)
        ax.scatter(float(blend[pat]), yi, s=44, c=C_BLEND, marker="x", linewidth=1.4, zorder=4)
        ax.scatter(float(sw[pat]), yi, s=40, c=C_SWITCH, edgecolor="black", linewidth=0.3, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace("Pat_", "") for p in order], fontsize=6.6)
    ax.set_xlabel("off-shaft distant-discovery AUC")
    ax.set_ylabel("patient")
    ax.set_xlim(0.12, 1.0)
    ax.text(0.02, 1.03, "c  switch preserves responders; blend (x) damages them",
            transform=ax.transAxes, va="bottom", fontsize=8.0)


def panel_d(ax, cal):
    lim = float(max(cal.p_pred_mean.max(), cal.soz_observed.max())) * 1.15 + 0.01
    ax.plot([0, lim], [0, lim], "-", color="0.7", lw=0.9, zorder=0)
    sizes = 18 + 60 * (cal.n / cal.n.max())
    ax.scatter(cal.p_pred_mean, cal.soz_observed, s=sizes, c=C_SWITCH,
               edgecolor="black", linewidth=0.4, zorder=3)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("predicted P(SOZ)")
    ax.set_ylabel("observed SOZ fraction")
    ax.text(0.03, 0.95, "on diagonal = calibrated", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.4, color="0.4")
    ax.text(0.02, 1.03, "d  leave-one-patient-out P(SOZ)", transform=ax.transAxes,
            va="bottom", fontsize=8.3)


def main():
    pp = pd.read_csv(CP / "compound_auc_per_patient.csv")
    coh = pd.read_csv(CP / "compound_auc_cohort.csv")
    cal = pd.read_csv(CP / "calibration_curve.csv")

    fig = plt.figure(figsize=(11.0, 8.0))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30)
    panel_a(fig.add_subplot(gs[0, 0]), pp)
    panel_b(fig.add_subplot(gs[0, 1]), coh)
    panel_c(fig.add_subplot(gs[1, 0]), pp)
    panel_d(fig.add_subplot(gs[1, 1]), cal)

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_AFF,
               markeredgecolor="black", markersize=8, label="propagator (affinity)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_STR,
               markeredgecolor="black", markersize=8, label="node strength"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_SWITCH,
               markeredgecolor="black", markersize=8, label="switch (seed-regime gate)"),
        Line2D([0], [0], marker="x", color=C_BLEND, markersize=8, linestyle="none",
               label="blend (cautionary — damages responders)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.015),
               ncol=4, frameon=False, fontsize=7.4)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / "fig_epi_compound_probability.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_40] -> {out.name}")
    print("  compound cohort lift@5:", dict(zip(coh.compound, coh.med_lift5.round(2))))
    print("  adaptive AUC per patient:",
          dict(zip(pp[pp.compound == 'adaptive'].patient,
                   pp[pp.compound == 'adaptive'].auc.round(3))))


if __name__ == "__main__":
    main()
