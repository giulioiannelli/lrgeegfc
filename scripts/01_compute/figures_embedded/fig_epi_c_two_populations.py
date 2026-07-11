#!/usr/bin/env python3
r"""fig:epi_c — the cohort splits into two seizure-network shapes (§3, R3.3).

CORE MESSAGE: the cohort mean hides two populations. In eight patients the seizure zone is a
co-diffusing COMMUNITY and the detector is strong (AUC 0.90-0.99 in the clearest). In the
other two -- both right-hemisphere implants -- the seizure zone IS a network HUB rather than
a community, and the propagator read-out fails (AUC 0.49, 0.57); there, plain coupling
strength predicts the seizure contacts the propagator misses (Pat_15: strength 0.93 vs
propagator affinity 0.21). Choosing per patient between the affinity and the strength read-
out for the hub-type implants recovers nine of ten (off-shaft AUC 0.716 -> 0.750); in the
tenth (Pat_10) both read-outs fail. Performance is reported per patient, not pooled -- the
mean would average a working detector with a known failure mode.

Panel a: per-patient detector AUC, split by population. Panel b: propagator affinity vs node
strength per patient -- community patients below the diagonal (affinity wins), hub patients
above it or at the floor (strength wins, or both fail); the switch takes the better regime.

Reads : data/audit/epi_propagator_detector/detector_lopo_per_patient.csv
        data/audit/epi_marker_compound/compound_auc_per_patient.csv
Writes: data/reports/results_section3/fig_epi_c_two_populations.pdf
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

DET = ROOT / "data/audit/epi_propagator_detector/detector_lopo_per_patient.csv"
CMP = ROOT / "data/audit/epi_marker_compound/compound_auc_per_patient.csv"
OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_c_two_populations.pdf"

HUB_AUC = 0.60                                   # detector AUC below -> hub-type implant
C_COMM, C_HUB = "#3a9a4f", "#d1352b"
COH_AFF, COH_SW = 0.716, 0.750                   # off-shaft affinity -> switch (audit_104)


def draw(target, letters=True):
    """Render the two-population panel (per-patient AUC + regime switch) onto ``target``.

    Exposed so the §3 compound can tile it; ``main`` calls it on a standalone figure.
    ``letters=False`` suppresses the internal per-axis a/b tags for compound use."""
    det = pd.read_csv(DET).set_index("patient")
    cmp = pd.read_csv(CMP)
    aff = cmp[cmp.compound == "affinity"].set_index("patient").auc
    stg = cmp[cmp.compound == "strength"].set_index("patient").auc
    sw = cmp[cmp.compound == "switch"].set_index("patient").auc

    det = det.sort_values("auc", ascending=True)
    is_hub = det.auc < HUB_AUC

    axa, axb = target.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.05, 1.0]))
    target.subplots_adjust(left=0.09, right=0.975, top=0.91, bottom=0.20, wspace=0.26)

    # -- panel a: per-patient detector AUC, two populations --------------------
    y = np.arange(len(det))
    axa.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)          # chance
    for yi, (pat, row) in zip(y, det.iterrows()):
        c = C_HUB if row.auc < HUB_AUC else C_COMM
        axa.plot([0.5, row.auc], [yi, yi], color=c, lw=2.6, zorder=2,
                 solid_capstyle="round")
        axa.scatter(row.auc, yi, s=135, color=c, edgecolor="white", linewidth=1.1,
                    zorder=4)
    axa.set_yticks(y)
    axa.set_yticklabels(det.index)
    for tick, pat in zip(axa.get_yticklabels(), det.index):
        tick.set_color(C_HUB if is_hub[pat] else C_COMM)
    axa.set_ylim(-0.7, len(det) - 0.3)
    axa.set_xlim(0.42, 1.0)
    axa.set_xlabel("detector AUC (leave-one-patient-out)", fontsize=12)
    axa.tick_params(axis="y", length=0)
    axa.spines[["top", "right"]].set_visible(False)
    if letters:
        axa.text(0.0, 1.02, r"$\mathbf{a}$", transform=axa.transAxes, fontsize=17,
                 va="bottom", fontweight="bold")
    n_comm, n_hub = int((~is_hub).sum()), int(is_hub.sum())
    axa.text(0.62, len(det) - 0.5, f"community  ({n_comm})", color=C_COMM,
             fontsize=11.5, fontweight="bold")
    axa.text(0.62, 0.55, f"hub  ({n_hub})", color=C_HUB, fontsize=11.5,
             fontweight="bold")

    # -- panel b: affinity -> switch per patient (the hub rescue) --------------
    # switch = take node strength for hub-like seeds, else the propagator affinity.
    # Community affinity already works (preserved); Pat_15 is rescued across chance;
    # Pat_10 is unrecoverable (both read-outs fail).
    order_b = sw.reindex(det.index).sort_values().index          # by switch AUC
    yb = np.arange(len(order_b))
    axb.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)      # chance
    for yi, pat in zip(yb, order_b):
        if pat not in aff.index:
            continue
        c = C_HUB if is_hub.get(pat, False) else C_COMM
        a, s = float(aff[pat]), float(sw[pat])
        jump = abs(s - a) > 0.03
        if jump:
            axb.annotate("", xy=(s, yi), xytext=(a, yi),
                         arrowprops=dict(arrowstyle="->", color=c, lw=2.0), zorder=3)
        else:
            axb.plot([a, s], [yi, yi], color=c, lw=2.0, alpha=0.5, zorder=2)
        axb.scatter(a, yi, s=70, facecolors="white", edgecolors=c, linewidths=1.6,
                    zorder=4)                                     # affinity (start)
        axb.scatter(s, yi, s=130, color=c, edgecolor="white", linewidth=1.0,
                    zorder=5)                                     # switch (end)
    axb.set_yticks(yb)
    axb.set_yticklabels(order_b)
    for tick, pat in zip(axb.get_yticklabels(), order_b):
        tick.set_color(C_HUB if is_hub.get(pat, False) else C_COMM)
    axb.set_ylim(-0.7, len(order_b) - 0.3)
    axb.set_xlim(0.15, 1.0)
    axb.set_xlabel("off-shaft AUC:  affinity  $\\rightarrow$  per-patient switch",
                   fontsize=12)
    axb.tick_params(axis="y", length=0)
    axb.spines[["top", "right"]].set_visible(False)
    if letters:
        axb.text(0.0, 1.02, r"$\mathbf{b}$", transform=axb.transAxes, fontsize=17,
                 va="bottom", fontweight="bold")
    axb.text(0.97, 0.06,
             rf"cohort: {COH_AFF:.2f} $\rightarrow$ {COH_SW:.2f},  $8/10 \rightarrow 9/10$"
             "\n" r"(Pat\_10 unrecoverable — both fail)",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=9.0,
             color="0.2",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.8", alpha=0.9))

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc=C_COMM, mec="white", ms=11,
               label="community population (propagator strong)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_HUB, mec="white", ms=11,
               label="hub population (propagator fails; strength predicts)"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.6, ms=9,
               label="affinity (b, start)"),
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="white", ms=11,
               label="per-patient switch (b, end)"),
    ]
    target.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
                  ncol=2, frameon=False, fontsize=10.0, handletextpad=0.5,
                  columnspacing=2.0)
    return det, aff, stg, sw


def main():
    fig = plt.figure(figsize=(12.4, 5.6))
    det, aff, stg, sw = draw(fig)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:epi_c — two seizure-network populations\n")
    for pat, row in det.iterrows():
        tag = "HUB" if row.auc < HUB_AUC else "community"
        print(f"  {pat} detector AUC={row.auc:.3f} [{tag}]  affinity={aff.get(pat,float('nan')):.3f}"
              f" strength={stg.get(pat,float('nan')):.3f} switch={sw.get(pat,float('nan')):.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
