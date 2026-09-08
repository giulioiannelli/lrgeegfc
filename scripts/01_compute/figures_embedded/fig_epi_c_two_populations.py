#!/usr/bin/env python3
r"""fig:epi_c — filtering the graph to its planar skeleton recovers the whole cohort (§3, R3.3).

CORE MESSAGE (revised 2026-07-15, TMFG marker): on the fully connected / strength-preserving
graph the SOZ detector splits the cohort -- eight co-diffusing-community patients work, two
right-hemisphere hub implants fail (AUC 0.49, 0.57). Filtering the graph to its planar
(triangulated maximally-filtered) backbone recovers BOTH hubs (0.49->0.75, 0.57->0.88): every
patient now ranks seizure above healthy (median AUC 0.87->0.91, floor 0.49->0.75, six above
0.90). There is no failure population on the planar backbone -- the recovery is specific to the
planar filter, since the strength-preserving backbone that carries the cognitive trace does not
achieve it. Performance is reported across all ten patients.

Panel a: per-patient detector AUC, dense/strength graph (open) -> planar backbone (filled),
the two rescued hub implants in accent. Panel b: before-vs-after AUC scatter (diagonal = no
change); the two hubs leap far above the diagonal, the already-strong patients hold near the
ceiling.

Reads : data/sparsified_arc/marker_detector_tmfg/detector_lopo_per_patient.csv   (after, TMFG; script 27)
        data/audit/epi_propagator_detector/detector_lopo_per_patient.csv          (before, dense)
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

DET = ROOT / "data/sparsified_arc/marker_detector_tmfg/detector_lopo_per_patient.csv"   # after (TMFG)
DET_BEFORE = ROOT / "data/audit/epi_propagator_detector/detector_lopo_per_patient.csv"  # before (dense)
OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_c_two_populations.pdf"

RESCUED = ["Pat_10", "Pat_15"]                   # right-hemisphere hub implants
C_MAIN, C_RESCUE = "#2b6cb0", "#dd6b20"


def draw(target, letters=True):
    """Render the cohort-recovery panel (before->after AUC) onto ``target``.

    Exposed so the §3 compound can tile it; ``main`` calls it on a standalone figure.
    ``letters=False`` suppresses the internal per-axis a/b tags for compound use."""
    aft = pd.read_csv(DET).set_index("patient").auc
    bef = pd.read_csv(DET_BEFORE).set_index("patient").auc
    pats = [p for p in aft.index if p in bef.index]
    df = pd.DataFrame({"before": bef.reindex(pats), "after": aft.reindex(pats)}).sort_values("after")

    axa, axb = target.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.05, 1.0]))
    target.subplots_adjust(left=0.09, right=0.975, top=0.91, bottom=0.20, wspace=0.30)

    # -- panel a: per-patient AUC, dense/strength (open) -> planar (filled) -----
    y = np.arange(len(df))
    axa.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)          # chance
    for yi, (pat, row) in zip(y, df.iterrows()):
        resc = pat in RESCUED
        c = C_RESCUE if resc else C_MAIN
        axa.annotate("", xy=(row.after, yi), xytext=(row.before, yi),
                     arrowprops=dict(arrowstyle="->", color=c, lw=2.4 if resc else 1.5,
                                     alpha=0.95 if resc else 0.5), zorder=2)
        axa.scatter(row.before, yi, s=70, facecolors="white", edgecolors=c, linewidths=1.6,
                    zorder=4)                                         # before = open
        axa.scatter(row.after, yi, s=150 if resc else 120, color=c, edgecolor="white",
                    linewidth=1.0, zorder=5)                          # after = filled
    axa.set_yticks(y)
    axa.set_yticklabels(df.index)
    for tick, pat in zip(axa.get_yticklabels(), df.index):
        tick.set_color(C_RESCUE if pat in RESCUED else "0.2")
    axa.set_ylim(-0.7, len(df) - 0.3)
    axa.set_xlim(0.42, 1.0)
    axa.set_xlabel("detector AUC (leave-one-patient-out)", fontsize=12)
    axa.tick_params(axis="y", length=0)
    axa.spines[["top", "right"]].set_visible(False)
    if letters:
        axa.text(0.0, 1.02, r"$\mathbf{a}$", transform=axa.transAxes, fontsize=17,
                 va="bottom", fontweight="bold")
    axa.text(0.44, len(df) - 0.35, r"dense / strength $\rightarrow$ planar backbone",
             fontsize=10.5, color="0.35", style="italic")

    # -- panel b: before-vs-after AUC scatter (diagonal = no change) ------------
    axb.plot([0.42, 1.0], [0.42, 1.0], color="0.6", lw=1.1, ls="--", zorder=1)   # y = x
    for pat, row in df.iterrows():
        resc = pat in RESCUED
        c = C_RESCUE if resc else C_MAIN
        axb.scatter(row.before, row.after, s=170 if resc else 90, color=c,
                    edgecolor="white", linewidth=1.0, zorder=4)
        if resc:
            axb.annotate(pat.replace("_", r"\_"), (row.before, row.after),
                         textcoords="offset points", xytext=(9, -3), fontsize=9.5,
                         color=c, fontweight="bold")
    axb.set_xlim(0.42, 1.0)
    axb.set_ylim(0.42, 1.0)
    axb.set_xlabel("AUC — dense / strength graph", fontsize=12)
    axb.set_ylabel("AUC — planar backbone", fontsize=12)
    axb.spines[["top", "right"]].set_visible(False)
    if letters:
        axb.text(0.0, 1.02, r"$\mathbf{b}$", transform=axb.transAxes, fontsize=17,
                 va="bottom", fontweight="bold")
    med_b, med_a = float(df.before.median()), float(df.after.median())
    axb.text(0.97, 0.06,
             rf"median $ {med_b:.2f}\!\rightarrow\!{med_a:.2f}$,  floor $0.49\!\rightarrow\!0.75$"
             "\n" r"both hub implants rescued",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=9.0, color="0.2",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.8", alpha=0.9))

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.6, ms=9,
               label="dense / strength graph (before)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_MAIN, mec="white", ms=11,
               label="planar backbone (after)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_RESCUE, mec="white", ms=11,
               label="right-hemisphere hub implant (rescued)"),
    ]
    target.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
                  ncol=3, frameon=False, fontsize=9.5, handletextpad=0.5,
                  columnspacing=1.6)
    return df


def main():
    fig = plt.figure(figsize=(12.4, 5.6))
    df = draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:epi_c — planar backbone recovers the whole cohort\n")
    for pat, row in df.sort_values("after", ascending=False).iterrows():
        tag = "  <- rescued hub" if pat in RESCUED else ""
        print(f"  {pat}  dense={row.before:.3f} -> planar={row.after:.3f}{tag}")
    print(f"\n  median {df.before.median():.3f} -> {df.after.median():.3f}  "
          f"| floor {df.before.min():.3f} -> {df.after.min():.3f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
