#!/usr/bin/env python3
r"""fig:epi_a — seizure-onset contacts form a co-organized diffusion group (§3, R3.1).

CORE MESSAGE: a contact-by-contact SOZ label does not exist (within a patient it collapses
into coupling strength; across patients nothing transfers). Read RELATIONALLY, the seizure
contacts diffuse together: a contact's heat-kernel affinity to the rest of the seizure set
separates seizure from healthy contacts OVER AND ABOVE node strength. With the strength
baseline removed (so strength alone = chance, 0.5), the multiscale affinity reaches (on the
sparsified mst@0.20 backbone, rest_post):
  delta 0.83 | low_gamma 0.82 | beta 0.745 | alpha 0.61 (marginal)
higher, in each band, than seizure sets of matched coupling strength reach by chance,
clearing the matched-strength null in 7/8/8/5 of 10 patients. Reading the affinity ACROSS
diffusion scales (multiscale mean) rather than at the single fine scale raises the ranking
in the majority of patients (delta 7/10, beta 8/10, low_gamma 7/10) — tau sharpens the
ranking (this panel), not the top-of-list precision (that is a band-fusion effect, panel c).

Single panel: strength-RESIDUAL AUC on x (0.5 = strength-alone / chance), one row per band.
Per patient: OPEN dot = single fine scale, FILLED dot = multiscale, connected; filled dot
coloured by whether the patient beats its OWN matched-strength fake-seizure null.

Reads : data/sparsified_arc/epi_arc_mst020/per_cell.csv
Writes: data/preprint/figures/results_section3/fig_epi_a_relational_marker.pdf
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

SRC = ROOT / "data/sparsified_arc/epi_arc_mst020/per_cell.csv"
OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_a_relational_marker.pdf"

BANDS = ["delta", "low_gamma", "beta", "alpha"]      # descending median multiscale AUC
C_BAND = {"delta": "#08519c", "low_gamma": "#2b8cbe",
          "beta": "#d1491c", "alpha": "#807dba"}
C_BEAT, C_MISS = "#3a9a4f", "#b9b9b9"                 # filled dot: beats / within own null
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
Y_OFF = dict(zip(COHORT, np.linspace(-0.26, 0.26, len(COHORT))))


def draw(target):
    """Render the relational-marker panel onto ``target`` (a Figure or SubFigure).

    Exposed so the §3 compound can tile it; ``main`` calls it on a standalone figure."""
    pc = pd.read_csv(SRC)

    yof = {b: len(BANDS) - 1 - i for i, b in enumerate(BANDS)}   # delta at top

    ax = target.subplots()
    target.subplots_adjust(left=0.11, right=0.97, top=0.93, bottom=0.28)

    ax.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)      # strength-alone / chance

    counts = {}
    for b in BANDS:
        y = yof[b]
        d = pc[pc.band == b].set_index("patient")
        n_beat = 0
        for pat in COHORT:
            if pat not in d.index:
                continue
            a_s = float(d.loc[pat, "auc_single"])
            a_m = float(d.loc[pat, "auc_multi"])
            beats = bool(d.loc[pat, "p_null_multi"] < 0.05)
            n_beat += int(beats)
            yy = y + Y_OFF[pat]
            # single -> multi connector (the tau-sweep effect for this patient)
            ax.plot([a_s, a_m], [yy, yy], color="0.75", lw=1.0, zorder=3,
                    solid_capstyle="round")
            ax.scatter(a_s, yy, s=34, facecolor="white",
                       edgecolor=C_BAND[b], linewidth=1.1, zorder=4)   # single = open
            ax.scatter(a_m, yy, s=96, color=C_BEAT if beats else C_MISS,
                       edgecolor="white", linewidth=0.8, zorder=5)      # multi = filled
        med = float(d["auc_multi"].median())
        counts[b] = n_beat
        ax.plot([0.5, med], [y, y], color=C_BAND[b], lw=3.4, zorder=6,
                solid_capstyle="round")
        ax.scatter(med, y, s=190, color=C_BAND[b], edgecolor="white", linewidth=1.3,
                   zorder=7)
        ax.text(med, y + 0.33, f"{n_beat}/10", color=C_BAND[b],
                fontsize=11, ha="center", va="bottom", fontweight="bold")

    ax.set_yticks([yof[b] for b in BANDS])
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    for tick, b in zip(ax.get_yticklabels(), BANDS):
        tick.set_color(C_BAND[b])
        tick.set_fontsize(19)
    ax.set_ylim(-0.7, len(BANDS) - 0.3)
    ax.set_xlim(0.30, 1.0)
    ax.set_xlabel(r"strength-residual affinity AUC   "
                  r"(seizure vs healthy; $0.5$ = strength alone / chance)", fontsize=12.5)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="#444", mew=1.1, ms=8,
               label="single fine scale ($\\tau_{\\min}$)"),
        Line2D([0], [0], marker="o", ls="", mfc=C_BEAT, mec="white", ms=11,
               label="multiscale, beats own null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_MISS, mec="white", ms=11,
               label="multiscale, within own null"),
        Line2D([0], [0], marker="o", ls="-", color="#444", mfc="#444", mec="white",
               ms=12, lw=3.2, label="cohort median AUC ($n$/10 beat own null)"),
    ]
    target.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08),
                  ncol=2, frameon=False, fontsize=10, handletextpad=0.5,
                  columnspacing=1.6)
    return dict(pc=pc, counts=counts)


def main():
    fig = plt.figure(figsize=(8.8, 5.2))
    stats = draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    pc, counts = stats["pc"], stats["counts"]
    print("fig:epi_a — relational SOZ marker (strength-residual, all contacts, mst@0.20)\n")
    for b in BANDS:
        x = pc[pc.band == b]
        print(f"  {b:11s} median AUC single={x.auc_single.median():.3f} "
              f"multi={x.auc_multi.median():.3f}  beats={counts[b]}/10  "
              f"multi>single={(x.auc_gain_multi_vs_single > 0).sum()}/10")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
