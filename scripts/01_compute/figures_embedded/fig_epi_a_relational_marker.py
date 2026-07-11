#!/usr/bin/env python3
r"""fig:epi_a — seizure-onset contacts form a co-organized diffusion group (§3, R3.1).

CORE MESSAGE: a contact-by-contact SOZ label does not exist (within a patient it collapses
into coupling strength; across patients nothing transfers). Read RELATIONALLY, the seizure
contacts diffuse together: a contact's heat-kernel affinity to the rest of the seizure set
separates seizure from healthy contacts OVER AND ABOVE node strength. With the strength
baseline removed (so strength alone = chance, 0.5), the affinity still reaches:
  delta 0.80 | low_gamma 0.74 | beta 0.69 | alpha 0.60 (marginal)
higher, in each band, than matched-coupling contacts reach by chance (no reshuffle of 300
met it), holding in 7/6/7/6 of 10 patients.

Single panel: strength-RESIDUAL AUC on x (0.5 = strength-alone / chance), one row per band,
per-patient dots coloured by whether each patient beats its OWN matched-strength null.

(The 'not proximity' point is carried in prose: |ImCoh| is zero-lag-immune to volume
conduction. The off-shaft distant-marker demonstration lives in prose / supplement, not
here -- it is a different, harder contact set and would read as a second delta value on a
shared axis.)

Reads : data/audit/epi_marker_allcontacts/{allcontacts_per_patient,allcontacts_stratnull_per_patient,allcontacts_nulls}.csv
Writes: data/reports/results_section3/fig_epi_a_relational_marker.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

AC = ROOT / "data/audit/epi_marker_allcontacts"
OUT = ROOT / "data/preprint/figures/results_section3/fig_epi_a_relational_marker.pdf"

MARKER = "heat_t5"                                   # slow heat-kernel (README headline)
BANDS = ["delta", "low_gamma", "beta", "alpha"]      # descending median AUC
C_BAND = {"delta": "#08519c", "low_gamma": "#2b8cbe",
          "beta": "#d1491c", "alpha": "#807dba"}
C_BEAT, C_MISS = "#3a9a4f", "#b9b9b9"                 # dot: beats / within own matched-strength null
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
Y_OFF = dict(zip(COHORT, np.linspace(-0.26, 0.26, len(COHORT))))


def draw(target):
    """Render the relational-marker panel onto ``target`` (a Figure or SubFigure).

    Exposed so the §3 compound can tile it; ``main`` calls it on a standalone figure."""
    pp = pd.read_csv(AC / "allcontacts_per_patient.csv")
    pp = pp[pp.marker == MARKER]
    sn = pd.read_csv(AC / "allcontacts_stratnull_per_patient.csv")
    nulls = pd.read_csv(AC / "allcontacts_nulls.csv").set_index("band")

    yof = {b: len(BANDS) - 1 - i for i, b in enumerate(BANDS)}   # delta at top

    ax = target.subplots()
    target.subplots_adjust(left=0.11, right=0.97, top=0.93, bottom=0.26)

    # matched-strength null envelope (label-shuffle): median ~0.49, p95 ~0.55-0.57
    nmed = float(nulls.null_median_mean.mean())
    np95 = float(nulls.null_median_p95.max())
    ax.axvspan(nmed, np95, color="0.82", alpha=0.6, zorder=0)
    ax.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)      # strength-alone / chance

    for b in BANDS:
        y = yof[b]
        dots = pp[pp.band == b].set_index("patient").auc_resid
        beats = sn[sn.band == b].set_index("patient").beats_p95
        for pat in COHORT:
            if pat not in dots.index:
                continue
            c = C_BEAT if bool(beats.get(pat, False)) else C_MISS
            ax.scatter(dots[pat], y + Y_OFF[pat], s=100, color=c, edgecolor="white",
                       linewidth=0.8, zorder=4)
        med = float(nulls.loc[b, "real_median_auc"])
        ax.plot([0.5, med], [y, y], color=C_BAND[b], lw=3.2, zorder=5,
                solid_capstyle="round")
        ax.scatter(med, y, s=185, color=C_BAND[b], edgecolor="white", linewidth=1.3,
                   zorder=6)
        ax.text(med, y + 0.32, f"{int(beats.sum())}/10", color=C_BAND[b],
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
        Line2D([0], [0], marker="o", ls="", mfc=C_BEAT, mec="white", ms=11,
               label="patient beats own matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc=C_MISS, mec="white", ms=11,
               label="within own null"),
        Line2D([0], [0], marker="o", ls="-", color="#444", mfc="#444", mec="white",
               ms=12, lw=3.2, label="cohort median AUC ($n$/10 beat own null)"),
        Patch(facecolor="0.82", alpha=0.6,
              label="matched-strength null (label-shuffle, $0/300$ reached obs.)"),
    ]
    target.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.06),
                  ncol=2, frameon=False, fontsize=10, handletextpad=0.5,
                  columnspacing=1.6)
    return dict(nulls=nulls, sn=sn)


def main():
    fig = plt.figure(figsize=(8.8, 5.2))
    stats = draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    nulls, sn = stats["nulls"], stats["sn"]
    print("fig:epi_a — relational SOZ marker (strength-residual, all contacts)\n")
    for b in BANDS:
        beats = int(sn[sn.band == b].beats_p95.sum())
        print(f"  {b:11s} median AUC={nulls.loc[b,'real_median_auc']:.3f} "
              f"null={nulls.loc[b,'null_median_mean']:.3f} p95={nulls.loc[b,'null_median_p95']:.3f} "
              f"beats={beats}/10 p_emp={nulls.loc[b,'p_empirical']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
