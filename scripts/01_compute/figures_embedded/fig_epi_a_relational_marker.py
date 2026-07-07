#!/usr/bin/env python3
r"""fig:epi_a — seizure-onset contacts form a co-organized diffusion group (§3, R3.1).

CORE MESSAGE: a contact-by-contact SOZ label does not exist (within a patient it collapses
into coupling strength; across patients nothing transfers). Read RELATIONALLY, the seizure
contacts diffuse together: a contact's heat-kernel affinity to the rest of the seizure set
separates seizure from healthy contacts OVER AND ABOVE node strength. With the strength
baseline removed (so strength alone = chance, 0.5), the affinity still reaches:
  delta 0.80 | low_gamma 0.74 | beta 0.69 | alpha 0.60 (marginal)
higher, in each band, than matched-coupling contacts reach by chance (no reshuffle of 300
met it), holding in 7/6/7/6 of 10 patients. And because |ImCoh| is zero-lag-immune, this is
not proximity: given seizure contacts on some electrodes as seeds, the same delta affinity
ranks the seizure contacts on OTHER electrodes above healthy ones (off-shaft AUC 0.72, 8/10,
label-shuffle collapses to chance) -- the case where distance to a known contact is useless.

Single panel, strength-RESIDUAL AUC on x (0.5 = strength-alone / chance): top group = the
four all-contacts bands; bottom (separated) = the off-shaft delta demonstration.

Reads : data/audit/epi_marker_allcontacts/{allcontacts_per_patient,allcontacts_stratnull_per_patient,allcontacts_nulls}.csv
        data/audit/epi_marker_library/marker_library_per_patient.csv
        data/audit/epi_marker_library/verify_nulls.csv
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
LIB = ROOT / "data/audit/epi_marker_library"
OUT = ROOT / "data/reports/results_section3/fig_epi_a_relational_marker.pdf"

MARKER = "heat_t5"                                   # slow heat-kernel (README headline)
BANDS = ["delta", "low_gamma", "beta", "alpha"]      # descending median AUC
C_BAND = {"delta": "#08519c", "low_gamma": "#2b8cbe",
          "beta": "#d1491c", "alpha": "#807dba"}
C_BEAT, C_MISS = "#3a9a4f", "#b9b9b9"                 # dot: beats / within own matched-strength null
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
Y_OFF = dict(zip(COHORT, np.linspace(-0.26, 0.26, len(COHORT))))


def main():
    pp = pd.read_csv(AC / "allcontacts_per_patient.csv")
    pp = pp[pp.marker == MARKER]
    sn = pd.read_csv(AC / "allcontacts_stratnull_per_patient.csv")
    nulls = pd.read_csv(AC / "allcontacts_nulls.csv").set_index("band")
    lib = pd.read_csv(LIB / "marker_library_per_patient.csv")
    lib = lib[(lib.marker == MARKER) & (lib.band == "delta")].set_index("patient")
    vnull = pd.read_csv(LIB / "verify_nulls.csv").set_index("marker").loc[MARKER]

    # y layout: 4 all-contacts bands (top) + 1 off-shaft delta row (bottom, separated)
    rows = list(BANDS) + ["__offshaft__"]
    yof = {"delta": 5.4, "low_gamma": 4.4, "beta": 3.4, "alpha": 2.4,
           "__offshaft__": 0.9}

    fig, ax = plt.subplots(figsize=(9.4, 6.6))

    # matched-strength null envelope (label-shuffle): median ~0.49, p95 ~0.55-0.57
    nmed = float(nulls.null_median_mean.mean())
    np95 = float(nulls.null_median_p95.max())
    ax.axvspan(nmed, np95, color="0.82", alpha=0.6, zorder=0)
    ax.axvline(0.5, color="0.5", lw=1.1, ls="--", zorder=1)      # strength-alone / chance

    def draw_row(y, dots, beats, med, band_c, count_txt):
        for pat in dots.index:
            c = C_BEAT if bool(beats.get(pat, False)) else C_MISS
            ax.scatter(dots[pat], y + Y_OFF.get(pat, 0.0), s=95, color=c,
                       edgecolor="white", linewidth=0.8, zorder=4)
        ax.plot([0.5, med], [y, y], color=band_c, lw=3.0, zorder=5,
                solid_capstyle="round")
        ax.scatter(med, y, s=175, color=band_c, edgecolor="white", linewidth=1.3,
                   zorder=6)
        ax.text(med, y + 0.30, count_txt, color=band_c, fontsize=10.5, ha="center",
                va="bottom", fontweight="bold")

    # all-contacts bands
    for b in BANDS:
        dots = pp[pp.band == b].set_index("patient").auc_resid
        beats = sn[sn.band == b].set_index("patient").beats_p95
        med = float(nulls.loc[b, "real_median_auc"])
        n = int(beats.sum())
        draw_row(yof[b], dots, beats, med, C_BAND[b], f"{n}/10")

    # off-shaft delta demonstration
    off = lib.auc_resid
    off_beats = (off > float(vnull.null_median_p95))
    draw_row(yof["__offshaft__"], off, off_beats, float(vnull.real_median_auc),
             C_BAND["delta"], f"{int(off_beats.sum())}/10")

    # group separator + labels
    ax.axhline(1.65, color="0.7", lw=0.8, ls=":", zorder=1)
    ax.text(0.5 - 0.005, 6.15, "all contacts", fontsize=11, style="italic",
            color="0.35", ha="right")
    ax.text(0.5 - 0.005, 1.45, r"off-shaft ($\delta$): distance to a known contact useless",
            fontsize=10, style="italic", color="0.35", ha="right")

    ax.set_yticks([yof[b] for b in BANDS] + [yof["__offshaft__"]])
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS] +
                       [BRAIN_BAND_TEX_DICT["delta"] + r" off-shaft"])
    for tick, b in zip(ax.get_yticklabels(), BANDS + ["delta"]):
        tick.set_color(C_BAND[b])
        tick.set_fontsize(16)
    ax.set_ylim(0.2, 6.4)
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
               ms=12, lw=3.0, label="cohort median AUC ($n$/10 beat own null)"),
        Patch(facecolor="0.82", alpha=0.6,
              label="matched-strength null (label-shuffle, $0/300$ reached obs.)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.035),
               ncol=2, frameon=False, fontsize=10, handletextpad=0.5,
               columnspacing=1.6)

    fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.20)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:epi_a — relational SOZ marker (strength-residual, all-contacts + off-shaft)\n")
    for b in BANDS:
        beats = int(sn[sn.band == b].beats_p95.sum())
        print(f"  {b:11s} median AUC={nulls.loc[b,'real_median_auc']:.3f} "
              f"null={nulls.loc[b,'null_median_mean']:.3f} p95={nulls.loc[b,'null_median_p95']:.3f} "
              f"beats={beats}/10 p_emp={nulls.loc[b,'p_empirical']}")
    print(f"  off-shaft delta: median AUC={float(vnull.real_median_auc):.3f} "
          f"null={float(vnull.null_median_mean):.3f} beats={int((off>float(vnull.null_median_p95)).sum())}/10")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
