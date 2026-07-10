#!/usr/bin/env python3
"""fig_drift_infspec_duration_matched — real vs duration-matched drift (audit_168).

Cophenetic arc functionals per band, real (task genuinely happened) vs a drift arc
built from 5 temporally-ordered slots of pure pre-task rest_pre, slot lengths
proportional to the real phase durations. Two panels:

  LEFT  whole-task T_test: real competes with the (pessimistic, single-recording)
        drift upper bound; the FAIR full-duration locked C2 null is what certifies
        beta whole-task drift-clean (0.0137) -- audit_168's over-aggressive drift
        only draws level, it does not overturn C2.
  RIGHT inference-specific T_infspec_pe: real sits BELOW drift in every band (drift
        reproduces and exceeds it) -> the N2 inference-specific decomposition is
        drift-unverified / drift-confounded. No fair full-duration null is
        constructible (5 phase durations sum to > rest_pre).

Reads data/audit/duration_matched_drift_infspec/per_patient_per_band.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

SRC = ROOT / "data/audit/duration_matched_drift_infspec/per_patient_per_band.csv"
OUT = FIGURES_ROOT / "drift_infspec_duration_matched.pdf"
BANDS = list(BRAIN_BANDS)


def main():
    use_lrg_style()
    df = pd.read_csv(SRC)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.9, bottom=0.14, wspace=0.24)

    for ax, (rcol, dcol, ttl) in zip(
        (axL, axR),
        (("coph_real_Ttest", "coph_drift_Ttest", r"whole-task  $T_{\mathrm{test}}$"),
         ("coph_real_ispe", "coph_drift_ispe", r"inference-specific  $T_{\mathrm{infspec\cdot e}}$"))):
        x = np.arange(len(BANDS))
        for i, b in enumerate(BANDS):
            s = df[df.band == b]
            col = band_color(b)
            # per-patient points (real filled, drift open) + median markers
            ax.scatter(np.full(len(s), i) - 0.16 + np.linspace(-0.05, 0.05, len(s)),
                       s[rcol], s=16, color=col, alpha=0.55, edgecolor="none", zorder=2)
            ax.scatter(np.full(len(s), i) + 0.16 + np.linspace(-0.05, 0.05, len(s)),
                       s[dcol], s=16, facecolor="none", edgecolor=col, linewidth=0.7, zorder=2)
            ax.plot([i - 0.28, i - 0.04], [s[rcol].median()] * 2, color=col, lw=3, zorder=3)
            ax.plot([i + 0.04, i + 0.28], [s[dcol].median()] * 2, color=col, lw=3,
                    alpha=0.5, zorder=3)
        ax.axhline(0, color="0.6", lw=0.7, ls=(0, (4, 4)))
        ax.set_xticks(x)
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS])
        ax.set_ylabel(r"cophenetic $\rho_{\mathrm{sym}}$ arc functional")
        ax.set_title(ttl, fontsize=10)

    # legend proxy: filled = real, open = drift
    from matplotlib.lines import Line2D
    proxies = [Line2D([0], [0], marker="o", color="0.3", ls="none", label="real (task)"),
               Line2D([0], [0], marker="o", markerfacecolor="none", color="0.3",
                      ls="none", label="duration-matched drift")]
    fig.legend(handles=proxies, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=2, frameon=False, fontsize=8.5)
    axL.text(0.5, 0.96, "real ~ drift here, but fair C2 certifies\n"
             r"$\beta$ whole-task drift-clean (0.0137)", transform=axL.transAxes,
             ha="center", va="top", fontsize=7, color="0.35")
    axR.text(0.5, 0.96, "real < drift every band\n(drift reproduces + exceeds it)",
             transform=axR.transAxes, ha="center", va="top", fontsize=7, color="#b2182b")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {OUT}")


if __name__ == "__main__":
    main()
