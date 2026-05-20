#!/usr/bin/env python3
"""§5.7 Figure 1 — Raw |ImCoh| substrate ρ_split under matched-strength null.

For every band (δ, θ, α, β, γ_l, γ_h) overlay per-patient observed
substrate ρ_split against the matched-strength surrogate 5–95% strip and
the surrogate median. Cohort ratio + paired Wilcoxon p in panel titles.

Inputs:
    data/audit/raw_fc_matched_strength/per_patient_per_band_all_bands.csv
    data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv

Output:
    data/reports/section_5_matched_strength_refinement/figures/
        fig_1_raw_fc_rho_split_matched_strength.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt

from lrg_eegfc.visuals.styles import use_lrg_style

from _shared_ms import (
    ALL_BANDS, BAND_LABELS, CLR_OBS, CLR_SURR, CLR_SURR_FILL, CLR_ZERO,
    COHORT, FIG_DIR, RAW_FC_DIR, cohort_ratio, patient_short,
)


def main() -> Path:
    use_lrg_style()

    per_pat = pd.read_csv(RAW_FC_DIR / "per_patient_per_band_all_bands.csv")
    cohort = pd.read_csv(RAW_FC_DIR / "cohort_summary_all_bands.csv")

    bands = ALL_BANDS

    fig, axes = plt.subplots(1, len(bands),
                             figsize=(3.0 * len(bands), 4.6),
                             sharey=False)
    plt.subplots_adjust(wspace=0.55)

    for ax, band in zip(axes, bands):
        sub = per_pat[per_pat.band == band].copy()
        sub = sub.sort_values("obs_rho", ascending=True).reset_index(drop=True)
        ys = np.arange(len(sub))
        for y, (_, r) in zip(ys, sub.iterrows()):
            ax.hlines(y, r.surr_p5, r.surr_p95,
                      color=CLR_SURR_FILL, lw=4.0, alpha=0.95)
            ax.plot(r.surr_p50, y, marker="|", color=CLR_SURR,
                    markersize=9, mew=1.2)
            ax.plot(r.obs_rho, y, marker="o", color=CLR_OBS,
                    markersize=6.5, mec="white", mew=0.7)
        ax.set_yticks(ys)
        ax.set_yticklabels([patient_short(p) for p in sub.patient],
                           fontsize=7)
        ax.axvline(0, color=CLR_ZERO, lw=0.6, ls="--", zorder=0)
        ax.set_xlabel(r"$\rho_{\mathrm{split}}^{|\mathrm{ImCoh}|}$")
        cohort_row = cohort[cohort.band == band].iloc[0]
        ratio = cohort_ratio(cohort_row.obs_median_rho,
                              cohort_row.surr_median_rho_median)
        wp = float(cohort_row.paired_wilcoxon_p)
        n_above = str(cohort_row.n_above_surrogate)
        ax.set_title(
            f"{BAND_LABELS[band]}\n"
            rf"ratio $={ratio:.1f}\times$" "\n"
            rf"$p={wp:.3f}$,  $n_{{>}}={n_above}$",
            fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=7)

    # Figure-level legend (rule 1: shared legend lives on the figure).
    proxies = [
        plt.Line2D([0], [0], marker="o", color="none", mfc=CLR_OBS,
                   mec="white", mew=0.6, markersize=7,
                   label="observed"),
        plt.Line2D([0], [0], color=CLR_SURR_FILL, lw=5,
                   label="surrogate 5–95%"),
        plt.Line2D([0], [0], marker="|", color=CLR_SURR,
                   linestyle="none", markersize=9, mew=1.2,
                   label="surrogate median"),
    ]
    fig.legend(handles=proxies, loc="lower center",
               bbox_to_anchor=(0.5, -0.04),
               ncol=len(proxies), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 1))

    out = FIG_DIR / "fig_1_raw_fc_rho_split_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
