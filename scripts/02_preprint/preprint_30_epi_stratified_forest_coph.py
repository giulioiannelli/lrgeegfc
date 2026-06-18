#!/usr/bin/env python3
"""Per-config cophenetic ρ_split forest — per-patient dots vs surrogate.

2×3 band grid. Within each band panel, the six configs on the x-axis; for each
config the per-patient observed ρ_split (filled dots, jittered), the per-patient
matched-strength surrogate median (grey ×, the null cloud), and the cohort
median (wide dash). Pat_13 (30 epi contacts, cohort-max) is drawn as a ringed
marker so its leverage on the epi configs is visible; the epi-only lane's wide
scatter is the visual proof of its underpower.

No in-axes numeric text (project rule). Reads ONLY the audit_77 per-patient cache.

Input  : data/audit/epi_stratified/cophenetic_per_patient.csv
Output : data/preprint/figures/all_bands/fig_epi_stratified_forest_coph.pdf
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

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CONFIG_ORDER = ["full", "nonepi_nonepi", "exclude_epi", "cross",
                "epi_epi", "epi_only"]
CONFIG_LABEL = {"full": "full", "nonepi_nonepi": "non–non",
                "exclude_epi": "excl-epi", "cross": "cross",
                "epi_epi": "epi–epi", "epi_only": "epi-only"}
OBS_C = "#1f3d6e"
SURR_C = "#999999"
MED_C = "#c0392b"


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_stratified" / "cophenetic_per_patient.csv"
    if not src.exists():
        raise SystemExit(f"[preprint_30] missing {src}; run audit_77 first")
    df = pd.read_csv(src)
    df = df[df.defined]
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(2, 3, figsize=(13, 6.4), sharey=True)
    for bi, band in enumerate(BANDS):
        ax = axes[bi // 3][bi % 3]
        for ci, cfg in enumerate(CONFIG_ORDER):
            sub = df[(df.band == band) & (df.config == cfg)]
            if sub.empty:
                continue
            jit = rng.uniform(-0.16, 0.16, size=len(sub))
            is13 = (sub.patient == "Pat_13").values
            # surrogate null centre (grey x)
            ax.scatter(ci + jit, sub.surr_p50.values, s=14, marker="x",
                       color=SURR_C, linewidths=0.8, zorder=2)
            # observed (filled), Pat_13 ringed
            ax.scatter((ci + jit)[~is13], sub.obs_stat.values[~is13], s=24,
                       color=OBS_C, edgecolor="white", linewidth=0.4, zorder=3)
            ax.scatter((ci + jit)[is13], sub.obs_stat.values[is13], s=52,
                       facecolor=OBS_C, edgecolor=MED_C, linewidth=1.6,
                       zorder=4)
            # cohort median (wide dash)
            med = float(np.median(sub.obs_stat.values))
            ax.plot([ci - 0.32, ci + 0.32], [med, med], color=MED_C, lw=2.2,
                    zorder=5)
        ax.axhline(0, color="0.5", lw=0.7, ls="--", zorder=0)
        ax.set_xticks(range(len(CONFIG_ORDER)))
        ax.set_xticklabels([CONFIG_LABEL[c] for c in CONFIG_ORDER],
                           rotation=30, ha="right", fontsize=8)
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        if bi % 3 == 0:
            ax.set_ylabel(r"$\rho_{\mathrm{split}}$")

    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=OBS_C,
               markeredgecolor="white", markersize=7, label="observed (patient)"),
        Line2D([], [], marker="o", color="none", markerfacecolor=OBS_C,
               markeredgecolor=MED_C, markersize=9, markeredgewidth=1.6,
               label="Pat_13 (30 epi)"),
        Line2D([], [], marker="x", color=SURR_C, linestyle="none",
               markersize=7, label="surrogate median"),
        Line2D([], [], color=MED_C, lw=2.2, label="cohort median"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_stratified_forest_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_30] -> {out}")


if __name__ == "__main__":
    main()
