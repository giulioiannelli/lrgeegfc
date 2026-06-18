#!/usr/bin/env python3
"""Trace-localization gradient across epi pair-classes — cophenetic.

2×3 band grid. Per band, the cohort-median ρ_split across the three pair
classes in localization order — non–non → cross → epi–epi — connected as a
gradient, against the full-graph reference (dashed) and the surrogate-median
line (grey). Marker filled when the cohort verdict clears (separated /
persist vs full), open otherwise. This is the cleanest "where does the
established trace live relative to epileptic tissue" read, and doubles as the
per-pair taxonomy stratified by epi class.

No in-axes numeric text (project rule). Reads ONLY the audit_77 cohort cache.

Input  : data/audit/epi_stratified/cophenetic_cohort.csv
Output : data/preprint/figures/all_bands/fig_epi_stratified_pairclass_gradient_coph.pdf
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
PAIRCLASS = ["nonepi_nonepi", "cross", "epi_epi"]
PC_LABEL = {"nonepi_nonepi": "non–non", "cross": "cross", "epi_epi": "epi–epi"}
OBS_C = "#1f3d6e"
FULL_C = "#444444"
SURR_C = "#aa8888"
CLEARED = {"separated", "persist"}


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_stratified" / "cophenetic_cohort.csv"
    if not src.exists():
        raise SystemExit(f"[preprint_31] missing {src}; run audit_77 first")
    df = pd.read_csv(src)

    def cell(band, cfg):
        r = df[(df.band == band) & (df.config == cfg)]
        return None if r.empty else r.iloc[0]

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.2), sharex=True)
    xs = np.arange(len(PAIRCLASS))
    for bi, band in enumerate(BANDS):
        ax = axes[bi // 3][bi % 3]
        full = cell(band, "full")
        if full is not None and np.isfinite(full.obs_median):
            ax.axhline(full.obs_median, color=FULL_C, lw=1.1, ls="--",
                       zorder=1, label="full")
        ys, ss, filled = [], [], []
        for cfg in PAIRCLASS:
            r = cell(band, cfg)
            if r is None or not np.isfinite(r.obs_median):
                ys.append(np.nan); ss.append(np.nan); filled.append(False)
            else:
                ys.append(r.obs_median); ss.append(r.surr_median_median)
                filled.append(str(r.verdict) in CLEARED
                              or str(r.sensitivity_flag) in CLEARED)
        ys = np.array(ys); ss = np.array(ss)
        ax.plot(xs, ss, color=SURR_C, lw=1.0, ls=":", marker="_",
                zorder=2)
        ax.plot(xs, ys, color=OBS_C, lw=1.6, zorder=3)
        for x, y, f in zip(xs, ys, filled):
            if not np.isfinite(y):
                continue
            ax.scatter(x, y, s=58, zorder=4,
                       facecolor=OBS_C if f else "white",
                       edgecolor=OBS_C, linewidth=1.6)
        ax.axhline(0, color="0.6", lw=0.6, ls="-", zorder=0)
        ax.set_xticks(xs)
        ax.set_xticklabels([PC_LABEL[c] for c in PAIRCLASS])
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        if bi % 3 == 0:
            ax.set_ylabel(r"cohort median $\rho_{\mathrm{split}}$")

    handles = [
        Line2D([], [], color=OBS_C, lw=1.6, marker="o", markerfacecolor=OBS_C,
               markeredgecolor=OBS_C, label="pair-class (filled = cleared)"),
        Line2D([], [], color=OBS_C, lw=0, marker="o", markerfacecolor="white",
               markeredgecolor=OBS_C, label="not cleared"),
        Line2D([], [], color=FULL_C, lw=1.1, ls="--", label="full reference"),
        Line2D([], [], color=SURR_C, lw=1.0, ls=":", label="surrogate median"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_stratified_pairclass_gradient_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_31] -> {out}")


if __name__ == "__main__":
    main()
