#!/usr/bin/env python3
"""Per-band epi-vs-non-epi node behaviour under surrogate-z trace features.

2×3 band grid. Within each band panel, the three per-node surrogate-z trace
features (ρ_split_node, coherence_node, Grassmann mode-load); for each feature
two strips — non-epi (blue) and epi (red) — with median crossbars. The
surrogate expectation is z=0 (dashed). Reading rule: *a band+feature where the
red (epi) strip sits above the blue (non-epi) strip AND above 0 is a marker
candidate; if epi tracks non-epi, there is no trace-based separation.* The
quantitative per-patient AUC vs strength/probe baselines lives in the audit_79
table — this figure is the visual companion.

No in-axes numeric text (project rule). Reads ONLY the audit_79 feature cache.

Input  : data/audit/epi_marker/node_features_per_band.csv
Output : data/preprint/figures/all_bands/fig_epi_node_behavior_per_band.pdf
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
FEATURES = [("z_rho_split", r"$z\,\rho_{\mathrm{split}}$"),
            ("z_coherence", r"$z\,c$"),
            ("z_grass_modeload", r"$z\,m_{\mathrm{load}}$")]
NON_C = "#2c5f8a"
EPI_C = "#c0392b"


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_marker" / "node_features_per_band.csv"
    if not src.exists():
        raise SystemExit(f"[preprint_32] missing {src}; run audit_79 first")
    df = pd.read_csv(src)
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(2, 3, figsize=(13, 6.6), sharey=True)
    for bi, band in enumerate(BANDS):
        ax = axes[bi // 3][bi % 3]
        fb = df[df.band == band]
        for fi, (col, _lab) in enumerate(FEATURES):
            for gi, (is_epi, color) in enumerate([(False, NON_C), (True, EPI_C)]):
                vals = fb[fb.is_epi == is_epi][col].values
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    continue
                x0 = fi + (0.78 if gi else 0.22)
                jit = rng.uniform(-0.07, 0.07, size=vals.size)
                ax.scatter(np.full(vals.size, x0) + jit, vals, s=7,
                           color=color, alpha=0.45, edgecolor="none", zorder=2)
                med = float(np.median(vals))
                ax.plot([x0 - 0.13, x0 + 0.13], [med, med], color=color,
                        lw=2.4, zorder=3)
        ax.axhline(0, color="0.5", lw=0.7, ls="--", zorder=0)
        ax.set_xticks(range(len(FEATURES)))
        ax.set_xticklabels([lab for _, lab in FEATURES])
        ax.set_xlim(-0.3, len(FEATURES) - 0.1)
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        if bi % 3 == 0:
            ax.set_ylabel("surrogate-z")

    # clip extreme y for readability (keep most mass)
    ymax = np.nanpercentile(np.abs(df[[c for c, _ in FEATURES]].values), 99)
    for ax in axes.flat:
        ax.set_ylim(-ymax, ymax)

    handles = [Line2D([], [], marker="o", linestyle="none", color=NON_C,
                      markersize=7, label="non-epi node"),
               Line2D([], [], marker="o", linestyle="none", color=EPI_C,
                      markersize=7, label="epi node")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_node_behavior_per_band.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_32] -> {out}")


if __name__ == "__main__":
    main()
