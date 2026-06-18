#!/usr/bin/env python3
"""Headline epi-stratified verdict matrix — cophenetic substrate.

6 bands (rows) × 6 configs (cols) categorical heatmap of the cophenetic
ρ_split trace sensitivity flag (config vs full-graph baseline). Column order
tells the localization story left→right:

    full → non–non → excl-epi → cross → epi–epi → epi-only

Color = sensitivity_flag (persist / weaken / emerge / absent / baseline /
undefined). No in-axes numeric text (project rule): the full numeric verdict
table is printed to stdout. Reads ONLY the audit_77 cohort cache.

Input  : data/audit/epi_stratified/cophenetic_cohort.csv
Output : data/preprint/figures/all_bands/fig_epi_stratified_verdict_matrix_coph.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

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
FLAG_COLOR = {"baseline": "#2c5f8a", "persist": "#1f7a1f", "weaken": "#c0392b",
              "emerge": "#e67e22", "absent": "#95a5a6", "undefined": "#ffffff",
              "unknown": "#d9d9d9"}
FLAG_ORDER = ["baseline", "persist", "weaken", "emerge", "absent", "undefined"]


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_stratified" / "cophenetic_cohort.csv"
    if not src.exists():
        raise SystemExit(f"[preprint_28] missing {src}; run audit_77 first")
    df = pd.read_csv(src)

    # index by (band, config) for fast lookup
    flag = {(r.band, r.config): str(r.sensitivity_flag) for r in df.itertuples()}
    obsm = {(r.band, r.config): r.obs_median for r in df.itertuples()}
    pval = {(r.band, r.config): r.paired_wilcoxon_p for r in df.itertuples()}

    nb, nc = len(BANDS), len(CONFIG_ORDER)
    fig, ax = plt.subplots(figsize=(1.05 * nc + 1.6, 0.85 * nb + 1.2))

    # stdout numeric table (the rule: numbers to stdout, color to the figure)
    print(f"{'band':<11}" + "".join(f"{CONFIG_LABEL[c]:>11}" for c in CONFIG_ORDER))
    for bi, band in enumerate(BANDS):
        line = f"{band:<11}"
        for ci, cfg in enumerate(CONFIG_ORDER):
            f = flag.get((band, cfg), "undefined")
            color = FLAG_COLOR.get(f, "#d9d9d9")
            y = nb - 1 - bi
            if f == "undefined":
                ax.add_patch(plt.Rectangle((ci, y), 1, 1, facecolor="white",
                                           edgecolor="#bbbbbb", hatch="////",
                                           linewidth=0.5))
            else:
                ax.add_patch(plt.Rectangle((ci, y), 1, 1, facecolor=color,
                                           edgecolor="white", linewidth=1.2))
            om = obsm.get((band, cfg), np.nan)
            pv = pval.get((band, cfg), np.nan)
            cell = "n/d" if f == "undefined" else f"{om:+.2f}/p{pv:.3f}"
            line += f"{cell:>11}"
        print(line)

    ax.set_xlim(0, nc); ax.set_ylim(0, nb)
    ax.set_xticks(np.arange(nc) + 0.5)
    ax.set_xticklabels([CONFIG_LABEL[c] for c in CONFIG_ORDER], rotation=30,
                       ha="right")
    ax.set_yticks(np.arange(nb) + 0.5)
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in reversed(BANDS)])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    handles = [Patch(facecolor=FLAG_COLOR[f], edgecolor="white", label=f)
               if f != "undefined" else
               Patch(facecolor="white", edgecolor="#bbbbbb", hatch="////",
                     label="undefined")
               for f in FLAG_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=len(FLAG_ORDER),
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_stratified_verdict_matrix_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_28] -> {out}")


if __name__ == "__main__":
    main()
