#!/usr/bin/env python3
"""Epi-stratified Grassmann verdict strips — k-resolved, per band.

2×3 grid of bands. Within each band panel, three horizontal lanes
(full / excl-epi / epi-only); each lane is a strip over spectral cutoff k
colored by the per-k sensitivity flag (config vs full baseline). The epi-only
lane visibly terminates early — the honest k-coverage collapse (needs n_epi>k;
≥8-patient frontier). k cells with <8 defined patients are drawn as
"underpowered" (not their unreliable flag).

No in-axes numeric text (project rule). Reads ONLY the audit_78 caches.

Input  : data/audit/epi_stratified/grassmann_cohort.csv
Output : data/preprint/figures/all_bands/fig_epi_stratified_verdict_matrix_grass.pdf
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
CONFIGS = ["full", "exclude_epi", "epi_only"]
CONFIG_LABEL = {"full": "full", "exclude_epi": "excl-epi", "epi_only": "epi-only"}
FLAG_COLOR = {"baseline": "#2c5f8a", "persist": "#1f7a1f", "weaken": "#c0392b",
              "emerge": "#e67e22", "absent": "#95a5a6"}
UNDERPOWERED = "#e0e0e0"
MIN_DEFINED = 8


def main() -> None:
    src = ROOT / "data" / "audit" / "epi_stratified" / "grassmann_cohort.csv"
    if not src.exists():
        raise SystemExit(f"[preprint_29] missing {src}; run audit_78 first")
    df = pd.read_csv(src)
    kmax = int(df.k.max())

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 5.2), sharex=True)
    for bi, band in enumerate(BANDS):
        ax = axes[bi // 3][bi % 3]
        for li, cfg in enumerate(CONFIGS):
            y = len(CONFIGS) - 1 - li
            sub = df[(df.band == band) & (df.config == cfg)].sort_values("k")
            for r in sub.itertuples():
                if r.n_defined >= MIN_DEFINED:
                    c = FLAG_COLOR.get(str(r.sensitivity_flag), "#d9d9d9")
                else:
                    c = UNDERPOWERED
                ax.add_patch(plt.Rectangle((r.k, y), 1, 0.82, facecolor=c,
                                           edgecolor="none"))
        ax.set_xlim(2, kmax + 1)
        ax.set_ylim(-0.1, len(CONFIGS))
        ax.set_yticks([len(CONFIGS) - 1 - i + 0.41 for i in range(len(CONFIGS))])
        ax.set_yticklabels([CONFIG_LABEL[c] for c in CONFIGS], fontsize=8)
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=11)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0, axis="y")
        if bi // 3 == 1:
            ax.set_xlabel(r"spectral cutoff $k$")

    handles = [Patch(facecolor=FLAG_COLOR[f], label=f)
               for f in ["baseline", "persist", "weaken", "emerge", "absent"]]
    handles.append(Patch(facecolor=UNDERPOWERED, label="underpowered (<8 pat)"))
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_epi_stratified_verdict_matrix_grass.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[preprint_29] -> {out}")


if __name__ == "__main__":
    main()
