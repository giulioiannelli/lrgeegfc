#!/usr/bin/env python3
"""H2d — per-patient × per-band Δρ heatmap (k-averaged over k ∈ [2, 49]).

Evidence-at-a-glance figure for the H2d band-heterogeneity result:
9 patients × 6 bands of Δρ values. The per-patient view behind the
bar chart — reveals individual variability and shows that θ's
lower mean Δρ is consistent across most patients, not driven by
outliers.

Reads:
    data/reports/imcoh_vi/h2d_persistence_raw.csv

Writes:
    data/reports/imcoh_vi/figures/h2d_unanimity.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT


def main() -> None:
    raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    dr_pat = (
        raw.groupby(["patient", "band"])["delta_rho"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    patients = list(dr_pat.index)
    values = dr_pat.to_numpy()

    vmin = min(-0.10, float(np.nanmin(values)))
    vmax = max(+0.60, float(np.nanmax(values)))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=150)
    im = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients, fontsize=10)

    n_neg = 0
    for i, pat in enumerate(patients):
        for j, b in enumerate(BANDS):
            v = values[i, j]
            if not np.isfinite(v):
                continue
            text_color = "white" if v > 0.35 else "black"
            ax.text(j, i, f"{v:+.2f}",
                    ha="center", va="center", fontsize=9, color=text_color)
            if v < 0:
                ax.add_patch(plt.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    fill=False, edgecolor="black", linewidth=2.0, zorder=4,
                ))
                n_neg += 1

    # Column-mean annotation under the heatmap
    band_means = np.nanmean(values, axis=0)
    for j, (b, m) in enumerate(zip(BANDS, band_means)):
        ax.text(j, len(patients) - 0.35,
                f"mean\n{m:+.2f}",
                ha="center", va="top", fontsize=9,
                transform=ax.transData,
                color="#333", fontweight="bold" if b == "theta" else "normal")

    cb = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label(r"$\Delta\rho = \rho_{\mathrm{task}} - \rho_{\mathrm{inert}}$",
                 fontsize=10)

    n_total = int(np.isfinite(values).sum())
    n_pos = n_total - n_neg

    ax.set_title(
        f"H2d — per-patient block-persistence Δρ (k-averaged, k ∈ [2, 49], n = 9)\n"
        f"{n_pos}/{n_total} patient × band cells positive;  column-mean θ lowest across the cohort",
        fontsize=11,
    )
    ax.set_xlabel("frequency band", fontsize=11)
    ax.set_ylabel("patient", fontsize=11)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"h2d_unanimity.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
