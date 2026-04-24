#!/usr/bin/env python3
"""H2c — per-patient × per-band Spearman ρ heatmap (target = task_test).

Evidence-at-a-glance figure for the H2c universal-direction result:
9 patients × 6 bands = 54 cells, 53 of which are positive. The single
negative cell (Pat_15 γ_h = −0.09) is visually boxed.

Reads:
    data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv

Writes:
    data/reports/imcoh_vi/figures/h2c_unanimity.{pdf,png}
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
    df = pd.read_csv(IN_DIR / "h2c_ultrametric_drift_raw.csv")
    df = df[df["has_rpre"] & df["has_rpost"] & df["has_ttest"]]

    mat = (
        df.pivot(index="patient", columns="band", values="rho_task")
        .reindex(columns=BANDS)
    )
    patients = list(mat.index)
    values = mat.to_numpy()

    # Diverging scale centred on 0, span informed by data
    vmin = min(-0.15, float(np.nanmin(values)))
    vmax = max(+1.00, float(np.nanmax(values)))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=150)
    im = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients, fontsize=10)

    # Annotate each cell with its ρ
    n_neg = 0
    for i, pat in enumerate(patients):
        for j, b in enumerate(BANDS):
            v = values[i, j]
            if not np.isfinite(v):
                continue
            text_color = "white" if v > 0.6 else "black"
            ax.text(j, i, f"{v:+.2f}",
                    ha="center", va="center", fontsize=9, color=text_color)
            if v < 0:
                # Frame the (very few) negative cells
                ax.add_patch(plt.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    fill=False, edgecolor="black", linewidth=2.0, zorder=4,
                ))
                n_neg += 1

    n_total = int(np.isfinite(values).sum())
    n_pos = n_total - n_neg

    cb = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label(r"Spearman $\rho(\Delta_{\mathrm{rest}},\, \Delta_{\mathrm{task}})$",
                 fontsize=10)

    ax.set_title(
        f"H2c — per-patient drift-direction correlation (target = task_test, n = 9)\n"
        f"{n_pos}/{n_total} patient × band cells positive;  all 6 bands q < 0.005 FDR",
        fontsize=11,
    )
    ax.set_xlabel("frequency band", fontsize=11)
    ax.set_ylabel("patient", fontsize=11)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"h2c_unanimity.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
