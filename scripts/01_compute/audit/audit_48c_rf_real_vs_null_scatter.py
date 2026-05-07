#!/usr/bin/env python3
"""Audit 48c -- per-patient real-vs-null scatter for hard RF (theta=0.70) and
soft RF (mean best-Jaccard). The diagonal y = x is "real = null"; trace zone
is real > null (above the diagonal). The figure shows 9/10 (or 10/10) above
the diagonal for cohort-passing bands.

Reads from existing audit_48 + audit_48b CSVs; no recompute. ~80 lines.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT

BANDS = list(BRAIN_BANDS_NAMES)
BASE = ROOT / "data" / "reports" / "section_5_lrg_trace" / "14_rf_clade_persistence"
TBL = BASE / "tables"
FIG = BASE / "figures"

HARD = pd.read_csv(TBL / "Td_per_patient_per_band.csv")
HARD = HARD[HARD["threshold"] == 0.70]
SOFT = pd.read_csv(TBL / "Td_per_patient_per_band_soft.csv")


def panel(ax, df, real_col, null_col, band, title_fmt):
    sub = df[df["band"] == band]
    real = sub[real_col].values.astype(float)
    null = sub[null_col].values.astype(float)
    pats = sub["patient"].tolist()
    n_above = int(np.sum(real > null))

    xy = np.concatenate([real, null])
    xy = xy[~np.isnan(xy)]
    if len(xy) == 0:
        return
    lim = max(abs(xy.min()), abs(xy.max())) * 1.15
    if lim < 1e-3:
        lim = 0.05
    ax.fill_between([-lim, lim], [-lim, lim], lim, color="#ffe8d6",
                    alpha=0.55, zorder=0, label="real > null (trace)")
    ax.plot([-lim, lim], [-lim, lim], color="#777", lw=0.7, ls="--",
            zorder=1)
    ax.axhline(0, color="0.55", lw=0.4, ls=":", zorder=1)
    ax.axvline(0, color="0.55", lw=0.4, ls=":", zorder=1)
    for p, x, y in zip(pats, null, real):
        c = "#d62728" if p == "Pat_03" else "#1a4f73"
        ax.scatter(x, y, s=24, color=c, edgecolor="white", linewidth=0.5,
                   zorder=3)
        ax.annotate(p[-2:], (x, y), fontsize=6, alpha=0.75,
                    xytext=(2.5, 2.5), textcoords="offset points", zorder=4)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title_fmt.format(band=BRAIN_BAND_TEX_DICT[band], n=n_above),
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


fig, axes = plt.subplots(2, 6, figsize=(15, 5.6))
for j, band in enumerate(BANDS):
    panel(axes[0, j], HARD, "T_RF_cohort", "T_RF_null_cohort", band,
          "{band}  hard $\\theta{{=}}0.70$  ({n}/10)")
    panel(axes[1, j], SOFT, "T_RFsoft_cohort", "T_RFsoft_null_cohort", band,
          "{band}  soft mean-J  ({n}/10)")
for i, lbl in enumerate(["hard RF\nreal $T_{RF}$", "soft RF\nreal $T_{RFsoft}$"]):
    axes[i, 0].set_ylabel(lbl, fontsize=9)
for j in range(6):
    axes[1, j].set_xlabel("null  $T$ (within-baseline)", fontsize=8)

handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1a4f73",
               markersize=7, label="patient (Pat_05–15)"),
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728",
               markersize=7, label="Pat_03 (1024 Hz outlier)"),
    plt.Line2D([0], [0], color="#777", ls="--", lw=1.0, label="real = null"),
    plt.Rectangle((0, 0), 1, 1, color="#ffe8d6", alpha=0.7,
                  label="trace zone (real > null)"),
]
fig.legend(handles=handles, loc="lower center", ncol=4,
           bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(FIG / "rf_real_vs_null_hard_vs_soft.pdf")
plt.close(fig)

print(f"\n=== panel-by-panel n_above_null counts ===")
for band in BANDS:
    h = HARD[HARD["band"] == band]
    s = SOFT[SOFT["band"] == band]
    n_h = int(np.sum(h["T_RF_cohort"].values > h["T_RF_null_cohort"].values))
    n_s = int(np.sum(s["T_RFsoft_cohort"].values > s["T_RFsoft_null_cohort"].values))
    print(f"  {band:>11}:  hard={n_h}/10  soft={n_s}/10")
print(f"\n[audit_48c] figure -> {FIG / 'rf_real_vs_null_hard_vs_soft.pdf'}")
