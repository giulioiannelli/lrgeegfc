#!/usr/bin/env python3
"""Three-control consolidation figure for the continuous-trace matrix.

Reads:  data/reports/imcoh_continuous_trace/controls_summary.csv
Writes: data/reports/imcoh_continuous_trace/figures/controls_overview.pdf

Layout
------
A 2-row figure.

Row 1 — Per-band cohort dot-plots, one panel per band. Each panel shows
        (left to right) shared / split / null_drift / cross-probe /
        same-probe ρ values per patient. Cohort median + IQR per
        column. Reference line at 0. Reading: surviving columns are
        clouds above 0; controlled signal sits in the *split* /
        *cross-probe* columns.

Row 2 — Per-patient pairwise comparison ρ_split vs ρ_null_drift, one
        scatter per band. Diagonal reference. Above diagonal = signal
        above drift floor (good). Annotated with paired Wilcoxon p.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z


CT_DIR = REPORTS_ROOT / "imcoh_continuous_trace"
OUT_DIR = CT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

COLS = [
    ("rho_shared",            "shared",            "#1f77b4"),
    ("rho_split",             "split",             "#d62728"),
    ("rho_null_drift",        "drift",             "#7e7e7e"),
    ("rho_split_cross_probe", "cross",             "#2ca02c"),
    ("rho_split_same_probe",  "same",              "#9467bd"),
]


def _band_dotplot(ax, sub: pd.DataFrame, title: str) -> None:
    rng = np.random.default_rng(0)
    for j, (col, _, color) in enumerate(COLS):
        vals = sub[col].dropna().to_numpy()
        if len(vals) == 0:
            continue
        xs = j + rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(xs, vals, s=34, color=color, alpha=0.85,
                   edgecolor="black", linewidth=0.3, zorder=3)
        med = float(np.median(vals))
        q1, q3 = float(np.quantile(vals, 0.25)), float(np.quantile(vals, 0.75))
        ax.vlines(j, q1, q3, color="black", lw=1.6, zorder=4, alpha=0.7)
        ax.hlines(med, j - 0.24, j + 0.24, color="black", lw=2.0, zorder=5)
        n_pos = int((vals > 0).sum())
        ax.text(j, 1.01, f"{n_pos}/{len(vals)}", ha="center",
                va="bottom", fontsize=8, color=color,
                fontweight="bold", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels([c[1] for c in COLS], fontsize=8.5, rotation=20)
    ax.set_xlim(-0.5, len(COLS) - 0.5)
    ax.set_ylim(-0.55, 1.0)
    ax.set_title(title, fontsize=11)


def _pair_scatter(ax, sub: pd.DataFrame, title: str) -> None:
    d = sub.dropna(subset=["rho_split", "rho_null_drift"])
    if len(d) == 0:
        ax.set_title(title + "\n(no data)", fontsize=10)
        return
    x = d["rho_null_drift"].to_numpy()
    y = d["rho_split"].to_numpy()
    diff = y - x
    z, p = wilcoxon_z(diff) if len(diff) >= 3 else (np.nan, np.nan)
    n_above = int((y > x).sum())

    lo = min(x.min(), y.min()) - 0.05
    hi = max(x.max(), y.max()) + 0.05
    ax.plot([lo, hi], [lo, hi], color="#444", lw=0.8, ls="--", zorder=1)
    ax.axhline(0, color="#aaa", lw=0.5)
    ax.axvline(0, color="#aaa", lw=0.5)
    above = y > x
    ax.scatter(x[above], y[above], s=46, color="#2ca02c", alpha=0.9,
               edgecolor="black", linewidth=0.3, zorder=3,
               label="above drift")
    ax.scatter(x[~above], y[~above], s=46, color="#d62728", alpha=0.9,
               edgecolor="black", linewidth=0.3, zorder=3,
               label="below drift")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\rho_{\mathrm{drift}}$", fontsize=9)
    ax.set_ylabel(r"$\rho_{\mathrm{split}}$", fontsize=9)
    ttl = (rf"{title}"
           f"\n{n_above}/{len(d)} above   "
           f"$z={z:+.2f}$, $p={p:.3f}$")
    ax.set_title(ttl, fontsize=9.5)


def render() -> None:
    df = pd.read_csv(CT_DIR / "controls_summary.csv")

    fig = plt.figure(figsize=(18.0, 8.5), dpi=150)
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.0],
                          hspace=0.42, wspace=0.32)

    # Row 1 — per-band 5-column dotplots, ρ axis shared
    axes_top: list[plt.Axes] = []
    for j, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[0, j])
        sub = df[df["band"] == band]
        _band_dotplot(ax, sub, f"{TEX[band]}")
        axes_top.append(ax)
    axes_top[0].set_ylabel(r"$\rho$", fontsize=11)
    for ax in axes_top[1:]:
        ax.set_yticklabels([])
    axes_top[0].annotate(
        "Row A — Per-patient ρ across 5 estimators per band.\n"
        "Numbers above each column: patients with ρ>0.\n"
        "Black bar = cohort median; tick = IQR.",
        xy=(0.0, 1.18), xycoords="axes fraction",
        fontsize=9, ha="left",
    )

    # Row 2 — split vs drift scatter per band
    axes_bot: list[plt.Axes] = []
    for j, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[1, j])
        sub = df[df["band"] == band]
        _pair_scatter(ax, sub, TEX[band])
        axes_bot.append(ax)
    axes_bot[0].annotate(
        "Row B — Per-patient ρ_split vs ρ_drift (same noise regime).\n"
        "Above dashed diagonal = task-induced trace exceeds drift floor.\n"
        "Annotated: count above + paired Wilcoxon (split > drift).",
        xy=(0.0, 1.18), xycoords="axes fraction",
        fontsize=9, ha="left",
    )

    out = OUT_DIR / "controls_overview.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    render()
