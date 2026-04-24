#!/usr/bin/env python3
"""Band typology scatter — ρ (H2c) vs Δρ (H2d) per band.

Visualizes the ρ vs Δρ dissociation: continuous drift direction (x) vs
block-level persistence (y) per band. Each band is one labeled point
with bootstrap 95 % CI crosshairs. θ highlighted as the ergodic-like
band. The fact that all points are in the top-right quadrant (both
positive) is the universality; the vertical spread is the
band-heterogeneity.

Reads:
    data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv
    data/reports/imcoh_vi/h2d_persistence_raw.csv

Writes:
    data/reports/imcoh_vi/figures/band_typology.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path
import sys as _sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

_sys.path.insert(0, str(Path(__file__).parent))
from _shared import boot_ci_mean  # noqa: E402


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

BAND_COLORS = {
    "delta":      "#7c4b9e",
    "theta":      "#c44e4e",
    "alpha":      "#4e6cc4",
    "beta":       "#2c8a5c",
    "low_gamma":  "#d58a2b",
    "high_gamma": "#8a5a2b",
}


def main() -> None:
    # H2c ρ per patient × band (target = task_test, n = 9)
    rho_df = pd.read_csv(IN_DIR / "h2c_ultrametric_drift_raw.csv")
    rho_df = rho_df[rho_df["has_rpre"] & rho_df["has_rpost"] & rho_df["has_ttest"]]
    rho_pat = (
        rho_df.pivot(index="patient", columns="band", values="rho_task")
        .reindex(columns=BANDS)
    )

    # H2d Δρ per patient × band (k-averaged)
    dr_raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    dr_pat = (
        dr_raw.groupby(["patient", "band"])["delta_rho"].mean()
        .unstack("band").reindex(columns=BANDS)
    )

    # Align patient sets
    common = rho_pat.index.intersection(dr_pat.index)
    rho_pat = rho_pat.loc[common]
    dr_pat = dr_pat.loc[common]

    # Bootstrap means + CIs
    x_mean, y_mean = [], []
    x_lo, x_hi = [], []
    y_lo, y_hi = [], []
    for b in BANDS:
        _, xl, xh = boot_ci_mean(rho_pat[b].to_numpy(), B=10000)
        _, yl, yh = boot_ci_mean(dr_pat[b].to_numpy(), B=10000)
        x_mean.append(rho_pat[b].mean())
        y_mean.append(dr_pat[b].mean())
        x_lo.append(xl); x_hi.append(xh)
        y_lo.append(yl); y_hi.append(yh)

    x_mean = np.array(x_mean); y_mean = np.array(y_mean)
    x_err = np.vstack([x_mean - np.array(x_lo), np.array(x_hi) - x_mean])
    y_err = np.vstack([y_mean - np.array(y_lo), np.array(y_hi) - y_mean])

    # ─── figure ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.4, 6.4), dpi=150)

    for i, b in enumerate(BANDS):
        color = BAND_COLORS[b]
        is_theta = (b == "theta")
        edge = "black" if is_theta else color
        size = 200 if is_theta else 130
        ax.errorbar(
            x_mean[i], y_mean[i],
            xerr=x_err[:, i:i+1], yerr=y_err[:, i:i+1],
            fmt="none", ecolor=color, elinewidth=1.0, capsize=3, alpha=0.55,
            zorder=2,
        )
        ax.scatter(
            x_mean[i], y_mean[i],
            s=size, color=color, edgecolor=edge, linewidth=1.5 if is_theta else 0.8,
            zorder=3,
        )
        ax.annotate(
            TEX[b],
            xy=(x_mean[i], y_mean[i]),
            xytext=(8, 6), textcoords="offset points",
            fontsize=14, fontweight="bold" if is_theta else "normal",
            color=edge,
        )

    # Inter-band medians as quadrant lines
    xm = np.median(x_mean); ym = np.median(y_mean)
    ax.axvline(xm, color="#bbb", linestyle=":", linewidth=1.0, zorder=1)
    ax.axhline(ym, color="#bbb", linestyle=":", linewidth=1.0, zorder=1)
    ax.text(xm + 0.002, ax.get_ylim()[1] - 0.01,
            f" median ρ = {xm:.2f}", fontsize=7.5, color="#888",
            ha="left", va="top")
    ax.text(ax.get_xlim()[1] - 0.002, ym + 0.002,
            f"median Δρ = {ym:.2f} ", fontsize=7.5, color="#888",
            ha="right", va="bottom")

    # Quadrant labels (TRACE / ERGODIC / DEAD)
    def _quad_label(x, y, text, color):
        ax.text(
            x, y, text, fontsize=9, color=color, alpha=0.7,
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor=color,
                      linewidth=0.7, alpha=0.85),
        )

    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()
    _quad_label((xhi + xm) / 2, (yhi + ym) / 2, "TRACE\n(ρ↑, Δρ↑)", "#2c8a5c")
    _quad_label((xlo + xm) / 2, (yhi + ym) / 2, "atypical\n(ρ↓, Δρ↑)", "#888888")
    _quad_label((xhi + xm) / 2, (ylo + ym) / 2, "ERGODIC\n(ρ↑, Δρ↓)", "#c44e4e")
    _quad_label((xlo + xm) / 2, (ylo + ym) / 2, "weak\n(ρ↓, Δρ↓)", "#888888")

    ax.set_xlabel(r"$\rho$  —  H2c continuous drift direction (Spearman)",
                  fontsize=11)
    ax.set_ylabel(r"$\Delta\rho$  —  H2d block-level persistence", fontsize=11)
    ax.set_title(
        "Band typology: continuous drift vs block persistence (n = 9)\n"
        "error bars = 95 % bootstrap CI; medians split the plane for visual classification",
        fontsize=11,
    )
    ax.grid(True, linestyle=":", alpha=0.4)

    # Caption
    caption = (
        "All bands lie in the upper-right quadrant (both ρ and Δρ positive = universal task-trace). "
        "θ stands apart with the smallest Δρ at comparable ρ (ergodic-like). "
        "α shows the strongest block persistence (Bonferroni θ-vs-α p = 0.010★)."
    )
    fig.text(0.5, -0.02, caption, ha="center", va="top",
             fontsize=8.5, color="#333", wrap=True)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"band_typology.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
