#!/usr/bin/env python3
"""H2e — drift noise-floor figure: ρ_cross vs ρ_null_drift, Δρ likewise.

Two-row figure conveying the H2e result at a glance:

  Top — per-band strip of per-patient ρ:
        * ρ_cross (H2c Spearman on ultrametric Δ vectors)  — blue dots
        * ρ_null_drift (within-session split-half drift correlation) — grey dots
        * ρ_within (split-half reliability ceiling) — horizontal band
        Wilcoxon q-value (FDR m=6) annotated per band.

  Bottom — per-band Δρ_cross vs Δρ_null bar pairs:
        * Δρ_cross (H2d k-averaged) — coloured bar
        * Δρ_null (drift-triple Δρ)  — grey bar
        θ highlighted; FDR q on the contrast annotated.

Reads:
    data/reports/imcoh_vi/h2e_split_half_rho_raw.csv
    data/reports/imcoh_vi/h2e_split_half_drho_raw.csv
    data/reports/imcoh_vi/h2e_split_half.csv

Writes:
    data/reports/imcoh_vi/figures/h2e_drift_floor.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

CROSS_COLOR = "#2a5d9f"   # deep blue — H2c / H2d cross-phase
NULL_COLOR  = "#9a9a9a"   # neutral grey — drift null
CEIL_COLOR  = "#d7e8f6"   # pale band — reliability ceiling
THETA_COLOR = "#c44e4e"


# ─────────────────────────── top panel: ρ strip ───────────────────────────


def _q_star(q: float) -> str:
    if not np.isfinite(q):
        return ""
    if q < 0.001: return "***"
    if q < 0.01:  return "**"
    if q < 0.05:  return "★"
    return ""


def draw_rho_panel(ax: plt.Axes, rho_raw: pd.DataFrame, stats: pd.DataFrame) -> None:
    rng = np.random.default_rng(0)
    xpos = np.arange(len(BANDS))
    for i, band in enumerate(BANDS):
        sub = rho_raw[rho_raw["band"] == band]
        if sub.empty:
            continue
        # Reliability ceiling: shaded band between min/max of ρ_within across phases
        ceil_vals = np.concatenate([
            sub["rho_within_rpre"].dropna().to_numpy(),
            sub["rho_within_rpost"].dropna().to_numpy(),
        ])
        if ceil_vals.size:
            lo, hi = np.quantile(ceil_vals, [0.1, 0.9])
            ax.fill_between([i - 0.38, i + 0.38], lo, hi,
                            color=CEIL_COLOR, zorder=1, alpha=0.9)
        # Strip with horizontal jitter
        jit = rng.uniform(-0.10, 0.10, size=len(sub))
        ax.scatter(np.full(len(sub), i - 0.18) + jit, sub["rho_cross"],
                   s=26, color=CROSS_COLOR, edgecolor="black", linewidth=0.4,
                   zorder=3, label="ρ_cross" if i == 0 else None)
        ax.scatter(np.full(len(sub), i + 0.18) + jit, sub["rho_null_drift"],
                   s=26, color=NULL_COLOR, edgecolor="black", linewidth=0.4,
                   zorder=3, label="ρ_null (drift)" if i == 0 else None)
        # Mean markers
        ax.scatter(i - 0.18, sub["rho_cross"].mean(), s=90, marker="_",
                   color="black", linewidth=2.0, zorder=4)
        ax.scatter(i + 0.18, sub["rho_null_drift"].mean(), s=90, marker="_",
                   color="black", linewidth=2.0, zorder=4)
        # FDR star
        q = stats[(stats["metric"] == "rho") & (stats["band"] == band)]["q"]
        if not q.empty and np.isfinite(q.iloc[0]):
            star = _q_star(float(q.iloc[0]))
            y_top = max(sub["rho_cross"].max(), 0.75)
            ax.text(i, y_top + 0.08, star, ha="center", va="bottom",
                    fontsize=14, fontweight="bold", color="#222")

    ax.axhline(0, color="black", lw=0.6, alpha=0.6, zorder=2)
    ax.set_xticks(xpos)
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=13)
    ax.set_ylabel(r"Spearman $\rho$ on $\Delta$ ultrametric matrices", fontsize=11)
    ax.set_ylim(-0.4, 1.0)
    ax.set_title(
        "H2c ρ (blue) vs within-session split-half drift null (grey); "
        "pale band = reliability ceiling (10–90 % of within-phase ρ). "
        "Stars: FDR q on paired cross>null.",
        loc="left", fontsize=10,
    )
    ax.legend(loc="lower right", frameon=False, fontsize=9)


# ─────────────────────────── bottom panel: Δρ bar pairs ───────────────────────────


def draw_drho_panel(ax: plt.Axes, drho_raw: pd.DataFrame, stats: pd.DataFrame) -> None:
    # Per-patient k-averaged means
    per = (drho_raw.groupby(["patient", "band"])[["drho_cross", "drho_null"]]
           .mean().reset_index())
    xpos = np.arange(len(BANDS))
    w = 0.36
    for i, band in enumerate(BANDS):
        sub = per[per["band"] == band]
        if sub.empty:
            continue
        m_cross = sub["drho_cross"].mean()
        m_null  = sub["drho_null"].mean()
        # 95 % bootstrap CI on cross and null
        ci_cross = _boot_ci(sub["drho_cross"].to_numpy())
        ci_null  = _boot_ci(sub["drho_null"].to_numpy())
        col_cross = THETA_COLOR if band == "theta" else CROSS_COLOR
        ax.bar(i - w / 2, m_cross, w, color=col_cross, edgecolor="black",
               linewidth=0.6, label="Δρ_cross" if i == 0 else None)
        ax.bar(i + w / 2, m_null, w, color=NULL_COLOR, edgecolor="black",
               linewidth=0.6, label="Δρ_null (drift)" if i == 0 else None)
        ax.errorbar(i - w / 2, m_cross,
                    yerr=[[m_cross - ci_cross[0]], [ci_cross[1] - m_cross]],
                    fmt="none", ecolor="black", capsize=3, linewidth=0.8)
        ax.errorbar(i + w / 2, m_null,
                    yerr=[[m_null - ci_null[0]], [ci_null[1] - m_null]],
                    fmt="none", ecolor="black", capsize=3, linewidth=0.8)
        q = stats[(stats["metric"] == "drho") & (stats["band"] == band)]["q"]
        if not q.empty and np.isfinite(q.iloc[0]):
            star = _q_star(float(q.iloc[0]))
            ax.text(i, max(m_cross, m_null) + 0.03, star,
                    ha="center", va="bottom", fontsize=13,
                    fontweight="bold", color="#222")

    ax.axhline(0, color="black", lw=0.6, alpha=0.6, zorder=2)
    ax.set_xticks(xpos)
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=13)
    ax.set_ylabel(r"mean $\Delta\rho$ (k-averaged, $k\in[2,49]$)", fontsize=11)
    ax.set_ylim(0.0, 0.40)
    ax.set_title(
        "H2d block persistence (coloured) vs drift-triple Δρ null (grey); "
        "θ in red. Stars: FDR q on paired cross>null.",
        loc="left", fontsize=10,
    )
    ax.legend(loc="upper right", frameon=False, fontsize=9)


def _boot_ci(x: np.ndarray, B: int = 10_000) -> tuple[float, float]:
    rng = np.random.default_rng(0)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return (float("nan"), float("nan"))
    idx = rng.integers(0, len(x), size=(B, len(x)))
    means = x[idx].mean(axis=1)
    return tuple(float(v) for v in np.quantile(means, [0.025, 0.975]))


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    rho_raw  = pd.read_csv(IN_DIR / "h2e_split_half_rho_raw.csv")
    drho_raw = pd.read_csv(IN_DIR / "h2e_split_half_drho_raw.csv")
    stats    = pd.read_csv(IN_DIR / "h2e_split_half.csv")

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.8), dpi=150,
                              gridspec_kw={"height_ratios": [1.0, 0.85]})
    draw_rho_panel(axes[0], rho_raw, stats)
    draw_drho_panel(axes[1], drho_raw, stats)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"h2e_drift_floor.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
