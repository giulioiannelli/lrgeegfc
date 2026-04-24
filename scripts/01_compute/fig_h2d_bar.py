#!/usr/bin/env python3
"""H2d — per-band Δρ bar chart with Bonferroni post-hoc annotations.

Headline figure conveying the H2d result at a glance:
 - Mean Δρ per band with 95 % bootstrap CI whiskers (all 6 bands > 0 FDR).
 - θ bar drawn in a distinct color; significantly lower than α (p=0.010)
   and δ (p=0.049) per Bonferroni-corrected paired Wilcoxon.
 - FDR q-values (all < 0.005) and ρ_task / ρ_inert numbers annotated.

Reads:
    data/reports/imcoh_vi/h2d_persistence_raw.csv
    data/reports/imcoh_vi/h2_band_selectivity.csv

Writes:
    data/reports/imcoh_vi/figures/h2d_band_persistence.pdf
    data/reports/imcoh_vi/figures/h2d_band_persistence.png
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

# Canonical bootstrap from the shared H2 stats helpers
import sys as _sys
_sys.path.insert(0, str(Path(__file__).parent))
from _shared import boot_ci_mean  # noqa: E402


IN_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT

THETA_COLOR = "#c44e4e"  # distinct for the ergodic band
TRACE_COLOR = "#4e6cc4"
EDGE_COLOR = "#1a1a1a"


def main() -> None:
    # ─── per-patient×band Δρ (k-averaged across k ∈ [2, 49]) ──────────
    raw = pd.read_csv(IN_DIR / "h2d_persistence_raw.csv")
    dr_pat = (
        raw.groupby(["patient", "band"])["delta_rho"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    rho_task = (
        raw.groupby(["patient", "band"])["rho_task"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    rho_inert = (
        raw.groupby(["patient", "band"])["rho_inert"].mean()
        .unstack("band").reindex(columns=BANDS)
    )

    # ─── band-level stats: mean, bootstrap 95 % CI ────────────────────
    mean_dr = np.array([dr_pat[b].mean() for b in BANDS])
    ci_lo = np.empty(len(BANDS))
    ci_hi = np.empty(len(BANDS))
    for i, b in enumerate(BANDS):
        _, lo, hi = boot_ci_mean(dr_pat[b].to_numpy(), B=10000)
        ci_lo[i], ci_hi[i] = lo, hi

    mean_task = np.array([rho_task[b].mean() for b in BANDS])
    mean_inert = np.array([rho_inert[b].mean() for b in BANDS])

    # ─── post-hoc θ-vs-others (Bonferroni-corrected) ──────────────────
    bs_df = pd.read_csv(IN_DIR / "h2_band_selectivity.csv")
    h2d_row = bs_df[bs_df["measure"] == "H2d Δρ (block persistence)"].iloc[0]
    theta_bonf = {
        "delta":      float(h2d_row["theta_vs_delta_p_bonf"]),
        "alpha":      float(h2d_row["theta_vs_alpha_p_bonf"]),
        "beta":       float(h2d_row["theta_vs_beta_p_bonf"]),
        "low_gamma":  float(h2d_row["theta_vs_low_gamma_p_bonf"]),
        "high_gamma": float(h2d_row["theta_vs_high_gamma_p_bonf"]),
    }
    friedman_p = float(h2d_row["friedman_p"])

    # ─── figure ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=150)

    xs = np.arange(len(BANDS))
    colors = [THETA_COLOR if b == "theta" else TRACE_COLOR for b in BANDS]

    # Bars with asymmetric CI
    yerr = np.vstack([mean_dr - ci_lo, ci_hi - mean_dr])
    bars = ax.bar(
        xs, mean_dr, yerr=yerr,
        color=colors, edgecolor=EDGE_COLOR, linewidth=1.1,
        capsize=5, error_kw={"elinewidth": 1.0, "ecolor": EDGE_COLOR},
        width=0.72, zorder=2,
    )

    # ρ_inert reference line per band (baseline persistence)
    for i, (m_t, m_i) in enumerate(zip(mean_task, mean_inert)):
        ax.hlines(
            m_i, i - 0.36, i + 0.36,
            linestyles=":", colors="#555", linewidth=1.3, zorder=3,
        )

    # Numeric labels on top of each bar
    for i, (m, lo, hi) in enumerate(zip(mean_dr, ci_lo, ci_hi)):
        ax.text(
            i, hi + 0.008, f"{m:+.3f}",
            ha="center", va="bottom", fontsize=9, color=EDGE_COLOR,
        )

    # ρ_task / ρ_inert inside each bar
    for i, (m_t, m_i) in enumerate(zip(mean_task, mean_inert)):
        ax.text(
            i, 0.01,
            f"ρ_task={m_t:.2f}\nρ_inert={m_i:.2f}",
            ha="center", va="bottom", fontsize=7.5, color="white",
            fontweight="bold",
        )

    # Bonferroni brackets for θ vs α (p=0.010) and θ vs δ (p=0.049)
    theta_idx = BANDS.index("theta")
    alpha_idx = BANDS.index("alpha")
    delta_idx = BANDS.index("delta")

    def _draw_bracket(x1: int, x2: int, y: float, label: str) -> None:
        ax.plot([x1, x1, x2, x2], [y, y + 0.012, y + 0.012, y],
                color=EDGE_COLOR, linewidth=1.1, zorder=4)
        ax.text((x1 + x2) / 2, y + 0.014, label,
                ha="center", va="bottom", fontsize=9, zorder=4)

    top = max(ci_hi) + 0.03
    _draw_bracket(theta_idx, alpha_idx, top + 0.05,
                  rf"θ vs α  p$_{{\mathrm{{Bonf}}}}$={theta_bonf['alpha']:.3f} ★")
    _draw_bracket(theta_idx, delta_idx, top + 0.01,
                  rf"θ vs δ  p$_{{\mathrm{{Bonf}}}}$={theta_bonf['delta']:.3f} ★")

    # Cosmetics
    ax.set_xticks(xs)
    ax.set_xticklabels([TEX[b] for b in BANDS], fontsize=12)
    ax.set_ylabel(r"$\Delta\rho = \rho_{\mathrm{task}} - \rho_{\mathrm{inert}}$",
                  fontsize=12)
    ax.set_title(
        "H2d — block-level task-trace persistence per band (n = 9)\n"
        "all q < 0.005 FDR;  dashed line inside each bar = baseline spontaneous co-clustering rate $\\rho_{\\mathrm{inert}}$",
        fontsize=11,
    )
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_ylim(-0.02, top + 0.10)
    ax.grid(axis="y", linestyle=":", alpha=0.5, zorder=0)

    # Friedman context as a subtitle caption below the x-axis
    caption = (
        f"Friedman χ²(5) = {float(h2d_row['friedman_chi2']):.2f}, p = {friedman_p:.2f} "
        "(omnibus n.s.); band-heterogeneity carried by Bonferroni-corrected θ-vs-α (p=0.010★) "
        "and θ-vs-δ (p=0.049★) post-hoc paired Wilcoxon."
    )
    fig.text(0.5, -0.02, caption, ha="center", va="top",
             fontsize=8.5, color="#333", wrap=True)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"h2d_band_persistence.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
