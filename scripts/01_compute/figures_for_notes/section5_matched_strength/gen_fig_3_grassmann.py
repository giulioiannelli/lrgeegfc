#!/usr/bin/env python3
"""§5.7 Figure 3 — Grassmann principal-angle decomposition refined by
matched-strength + epi-zone-exclusion controls.

Re-uses the audit_46 cohort-median Δθ_i(k) triangular heatmap (the only
representation that exposes which spectral modes carry the trace
direction at each spectral cutoff k) and overlays the new
matched-strength + epi-X cohort verdicts on the bottom strip.

Layout (per band, 2×3 grid for δ θ α β γ_l γ_h):
  - top (triangular heatmap): cohort-median Δθ_i(k) at mode index i
    (rows) and spectral cutoff k (columns). Red = trace direction
    (post-task subspace is closer to task than pre-task is).
  - bottom strip: two lines —
      • black: `n_below_own_surrogate / 10` under matched-strength
        (audit_66) per k
      • rust:  `n_below_own_surrogate / 10` under epi-zone exclusion
        (audit_67) per k
    Cells where the audit_66 cohort-paired Wilcoxon `p < 0.05` are
    tick-marked in green at the top of the strip; the manuscript
    matched-strength-surviving windows (β k=27..55, γ_l k=12..23,
    γ_h k=19..27) are shaded green on both the heatmap and the strip.

Inputs:
    data/audit/section5_v2_round3_redo/tables/grassmann_principal_angles_cohort.csv
    data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv
    data/audit/grassmann_epi_exclusion/cohort_summary.csv

Output:
    data/reports/section_5_matched_strength_refinement/figures/
        fig_3_grassmann_matched_strength.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from lrg_eegfc.visuals.styles import use_lrg_style

from _shared_ms import (
    ALL_BANDS, BAND_LABELS, CLR_EPIX, CLR_HIGHLIGHT, FIG_DIR,
    GRASSMANN_DIR, GRASSMANN_EPIX_DIR,
)


# Cohort principal-angle CSV (audit_46), already at COHORT_N10 / imcoh_abs.
PRINCIPAL_CSV = (ROOT / "data" / "audit" / "section5_v2_round3_redo"
                 / "tables" / "grassmann_principal_angles_cohort.csv")

# k-grid the heatmap covers (matches audit_46).
K_GRID = list(range(2, 81))
K_MIN, K_MAX = K_GRID[0], K_GRID[-1]

# Manuscript matched-strength-surviving windows from the errata corrige.
MANUSCRIPT_WINDOW: dict[str, tuple[int, int] | None] = {
    "delta": None,
    "theta": None,
    "alpha": None,
    "beta": (27, 55),
    "low_gamma": (12, 23),
    "high_gamma": (19, 27),
}


def _build_heatmap(cohort_pa: pd.DataFrame, band: str) -> np.ndarray:
    H = np.full((K_MAX, len(K_GRID)), np.nan)
    sub = cohort_pa[cohort_pa.band == band]
    for _, row in sub.iterrows():
        k = int(row.k)
        i = int(row.mode_index) - 1
        if k < K_MIN or k > K_MAX or i >= K_MAX:
            continue
        j = K_GRID.index(k)
        H[i, j] = row.delta_theta_median_radians
    return H


def _ms_series(cohort_ms: pd.DataFrame, band: str
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (k, n_below/10, p) clipped to K_GRID."""
    sub = cohort_ms[(cohort_ms.band == band)
                    & (cohort_ms.k >= K_MIN)
                    & (cohort_ms.k <= K_MAX)].sort_values("k")
    ks = sub.k.values
    n = sub.n_patients_below_own_surrogate.values / 10.0
    p = sub.paired_wilcoxon_p.values
    return ks, n, p


def main() -> Path:
    use_lrg_style()

    pa_cohort = pd.read_csv(PRINCIPAL_CSV)
    ms_cohort = pd.read_csv(GRASSMANN_DIR / "cohort_summary.csv")
    epiX_cohort = pd.read_csv(GRASSMANN_EPIX_DIR / "cohort_summary.csv")

    # Global colour scale: 95th percentile of |Δθ| over all bands+cells.
    vmax = float(np.percentile(
        pa_cohort.delta_theta_median_radians.abs().values, 95))
    if vmax < 0.05:
        vmax = 0.05

    fig = plt.figure(figsize=(14.4, 10.5))
    outer = GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.30)

    panel_order = [
        ("delta", 0, 0), ("theta", 0, 1), ("alpha", 0, 2),
        ("beta", 1, 0), ("low_gamma", 1, 1), ("high_gamma", 1, 2),
    ]

    im_handle = None
    for band, r, c in panel_order:
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3.4, 1.15],
                                        hspace=0.16)
        ax_h = fig.add_subplot(inner[0, 0])
        ax_s = fig.add_subplot(inner[1, 0], sharex=ax_h)

        win = MANUSCRIPT_WINDOW[band]
        if win is not None:
            ax_h.set_facecolor("#fff7e6")
            ax_s.set_facecolor("#fff7e6")

        H = _build_heatmap(pa_cohort, band)
        im = ax_h.imshow(
            H, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            origin="upper", interpolation="nearest",
            extent=(K_MIN - 0.5, K_MAX + 0.5, K_MAX + 0.5, 0.5),
        )
        im_handle = im

        # Diagonal (mode index = k)
        ax_h.plot([K_MIN - 0.5, K_MAX + 0.5],
                  [K_MIN - 0.5, K_MAX + 0.5],
                  color="0.55", lw=0.5, ls="--", zorder=4)

        # Shade manuscript window on heatmap (vertical band).
        if win is not None:
            ax_h.axvspan(win[0] - 0.5, win[1] + 0.5,
                         facecolor=CLR_HIGHLIGHT, alpha=0.13,
                         zorder=3)
            ax_s.axvspan(win[0] - 0.5, win[1] + 0.5,
                         facecolor=CLR_HIGHLIGHT, alpha=0.13,
                         zorder=0)

        # Heatmap ticks
        x_ticks = [2, 10, 20, 30, 40, 50, 60, 70, 80]
        ax_h.set_xticks(x_ticks)
        ax_h.set_xticklabels([str(k) for k in x_ticks], fontsize=8)
        y_ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80]
        ax_h.set_yticks(y_ticks)
        ax_h.set_yticklabels([str(t) for t in y_ticks], fontsize=8)
        ax_h.set_ylabel(r"mode index $i$", fontsize=9)
        plt.setp(ax_h.get_xticklabels(), visible=False)

        suffix = (rf"  (matched-strength window $k=[{win[0]},{win[1]}]$)"
                  if win is not None else "")
        ax_h.set_title(f"{BAND_LABELS[band]}{suffix}",
                       fontsize=11, pad=4)
        for spine in ax_h.spines.values():
            spine.set_linewidth(0.6)

        # ----- bottom strip -----
        ks_ms, n_ms, p_ms = _ms_series(ms_cohort, band)
        ks_ex, n_ex, _ = _ms_series(epiX_cohort, band)

        ax_s.plot(ks_ms, n_ms, color="#1d1d1d", lw=1.2,
                  label="matched-strength")
        if ks_ex.size:
            ax_s.plot(ks_ex, n_ex, color=CLR_EPIX, lw=1.2,
                      label="epi-excluded")
        ax_s.axhline(0.7, color="0.5", lw=0.7, ls="--", zorder=1)
        ax_s.set_xlim(K_MIN - 0.5, K_MAX + 0.5)
        ax_s.set_ylim(-0.06, 1.08)
        ax_s.set_yticks([0.0, 0.5, 1.0])
        ax_s.set_yticklabels(["0", "5/10", "10/10"], fontsize=8)
        ax_s.set_xticks(x_ticks)
        ax_s.set_xticklabels([str(k) for k in x_ticks], fontsize=8)
        ax_s.set_xlabel(r"spectral cutoff $k$", fontsize=9)
        ax_s.set_ylabel(r"$n_{<\,\mathrm{surr}}$", fontsize=9)
        for spine in ax_s.spines.values():
            spine.set_linewidth(0.6)

        # Tick markers at top of strip wherever audit_66 cohort p<0.05
        sig = p_ms < 0.05
        if sig.any():
            ax_s.scatter(ks_ms[sig], np.full(sig.sum(), 1.04),
                         marker="|", color=CLR_HIGHLIGHT, s=18,
                         zorder=3,
                         label=r"$p_{\mathrm{cohort}}<0.05$")

        if band == "delta":
            ax_s.legend(loc="lower right", frameon=False, fontsize=7)

    # Shared colorbar — right edge of figure
    cax = fig.add_axes([0.94, 0.22, 0.012, 0.55])
    cb = fig.colorbar(im_handle, cax=cax)
    cb.set_label(
        r"cohort-median $\Delta\theta_i$ (rad) — "
        r"red $=$ trace direction (post closer to task than pre is)",
        fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.06)

    out = FIG_DIR / "fig_3_grassmann_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
