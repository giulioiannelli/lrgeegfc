#!/usr/bin/env python3
"""Six-band cohort coherence composite at the per-pair cophenetic layer.

Each band's panel is split vertically into two strips that share the
x-axis (decile of within-patient rank Δ_task):

    Top strip  — COHORT CONSENSUS at high resolution (Q = 30 quantiles).
                 Single row, cohort-averaged signed trace-concordance per
                 quantile. Tight color saturation (±0.04) so β's modest
                 cohort-level deviations saturate to deep green while
                 null bands stay pale.

    Bottom grid — PER-PATIENT detail (10 rows × 10 deciles). Each cell
                 colored by signed trace-concordance at normal saturation
                 (±0.10) so per-patient direction and cohort consistency
                 are both visible.

The visual diagnostic separates two regimes:
    - Strong trace (β): top strip is solid deep green across all 30
      quantiles; bottom grid has most patient rows green-dominated.
    - Null (θ): top strip is pale / mottled across quantiles; bottom
      grid is a near-random checkerboard.

The intermediate bands (α, γ_l, γ_h, δ) sit between these regimes —
top strip lighter or mixed, bottom grid more random.

Layout: 2 rows × 3 cols, canonical band order:
    δ θ α
    β γ_l γ_h
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import rankdata

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

TRACE_CMAP = LinearSegmentedColormap.from_list(
    "trace_rwg",
    [
        "#5b1212",
        "#a02525",
        "#d8584c",
        "#f0a89e",
        "#fdf2ef",
        "#bfddc7",
        "#5fac7e",
        "#1a7c3e",
        "#0a3a1d",
    ],
    N=512,
)

N_PATIENT_DECILES = 10
N_COHORT_QUANTILES = 30
COHORT_VABS = 0.03
PATIENT_VABS = 0.10


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
    )


def compute_coherence(band: str) -> dict:
    """Per-patient × per-decile and cohort consensus at high resolution."""
    pat_edges = np.linspace(0.0, 1.0, N_PATIENT_DECILES + 1)
    pat_centers = (pat_edges[:-1] + pat_edges[1:]) / 2
    sign_pat = np.sign(pat_centers - 0.5)

    coh_edges = np.linspace(0.0, 1.0, N_COHORT_QUANTILES + 1)
    coh_centers = (coh_edges[:-1] + coh_edges[1:]) / 2
    sign_coh = np.sign(coh_centers - 0.5)

    per_pat = []
    per_pat_cohort = []
    labels = []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        t_rank = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        r_rank = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        bin_pat = np.minimum((t_rank * N_PATIENT_DECILES).astype(int),
                             N_PATIENT_DECILES - 1)
        bin_coh = np.minimum((t_rank * N_COHORT_QUANTILES).astype(int),
                             N_COHORT_QUANTILES - 1)
        row_pat = np.full(N_PATIENT_DECILES, np.nan)
        for k in range(N_PATIENT_DECILES):
            m = bin_pat == k
            if m.any():
                row_pat[k] = r_rank[m].mean()
        row_coh = np.full(N_COHORT_QUANTILES, np.nan)
        for k in range(N_COHORT_QUANTILES):
            m = bin_coh == k
            if m.any():
                row_coh[k] = r_rank[m].mean()
        per_pat.append((row_pat - 0.5) * sign_pat)
        per_pat_cohort.append((row_coh - 0.5) * sign_coh)
        labels.append(pat.replace("Pat_", "P"))

    return dict(
        patient=np.asarray(per_pat),
        cohort_strip=np.nanmean(np.asarray(per_pat_cohort), axis=0),
        labels=labels,
    )


def plot_band_composite(fig, outer_spec, band: str, payload: dict,
                        cohort_row: pd.Series,
                        show_y: bool, show_x: bool) -> tuple:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    sub = outer_spec.subgridspec(2, 1, height_ratios=[1, 4.5], hspace=0.07)
    ax_top = fig.add_subplot(sub[0])
    ax_bot = fig.add_subplot(sub[1])

    im_top = ax_top.imshow(
        payload["cohort_strip"][np.newaxis, :],
        origin="lower",
        cmap=TRACE_CMAP,
        vmin=-COHORT_VABS,
        vmax=COHORT_VABS,
        aspect="auto",
        interpolation="nearest",
    )
    ax_top.set_yticks([0])
    ax_top.set_yticklabels(["COHORT"], fontsize=8.0, fontweight="bold")
    ax_top.set_xticks([])
    for spine_name, spine in ax_top.spines.items():
        spine.set_edgecolor("0.10")
        spine.set_linewidth(1.4)
    for k in range(1, N_COHORT_QUANTILES):
        ax_top.axvline(k - 0.5, color="white", lw=0.25, alpha=0.55)
    ax_top.axvline(N_COHORT_QUANTILES / 2 - 0.5, color="0.05",
                   lw=1.4, alpha=0.95)

    coherence_p = payload["patient"]
    P = coherence_p.shape[0]
    im_bot = ax_bot.imshow(
        coherence_p,
        origin="lower",
        cmap=TRACE_CMAP,
        vmin=-PATIENT_VABS,
        vmax=PATIENT_VABS,
        aspect="auto",
        interpolation="nearest",
    )
    ax_bot.set_yticks(np.arange(P))
    if show_y:
        ax_bot.set_yticklabels(payload["labels"], fontsize=8.5)
        ax_bot.set_ylabel("patient", fontsize=10)
    else:
        ax_bot.set_yticklabels(payload["labels"], fontsize=8.5)
    ax_bot.set_xticks(np.arange(N_PATIENT_DECILES))
    ax_bot.set_xticklabels([str(k + 1) for k in range(N_PATIENT_DECILES)],
                           fontsize=8.5)
    if show_x:
        ax_bot.set_xlabel(r"within-patient decile of rank $\Delta_{\mathrm{task}}(i,j)$",
                         fontsize=10)
    for k in range(1, N_PATIENT_DECILES):
        ax_bot.axvline(k - 0.5, color="white", lw=0.55, alpha=0.7)
    for p in range(1, P):
        ax_bot.axhline(p - 0.5, color="white", lw=0.55, alpha=0.7)
    ax_bot.axvline(N_PATIENT_DECILES / 2 - 0.5, color="0.10",
                   lw=1.2, alpha=0.75)

    rho = float(cohort_row["obs_median_rho"])
    p_val = float(cohort_row["paired_wilcoxon_p"])
    coh_strip = payload["cohort_strip"]
    cohort_strip_mean = float(np.nanmean(coh_strip))
    n_strip_green = int((coh_strip > 0).sum())
    n_grid_green = int((coherence_p > 0).sum())
    is_sig = p_val < 0.05

    title_color = "#0a3a1d" if is_sig else "0.25"
    ax_top.text(0.5, 1.55, band_tex, transform=ax_top.transAxes,
                ha="center", va="bottom",
                fontsize=18, fontweight="bold", color=title_color)
    ax_top.text(0.5, 1.10,
                rf"$\rho_{{\mathrm{{split}}}}^{{\mathrm{{coph}}}} = {rho:+.3f}$  |  "
                rf"$p = {p_val:.3f}$  |  "
                rf"COHORT$\langle c \rangle = {cohort_strip_mean:+.3f}$  |  "
                rf"strip {n_strip_green}/{N_COHORT_QUANTILES}, "
                rf"grid {n_grid_green}/{coherence_p.size}",
                transform=ax_top.transAxes, ha="center", va="bottom",
                fontsize=9.0, color="0.20")

    if is_sig:
        for spine in ax_bot.spines.values():
            spine.set_edgecolor("#0a3a1d")
            spine.set_linewidth(2.4)
    else:
        for spine in ax_bot.spines.values():
            spine.set_edgecolor("0.60")
            spine.set_linewidth(0.8)

    return im_top, im_bot


def main() -> Path:
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    payloads = {}
    for band in BAND_ORDER:
        payloads[band] = compute_coherence(band)
        c = payloads[band]["patient"]
        cs = payloads[band]["cohort_strip"]
        print(f"  {band}: cohort-strip <c>={np.nanmean(cs):+.3f}  "
              f"max|c|={np.max(np.abs(cs)):.3f}  "
              f"grid {int((c>0).sum())}/{c.size}")

    fig = plt.figure(figsize=(16.5, 11.5))
    outer = fig.add_gridspec(
        2, 3, wspace=0.20, hspace=0.62,
        left=0.06, right=0.97, top=0.86, bottom=0.13,
    )

    last_im_bot = None
    last_im_top = None
    for i, band in enumerate(BAND_ORDER):
        r, c = i // 3, i % 3
        cohort_row = cohort[cohort.band == band].iloc[0]
        im_top, im_bot = plot_band_composite(
            fig, outer[r, c], band, payloads[band], cohort_row,
            show_y=(c == 0), show_x=(r == 1),
        )
        last_im_bot = im_bot
        last_im_top = im_top

    cb_strip_ax = fig.add_axes([0.08, 0.063, 0.36, 0.020])
    cb_strip = fig.colorbar(last_im_top, cax=cb_strip_ax,
                            orientation="horizontal", extend="both")
    cb_strip.set_label(
        rf"cohort strip $c$ (tight scale, $\pm{COHORT_VABS:.2f}$)",
        fontsize=9.0,
    )
    cb_strip.ax.tick_params(labelsize=8)

    cb_grid_ax = fig.add_axes([0.56, 0.063, 0.36, 0.020])
    cb_grid = fig.colorbar(last_im_bot, cax=cb_grid_ax,
                           orientation="horizontal", extend="both")
    cb_grid.set_label(
        rf"per-patient grid $c$ (normal scale, $\pm{PATIENT_VABS:.2f}$)",
        fontsize=9.0,
    )
    cb_grid.ax.tick_params(labelsize=8)

    fig.text(0.5, 0.955,
             "Cohort coherence composite at the per-pair cophenetic "
             "communication-distance layer  "
             r"($n = 10$, ImCoh|·|, $R = 200$ matched-strength)",
             ha="center", va="bottom", fontsize=11)
    fig.text(0.5, 0.935,
             rf"top strip = cohort consensus across $Q = {N_COHORT_QUANTILES}$ "
             rf"quantiles (tight $\pm {COHORT_VABS:.2f}$); "
             rf"bottom grid = patient $\times$ decile (normal $\pm {PATIENT_VABS:.2f}$); "
             r"green = trace direction, red = anti direction",
             ha="center", va="bottom", fontsize=9, color="0.30")

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands" / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_cohort_coherence_map.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
