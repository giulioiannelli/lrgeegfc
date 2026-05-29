#!/usr/bin/env python3
"""Single-band cohort joint-rank density (standalone β panel).

Standalone version of one panel from
``preprint_09_bands_joint_density``.  Produces the β (or other band)
joint-rank density figure that the user accepts as Fig 1 candidate B:

* Main panel: signed enrichment-over-noise field
  ``signed_excess = sign(ρ_obs - 0) ·
       max(|copula(ρ_obs) - 1| - |copula(ρ_surr_p95) - 1|, 0)``
  rendered as a green-on-cream diagonal (trace) / red-on-cream
  anti-diagonal (anti) heatmap on the unit-square rank grid, with
  68% / 95% HDR contours of the joint density overlaid.
* Top marginal: cohort distribution of ``u = rank Δ_task`` split by
  ``v = rank Δ_rest`` half (green = top, red = bottom).
* Right marginal: symmetric — distribution of v split by u-half.

Trace bands → conditional histograms visibly separate; null bands →
they overlap.

Usage
-----
    python preprint_17_band_joint_density_single.py             # beta
    python preprint_17_band_joint_density_single.py --band alpha
    python preprint_17_band_joint_density_single.py --band high_gamma

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
data/reports/imcoh_continuous_trace/per_pair_split/<Pat>_<band>.npz

Output (PDF only)
-----------------
data/preprint/figures/<band>/per_pair_trace/fig_<band>_joint_density.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde, norm, rankdata

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

TRACE_CMAP = LinearSegmentedColormap.from_list(
    "trace_rwg",
    [
        "#4a0e0e", "#a02525", "#d8584c", "#f0a89e",
        "#fdf2ef",
        "#bfddc7", "#5fac7e", "#1a7c3e", "#0a3a1d",
    ],
    N=512,
)

GRID_N = 120
SHARED_VABS = 0.025
KDE_BW = 0.22

SOURCES = ("copula", "empirical")
SOURCE_LABEL = {
    "copula":    r"Gaussian copula at $\rho_{\mathrm{obs}}$",
    "empirical": r"cohort-mean empirical KDE",
}


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return dict(
        dD_task=np.asarray(d["dD_task"]),
        dD_rest=np.asarray(d["dD_rest"]),
    )


def pool_ranks(band: str) -> tuple[np.ndarray, np.ndarray]:
    t_all, r_all = [], []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        t = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        r = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        t_all.append(t); r_all.append(r)
    return np.concatenate(t_all), np.concatenate(r_all)


def cohort_average_empirical_density(band: str, grid_n: int,
                                      bw: float) -> np.ndarray:
    """Per-patient KDE on the within-patient rank pairs, then cohort
    mean.  Returns a ``(grid_n, grid_n)`` joint-density array on the
    same ``[0, 1]²`` mesh as the analytic Gaussian copula, so the two
    can be plugged into the same ``signed_excess`` recipe.
    """
    xy = np.linspace(0.0, 1.0, grid_n)
    xx, yy = np.meshgrid(xy, xy)
    grid_pts = np.vstack([xx.ravel(), yy.ravel()])
    per_pat = []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        t = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        r = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        k = gaussian_kde(np.vstack([t, r]), bw_method=bw)
        per_pat.append(k(grid_pts).reshape(grid_n, grid_n))
    return np.mean(per_pat, axis=0)


def gaussian_copula_density(rho: float, grid_n: int) -> np.ndarray:
    """Analytic bivariate Gaussian copula density on ``[0, 1]²``."""
    eps = 1.5e-3
    u = np.linspace(eps, 1.0 - eps, grid_n)
    uu, vv = np.meshgrid(u, u)
    nu = norm.ppf(uu); nv = norm.ppf(vv)
    rho2 = rho * rho
    norm_factor = 1.0 / np.sqrt(1.0 - rho2)
    expo = -0.5 / (1.0 - rho2) * (rho2 * (nu**2 + nv**2)
                                    - 2.0 * rho * nu * nv)
    return norm_factor * np.exp(expo)


def render(band: str, source: str = "copula") -> Path:
    if source not in SOURCES:
        raise ValueError(f"source must be in {SOURCES}, got {source!r}")
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    cohort_row = cohort[cohort.band == band].iloc[0]

    t, r = pool_ranks(band)
    n_pairs = t.size
    print(f"  {band}: pooled {n_pairs} pairs across {len(COHORT)} patients "
          f"— source = {source}")

    rho_cohort = float(cohort_row["obs_median_rho"])
    rho_surr_p95 = float(cohort_row["surr_median_rho_p95"])
    p_val = float(cohort_row["paired_wilcoxon_p"])
    is_sig = p_val < 0.05

    # Density-source switch.
    # - copula:    analytic Gaussian copula at ρ_obs (model).
    # - empirical: per-patient KDE on rank pairs, cohort-mean (data).
    # Both run through the SAME ``signed_excess`` recipe against the
    # analytic matched-strength-p95 noise floor, so the visual scale
    # (cmap range, contour HDR levels) is interchangeable between modes.
    if source == "copula":
        z = gaussian_copula_density(rho_cohort, GRID_N)
    else:
        z = cohort_average_empirical_density(band, GRID_N, KDE_BW)
    z_surr = gaussian_copula_density(rho_surr_p95, GRID_N)
    dev_obs = z - 1.0
    dev_surr = z_surr - 1.0
    abs_excess = np.maximum(np.abs(dev_obs) - np.abs(dev_surr), 0.0)
    signed_excess = np.sign(dev_obs) * abs_excess
    xy = np.linspace(0.0, 1.0, GRID_N)
    xx, yy = np.meshgrid(xy, xy)

    # ---------- Figure scaffold (single panel + 2 marginals) ----------
    fig = plt.figure(figsize=(7.6, 6.9))
    sub = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 3.6], width_ratios=[3.6, 1.0],
        left=0.11, right=0.86, top=0.97, bottom=0.10,
        hspace=0.04, wspace=0.04,
    )
    ax_main = fig.add_subplot(sub[1, 0])
    ax_top = fig.add_subplot(sub[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(sub[1, 1], sharey=ax_main)
    # Band label lives in the otherwise-empty top-right corner cell
    # — between the top and right marginal axes, above the main panel
    # — so it adds no extra figure height.
    ax_label = fig.add_subplot(sub[0, 1])
    ax_label.set_xticks([]); ax_label.set_yticks([])
    for s in ax_label.spines.values():
        s.set_visible(False)
    ax_label.text(0.5, 0.5, band_tex,
                  transform=ax_label.transAxes,
                  ha="center", va="center",
                  fontsize=24, fontweight="bold", color="0.18")

    for a in (ax_top, ax_right, ax_main):
        a.set_autoscalex_on(False); a.set_autoscaley_on(False)
    ax_main.set_box_aspect(1.0)
    ax_top.set_box_aspect(1.0 / 3.6)
    ax_right.set_box_aspect(3.6)
    # No ax_main.set_facecolor — figures default to transparent
    # background ([[feedback-default-transparent-figures]]).  The
    # data limits are locked at exactly (0, 1) so there is no
    # padding between the imshow extent and the panel edge.

    im = ax_main.imshow(
        signed_excess, origin="lower", extent=[0.0, 1.0, 0.0, 1.0],
        cmap=TRACE_CMAP, vmin=-SHARED_VABS, vmax=SHARED_VABS,
        aspect="auto", interpolation="bilinear",
    )

    ax_main.plot([0, 1], [0, 1], color="white", lw=2.6,
                 alpha=0.85, zorder=4)
    ax_main.plot([0, 1], [1, 0], color="white", lw=1.0, ls="--",
                 alpha=0.40, zorder=4)

    # HDR contours at 68% / 95% mass coverage of the same density
    # surface the heatmap is built on (analytic copula OR empirical KDE,
    # depending on --source).  Dark halo + bright white core so the
    # contours read on any patch of the green/red/cream cmap.
    dA = 1.0 / (GRID_N * GRID_N)
    def _hdr_level(zf: np.ndarray, coverage: float) -> float:
        zs = np.sort(zf.ravel())[::-1]
        cm = np.cumsum(zs) * dA
        idx = min(int(np.searchsorted(cm, coverage)), zs.size - 1)
        return float(zs[idx])
    lvl_68 = _hdr_level(z, 0.68)
    lvl_95 = _hdr_level(z, 0.95)
    ax_main.contour(xx, yy, z, levels=[lvl_95, lvl_68],
                    colors=["0.20", "0.20"],
                    linewidths=[3.6, 2.4],
                    alpha=[0.45, 0.55], zorder=4.5)
    ax_main.contour(xx, yy, z, levels=[lvl_95, lvl_68],
                    colors=["white", "white"],
                    linewidths=[1.4, 0.9],
                    alpha=[0.85, 0.85], zorder=5)

    ax_main.plot(0.5, 0.5, "o", markersize=5, color="white",
                 markeredgecolor="0.10", markeredgewidth=1.0, zorder=6)

    ax_main.set_xticks([0.0, 0.5, 1.0])
    ax_main.set_yticks([0.0, 0.5, 1.0])
    ax_main.set_xticklabels(["0", "½", "1"])
    ax_main.set_yticklabels(["0", "½", "1"])
    ax_main.set_xlabel(r"$u$ = within-patient rank of "
                       r"$\Delta_{\mathrm{task}}(i,j)$",
                       color="0.20", fontsize=10.5)
    ax_main.set_ylabel(r"$v$ = within-patient rank of "
                       r"$\Delta_{\mathrm{rest}}(i,j)$",
                       color="0.20", fontsize=10.5)

    # Marginals: conditional split histograms
    n_bins = 28
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    mask_high_r = r > 0.5
    hist_top_hi, _ = np.histogram(t[mask_high_r], bins=bin_edges, density=True)
    hist_top_lo, _ = np.histogram(t[~mask_high_r], bins=bin_edges, density=True)
    yn_top = max(hist_top_hi.max(), hist_top_lo.max(), 1e-9)
    hist_top_hi_n = 0.85 * hist_top_hi / yn_top
    hist_top_lo_n = 0.85 * hist_top_lo / yn_top
    mask_high_t = t > 0.5
    hist_rt_hi, _ = np.histogram(r[mask_high_t], bins=bin_edges, density=True)
    hist_rt_lo, _ = np.histogram(r[~mask_high_t], bins=bin_edges, density=True)
    yn_rt = max(hist_rt_hi.max(), hist_rt_lo.max(), 1e-9)
    hist_rt_hi_n = 0.85 * hist_rt_hi / yn_rt
    hist_rt_lo_n = 0.85 * hist_rt_lo / yn_rt

    ax_top.stairs(hist_top_lo_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#c0392b", alpha=0.40,
                  edgecolor="#7f1d1d", linewidth=0.9, zorder=2)
    ax_top.stairs(hist_top_hi_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#1a7c3e", alpha=0.45,
                  edgecolor="#0a3a1d", linewidth=0.9, zorder=2)
    ax_top.axvline(0.5, color="0.35", lw=0.7, alpha=0.6, zorder=1)
    ax_top.set_ylim(0, 1)
    ax_top.tick_params(axis="both", which="both",
                       bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
    ax_top.spines[["top", "right", "left"]].set_visible(False)
    ax_top.spines["bottom"].set_color("0.70")
    ax_top.text(0.03, 0.10, r"$P(u)$ split by $v$ half",
                transform=ax_top.transAxes, ha="left", va="bottom",
                color="0.20", fontstyle="italic",
                fontsize=9.5,
                bbox=dict(facecolor="white", alpha=0.65,
                          edgecolor="none", boxstyle="round,pad=0.25"),
                zorder=10)

    ax_right.stairs(hist_rt_lo_n, bin_edges, baseline=0.0, fill=True,
                    orientation="horizontal",
                    facecolor="#c0392b", alpha=0.40,
                    edgecolor="#7f1d1d", linewidth=0.9, zorder=2)
    ax_right.stairs(hist_rt_hi_n, bin_edges, baseline=0.0, fill=True,
                    orientation="horizontal",
                    facecolor="#1a7c3e", alpha=0.45,
                    edgecolor="#0a3a1d", linewidth=0.9, zorder=2)
    ax_right.axhline(0.5, color="0.35", lw=0.7, alpha=0.6, zorder=1)
    ax_right.set_xlim(0, 1)
    ax_right.tick_params(axis="both", which="both",
                         bottom=False, top=False, left=False, right=False,
                         labelbottom=False, labeltop=False,
                         labelleft=False, labelright=False)
    ax_right.spines[["top", "right", "bottom"]].set_visible(False)
    ax_right.spines["left"].set_color("0.70")
    ax_right.text(0.22, 0.97, r"$P(v)$ split by $u$ half",
                  transform=ax_right.transAxes,
                  ha="left", va="top",
                  rotation=270, rotation_mode="anchor",
                  color="0.20", fontstyle="italic", fontsize=9.5,
                  bbox=dict(facecolor="white", alpha=0.65,
                            edgecolor="none", boxstyle="round,pad=0.15"),
                  clip_on=False, zorder=10)

    # ρ_obs / ρ_surr / source are reported to stdout only; no text on
    # the figure — band label is the only annotation, lives in the
    # otherwise-empty top-right corner cell.
    print(f"    rho_obs={rho_cohort:+.3f}  "
          f"rho_surr_p95={rho_surr_p95:+.3f}  "
          f"wilcoxon_p={p_val:.4f}  source={source}")

    for spine in ax_main.spines.values():
        spine.set_edgecolor("0.55")
        spine.set_linewidth(0.9)

    # Colorbar — sized to match the main panel's vertical extent
    # (main row spans y ∈ [0.10, 0.781] inside the gridspec).
    cbar_ax = fig.add_axes([0.89, 0.10, 0.035, 0.681])
    cb = fig.colorbar(im, cax=cbar_ax, orientation="vertical", extend="both")
    # Δc⋆ = sgn(c − 1) · [|c − 1| − |c_surr^p95 − 1|]_+
    # where c is the joint copula (or empirical KDE) density on the
    # rank grid and c_surr^p95 is the matched-strength noise floor.
    # Full expansion lives in the manuscript caption / companion .md.
    cb.set_label(
        r"$\Delta c^{\!\star}$",
        rotation=270, labelpad=22, fontsize=14,
    )
    cb.ax.tick_params(labelsize=10.5)

    ax_main.set_xlim(0.0, 1.0, auto=False)
    ax_main.set_ylim(0.0, 1.0, auto=False)

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_{band}_joint_density_{source}.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def main() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--source", default="copula",
                        choices=list(SOURCES),
                        help="Density surface that fills the heatmap. "
                             "'copula' (default): analytic Gaussian "
                             "copula at ρ_obs (the model). "
                             "'empirical': per-patient KDE on rank "
                             "pairs, cohort-mean (the data).")
    args = parser.parse_args()
    return render(args.band, args.source)


if __name__ == "__main__":
    main()
