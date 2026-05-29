#!/usr/bin/env python3
"""Six-band cohort joint density at the per-pair cophenetic layer.

For each band, plot the smooth 2D joint density of the cohort-pooled
within-patient rank pair (rank Δ_task(i,j), rank Δ_rest(i,j)) as a
heatmap of the density residual (KDE − uniform), with an overlaid
principal-axis confidence ellipse.

A trace band reads as a smooth diagonal ridge of excess density along
y = x and a *tilted, elongated* ellipse. A null band reads as a near-
uniform colorfield and a *near-circular, untilted* ellipse. The
geometric primitive (ellipse axes ratio + tilt angle) is what makes
β stand out against the rest instantly.

The marginal sub-panels show the CONDITIONAL split — distribution of
rank Δ_rest for pairs in the top vs bottom half of rank Δ_task (and
vice versa). Trace bands separate these conditional distributions
visibly; null bands let them overlap.

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
from matplotlib.patches import Ellipse
from scipy.stats import rankdata, gaussian_kde

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

# Per-script font bump. All non-band-label text inherits these — no
# `fontsize=` hardcoded line-by-line below. Band labels stay at the
# fixed 22 pt because they're the dominant typographic element.
matplotlib.rcParams.update({
    "font.size": 12.5,
    "axes.labelsize": 12.5,
    "axes.titlesize": 12.5,
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "legend.fontsize": 11.5,
})


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
        "#4a0e0e",
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

GRID_N = 120
BW = 0.22
SHARED_VABS = 0.025

# Toggle the joint-density envelope drawn on top of the heatmap.
# - "hdr"     : highest-density-region contours of z at 68% / 95% mass coverage
#               (closed curves; always inside [0,1]² by construction).
# - "ellipse" : analytic 2σ ellipse of the joint-rank covariance
#               (extends past [0,1] when |ρ| > 0; geometric artefact of the
#               Uniform marginal scale).
# - "both"    : draw both.
# - "none"    : draw neither.
ENVELOPE_STYLE = "hdr"


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
        t_rank = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        r_rank = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        t_all.append(t_rank)
        r_all.append(r_rank)
    return np.concatenate(t_all), np.concatenate(r_all)


def per_patient_rhos(band: str) -> tuple[np.ndarray, list[str]]:
    """Per-patient Spearman rho_p of (Δ_task, Δ_rest) for this band."""
    from scipy.stats import spearmanr
    rhos, labels = [], []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        if len(pair["dD_task"]) < 2:
            continue
        r, _ = spearmanr(pair["dD_task"], pair["dD_rest"])
        rhos.append(float(r))
        labels.append(pat.replace("Pat_", "P"))
    return np.asarray(rhos), labels


def cohort_average_density(band: str, grid_n: int, bw: float) -> np.ndarray:
    """Per-patient KDE, then cohort-mean. Returns density / uniform (centered at 1)."""
    xy = np.linspace(0.0, 1.0, grid_n)
    xx, yy = np.meshgrid(xy, xy)
    grid_points = np.vstack([xx.ravel(), yy.ravel()])
    per_pat = []
    for pat in COHORT:
        pair = load_pair_data(pat, band)
        n = len(pair["dD_task"])
        if n < 2:
            continue
        t = (rankdata(pair["dD_task"], method="ordinal") - 0.5) / n
        r = (rankdata(pair["dD_rest"], method="ordinal") - 0.5) / n
        k = gaussian_kde(np.vstack([t, r]), bw_method=bw)
        per_pat.append(k(grid_points).reshape(grid_n, grid_n))
    return np.mean(per_pat, axis=0)


def gaussian_copula_density(rho: float, grid_n: int) -> np.ndarray:
    """Analytic bivariate Gaussian copula density on [0, 1]² grid.

    For rank pairs marginally uniform on [0, 1] and joint Gaussian copula
    with correlation rho, returns the joint density on a grid_n × grid_n
    mesh. Density = 1 under independence (rho = 0); > 1 along the y = x
    diagonal under positive rho; > 1 along y = 1 - x under negative rho.
    """
    from scipy.stats import norm
    eps = 1.5e-3
    u = np.linspace(eps, 1.0 - eps, grid_n)
    uu, vv = np.meshgrid(u, u)
    nu = norm.ppf(uu)
    nv = norm.ppf(vv)
    rho2 = rho * rho
    norm_factor = 1.0 / np.sqrt(1.0 - rho2)
    expo = -0.5 / (1.0 - rho2) * (rho2 * (nu**2 + nv**2) - 2.0 * rho * nu * nv)
    return norm_factor * np.exp(expo)


def plot_band_density(fig, outer_spec, band: str, t: np.ndarray, r: np.ndarray,
                      cohort_row: pd.Series) -> "matplotlib.image.AxesImage":
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    sub = outer_spec.subgridspec(
        2, 2, height_ratios=[1.0, 3.6], width_ratios=[3.6, 1.0],
        hspace=0.04, wspace=0.04,
    )
    ax_main = fig.add_subplot(sub[1, 0])
    ax_top = fig.add_subplot(sub[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(sub[1, 1], sharey=ax_main)
    # Disable autoscale on the SHARED direction of the marginal axes
    # immediately. Otherwise ax_top.stairs() (data in x ∈ [0, 1]) and
    # ax_right.stairs() (data in y ∈ [0, 1]) will trigger autoscale_view
    # on the shared X / Y axis at render time, which silently overrides
    # ax_main.set_xlim(-0.12, 1.12) and clips the 2σ ellipse.
    ax_top.set_autoscalex_on(False)
    ax_right.set_autoscaley_on(False)
    ax_main.set_autoscalex_on(False)
    ax_main.set_autoscaley_on(False)
    # Lock box aspects so the three axes' boxes track the subgrid ratios:
    # main is square, top is (1 : 3.6) tall:wide (same width as main),
    # right is (3.6 : 1) tall:wide (same height as main). This replaces
    # the older aspect="equal" on imshow, which shrank the main box but
    # left the marginal boxes at the gridspec slot width (visible mismatch).
    ax_main.set_box_aspect(1.0)
    ax_top.set_box_aspect(1.0 / 3.6)
    ax_right.set_box_aspect(3.6)
    # The main panel's data limits extend past the unit square so the 2σ
    # ellipse is fully shown, but the Gaussian-copula density is only
    # defined on [0, 1]². Paint the axes background with the colormap's
    # centre colour (#fdf2ef, the "z = 0" cream) so the padding outside
    # the imshow merges seamlessly with the heatmap's neutral regions —
    # no white frame between the unit-square density and the panel edge.
    ax_main.set_facecolor("#fdf2ef")

    rho_cohort = float(cohort_row["obs_median_rho"])
    rho_surr_p95 = float(cohort_row["surr_median_rho_p95"])
    z = gaussian_copula_density(rho_cohort, GRID_N)
    z_surr = gaussian_copula_density(rho_surr_p95, GRID_N)
    z_centered = z - 1.0
    # Enrichment-over-noise field: |obs - independence| minus
    # |surrogate-p95 - independence|, clipped at zero, signed by the
    # observed deviation direction. Above noise floor → green diagonal +
    # red anti-diagonal (trace pattern). At-or-below noise floor →
    # everywhere zero → cream. This is what visually separates α (above
    # noise) from γ_l (below noise) despite their similar cohort ρ.
    dev_obs = z_centered
    dev_surr = z_surr - 1.0
    abs_excess = np.maximum(np.abs(dev_obs) - np.abs(dev_surr), 0.0)
    signed_excess = np.sign(dev_obs) * abs_excess
    xy = np.linspace(0.0, 1.0, GRID_N)
    xx, yy = np.meshgrid(xy, xy)

    im = ax_main.imshow(
        signed_excess, origin="lower", extent=[0.0, 1.0, 0.0, 1.0],
        cmap=TRACE_CMAP, vmin=-SHARED_VABS, vmax=SHARED_VABS,
        aspect="auto", interpolation="bilinear",
    )

    ax_main.plot([0, 1], [0, 1], color="white", lw=2.6, alpha=0.85, zorder=4)
    ax_main.plot([0, 1], [1, 0], color="white", lw=1.0, ls="--",
                 alpha=0.40, zorder=4)

    if ENVELOPE_STYLE in ("hdr", "both"):
        # Highest-density-region contours at fixed quantile coverage of the
        # joint density z (not z_centered). Always closed inside [0,1]².
        dA = 1.0 / (GRID_N * GRID_N)
        def _hdr_level(z_field: np.ndarray, coverage: float) -> float:
            z_sorted = np.sort(z_field.ravel())[::-1]
            cum_mass = np.cumsum(z_sorted) * dA
            idx = int(np.searchsorted(cum_mass, coverage))
            idx = min(idx, z_sorted.size - 1)
            return float(z_sorted[idx])
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

    # Analytic 2σ ellipse geometry from the joint-rank covariance is computed
    # unconditionally so axis_ratio + tilt can annotate the panel even when
    # the ellipse itself is not drawn. Uniform marginal variance is 1/12, so
    # the 2σ tips reach 0.5 ± 0.577 — extends past [0,1] for any |ρ|.
    var_uniform = 1.0 / 12.0
    cov_off = rho_cohort * var_uniform
    cohort_paired_cov = np.array([[var_uniform, cov_off],
                                  [cov_off, var_uniform]])
    eigvals, eigvecs = np.linalg.eigh(cohort_paired_cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    angle = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    width = 2 * 2 * float(np.sqrt(eigvals[0]))
    height = 2 * 2 * float(np.sqrt(eigvals[1]))
    axis_ratio = width / max(height, 1e-12)

    if ENVELOPE_STYLE in ("ellipse", "both"):
        ell = Ellipse(
            xy=(0.5, 0.5), width=width, height=height, angle=angle,
            fill=False, edgecolor="white", lw=2.4, alpha=0.85, zorder=5,
        )
        ax_main.add_patch(ell)
        ell_outer = Ellipse(
            xy=(0.5, 0.5), width=width, height=height, angle=angle,
            fill=False, edgecolor="0.20", lw=4.6, alpha=0.40, zorder=4.5,
        )
        ax_main.add_patch(ell_outer)

    ax_main.plot(0.5, 0.5, "o", markersize=5, color="white",
                 markeredgecolor="0.10", markeredgewidth=1.0, zorder=6)

    ax_main.set_xticks([0.0, 0.5, 1.0])
    ax_main.set_yticks([0.0, 0.5, 1.0])
    ax_main.set_xticklabels(["0", "½", "1"])
    ax_main.set_yticklabels(["0", "½", "1"])
    if band == BAND_ORDER[0]:
        ax_main.set_xlabel(r"$u$ = within-patient rank of "
                           r"$\Delta_{\mathrm{task}}(i,j)$",
                           color="0.20")
        ax_main.set_ylabel(r"$v$ = within-patient rank of "
                           r"$\Delta_{\mathrm{rest}}(i,j)$",
                           color="0.20")
    else:
        ax_main.set_xlabel(r"$u$", color="0.20")
        ax_main.set_ylabel(r"$v$", color="0.20")
    # set_xlim/set_ylim deliberately moved to the END of this function
    # (just before return). Reason: ax_top.stairs() and ax_right.stairs()
    # below add artists to the shared X/Y axes, which triggers matplotlib's
    # autoscale_view and resets the main-panel limits. Calling set_xlim
    # AFTER all marginal artists, with auto=False, locks the limits.

    # Conditional-split histograms — the marginals that ACTUALLY combine
    # with the main joint density.
    #   Top marginal: cohort-pooled distribution of rank Δ_task (= u),
    #     split by whether each pair sits in the TOP half (green) or
    #     BOTTOM half (red) of rank Δ_rest (= v).
    #   Right marginal: symmetric view — distribution of v, split by
    #     top/bottom half of u.
    # Trace bands → the two conditional histograms visibly SEPARATE
    # (top-v pairs cluster at high u, bottom-v pairs at low u). Null
    # bands → they OVERLAP. The visible separation IS the empirical
    # 1-D projection of the joint rank correlation shown in the main
    # panel.
    n_bins = 28
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)

    mask_high_r = r > 0.5
    hist_top_hi, _ = np.histogram(t[mask_high_r], bins=bin_edges, density=True)
    hist_top_lo, _ = np.histogram(t[~mask_high_r], bins=bin_edges, density=True)
    y_norm_top = max(hist_top_hi.max(), hist_top_lo.max(), 1e-9)
    hist_top_hi_n = 0.85 * hist_top_hi / y_norm_top
    hist_top_lo_n = 0.85 * hist_top_lo / y_norm_top

    mask_high_t = t > 0.5
    hist_rt_hi, _ = np.histogram(r[mask_high_t], bins=bin_edges, density=True)
    hist_rt_lo, _ = np.histogram(r[~mask_high_t], bins=bin_edges, density=True)
    y_norm_rt = max(hist_rt_hi.max(), hist_rt_lo.max(), 1e-9)
    hist_rt_hi_n = 0.85 * hist_rt_hi / y_norm_rt
    hist_rt_lo_n = 0.85 * hist_rt_lo / y_norm_rt

    # Top marginal: u histograms split by v
    ax_top.stairs(hist_top_lo_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#c0392b", alpha=0.40,
                  edgecolor="#7f1d1d", linewidth=0.9, zorder=2)
    ax_top.stairs(hist_top_hi_n, bin_edges, baseline=0.0, fill=True,
                  facecolor="#1a7c3e", alpha=0.45,
                  edgecolor="#0a3a1d", linewidth=0.9, zorder=2)
    ax_top.axvline(0.5, color="0.35", lw=0.7, alpha=0.6, zorder=1)
    ax_top.set_ylim(0, 1)
    # Hide tick rendering only (do NOT call set_xticks([]) — that clears
    # the FixedLocator on the shared X axis and wipes main's [0, ½, 1]
    # ticks). tick_params is per-axes and leaves the shared locator alone.
    ax_top.tick_params(axis="both", which="both",
                       bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
    ax_top.spines[["top", "right", "left"]].set_visible(False)
    ax_top.spines["bottom"].set_color("0.70")
    if band == BAND_ORDER[0]:
        ax_top.text(0.03, 0.10,
                    r"$P(u)$ split by $v$ half",
                    transform=ax_top.transAxes, ha="left", va="bottom",
                    color="0.20", fontstyle="italic",
                    bbox=dict(facecolor="white", alpha=0.65,
                              edgecolor="none",
                              boxstyle="round,pad=0.25"),
                    zorder=10)

    # Right marginal: v histograms split by u
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
    if band == BAND_ORDER[0]:
        ax_right.text(0.22, 0.97,
                      r"$P(v)$ split by $u$ half",
                      transform=ax_right.transAxes,
                      ha="left", va="top",
                      rotation=270, rotation_mode="anchor",
                      color="0.20", fontstyle="italic",
                      bbox=dict(facecolor="white", alpha=0.65,
                                edgecolor="none",
                                boxstyle="round,pad=0.15"),
                      clip_on=False, zorder=10)

    rho = float(cohort_row["obs_median_rho"])
    p_val = float(cohort_row["paired_wilcoxon_p"])
    is_sig = p_val < 0.05
    title_color = "#0a3a1d" if is_sig else "0.25"
    ax_top.text(0.5, 1.42, band_tex, transform=ax_top.transAxes,
                ha="center", va="bottom",
                fontsize=22, fontweight="bold", color=title_color)
    ax_top.text(0.5, 1.10,
                rf"$\rho_{{\mathrm{{obs}}}} = {rho:+.3f}$    "
                rf"$\rho_{{\mathrm{{surr}}}}^{{p95}} = {rho_surr_p95:+.3f}$",
                transform=ax_top.transAxes, ha="center", va="bottom",
                color="0.25")

    if is_sig:
        for spine in ax_main.spines.values():
            spine.set_edgecolor("#0a3a1d")
            spine.set_linewidth(2.4)
    else:
        for spine in ax_main.spines.values():
            spine.set_edgecolor("0.55")
            spine.set_linewidth(0.9)

    # Lock the main-panel limits AFTER every marginal artist has been
    # added (otherwise the marginals' autoscale propagates through the
    # sharex / sharey link and resets the limits, clipping the heatmap
    # and ellipse). Limits set EXACTLY to the imshow extent above —
    # zero gap between the heatmap and the panel border.
    ax_main.set_xlim(0.0, 1.0, auto=False)
    ax_main.set_ylim(0.0, 1.0, auto=False)

    return im


def main() -> Path:
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    pooled = {}
    for band in BAND_ORDER:
        pooled[band] = pool_ranks(band)
        print(f"  {band}: pooled {len(pooled[band][0])} pairs")

    fig = plt.figure(figsize=(17.0, 11.0))
    outer = fig.add_gridspec(
        2, 3, wspace=0.20, hspace=0.34,
        left=0.05, right=0.92, top=0.94, bottom=0.07,
    )

    last_im = None
    for i, band in enumerate(BAND_ORDER):
        r, c = i // 3, i % 3
        cohort_row = cohort[cohort.band == band].iloc[0]
        im = plot_band_density(
            fig, outer[r, c], band, pooled[band][0], pooled[band][1], cohort_row,
        )
        last_im = im

    cbar_ax = fig.add_axes([0.94, 0.16, 0.015, 0.62])
    cb = fig.colorbar(last_im, cax=cbar_ax, orientation="vertical",
                      extend="both")
    cb.set_label(
        r"signed enrichment over matched-strength noise floor",
        rotation=270, labelpad=18,
    )

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands" / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_joint_density.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
