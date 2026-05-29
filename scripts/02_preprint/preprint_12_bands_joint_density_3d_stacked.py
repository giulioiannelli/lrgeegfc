#!/usr/bin/env python3
"""3D stacked joint-density figure — raw ρ vs matched-strength-corrected.

Two stacks side by side, one big 3D axis per column. Six surfaces
stacked vertically per axis, ordered top → bottom by signed cohort ρ
(β at top, θ at bottom). LEFT column shows the raw centered copula
density; RIGHT column shows the matched-strength-corrected signed
enrichment. The contrast between the two columns IS the visual proof
that the matched-strength surrogate isolates the genuine memory
effect from a same-strength rewire's artifact:

    γ_l carries a saddle in the left column (its raw ρ ≈ +0.083 is
    comparable to α's +0.105) but goes dead flat in the right column
    because its surrogate p95 (+0.109) sits above its observed ρ.

LEFT column — raw centered copula density:
    h_raw(u, v; b) = Z_CHAR · tanh((c(u, v; rho_obs) − 1) / Z_CHAR)
RIGHT column — matched-strength-corrected signed enrichment:
    h_corr(u, v; b) = Z_CHAR · tanh(signed_excess / Z_CHAR)
    signed_excess = sign(c_obs − 1)
                    · max(|c_obs − 1| − |c_surr_p95 − 1|, 0)

No overlays — no 2σ ellipse, no floor heatmap, no diagonals.

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv

Output (PDF only)
-----------------
data/preprint/figures/all_bands/per_pair_trace/fig_bands_joint_density_3d_stacked.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import norm as sp_norm

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

matplotlib.rcParams.update({
    "font.size": 12.5,
    "axes.labelsize": 12.0,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
})


BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"

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

GRID_N = 110
RANK_LO = 0.02
RANK_HI = 0.98
Z_CHAR = 0.15

DZ_LAYER = 1.0
HEIGHT_FRAC = 0.42
HEIGHT_SCALE = HEIGHT_FRAC * DZ_LAYER / Z_CHAR


def _copula_centered(rho: float, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    nu = sp_norm.ppf(U)
    nv = sp_norm.ppf(V)
    rho2 = rho * rho
    norm_factor = 1.0 / np.sqrt(1.0 - rho2 + 1e-18)
    expo = -0.5 / (1.0 - rho2 + 1e-18) * (
        rho2 * (nu**2 + nv**2) - 2.0 * rho * nu * nv
    )
    return norm_factor * np.exp(expo) - 1.0


def raw_height(rho_obs: float, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    dev_obs = _copula_centered(rho_obs, U, V)
    return Z_CHAR * np.tanh(dev_obs / Z_CHAR)


def corrected_height(rho_obs: float, rho_surr_p95: float,
                     U: np.ndarray, V: np.ndarray) -> np.ndarray:
    dev_obs = _copula_centered(rho_obs, U, V)
    dev_surr = _copula_centered(rho_surr_p95, U, V)
    abs_excess = np.maximum(np.abs(dev_obs) - np.abs(dev_surr), 0.0)
    signed_excess = np.sign(dev_obs) * abs_excess
    return Z_CHAR * np.tanh(signed_excess / Z_CHAR)


def _style_axes(ax, n_bands: int) -> None:
    ax.set_xlim(RANK_LO, RANK_HI)
    ax.set_ylim(RANK_LO, RANK_HI)
    z_lo = -HEIGHT_FRAC * DZ_LAYER * 1.0
    z_hi = (n_bands - 1) * DZ_LAYER + HEIGHT_FRAC * DZ_LAYER * 1.0
    ax.set_zlim(z_lo, z_hi)
    ax.view_init(elev=14, azim=-62)
    ax.set_box_aspect((1.0, 1.0, 2.3))

    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", "½", "1"])
    ax.set_yticklabels(["0", "½", "1"])
    ax.set_zticks([])
    ax.set_xlabel(r"$u$ rank of $\Delta_{\mathrm{task}}$", color="0.20", labelpad=-4)
    ax.set_ylabel(r"$v$ rank of $\Delta_{\mathrm{rest}}$", color="0.20", labelpad=-4)
    ax.set_zlabel("")
    ax.tick_params(axis="x", pad=-3)
    ax.tick_params(axis="y", pad=-3)

    pane_color = (0.99, 0.99, 0.99, 0.4)
    ax.xaxis.set_pane_color(pane_color)
    ax.yaxis.set_pane_color(pane_color)
    ax.zaxis.set_pane_color(pane_color)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = (0.88, 0.88, 0.88, 0.45)
        axis._axinfo["grid"]["linewidth"] = 0.3


def _plot_stack(ax, U, V, bands_sorted, h_field_per_band, norm_z) -> None:
    n_bands = len(bands_sorted)
    for i, band in enumerate(bands_sorted):
        z_offset = (n_bands - 1 - i) * DZ_LAYER
        h_field = h_field_per_band[band]
        Z_surface = HEIGHT_SCALE * h_field + z_offset
        facecolors = TRACE_CMAP(norm_z(h_field))
        ax.plot_surface(
            U, V, Z_surface,
            facecolors=facecolors,
            rcount=GRID_N, ccount=GRID_N,
            antialiased=True,
            shade=True,
            linewidth=0,
        )


def main() -> Path:
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    rhos: dict[str, float] = {}
    rho_surrs: dict[str, float] = {}
    p_vals: dict[str, float] = {}
    for band in BAND_ORDER:
        row = cohort[cohort.band == band].iloc[0]
        rhos[band] = float(row["obs_median_rho"])
        rho_surrs[band] = float(row["surr_median_rho_p95"])
        p_vals[band] = float(row["paired_wilcoxon_p"])

    bands_sorted = sorted(BAND_ORDER, key=lambda b: -rhos[b])
    print("  stack order (top → bottom, descending ρ_obs):")
    for b in bands_sorted:
        print(f"    {b}: ρ_obs = {rhos[b]:+.3f}, "
              f"ρ_surr_p95 = {rho_surrs[b]:+.3f}, p = {p_vals[b]:.3f}")

    n_bands = len(bands_sorted)
    u = np.linspace(RANK_LO, RANK_HI, GRID_N)
    U, V = np.meshgrid(u, u)
    norm_z = Normalize(vmin=-Z_CHAR, vmax=Z_CHAR)

    raw_fields = {b: raw_height(rhos[b], U, V) for b in bands_sorted}
    corr_fields = {b: corrected_height(rhos[b], rho_surrs[b], U, V) for b in bands_sorted}

    # Tall vertical figure with two 3D axes side by side in a 1x3
    # gridspec (left ax / label gutter / right ax). The figure is
    # cropped to the 3D content via savefig bbox_inches='tight', which
    # removes the empty matplotlib padding margins from the saved PDF.
    fig = plt.figure(figsize=(11.5, 10.5))

    # ax bounding boxes deliberately extend ABOVE the figure top and
    # BELOW the figure bottom so matplotlib's invisible 3D padding ends
    # up off-screen instead of eating visible figure space. The visible
    # 3D box ends up filling nearly all of the figure height.
    ratios = (1.0, 0.28, 1.0)
    total = sum(ratios)
    fig_left, fig_right = 0.005, 0.995
    span = fig_right - fig_left
    left_ax_x = fig_left
    left_ax_w = (ratios[0] / total) * span
    right_ax_x = fig_left + ((ratios[0] + ratios[1]) / total) * span
    right_ax_w = (ratios[2] / total) * span

    ax_y = -0.25
    ax_h = 1.50

    ax_raw = fig.add_axes([left_ax_x, ax_y, left_ax_w, ax_h], projection="3d")
    ax_corr = fig.add_axes([right_ax_x, ax_y, right_ax_w, ax_h], projection="3d")

    for ax in (ax_raw, ax_corr):
        _style_axes(ax, n_bands)

    _plot_stack(ax_raw, U, V, bands_sorted, raw_fields, norm_z)
    _plot_stack(ax_corr, U, V, bands_sorted, corr_fields, norm_z)

    # Band labels — in the middle gutter between the two stacks.
    # Y-positions hand-tuned to match where each layer projects.
    mid_left = fig_left + (ratios[0] / total) * span
    mid_right = fig_left + ((ratios[0] + ratios[1]) / total) * span
    mid_x = 0.5 * (mid_left + mid_right)

    band_y_positions = np.linspace(0.85, 0.13, n_bands)

    for i, band in enumerate(bands_sorted):
        y_band = band_y_positions[i]
        is_sig = p_vals[band] < 0.05
        label_color = "#0a3a1d" if is_sig else "0.35"
        band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
        fig.text(mid_x, y_band + 0.022, rf"{band_tex}",
                 color=label_color, fontsize=26, fontweight="bold",
                 ha="center", va="center")
        fig.text(mid_x, y_band - 0.013,
                 rf"$\rho_{{\mathrm{{obs}}}}={rhos[band]:+.3f}$",
                 color=label_color, fontsize=11.5,
                 ha="center", va="center")
        fig.text(mid_x, y_band - 0.031,
                 rf"$\rho_{{\mathrm{{surr}}}}^{{p95}}={rho_surrs[band]:+.3f}$,  "
                 rf"$p={p_vals[band]:.3f}$",
                 color=label_color, fontsize=9.5,
                 ha="center", va="center")

    # Column titles — centered horizontally over each ax's slot.
    left_ax_cx = left_ax_x + left_ax_w / 2
    right_ax_cx = right_ax_x + right_ax_w / 2

    fig.text(left_ax_cx, 0.985,
             r"raw  $c(u,v;\rho_{\mathrm{obs}})-1$",
             ha="center", va="top", fontsize=14.5,
             color="0.20", fontweight="bold")
    fig.text(left_ax_cx, 0.962,
             r"(magnitude of the rank dependence, no null subtraction)",
             ha="center", va="top", fontsize=10.0, color="0.40",
             fontstyle="italic")

    fig.text(right_ax_cx, 0.985,
             r"matched-strength noise floor subtracted",
             ha="center", va="top", fontsize=14.5,
             color="0.20", fontweight="bold")
    fig.text(right_ax_cx, 0.962,
             r"sign$(c_{\mathrm{obs}}-1) \cdot \max(|c_{\mathrm{obs}}-1| - "
             r"|c_{\mathrm{surr}}^{p95}-1|,\,0)$",
             ha="center", va="top", fontsize=10.0, color="0.40",
             fontstyle="italic")

    # Shared colorbar at the bottom — pulled up close to the bottom-most
    # surface since matplotlib's 3D padding pushes its bbox off-screen.
    cb_ax = fig.add_axes([0.32, 0.025, 0.36, 0.010])
    sm = matplotlib.cm.ScalarMappable(
        norm=Normalize(vmin=-Z_CHAR, vmax=Z_CHAR),
        cmap=TRACE_CMAP,
    )
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cb_ax, orientation="horizontal", extend="both")
    cb.set_label(
        r"tanh-compressed height field "
        r"(green = trace ridge, red = anti)",
        labelpad=4,
    )
    cb.set_ticks([-Z_CHAR, 0.0, Z_CHAR])
    cb.set_ticklabels([r"$-Z_c$", "0", r"$+Z_c$"])

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands" / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_joint_density_3d_stacked.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
