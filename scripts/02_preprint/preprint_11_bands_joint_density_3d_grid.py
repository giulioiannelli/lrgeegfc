#!/usr/bin/env python3
"""3D version of the bands joint-density figure — 2x3 panel grid.

Each band gets a clean 3D surface of the **matched-strength noise-
floor-subtracted** signed enrichment, lifted into height. Trace bands
saddle UP at the y = x corners; anti-trace bands at the y = 1 − x
corners; bands below the matched-strength noise floor (e.g. γ_l with
surr p95 > obs ρ) sit dead-flat.

No overlays — no 2σ ellipse, no floor projection, no diagonals. Just
the surface.

Height + color (per panel)
--------------------------
    signed_excess(u, v; b) =
        sign(c_obs − 1)
        × max(|c_obs − 1| − |c_surr_p95 − 1|, 0)

    h(u, v; b) = Z_CHAR · tanh(signed_excess / Z_CHAR)

    Z_CHAR = 0.15 — tanh inflection scale. Trace bands (β) get a
    visibly tall saddle that saturates near (but not at) the corners;
    α has a moderate saddle; γ_l with surr p95 ABOVE obs ρ goes to
    zero everywhere → completely flat.

Significance routing: bands with cohort matched-strength p < 0.05 get
a dark-green title; non-significant bands a gray title.

Inputs
------
data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv

Output (PDF only)
-----------------
data/preprint/figures/all_bands/per_pair_trace/fig_bands_joint_density_3d_grid.pdf
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
    "axes.labelsize": 11.5,
    "axes.titlesize": 13.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
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
Z_CHAR = 0.15  # tanh inflection scale (also color vmin/vmax magnitude)


def _copula_centered(rho: float, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Centered Gaussian copula density c(u, v; rho) − 1 on the (U, V) grid."""
    nu = sp_norm.ppf(U)
    nv = sp_norm.ppf(V)
    rho2 = rho * rho
    norm_factor = 1.0 / np.sqrt(1.0 - rho2 + 1e-18)
    expo = -0.5 / (1.0 - rho2 + 1e-18) * (
        rho2 * (nu**2 + nv**2) - 2.0 * rho * nu * nv
    )
    return norm_factor * np.exp(expo) - 1.0


def signed_enrichment_field(rho_obs: float, rho_surr_p95: float, grid_n: int):
    """Return (U, V, h) where h is tanh-compressed signed enrichment."""
    u = np.linspace(RANK_LO, RANK_HI, grid_n)
    U, V = np.meshgrid(u, u)
    dev_obs = _copula_centered(rho_obs, U, V)
    dev_surr = _copula_centered(rho_surr_p95, U, V)
    abs_excess = np.maximum(np.abs(dev_obs) - np.abs(dev_surr), 0.0)
    signed_excess = np.sign(dev_obs) * abs_excess
    h = Z_CHAR * np.tanh(signed_excess / Z_CHAR)
    return U, V, h


def plot_band_surface(ax, band: str, cohort_row: pd.Series):
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    rho = float(cohort_row["obs_median_rho"])
    rho_surr_p95 = float(cohort_row["surr_median_rho_p95"])
    p_val = float(cohort_row["paired_wilcoxon_p"])
    is_sig = p_val < 0.05

    U, V, h = signed_enrichment_field(rho, rho_surr_p95, GRID_N)

    ax.plot_surface(
        U, V, h,
        cmap=TRACE_CMAP, vmin=-Z_CHAR, vmax=Z_CHAR,
        rcount=GRID_N, ccount=GRID_N,
        antialiased=True,
        shade=True,
        linewidth=0,
    )

    ax.set_xlim(RANK_LO, RANK_HI)
    ax.set_ylim(RANK_LO, RANK_HI)
    ax.set_zlim(-Z_CHAR * 1.05, Z_CHAR * 1.05)

    ax.view_init(elev=28, azim=-58)
    ax.set_box_aspect((1.0, 1.0, 0.55))

    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", "½", "1"])
    ax.set_yticklabels(["0", "½", "1"])
    ax.set_zticks([])

    ax.set_xlabel(r"$u$", labelpad=-7, color="0.20")
    ax.set_ylabel(r"$v$", labelpad=-7, color="0.20")
    ax.set_zlabel("")

    ax.tick_params(axis="x", pad=-3)
    ax.tick_params(axis="y", pad=-3)

    pane_color = (0.99, 0.99, 0.99, 0.45)
    ax.xaxis.set_pane_color(pane_color)
    ax.yaxis.set_pane_color(pane_color)
    ax.zaxis.set_pane_color(pane_color)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = (0.88, 0.88, 0.88, 0.45)
        axis._axinfo["grid"]["linewidth"] = 0.3

    title_color = "#0a3a1d" if is_sig else "0.35"
    ax.set_title(
        rf"{band_tex}    $\rho_{{\mathrm{{obs}}}} = {rho:+.3f}$,  "
        rf"$\rho_{{\mathrm{{surr}}}}^{{p95}} = {rho_surr_p95:+.3f}$,  "
        rf"$p = {p_val:.3f}$",
        color=title_color, pad=4, fontweight="bold",
    )


def main() -> Path:
    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    fig = plt.figure(figsize=(17.5, 11.0))

    for i, band in enumerate(BAND_ORDER):
        ax = fig.add_subplot(2, 3, i + 1, projection="3d")
        cohort_row = cohort[cohort.band == band].iloc[0]
        plot_band_surface(ax, band, cohort_row)
        print(f"  {band}: rho = {cohort_row['obs_median_rho']:+.3f}, "
              f"p = {cohort_row['paired_wilcoxon_p']:.3f}")

    cb_ax = fig.add_axes([0.945, 0.25, 0.012, 0.50])
    sm = matplotlib.cm.ScalarMappable(
        norm=Normalize(vmin=-Z_CHAR, vmax=Z_CHAR),
        cmap=TRACE_CMAP,
    )
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cb_ax, orientation="vertical", extend="both")
    cb.set_label(
        r"signed enrichment over matched-strength noise floor "
        r"(tanh-compressed; green = trace, red = anti)",
        rotation=270, labelpad=20,
    )
    cb.set_ticks([-Z_CHAR, 0.0, Z_CHAR])
    cb.set_ticklabels([r"$-Z_c$", "0", r"$+Z_c$"])

    plt.subplots_adjust(
        left=0.02, right=0.92, top=0.96, bottom=0.03,
        wspace=0.04, hspace=0.12,
    )

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands" / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_joint_density_3d_grid.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
