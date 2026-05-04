#!/usr/bin/env python3
"""MRL landscape — three-panel headline figure.

Panel A: ``M̄_step(b, h_rel; J_min = 0.9)`` — strict step landscape.
Panel B: ``M̄_smooth(b, h_rel)``           — continuous co-primary
                                              (positive-half mean ΔJ).
Panel C: ``Π(b, h_rel; J_min = 0.9)``      — cohort proportion of
                                              patients with ≥ 1 retained
                                              subtree in the cell.

Annotations on every panel:
- Hatched overlay where ``h_rel > 0.7`` ("probe-bias suspect" range,
  see ``.agents/guides/02_methods/probe-bias-guide.md``).
- Red-bordered cells: ``Π ≥ 7/9`` (cohort-wide retention).
- Grey-bordered cells: ``n_pat < 7`` (underpowered for cohort claim).

Reads (from ``data/reports/imcoh_mrl/``):
    mrl_cohort_step.csv, mrl_cohort_smooth.csv, mrl_per_node.csv

Writes:
    data/reports/imcoh_mrl/figures/mrl_landscape.{pdf,png}
    data/reports/imcoh_mrl/figures/mrl_landscape_supplementary.{pdf,png}  (J_min sweep)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_mrl"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT
PROBE_BIAS_H = 0.7    # h_rel above this is probe-bias suspect
COHORT_PI = 7 / 9     # cohort-wide threshold


def pivot_band_bin(df: pd.DataFrame, value_col: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pivot long df → (n_bands, n_bins) array.

    Returns:
      values  — (|BANDS|, K) float, NaN where missing
      h_lo    — (K,) bin lower edges, ascending
      h_hi    — (K,) bin upper edges
    """
    p = df.pivot_table(index="band", columns="bin", values=value_col,
                       aggfunc="first")
    p = p.reindex(BANDS)
    bin_idx = sorted(p.columns)
    p = p[bin_idx]
    edges_lo = (df.groupby("bin")["h_lo"].first()
                .reindex(bin_idx).to_numpy())
    edges_hi = (df.groupby("bin")["h_hi"].first()
                .reindex(bin_idx).to_numpy())
    return p.to_numpy(dtype=float), edges_lo, edges_hi


def add_cell_borders(ax: plt.Axes, mask_red: np.ndarray, mask_grey: np.ndarray) -> None:
    """Draw red borders where mask_red, grey borders where mask_grey & not red."""
    n_rows, n_cols = mask_red.shape
    for i in range(n_rows):
        for j in range(n_cols):
            if mask_red[i, j]:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       fill=False, edgecolor="#c00000",
                                       linewidth=1.6, zorder=4))
            elif mask_grey[i, j]:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       fill=False, edgecolor="#888888",
                                       linewidth=0.8, linestyle=":", zorder=3))


def add_probe_bias_hatch(ax: plt.Axes, h_lo: np.ndarray, h_hi: np.ndarray,
                         n_rows: int) -> None:
    """Hatched overlay on bins whose mid-edge sits above PROBE_BIAS_H."""
    h_mid = np.exp(0.5 * (np.log(h_lo) + np.log(h_hi)))
    for j, hm in enumerate(h_mid):
        if hm > PROBE_BIAS_H:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, n_rows,
                                   fill=False, hatch="///",
                                   edgecolor="#aaaaaa", linewidth=0,
                                   zorder=2))


def format_axes(ax: plt.Axes, h_lo: np.ndarray, h_hi: np.ndarray) -> None:
    n_bins = len(h_lo)
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([TEX[b] for b in BANDS], fontsize=11)
    # Show every other bin label to avoid crowding
    show = list(range(0, n_bins, 2))
    if (n_bins - 1) not in show:
        show.append(n_bins - 1)
    ax.set_xticks(show)
    ax.set_xticklabels([f"{h_lo[k]:.2f}" for k in show], fontsize=9, rotation=0)
    ax.set_xlabel(r"fractional dendrogram height $h_{\mathrm{rel}}$"
                  r" (lower edge, log-spaced)", fontsize=10)
    ax.set_ylabel("band", fontsize=10)


def headline_figure(step_df: pd.DataFrame, smooth_df: pd.DataFrame,
                    primary_jmin: float = 0.9) -> None:
    """Three-panel headline: M_step / M_smooth / Π at J_min = primary_jmin."""
    sub = step_df[step_df["J_min"] == primary_jmin].copy()
    M_step,  h_lo, h_hi = pivot_band_bin(sub, "M_step")
    Pi,      _,    _    = pivot_band_bin(sub, "Pi")
    M_smooth, _,  _     = pivot_band_bin(smooth_df, "M_smooth")
    npat,    _,    _    = pivot_band_bin(sub, "n_pat")
    cohort_red = (Pi >= COHORT_PI) & np.isfinite(Pi)
    underpow = (npat < 7) & np.isfinite(npat)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4), dpi=150)
    ax_a, ax_b, ax_c = axes

    n_bands, n_bins = M_step.shape

    im_a = ax_a.imshow(M_step, cmap="viridis", aspect="auto",
                       vmin=0, vmax=max(0.05, float(np.nanmax(M_step) or 0)))
    add_probe_bias_hatch(ax_a, h_lo, h_hi, n_bands)
    add_cell_borders(ax_a, cohort_red, underpow)
    format_axes(ax_a, h_lo, h_hi)
    cb_a = fig.colorbar(im_a, ax=ax_a, pad=0.02, shrink=0.9)
    cb_a.set_label(r"$\bar{M}_{\mathrm{step}}$ (retained fraction)", fontsize=9)
    ax_a.set_title(rf"A — $\bar{{M}}_{{\mathrm{{step}}}}$  ($J_{{\min}}={primary_jmin}$)",
                   fontsize=11)

    im_b = ax_b.imshow(M_smooth, cmap="magma", aspect="auto",
                       vmin=0, vmax=max(0.05, float(np.nanmax(M_smooth) or 0)))
    add_probe_bias_hatch(ax_b, h_lo, h_hi, n_bands)
    add_cell_borders(ax_b, cohort_red, underpow)
    format_axes(ax_b, h_lo, h_hi)
    cb_b = fig.colorbar(im_b, ax=ax_b, pad=0.02, shrink=0.9)
    cb_b.set_label(r"$\bar{M}_{\mathrm{smooth}} = \langle\max(\Delta J,0)\rangle$",
                   fontsize=9)
    ax_b.set_title(r"B — $\bar{M}_{\mathrm{smooth}}$  (continuous; no $J_{\min}$)",
                   fontsize=11)

    im_c = ax_c.imshow(Pi, cmap="cividis", aspect="auto", vmin=0, vmax=1)
    add_probe_bias_hatch(ax_c, h_lo, h_hi, n_bands)
    add_cell_borders(ax_c, cohort_red, underpow)
    format_axes(ax_c, h_lo, h_hi)
    cb_c = fig.colorbar(im_c, ax=ax_c, pad=0.02, shrink=0.9)
    cb_c.set_label(r"$\Pi$ (fraction patients with $\geq 1$ retained)", fontsize=9)
    ax_c.set_title(rf"C — $\Pi$  ($J_{{\min}}={primary_jmin}$;  red border = $\Pi\geq 7/9$)",
                   fontsize=11)

    for im in (im_a, im_b, im_c):
        im.set_rasterized(True)
    fig.tight_layout()
    out = OUT_DIR / "mrl_landscape.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def supplementary_jmin_sweep(step_df: pd.DataFrame) -> None:
    """One row × N_J columns of M_step heatmaps across the J_min sweep."""
    j_vals = sorted(step_df["J_min"].unique())
    fig, axes = plt.subplots(1, len(j_vals), figsize=(5 * len(j_vals), 4.4),
                             dpi=150, squeeze=False)
    axes = axes[0]

    for ax, J in zip(axes, j_vals):
        sub = step_df[step_df["J_min"] == J].copy()
        M_step, h_lo, h_hi = pivot_band_bin(sub, "M_step")
        Pi,     _,    _    = pivot_band_bin(sub, "Pi")
        npat,   _,    _    = pivot_band_bin(sub, "n_pat")
        cohort_red = (Pi >= COHORT_PI) & np.isfinite(Pi)
        underpow = (npat < 7) & np.isfinite(npat)

        im = ax.imshow(M_step, cmap="viridis", aspect="auto",
                       vmin=0, vmax=max(0.05, float(np.nanmax(M_step) or 0)),
                       rasterized=True)
        add_probe_bias_hatch(ax, h_lo, h_hi, len(BANDS))
        add_cell_borders(ax, cohort_red, underpow)
        format_axes(ax, h_lo, h_hi)
        fig.colorbar(im, ax=ax, pad=0.02, shrink=0.9)
        ax.set_title(rf"$\bar{{M}}_{{\mathrm{{step}}}}$ at $J_{{\min}}={J}$",
                     fontsize=11)

    fig.tight_layout()
    out = OUT_DIR / "mrl_landscape_supplementary.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def per_patient_strips(per_node_df: pd.DataFrame, step_df: pd.DataFrame,
                       primary_jmin: float = 0.9) -> None:
    """6-panel per-patient strip of m_p_step(b, bin) at primary J_min."""
    edges_df = step_df[step_df["J_min"] == primary_jmin].drop_duplicates("bin")
    edges_df = edges_df.sort_values("bin")
    h_lo = edges_df["h_lo"].to_numpy()
    h_hi = edges_df["h_hi"].to_numpy()
    K = len(h_lo)

    df = per_node_df.copy()
    # Recompute bin assignment in the figure script using the saved edges.
    df["bin"] = -1
    for k in range(K):
        if k < K - 1:
            mask = (df["h_rel"] >= h_lo[k]) & (df["h_rel"] < h_hi[k])
        else:
            mask = (df["h_rel"] >= h_lo[k]) & (df["h_rel"] <= h_hi[k])
        df.loc[mask, "bin"] = k
    df = df[df["bin"] >= 0]
    df["retained"] = ((df["J_pre"] < primary_jmin)
                      & (df["J_post"] >= primary_jmin))

    patients = sorted(df["patient"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 6.5), dpi=150,
                             sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, band in zip(axes, BANDS):
        cell_grid = np.full((len(patients), K), np.nan)
        for i, pat in enumerate(patients):
            sub = df[(df["patient"] == pat) & (df["band"] == band)]
            if sub.empty:
                continue
            for k in range(K):
                in_bin = sub[sub["bin"] == k]
                if in_bin.empty:
                    continue
                cell_grid[i, k] = float(in_bin["retained"].mean())

        im = ax.imshow(cell_grid, cmap="viridis", aspect="auto",
                       vmin=0, vmax=max(0.1, float(np.nanmax(cell_grid) or 0)),
                       rasterized=True)
        ax.set_yticks(range(len(patients)))
        ax.set_yticklabels(patients, fontsize=8)
        show = list(range(0, K, 2))
        if (K - 1) not in show:
            show.append(K - 1)
        ax.set_xticks(show)
        ax.set_xticklabels([f"{h_lo[k]:.2f}" for k in show],
                           fontsize=8, rotation=0)
        ax.set_title(TEX[band], fontsize=11)
        fig.colorbar(im, ax=ax, pad=0.02, shrink=0.9)

    for ax in axes[3:]:
        ax.set_xlabel(r"$h_{\mathrm{rel}}$", fontsize=10)

    fig.tight_layout()
    out = OUT_DIR / "mrl_per_patient_strips.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def main() -> None:
    step_df = pd.read_csv(IN_DIR / "mrl_cohort_step.csv")
    smooth_df = pd.read_csv(IN_DIR / "mrl_cohort_smooth.csv")
    per_node_df = pd.read_csv(IN_DIR / "mrl_per_node.csv")

    headline_figure(step_df, smooth_df, primary_jmin=0.9)
    supplementary_jmin_sweep(step_df)
    per_patient_strips(per_node_df, step_df, primary_jmin=0.9)


if __name__ == "__main__":
    main()
