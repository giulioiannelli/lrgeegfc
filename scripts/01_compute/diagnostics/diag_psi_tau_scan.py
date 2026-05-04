#!/usr/bin/env python
"""Psi(n;tau) and Psi_L(tau) diagnostic across the full LRG tau-range.

Empirical sanity check on whether argmax_n Psi gives non-trivial partitions on
imcoh_abs FC dendrograms when tau is scanned, or whether it is pinned at the
first split (probe edge / single big jump). Also reports Psi_L statistics so
we know whether cluster-level stability has any usable signal beyond Psi.

Outputs
-------
data/audit/psi_tau_scan/psi_grid.csv  -- one row per (patient, band, phase, tau)
data/audit/psi_tau_scan/psi_diagnostic.pdf -- multi-page diagnostic
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.visuals.lrg import compute_partition_stability_index
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrgsglib.nx_patches.funcs import get_giant_component

PATIENTS = ["Pat_02", "Pat_06", "Pat_03"]
PHASES = ["rest_pre", "task_test", "rest_post"]
BANDS = list(BRAIN_BANDS_NAMES)
N_TAU = 8
FC_METHOD = "imcoh_abs"

OUT_DIR = DATA_ROOT / "audit" / "psi_tau_scan"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_PATH = OUT_DIR / "psi_grid.csv"
PDF_PATH = OUT_DIR / "psi_diagnostic.pdf"


def rho_tau(eigvals: np.ndarray, eigvecs: np.ndarray, tau: float) -> np.ndarray:
    w = np.exp(-tau * eigvals)
    Z = w.sum()
    return (eigvecs * (w / Z)) @ eigvecs.T


def communication_distance(rho: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    finite = np.isfinite(D) & (D > 0)
    if finite.any():
        cap = D[finite].max() * 1e3
    else:
        cap = 1.0
    D[~np.isfinite(D)] = cap
    D[D < 0] = cap
    return D


def psi_local_per_node(Z: np.ndarray) -> np.ndarray:
    n_internal = Z.shape[0]
    n_leaves = n_internal + 1
    heights = Z[:, 2]
    parent = np.full(n_internal, -1, dtype=int)
    for j in range(n_internal):
        c0, c1 = int(Z[j, 0]), int(Z[j, 1])
        if c0 >= n_leaves:
            parent[c0 - n_leaves] = j
        if c1 >= n_leaves:
            parent[c1 - n_leaves] = j
    eps = 1e-12
    log_self = np.log10(np.clip(heights, eps, None))
    psi_L = np.full(n_internal, np.nan)
    has_parent = parent >= 0
    log_parent = np.log10(np.clip(heights[parent[has_parent]], eps, None))
    psi_L[has_parent] = log_parent - log_self[has_parent]
    return psi_L


def scan_one(
    fc_mat: np.ndarray, n_tau: int = N_TAU
) -> Tuple[pd.DataFrame, List[Tuple[np.ndarray, np.ndarray]], List[np.ndarray]]:
    A = np.abs(fc_mat).astype(float)
    np.fill_diagonal(A, 0.0)
    G = nx.from_numpy_array(A)
    G = get_giant_component(G)
    L = nx.laplacian_matrix(G, weight="weight").toarray().astype(float)
    eigvals, eigvecs = np.linalg.eigh(L)
    eigvals = np.clip(eigvals, 0.0, None)
    nonzero = eigvals > eigvals.max() * 1e-10
    lam_max = float(eigvals.max())
    lam_gap = float(eigvals[nonzero].min())
    tau_min = 1.0 / lam_max
    tau_max = 1.0 / lam_gap
    taus = np.geomspace(tau_min, tau_max, n_tau)

    rows: List[dict] = []
    psi_grid: List[Tuple[np.ndarray, np.ndarray]] = []
    psi_L_pool: List[np.ndarray] = []
    for ti, tau in enumerate(taus):
        rho = rho_tau(eigvals, eigvecs, tau)
        D = communication_distance(rho)
        dists = squareform(D, checks=False)
        Z = linkage(dists, method="average")
        psi_vals, n_comms = compute_partition_stability_index(Z)
        if len(psi_vals) == 0:
            continue
        sorted_psi = np.sort(psi_vals)[::-1]
        top1_top2 = float(sorted_psi[0] / max(sorted_psi[1], 1e-12)) if len(psi_vals) >= 2 else np.nan
        argmax_idx = int(np.argmax(psi_vals))
        psi_L = psi_local_per_node(Z)
        rows.append(
            dict(
                tau_idx=ti,
                tau=float(tau),
                tau_norm=float(tau * lam_max),
                lam_max=lam_max,
                lam_gap=lam_gap,
                argmax_n=int(n_comms[argmax_idx]),
                psi_max=float(psi_vals[argmax_idx]),
                psi_top1_over_top2=top1_top2,
                psi_L_median=float(np.nanmedian(psi_L)),
                psi_L_p90=float(np.nanpercentile(psi_L, 90)),
                psi_L_p95=float(np.nanpercentile(psi_L, 95)),
                n_psi_L_above_05=int(np.nansum(psi_L > 0.5)),
                n_psi_L_above_10=int(np.nansum(psi_L > 1.0)),
                n_internal=int(len(psi_L)),
            )
        )
        psi_grid.append((n_comms, psi_vals))
        psi_L_pool.append(psi_L)
    return pd.DataFrame(rows), psi_grid, psi_L_pool


def collect() -> Tuple[pd.DataFrame, dict]:
    all_rows: List[pd.DataFrame] = []
    psi_grids: dict = {}
    psi_L_pools: dict = {}
    for pat in PATIENTS:
        for band in BANDS:
            for phase in PHASES:
                fc = load_fc_matrix(pat, phase, band, FC_METHOD)
                if fc is None:
                    print(f"[skip] {pat} {phase} {band}: no FC cache")
                    continue
                try:
                    df, grid, pool = scan_one(fc)
                except Exception as exc:  # noqa: BLE001
                    print(f"[error] {pat} {phase} {band}: {exc}")
                    continue
                df["patient"] = pat
                df["band"] = band
                df["phase"] = phase
                all_rows.append(df)
                psi_grids[(pat, band, phase)] = grid
                psi_L_pools[(pat, band, phase)] = pool
                print(f"[ok] {pat} {phase} {band}: argmax_n={list(df.argmax_n)}")
    return pd.concat(all_rows, ignore_index=True), {"grids": psi_grids, "pools": psi_L_pools}


def fig_argmax_summary(df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(
        nrows=len(PATIENTS),
        ncols=len(BANDS),
        figsize=(3.0 * len(BANDS), 2.4 * len(PATIENTS)),
        sharey=True,
    )
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    phase_color = {"rest_pre": "#4477aa", "task_test": "#ee6677", "rest_post": "#228833"}
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            for phase in PHASES:
                sub = df[(df.patient == pat) & (df.band == band) & (df.phase == phase)].sort_values("tau_norm")
                if sub.empty:
                    continue
                ax.plot(
                    sub.tau_norm,
                    sub.argmax_n,
                    "-o",
                    color=phase_color[phase],
                    ms=4,
                    lw=1.2,
                    label=phase if (i == 0 and j == 0) else None,
                )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim(1.5, 200)
            ax.axhline(2, color="grey", ls=":", lw=0.6)
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == len(PATIENTS) - 1:
                ax.set_xlabel(r"$\tau \cdot \lambda_{\max}$")
            if j == 0:
                ax.set_ylabel(r"argmax$_n \,\Psi(n;\tau)$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return fig


def fig_top1_over_top2(df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(
        nrows=len(PATIENTS),
        ncols=len(BANDS),
        figsize=(3.0 * len(BANDS), 2.4 * len(PATIENTS)),
        sharey=True,
    )
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    phase_color = {"rest_pre": "#4477aa", "task_test": "#ee6677", "rest_post": "#228833"}
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            for phase in PHASES:
                sub = df[(df.patient == pat) & (df.band == band) & (df.phase == phase)].sort_values("tau_norm")
                if sub.empty:
                    continue
                ax.plot(
                    sub.tau_norm,
                    sub.psi_top1_over_top2,
                    "-o",
                    color=phase_color[phase],
                    ms=4,
                    lw=1.2,
                )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.axhline(1.0, color="grey", ls=":", lw=0.6)
            ax.axhline(2.0, color="grey", ls="--", lw=0.6)
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == len(PATIENTS) - 1:
                ax.set_xlabel(r"$\tau \cdot \lambda_{\max}$")
            if j == 0:
                ax.set_ylabel(r"$\Psi_{(1)}/\Psi_{(2)}$")
    fig.tight_layout()
    return fig


def fig_psi_L_summary(df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(
        nrows=len(PATIENTS),
        ncols=len(BANDS),
        figsize=(3.0 * len(BANDS), 2.4 * len(PATIENTS)),
        sharey=True,
    )
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    phase_color = {"rest_pre": "#4477aa", "task_test": "#ee6677", "rest_post": "#228833"}
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            for phase in PHASES:
                sub = df[(df.patient == pat) & (df.band == band) & (df.phase == phase)].sort_values("tau_norm")
                if sub.empty:
                    continue
                ax.plot(
                    sub.tau_norm,
                    sub.n_psi_L_above_05,
                    "-o",
                    color=phase_color[phase],
                    ms=4,
                    lw=1.2,
                )
            ax.set_xscale("log")
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == len(PATIENTS) - 1:
                ax.set_xlabel(r"$\tau \cdot \lambda_{\max}$")
            if j == 0:
                ax.set_ylabel(r"$\#\{\Psi_L > 0.5\}$")
    fig.tight_layout()
    return fig


def fig_psi_heatmaps(grids: dict, pat: str = "Pat_02") -> plt.Figure:
    fig, axes = plt.subplots(
        nrows=len(PHASES),
        ncols=len(BANDS),
        figsize=(2.4 * len(BANDS), 2.0 * len(PHASES)),
    )
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    for i, phase in enumerate(PHASES):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            grid = grids.get((pat, band, phase))
            if grid is None or len(grid) == 0:
                ax.set_facecolor("#eeeeee")
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            n_tau = len(grid)
            n_max = max(len(p) for _, p in grid)
            mat = np.full((n_max, n_tau), np.nan)
            for ti, (n_c, p) in enumerate(grid):
                mat[: len(p), ti] = p
            im = ax.imshow(
                mat,
                aspect="auto",
                origin="lower",
                cmap="viridis",
                extent=[-0.5, n_tau - 0.5, 1.5, n_max + 1.5],
            )
            im.set_rasterized(True)
            if i == 0:
                ax.set_title(band, fontsize=8)
            if j == 0:
                ax.set_ylabel(f"{phase}\nn-clusters", fontsize=7)
            if i == len(PHASES) - 1:
                ax.set_xlabel(r"$\tau$ idx", fontsize=7)
            ax.set_yscale("log")
    fig.tight_layout()
    return fig


def main() -> None:
    df, payload = collect()
    df.to_csv(CSV_PATH, index=False)
    print(f"[wrote] {CSV_PATH} ({len(df)} rows)")
    with PdfPages(PDF_PATH) as pdf:
        for fig in (
            fig_argmax_summary(df),
            fig_top1_over_top2(df),
            fig_psi_L_summary(df),
            fig_psi_heatmaps(payload["grids"], pat="Pat_02"),
        ):
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"[wrote] {PDF_PATH}")


if __name__ == "__main__":
    main()
