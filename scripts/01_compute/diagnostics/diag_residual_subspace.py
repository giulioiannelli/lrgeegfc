#!/usr/bin/env python
"""Residual-subspace alignment diagnostic.

For each (patient, band, τ̃), eigendecompose the cross-phase residual
ultrametric matrices S = D^task - D^pre and R = D^post - D^pre, then test
whether their leading eigenvector subspaces align (group-wise task-trace
signature) and whether each residual has dominant-mode energy concentration
(group-existence test).

α_k(τ̃) = mean canonical correlation between span(top-k u^S) and span(top-k u^R).
β_M(τ̃) = leading-mode-eigenvalue² / sum(eigenvalue²) — energy concentration.

Phase-label permutation null on α_1 included as a coarse null.

Outputs
-------
data/audit/residual_subspace/alpha_beta_grid.csv
data/audit/residual_subspace/loadings.npz   -- u1_S, u1_R per (pat, band, τ̃)
data/audit/residual_subspace/diagnostic.pdf -- 4 pages
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrgsglib.nx_patches.funcs import get_giant_component

PATIENTS = ["Pat_02", "Pat_06", "Pat_03"]
PHASES = ("rest_pre", "task_test", "rest_post")
BANDS = list(BRAIN_BANDS_NAMES)
N_TAU = 12
FC_METHOD = "imcoh_abs"
N_NULL_PERMS = 5
RNG_SEED = 20260428

OUT_DIR = DATA_ROOT / "audit" / "residual_subspace"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_PATH = OUT_DIR / "alpha_beta_grid.csv"
NPZ_PATH = OUT_DIR / "loadings.npz"
PDF_PATH = OUT_DIR / "diagnostic.pdf"


def laplacian_eig(fc: np.ndarray):
    A = np.abs(fc).astype(float)
    np.fill_diagonal(A, 0.0)
    G = nx.from_numpy_array(A)
    G_giant = get_giant_component(G)
    nodes_kept = sorted(G_giant.nodes())
    L = nx.laplacian_matrix(G_giant, nodelist=nodes_kept, weight="weight").toarray().astype(float)
    eigvals, eigvecs = np.linalg.eigh(L)
    eigvals = np.clip(eigvals, 0.0, None)
    lam_max = float(eigvals.max())
    nonzero = eigvals > lam_max * 1e-10
    lam_gap = float(eigvals[nonzero].min()) if nonzero.any() else lam_max
    return eigvals, eigvecs, lam_max, lam_gap, nodes_kept


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
    cap = D[finite].max() * 1e3 if finite.any() else 1.0
    D[~np.isfinite(D)] = cap
    D[D < 0] = cap
    return D


def subspace_align(U_a: np.ndarray, U_b: np.ndarray, k: int) -> float:
    Ua = U_a[:, :k]
    Ub = U_b[:, :k]
    sigmas = np.linalg.svd(Ua.T @ Ub, compute_uv=False)
    return float(np.mean(np.clip(sigmas, 0.0, 1.0)))


def eig_residual(M: np.ndarray):
    M = 0.5 * (M + M.T)
    eigvals, eigvecs = np.linalg.eigh(M)
    order = np.argsort(np.abs(eigvals))[::-1]
    return eigvals[order], eigvecs[:, order]


def scan_patient_band(pat: str, band: str, rng: np.random.Generator) -> Tuple[List[dict], dict]:
    eigs: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    lam_max: Dict[str, float] = {}
    lam_gap: Dict[str, float] = {}
    nodes_kept: Dict[str, list] = {}
    for phase in PHASES:
        fc = load_fc_matrix(pat, phase, band, FC_METHOD)
        if fc is None:
            print(f"[skip] {pat} {phase} {band}: no FC cache")
            return [], {}
        ev, eu, lmax, lgap, nks = laplacian_eig(fc)
        eigs[phase] = (ev, eu)
        lam_max[phase] = lmax
        lam_gap[phase] = lgap
        nodes_kept[phase] = list(nks)

    common = sorted(set.intersection(*[set(nodes_kept[p]) for p in PHASES]))
    if len(common) < 10:
        print(f"[skip] {pat} {band}: only {len(common)} common nodes")
        return [], {}

    common_idx_per_phase: Dict[str, np.ndarray] = {}
    for phase in PHASES:
        position = {n: i for i, n in enumerate(nodes_kept[phase])}
        common_idx_per_phase[phase] = np.array([position[c] for c in common])

    tau_tilde_max = min(lam_max[p] / lam_gap[p] for p in PHASES)
    tau_tilde_grid = np.geomspace(1.0, max(tau_tilde_max, 2.0), N_TAU)

    rows: List[dict] = []
    loadings_per_tau: Dict[int, dict] = {}
    for ti, tau_tilde in enumerate(tau_tilde_grid):
        D_phase: Dict[str, np.ndarray] = {}
        for phase in PHASES:
            ev, eu = eigs[phase]
            tau_phi = tau_tilde / lam_max[phase]
            rho_full = rho_tau(ev, eu, tau_phi)
            D_full = communication_distance(rho_full)
            idx = common_idx_per_phase[phase]
            D_phase[phase] = D_full[np.ix_(idx, idx)]

        S = D_phase["task_test"] - D_phase["rest_pre"]
        R = D_phase["rest_post"] - D_phase["rest_pre"]
        sv, su = eig_residual(S)
        rv, ru = eig_residual(R)
        beta_S = float(sv[0] ** 2 / (sv ** 2).sum())
        beta_R = float(rv[0] ** 2 / (rv ** 2).sum())
        alpha_1 = subspace_align(su, ru, 1)
        alpha_3 = subspace_align(su, ru, 3)
        alpha_5 = subspace_align(su, ru, 5)

        # Phase-label permutation null on alpha_1
        null_alphas = []
        D_list = [D_phase["rest_pre"], D_phase["task_test"], D_phase["rest_post"]]
        for _ in range(N_NULL_PERMS):
            perm = list(rng.permutation(3))
            ref = D_list[perm[0]]
            S_n = D_list[perm[1]] - ref
            R_n = D_list[perm[2]] - ref
            _, su_n = eig_residual(S_n)
            _, ru_n = eig_residual(R_n)
            null_alphas.append(subspace_align(su_n, ru_n, 1))
        null_alpha_1_median = float(np.median(null_alphas))
        null_alpha_1_p95 = float(np.percentile(null_alphas, 95))

        rows.append(
            dict(
                patient=pat, band=band,
                tau_idx=ti, tau_tilde=float(tau_tilde),
                n_common=len(common),
                lam_max_pre=lam_max["rest_pre"],
                lam_max_task=lam_max["task_test"],
                lam_max_post=lam_max["rest_post"],
                alpha_1=alpha_1, alpha_3=alpha_3, alpha_5=alpha_5,
                beta_S=beta_S, beta_R=beta_R,
                null_alpha_1_median=null_alpha_1_median,
                null_alpha_1_p95=null_alpha_1_p95,
                sv_top=float(sv[0]), rv_top=float(rv[0]),
                sv_top2=float(sv[1]) if len(sv) > 1 else float("nan"),
                rv_top2=float(rv[1]) if len(rv) > 1 else float("nan"),
            )
        )
        loadings_per_tau[ti] = dict(
            u1_S=su[:, 0].astype(np.float32),
            u1_R=ru[:, 0].astype(np.float32),
        )

    return rows, dict(
        common_nodes=np.array(common),
        loadings=loadings_per_tau,
        tau_tilde_grid=tau_tilde_grid,
    )


def collect():
    rng = np.random.default_rng(RNG_SEED)
    all_rows: List[pd.DataFrame] = []
    npz_payload: Dict[str, np.ndarray] = {}
    for pat in PATIENTS:
        for band in BANDS:
            print(f"[scan] {pat} {band}")
            rows, payload = scan_patient_band(pat, band, rng)
            if not rows:
                continue
            all_rows.append(pd.DataFrame(rows))
            npz_payload[f"{pat}__{band}__common_nodes"] = payload["common_nodes"]
            npz_payload[f"{pat}__{band}__tau_tilde_grid"] = payload["tau_tilde_grid"]
            for ti, ld in payload["loadings"].items():
                npz_payload[f"{pat}__{band}__u1_S__t{ti:02d}"] = ld["u1_S"]
                npz_payload[f"{pat}__{band}__u1_R__t{ti:02d}"] = ld["u1_R"]
    df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    return df, npz_payload


PHASE_COLOR = {"rest_pre": "#4477aa", "task_test": "#ee6677", "rest_post": "#228833"}
PATIENT_LS = {"Pat_02": "-", "Pat_06": "--", "Pat_03": ":"}


def fig_alpha_vs_null(df: pd.DataFrame) -> plt.Figure:
    n_p, n_b = len(PATIENTS), len(BANDS)
    fig, axes = plt.subplots(n_p, n_b, figsize=(2.6 * n_b, 2.2 * n_p), sharey=True)
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = df[(df.patient == pat) & (df.band == band)].sort_values("tau_tilde")
            if sub.empty:
                ax.set_facecolor("#f0f0f0"); ax.set_xticks([]); ax.set_yticks([]); continue
            ax.plot(sub.tau_tilde, sub.alpha_1, "-o", color="#222288", ms=4, lw=1.4, label="α₁ obs")
            ax.fill_between(sub.tau_tilde, 0, sub.null_alpha_1_p95, color="#cccccc", alpha=0.55, label="null p95")
            ax.plot(sub.tau_tilde, sub.null_alpha_1_median, "--", color="#888888", lw=1.0, label="null med")
            ax.set_xscale("log"); ax.set_ylim(0, 1.05)
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == n_p - 1:
                ax.set_xlabel(r"$\tilde{\tau} = \tau \cdot \lambda_{\max}$", fontsize=8)
            if j == 0:
                ax.set_ylabel(r"$\alpha_1(\tilde{\tau})$", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return fig


def fig_beta(df: pd.DataFrame) -> plt.Figure:
    n_p, n_b = len(PATIENTS), len(BANDS)
    fig, axes = plt.subplots(n_p, n_b, figsize=(2.6 * n_b, 2.2 * n_p), sharey=True)
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = df[(df.patient == pat) & (df.band == band)].sort_values("tau_tilde")
            if sub.empty:
                ax.set_facecolor("#f0f0f0"); ax.set_xticks([]); ax.set_yticks([]); continue
            ax.plot(sub.tau_tilde, sub.beta_S, "-o", color="#ee6677", ms=4, lw=1.2, label="β_S task−pre")
            ax.plot(sub.tau_tilde, sub.beta_R, "--s", color="#228833", ms=4, lw=1.2, label="β_R post−pre")
            ax.axhline(0.1, color="grey", ls=":", lw=0.5)
            ax.set_xscale("log"); ax.set_ylim(0, 1.05)
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == n_p - 1:
                ax.set_xlabel(r"$\tilde{\tau}$", fontsize=8)
            if j == 0:
                ax.set_ylabel(r"$\beta_M$", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return fig


def fig_alpha_k_compare(df: pd.DataFrame) -> plt.Figure:
    n_p, n_b = len(PATIENTS), len(BANDS)
    fig, axes = plt.subplots(n_p, n_b, figsize=(2.6 * n_b, 2.2 * n_p), sharey=True)
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = df[(df.patient == pat) & (df.band == band)].sort_values("tau_tilde")
            if sub.empty:
                ax.set_facecolor("#f0f0f0"); ax.set_xticks([]); ax.set_yticks([]); continue
            ax.plot(sub.tau_tilde, sub.alpha_1, "-o", color="#222288", ms=4, lw=1.2, label="α₁")
            ax.plot(sub.tau_tilde, sub.alpha_3, "--", color="#aa5599", lw=1.2, label="α₃")
            ax.plot(sub.tau_tilde, sub.alpha_5, ":", color="#cc8855", lw=1.2, label="α₅")
            ax.set_xscale("log"); ax.set_ylim(0, 1.05)
            ax.set_title(f"{pat} {band}", fontsize=8)
            if i == n_p - 1:
                ax.set_xlabel(r"$\tilde{\tau}$", fontsize=8)
            if j == 0:
                ax.set_ylabel(r"$\alpha_k$", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return fig


def fig_loadings_at_peak(df: pd.DataFrame, npz: dict) -> plt.Figure:
    n_p, n_b = len(PATIENTS), len(BANDS)
    fig, axes = plt.subplots(n_p, n_b, figsize=(2.6 * n_b, 2.2 * n_p))
    if axes.ndim == 1:
        axes = axes.reshape(1, -1)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = df[(df.patient == pat) & (df.band == band)].sort_values("tau_idx")
            if sub.empty:
                ax.set_facecolor("#f0f0f0"); ax.set_xticks([]); ax.set_yticks([]); continue
            best = sub.loc[(sub.alpha_1 - sub.null_alpha_1_p95).idxmax()]
            ti = int(best.tau_idx)
            key_R = f"{pat}__{band}__u1_R__t{ti:02d}"
            key_S = f"{pat}__{band}__u1_S__t{ti:02d}"
            if key_R not in npz:
                ax.set_facecolor("#f0f0f0"); ax.set_xticks([]); ax.set_yticks([]); continue
            uR = np.abs(npz[key_R]); uS = np.abs(npz[key_S])
            order = np.argsort(uR)[::-1]
            ax.plot(np.arange(len(uR)), uR[order], "-", color="#228833", lw=1.2, label="|u₁ᴿ|")
            ax.plot(np.arange(len(uS)), uS[order], "--", color="#ee6677", lw=0.9, alpha=0.8, label="|u₁ˢ| (R-order)")
            ax.set_yscale("log")
            ax.set_title(f"{pat} {band}\nτ̃-idx={ti} α₁={best.alpha_1:.2f}", fontsize=7)
            if i == n_p - 1:
                ax.set_xlabel("rank", fontsize=8)
            if j == 0:
                ax.set_ylabel(r"$|u_1|$", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    return fig


def main():
    df, npz_payload = collect()
    if df.empty:
        print("[!] No rows produced; aborting.")
        return
    df.to_csv(CSV_PATH, index=False)
    print(f"[wrote] {CSV_PATH} ({len(df)} rows)")
    np.savez_compressed(NPZ_PATH, **npz_payload)
    print(f"[wrote] {NPZ_PATH}")
    with PdfPages(PDF_PATH) as pdf:
        for fig in (
            fig_alpha_vs_null(df),
            fig_beta(df),
            fig_alpha_k_compare(df),
            fig_loadings_at_peak(df, npz_payload),
        ):
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"[wrote] {PDF_PATH}")


if __name__ == "__main__":
    main()
