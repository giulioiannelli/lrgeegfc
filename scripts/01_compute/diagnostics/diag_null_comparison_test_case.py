#!/usr/bin/env python
"""Specific-case dual-null comparison: matched-strength vs weight-permutation.

Question (one cell, Pat_02 β rest_pre↔rest_post):

  Does the cross-phase Grassmann probe give the same statistical picture
  under (a) the current matched-strength null (preserves strengths,
  drifts weight distribution) and (b) a pure weight-permutation null
  (preserves weight multiset exactly, drifts strengths)?

If both nulls reject observed at similar p-values, the matched-strength
weight-distribution drift documented in
``data/audit/matched_strength_weight_drift/`` does NOT materially affect
the verdict — and we can skip the (expensive) joint-null implementation.
If the two nulls diverge, the joint null is needed.

This script does NOT touch any cached matched-strength surrogate output.
It generates fresh ensembles for both nulls in parallel, runs them on a
single cell, and writes a side-by-side comparison.

Outputs:
- ``data/audit/null_comparison_test_case/per_k_pvalues.csv`` — one row
  per k with observed Grassmann distance and lower-tail p under each null.
- ``data/audit/null_comparison_test_case/grassmann_distance_comparison.pdf``
  — for each k, surrogate-distance histograms (matched-strength and
  weight-permutation) overlaid with observed distance.
- ``data/audit/null_comparison_test_case/README.md`` — summary +
  decision (verdict robust / divergent).

Runtime: ~10 minutes at R=200, K=6 k-values, on one phase pair.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.surrogate import strength_preserving_shuffle, verify_strengths
from lrg_eegfc.workflow.fc import load_fc_matrix


OUT_ROOT = Path("data/audit/null_comparison_test_case")
DEFAULT_PATIENT = "Pat_02"
DEFAULT_BAND = "beta"
DEFAULT_PHASE_A = "rest_pre"
DEFAULT_PHASE_B = "rest_post"
DEFAULT_R = 200
DEFAULT_SWAP_FACTOR = 20
DEFAULT_SEED = 20260528
K_LIST = (2, 4, 8, 16, 32, 48)


def weight_permutation_shuffle(W: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Random permutation of upper-triangular edge weights.

    Preserves the edge-weight multiset exactly. Per-node strengths are
    randomized (no constraint).
    """
    iu = np.triu_indices_from(W, k=1)
    weights = W[iu].copy()
    rng.shuffle(weights)
    W_perm = np.zeros_like(W)
    W_perm[iu] = weights
    W_perm = W_perm + W_perm.T
    return W_perm


def laplacian_eigs(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Combinatorial Laplacian eigendecomposition. Returns (eigvals_asc, eigvecs)."""
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return np.linalg.eigh(L)


def grassmann_chordal_distance(U1: np.ndarray, U2: np.ndarray) -> float:
    """Chordal (sin-of-principal-angles) Grassmann distance.

    For U1, U2 (N, k) orthonormal:
        d = sqrt( k - ||U1^T U2||_F^2 )  (equivalent to ||sin θ||_F).
    """
    M = U1.T @ U2
    inner_sq = float(np.sum(M * M))
    return float(np.sqrt(max(U1.shape[1] - inner_sq, 0.0)))


def _generate_paired_surrogates(
    W_A: np.ndarray,
    W_B: np.ndarray,
    R: int,
    null: str,
    swap_factor: int,
    seed: int,
    verbose: bool = True,
) -> tuple[list[np.ndarray], list[np.ndarray], int]:
    """Generate R paired surrogates for two phases under one null algorithm.

    Returns lists of length R: (eigvecs_A_list, eigvecs_B_list, n_ok).
    """
    N = W_A.shape[0]
    n_swaps = swap_factor * (N * (N - 1)) // 2
    rng = np.random.default_rng(seed)

    U_A_list: list[np.ndarray] = []
    U_B_list: list[np.ndarray] = []
    n_ok = 0
    for r in range(R):
        if null == "matched_strength":
            W_A_s = strength_preserving_shuffle(W_A, n_swaps, rng)
            W_B_s = strength_preserving_shuffle(W_B, n_swaps, rng)
            ok = (verify_strengths(W_A, W_A_s, tol=1e-4)
                  and verify_strengths(W_B, W_B_s, tol=1e-4))
        elif null == "weight_permutation":
            W_A_s = weight_permutation_shuffle(W_A, rng)
            W_B_s = weight_permutation_shuffle(W_B, rng)
            ok = True
        else:
            raise ValueError(f"unknown null: {null}")
        if not ok:
            continue
        _, U_A_s = laplacian_eigs(W_A_s)
        _, U_B_s = laplacian_eigs(W_B_s)
        U_A_list.append(U_A_s)
        U_B_list.append(U_B_s)
        n_ok += 1
        if verbose and (r + 1) % 25 == 0:
            print(f"  [{null}] {r + 1}/{R}", flush=True)
    return U_A_list, U_B_list, n_ok


def _plot_comparison(
    obs_dG: dict, ms_dG: dict, wp_dG: dict, out_path: Path
) -> None:
    import matplotlib.pyplot as plt

    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()

    n_panels = len(K_LIST)
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.0))
    axes = axes.flatten()
    for i, k in enumerate(K_LIST):
        ax = axes[i]
        ms_arr = np.array(ms_dG[k])
        wp_arr = np.array(wp_dG[k])
        lo = min(ms_arr.min(), wp_arr.min(), obs_dG[k]) * 0.95
        hi = max(ms_arr.max(), wp_arr.max(), obs_dG[k]) * 1.05
        bins = np.linspace(lo, hi, 30)
        ax.hist(ms_arr, bins=bins, density=True, alpha=0.55,
                color="#d62728", label="matched-strength" if i == 0 else None)
        ax.hist(wp_arr, bins=bins, density=True, alpha=0.55,
                color="#1f77b4",
                label="weight-permutation" if i == 0 else None)
        ax.axvline(obs_dG[k], color="black", linestyle="--", linewidth=1.5,
                   label="observed" if i == 0 else None)
        ax.set_title(f"k = {k}")
        ax.set_xlabel(r"$d_G(U_k^{\rm A}, U_k^{\rm B})$")
        if i % 3 == 0:
            ax.set_ylabel("density")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def _write_summary(
    patient: str, band: str, phase_A: str, phase_B: str, R: int,
    obs_dG: dict, ms_dG: dict, wp_dG: dict, df: pd.DataFrame, out_path: Path,
) -> None:
    # Decision rule on p-values
    ms_ps = df["p_lower_matched_strength"].to_numpy()
    wp_ps = df["p_lower_weight_permutation"].to_numpy()
    same_rejection = np.all((ms_ps < 0.05) == (wp_ps < 0.05))
    p_corr = float(np.corrcoef(ms_ps, wp_ps)[0, 1]) if len(ms_ps) > 1 else float("nan")
    max_diff = float(np.max(np.abs(ms_ps - wp_ps)))

    if same_rejection and p_corr > 0.7 and max_diff < 0.1:
        verdict = ("STATISTICAL EQUIVALENCE: both nulls give the same"
                   " rejection pattern with similar p-values. The"
                   " matched-strength weight-distribution drift does NOT"
                   " materially affect the Grassmann verdict on this case."
                   " Joint-null implementation can be SKIPPED.")
    elif same_rejection:
        verdict = ("DIRECTIONALLY CONSISTENT but quantitatively divergent:"
                   " both nulls reject the same k cells, but p-value"
                   " magnitudes differ. Joint null would TIGHTEN the"
                   " claim but is not strictly required for the verdict.")
    else:
        verdict = ("DIVERGENT NULLS: rejection patterns differ between"
                   " the two nulls. The matched-strength drift is"
                   " plausibly contributing to the rejection. Joint null"
                   " IS needed.")

    body = f"""---
name: null-comparison-test-case
era: IMCOH_ABS_COHORT_N10
status: current
kind: diagnostic-report
date: 2026-05-28
patient: {patient}
band: {band}
phase_A: {phase_A}
phase_B: {phase_B}
R: {R}
fc_method: imcoh_abs
---

# Dual-null specific-case comparison — {patient} {band} {phase_A}↔{phase_B}

Companion to
`data/audit/matched_strength_weight_drift/` (which showed material drift
of the weight distribution under the current matched-strength null).
This script asks whether the drift affects the actual verdict statistic.

## Per-k results

| k | observed $d_G$ | $p_{{\\rm lower}}^{{\\rm MS}}$ | $p_{{\\rm lower}}^{{\\rm WP}}$ | reject MS | reject WP |
|---|---|---|---|---|---|
"""
    for _, row in df.iterrows():
        rk = int(row["k"])
        ms_p = row["p_lower_matched_strength"]
        wp_p = row["p_lower_weight_permutation"]
        body += (f"| {rk} | {row['obs_dG']:.3f} | "
                 f"{ms_p:.4f} | {wp_p:.4f} | "
                 f"{'yes' if ms_p < 0.05 else 'no'} | "
                 f"{'yes' if wp_p < 0.05 else 'no'} |\n")
    body += f"""
- `p_lower_*` is the lower-tail p-value of observed vs surrogate
  distribution (small p = observed cross-phase distance is smaller than
  surrogates' = subspace preservation = trace).
- "reject" means `p_lower < 0.05`.

## Pearson correlation of p-values across k: {p_corr:.3f}
## Max |Δp| across k: {max_diff:.3f}

## Verdict

{verdict}

## How to read the figure

`grassmann_distance_comparison.pdf` overlays, for each k, the surrogate
distribution of cross-phase Grassmann distance under matched-strength
(red) and weight-permutation (blue). The observed distance is the
vertical dashed line.

- If observed sits to the LEFT of both distributions → trace direction
  detected by both nulls. If the distributions overlap, the two nulls
  agree statistically.
- If observed sits to the LEFT of red but NOT blue → matched-strength
  rejects more strongly than weight-permutation, suggesting weight-drift
  inflation.
- If both distributions are similar and observed is at their lower
  tail → joint-null implementation would not change the verdict.

## Algorithms compared

- **matched-strength** (`strength_preserving_shuffle`): 4-cycle ±δ
  perturbations. Preserves per-node strength exactly. Weight
  distribution drifts (KS ≈ 0.36 cohort median per audit_diag).
- **weight-permutation** (`weight_permutation_shuffle`): random
  permutation of upper-triangular edge weights. Preserves weight
  multiset exactly. Per-node strengths drift freely.

A joint null preserving both would require accept-reject MCMC on weight
permutations with a strength tolerance; ~1–2 days of work. This
diagnostic decides whether that work is required.
"""
    out_path.write_text(body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--patient", default=DEFAULT_PATIENT)
    ap.add_argument("--band", default=DEFAULT_BAND)
    ap.add_argument("--phase-A", default=DEFAULT_PHASE_A)
    ap.add_argument("--phase-B", default=DEFAULT_PHASE_B)
    ap.add_argument("--R", type=int, default=DEFAULT_R)
    ap.add_argument("--swap-factor", type=int, default=DEFAULT_SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    args = ap.parse_args()

    print(f"[diag] loading observed FC: {args.patient} {args.band} "
          f"{args.phase_A} + {args.phase_B}")
    W_A = load_fc_matrix(args.patient, args.phase_A, args.band, "imcoh_abs")
    W_B = load_fc_matrix(args.patient, args.phase_B, args.band, "imcoh_abs")
    N = W_A.shape[0]
    print(f"[diag] N = {N}; n_swaps = {args.swap_factor * N * (N - 1) // 2}")

    _, U_A_obs = laplacian_eigs(W_A)
    _, U_B_obs = laplacian_eigs(W_B)
    obs_dG = {k: grassmann_chordal_distance(U_A_obs[:, :k], U_B_obs[:, :k])
              for k in K_LIST}
    print("[diag] observed cross-phase Grassmann distances:")
    for k, d in obs_dG.items():
        print(f"  k = {k:>3d}  d_G = {d:.4f}")

    print(f"\n[diag] matched-strength ensemble (R={args.R}) ...")
    U_A_ms, U_B_ms, n_ok_ms = _generate_paired_surrogates(
        W_A, W_B, args.R, "matched_strength", args.swap_factor, args.seed,
    )
    print(f"  ok: {n_ok_ms}/{args.R}")

    print(f"[diag] weight-permutation ensemble (R={args.R}) ...")
    U_A_wp, U_B_wp, n_ok_wp = _generate_paired_surrogates(
        W_A, W_B, args.R, "weight_permutation", args.swap_factor,
        args.seed + 1,
    )
    print(f"  ok: {n_ok_wp}/{args.R}")

    ms_dG = {k: [] for k in K_LIST}
    wp_dG = {k: [] for k in K_LIST}
    for U_A_s, U_B_s in zip(U_A_ms, U_B_ms):
        for k in K_LIST:
            ms_dG[k].append(
                grassmann_chordal_distance(U_A_s[:, :k], U_B_s[:, :k])
            )
    for U_A_s, U_B_s in zip(U_A_wp, U_B_wp):
        for k in K_LIST:
            wp_dG[k].append(
                grassmann_chordal_distance(U_A_s[:, :k], U_B_s[:, :k])
            )

    rows = []
    print("\n[diag] per-k lower-tail p-values:")
    print(f"{'k':>4s}  {'obs':>8s}  {'MS p':>8s}  {'WP p':>8s}")
    for k in K_LIST:
        ms_arr = np.array(ms_dG[k])
        wp_arr = np.array(wp_dG[k])
        p_ms = float((1 + (ms_arr <= obs_dG[k]).sum()) / (n_ok_ms + 1))
        p_wp = float((1 + (wp_arr <= obs_dG[k]).sum()) / (n_ok_wp + 1))
        print(f"  {k:>4d}  {obs_dG[k]:>8.4f}  {p_ms:>8.4f}  {p_wp:>8.4f}")
        rows.append({
            "k": k,
            "obs_dG": obs_dG[k],
            "p_lower_matched_strength": p_ms,
            "p_lower_weight_permutation": p_wp,
            "ms_median": float(np.median(ms_arr)),
            "wp_median": float(np.median(wp_arr)),
            "ms_p05": float(np.percentile(ms_arr, 5)),
            "wp_p05": float(np.percentile(wp_arr, 5)),
        })
    df = pd.DataFrame(rows)
    args.out_root.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_root / "per_k_pvalues.csv", index=False)

    _plot_comparison(obs_dG, ms_dG, wp_dG,
                     args.out_root / "grassmann_distance_comparison.pdf")
    _write_summary(args.patient, args.band, args.phase_A, args.phase_B,
                   args.R, obs_dG, ms_dG, wp_dG, df,
                   args.out_root / "README.md")
    print(f"\n[diag] outputs → {args.out_root}/")


if __name__ == "__main__":
    main()
