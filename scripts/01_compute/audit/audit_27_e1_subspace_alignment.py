#!/usr/bin/env python3
"""Audit 27 — E1 spectral subspace alignment (eigenvector-direct pivot, rung 1).

Implements the scope at
``.agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md``.

For each (patient, band, k), computes the Grassmann chordal distance
between the top-k non-trivial eigenspaces of L̂^pre, L̂^test, L̂^post,
and the analogous within-rest null from the halves cache::

    Δ^E1_k       = d_chord(V^test_k, V^post_k) − d_chord(V^pre_k, V^post_k)
    Δ^E1_null,k  = d_chord(V^rpreB_k, V^rpostA_k) − d_chord(V^rpreA_k, V^rpostA_k)

Per-cell positivity: ``Δ^E1_k < 0  AND  |Δ^E1_k| > |Δ^E1_null,k|``.

Usage
-----
    python scripts/01_compute/audit/audit_27_e1_subspace_alignment.py [--sanity-only] [-v]

Outputs
-------
    data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv
    data/audit/spectral_subspace/e1_sbm_sanity.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_LIST
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
OUT_DIR = ROOT / "data" / "audit" / "spectral_subspace"
OUT_CSV = OUT_DIR / "e1_subspace_alignment_n10_imcoh_abs.csv"
OUT_SANITY = OUT_DIR / "e1_sbm_sanity.json"

K_GRID = [2, 3, 5, 8, 13, 21]


# ─────────────────────────── core math ───────────────────────────


def chordal_distance(V_a: np.ndarray, V_b: np.ndarray) -> float:
    """Grassmann chordal distance between two N×k orthonormal bases.

    `d_chord(V_a, V_b) = sqrt(k - sum(σ_i^2))` where σ_i are singular
    values of `V_a^T @ V_b`. Range `[0, sqrt(k)]`. Subspace-invariant
    (independent of basis sign / orientation).
    """
    if V_a.shape != V_b.shape:
        raise ValueError(f"shape mismatch: {V_a.shape} vs {V_b.shape}")
    k = V_a.shape[1]
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    inner = float(np.sum(sigma ** 2))
    # Numerical floor: σ_i ≤ 1 + ε for orthonormal V; clip to k.
    inner = min(inner, float(k))
    return float(np.sqrt(max(0.0, k - inner)))


def topk_basis(eigenvectors: np.ndarray, k: int) -> np.ndarray:
    """Return columns 1..k of an ascending-eigenvalue eigenvector matrix.

    Column 0 is the trivial mode (λ_1 = 0); always dropped. The next k
    columns are v_2…v_{k+1}.
    """
    if eigenvectors is None:
        raise ValueError("eigenvectors is None — refresh LRG cache (Phase 0 step 0a)")
    if k + 1 > eigenvectors.shape[1]:
        raise ValueError(f"k={k} too large for N={eigenvectors.shape[1]}")
    return np.ascontiguousarray(eigenvectors[:, 1 : k + 1])


# ─────────────────────────── SBM sanity gate ───────────────────────────


def _sbm_realisation(rng: np.random.Generator, n_blocks: int = 4,
                     block_size: int = 25, w_intra: float = 0.8,
                     w_inter: float = 0.1, sigma: float = 0.05,
                     permute: bool = False) -> np.ndarray:
    N = n_blocks * block_size
    A = np.full((N, N), w_inter)
    for b in range(n_blocks):
        s = slice(b * block_size, (b + 1) * block_size)
        A[s, s] = w_intra
    A = A + rng.normal(0.0, sigma, size=(N, N))
    A = 0.5 * (A + A.T)            # symmetric
    np.fill_diagonal(A, 0.0)
    A = np.maximum(A, 0.0)         # non-negative for Laplacian
    if permute:
        order = rng.permutation(N)
        A = A[order][:, order]
    return A


def _laplacian_eigvecs(A: np.ndarray) -> np.ndarray:
    deg = A.sum(axis=1)
    L = np.diag(deg) - A
    eigvals, eigvecs = np.linalg.eigh(L)
    # Sanity: eigvals ascending, eigvals[0] ≈ 0 (constant mode)
    return eigvecs


def run_sbm_sanity(seed: int = 0, n_blocks: int = 4) -> dict:
    """Three-way SBM gate for the chordal_distance implementation.

    Tests at ``k = n_blocks - 1`` (the natural number of block-indicator
    modes after dropping the trivial constant mode v_1). Going higher
    (k ≥ n_blocks) crosses the spectral gap into intra-block noise modes
    where two independent SBM realisations have genuinely different
    subspaces — a correct algorithm WILL show large d_chord there.

    Pass criteria (all required):
        within-realisation chordal at k=n_blocks-1 < 0.3
        SBM vs shuffled-block chordal at k=n_blocks-1 > 0.8
        reflexive d_chord(V, V) < 1e-5  (sqrt of float64 ε at k≈few)
    """
    k = n_blocks - 1
    rng = np.random.default_rng(seed)
    A1 = _sbm_realisation(rng, n_blocks=n_blocks, permute=False)
    A2 = _sbm_realisation(rng, n_blocks=n_blocks, permute=False)
    A_shuf = _sbm_realisation(rng, n_blocks=n_blocks, permute=True)

    V1 = _laplacian_eigvecs(A1)
    V2 = _laplacian_eigvecs(A2)
    Vs = _laplacian_eigvecs(A_shuf)

    Vk1 = topk_basis(V1, k)
    Vk2 = topk_basis(V2, k)
    Vks = topk_basis(Vs, k)

    d_within = chordal_distance(Vk1, Vk2)
    d_shuffle = chordal_distance(Vk1, Vks)
    d_reflex = chordal_distance(Vk1, Vk1)

    pass_within = d_within < 0.3
    pass_shuffle = d_shuffle > 0.8
    pass_reflex = d_reflex < 1e-5
    return {
        "k": k, "n_blocks": n_blocks, "seed": seed,
        "d_within_realisation": d_within,
        "d_shuffled_block":     d_shuffle,
        "d_reflexive":          d_reflex,
        "pass_within":  bool(pass_within),
        "pass_shuffle": bool(pass_shuffle),
        "pass_reflex":  bool(pass_reflex),
        "pass_all":     bool(pass_within and pass_shuffle and pass_reflex),
    }


# ─────────────────────────── per-patient cohort run ───────────────────────────


def _load_eigvecs_full(pat: str, band: str, phase: str) -> np.ndarray | None:
    r = load_lrg_result(pat, phase, band, "imcoh_abs")
    if r is None or r.eigenvectors is None:
        return None
    return r.eigenvectors


def _load_eigvecs_half(pat: str, band: str, phase_tag: str) -> np.ndarray | None:
    r = load_lrg_result(pat, phase_tag, band, "imcoh_abs",
                        cache_root=HALVES_CACHE)
    if r is None or r.eigenvectors is None:
        return None
    return r.eigenvectors


def cohort_run(verbose: bool = False) -> pd.DataFrame:
    rows: list[dict] = []
    for pat in PATIENTS_LIST:
        if verbose:
            print(f"=== {pat} ===")
        for band in BRAIN_BANDS_NAMES:
            V_pre  = _load_eigvecs_full(pat, band, "rest_pre")
            V_test = _load_eigvecs_full(pat, band, "task_test")
            V_post = _load_eigvecs_full(pat, band, "rest_post")
            V_rpreA  = _load_eigvecs_half(pat, band, "rest_pre_A")
            V_rpreB  = _load_eigvecs_half(pat, band, "rest_pre_B")
            V_rpostA = _load_eigvecs_half(pat, band, "rest_post_A")

            full_ok = all(v is not None for v in (V_pre, V_test, V_post))
            null_ok = all(v is not None for v in (V_rpreA, V_rpreB, V_rpostA))
            # Patients can have different N across bands only if the
            # giant-component differs; chordal needs N_a == N_b. We
            # require all three full phases (and all three halves) share
            # N within (patient, band) — which holds for n=10 imcoh_abs
            # by construction (giant component is full, N_pat is fixed).
            if full_ok:
                Ns = {V_pre.shape[0], V_test.shape[0], V_post.shape[0]}
                if len(Ns) != 1:
                    full_ok = False
            if null_ok:
                Ns = {V_rpreA.shape[0], V_rpreB.shape[0], V_rpostA.shape[0]}
                if len(Ns) != 1:
                    null_ok = False

            for k in K_GRID:
                row: dict = dict(patient=pat, band=band, k=k)
                if full_ok and (k + 1 <= V_pre.shape[1]):
                    Vp = topk_basis(V_pre, k)
                    Vt = topk_basis(V_test, k)
                    Vq = topk_basis(V_post, k)
                    d_pre_post  = chordal_distance(Vp, Vq)
                    d_test_post = chordal_distance(Vt, Vq)
                    delta = d_test_post - d_pre_post
                    row.update(d_pre_post=d_pre_post,
                               d_test_post=d_test_post,
                               delta_e1=delta)
                else:
                    row.update(d_pre_post=np.nan,
                               d_test_post=np.nan,
                               delta_e1=np.nan)
                if null_ok and (k + 1 <= V_rpreA.shape[1]):
                    Va = topk_basis(V_rpreA, k)
                    Vb = topk_basis(V_rpreB, k)
                    Vc = topk_basis(V_rpostA, k)
                    d_A_C = chordal_distance(Va, Vc)
                    d_B_C = chordal_distance(Vb, Vc)
                    delta_n = d_B_C - d_A_C
                    row.update(d_null_rpreA_rpostA=d_A_C,
                               d_null_rpreB_rpostA=d_B_C,
                               delta_e1_null=delta_n)
                else:
                    row.update(d_null_rpreA_rpostA=np.nan,
                               d_null_rpreB_rpostA=np.nan,
                               delta_e1_null=np.nan)

                de = row["delta_e1"]
                dn = row["delta_e1_null"]
                if np.isnan(de) or np.isnan(dn):
                    row.update(passes_null=False, sign_negative=False,
                               positive=False)
                else:
                    sign_neg = de < 0.0
                    passes = abs(de) > abs(dn)
                    row.update(passes_null=passes,
                               sign_negative=bool(sign_neg),
                               positive=bool(sign_neg and passes))
                rows.append(row)
    df = pd.DataFrame(rows)
    return df


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sanity-only", action="store_true",
                   help="Run SBM gate only, skip cohort.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("== SBM sanity gate ==")
    sanity = run_sbm_sanity(seed=0, n_blocks=4)
    OUT_SANITY.write_text(json.dumps(sanity, indent=2))
    for key, val in sanity.items():
        print(f"  {key:24s} = {val}")
    if not sanity["pass_all"]:
        print("\nSBM sanity gate FAILED — fix implementation before cohort run.")
        return 2

    if args.sanity_only:
        return 0

    print("\n== Cohort run ==")
    df = cohort_run(verbose=args.verbose)
    df.to_csv(OUT_CSV, index=False)
    n = len(df)
    n_pos = int(df["positive"].sum())
    print(f"\nWrote {OUT_CSV}")
    print(f"  {n_pos}/{n} cells positive")
    return 0


if __name__ == "__main__":
    sys.exit(main())
