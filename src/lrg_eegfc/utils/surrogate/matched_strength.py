"""Matched-strength Laplacian surrogates for FC adjacency matrices.

Strength-preserving 4-cycle ±δ rewiring (Maslov-Sneppen / Squartini family
adapted for weighted graphs) plus a disk cache of full Laplacian
eigendecompositions per surrogate, so any downstream statistic (KC tree
distance, Grassmann subspace, ρ propagator, eigenmode embedding, …) on the
*same* surrogate ensemble can be computed without regenerating the
surrogates.

Algorithm (4-cycle ±δ)
----------------------
For each swap, pick four distinct nodes (a, b, c, d) and a random δ inside
the feasible interval [lo, hi] that keeps every edge weight in [0, w_max]:

    W[a, b] += δ ;   W[c, d] += δ
    W[a, d] -= δ ;   W[c, b] -= δ

Per-node strengths are exactly preserved (each affected node sees +δ on one
edge and −δ on another). The default `n_swaps = SWAP_FACTOR · N(N−1)/2`
with `SWAP_FACTOR = 20` — empirically sufficient to randomize the
adjacency given strengths in N ≈ 100–125 graphs while keeping per-surrogate
runtime ~50 ms.

Reference: this is the algorithm used by audit_62 (shared-baseline ρ),
audit_63 (split-baseline ρ_split), audit_65 (KC tree distance), and
audit_66 (Grassmann subspace).

Cache layout
------------
``data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R{R}_swap{SF}_seed{S}_imcoh_abs.npz``

Each `.npz` contains:

- ``eigvals`` — float64 array of shape ``(R, N)``; ascending eigenvalues
  of L = D − W for each of the R surrogates.
- ``eigvecs`` — float64 array of shape ``(R, N, N)``; corresponding
  orthonormal eigenvectors. The trivial zero-mode is at column 0.

Filename encodes (R, swap_factor, seed) so different ensembles do not
clobber each other. Canonical seed for the 2026-05-11 R=200 ensemble is
``20260511``.

Disk footprint at R=200, N≈120: ~24 MB per file (uncompressed savez),
~2 GB for the full cohort × 3 trace bands × 3 phases (90 files).
Uncompressed for write/read speed; compressed `savez_compressed` reduces
size by ~40% but adds ~10 s per file.

Reuse pattern
-------------
::

    from lrg_eegfc.utils.surrogate import (
        load_or_compute_surrogate_eigs,
        strength_preserving_shuffle,
    )

    rng = np.random.default_rng(seed)
    for phase in PHASES:
        eigvals, eigvecs = load_or_compute_surrogate_eigs(
            pat, band, phase, W_phase,
            n_surr=200, swap_factor=20, seed=20260511, rng=rng,
        )
        # downstream: compute KC, Grassmann, ρ, ... from eigvals / eigvecs
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from lrg_eegfc.config.paths import CACHE_ROOT


SURROGATE_LRG_CACHE = CACHE_ROOT / "matched_strength_surrogate_lrg"


def strength_preserving_shuffle(
    W: np.ndarray,
    n_swaps: int,
    rng: np.random.Generator,
    w_max: float = 1.0,
) -> np.ndarray:
    """4-cycle ±δ strength-preserving rewiring.

    Returns a copy of ``W`` with `n_swaps` independent ±δ perturbations
    applied along randomly chosen 4-cycles. Per-node strengths are
    preserved exactly (modulo float64 numerical noise; verify with
    `verify_strengths`).

    Parameters
    ----------
    W : (N, N) ndarray
        Symmetric weighted adjacency, weights in ``[0, w_max]``.
    n_swaps : int
        Number of ±δ perturbations. ~20 · N(N−1)/2 is empirically
        sufficient for our N ≈ 100–125 FC matrices.
    rng : numpy random Generator
        Source of randomness for node sampling and δ draws.
    w_max : float, default 1.0
        Upper bound on edge weights (lower bound is 0).
    """
    W = W.copy()
    N = W.shape[0]
    samples = rng.integers(0, N, size=(n_swaps, 4))
    fracs = rng.uniform(0.0, 1.0, size=n_swaps)
    for i in range(n_swaps):
        a, b, c, d = samples[i]
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]; w2 = W[c, d]; w3 = W[a, d]; w4 = W[c, b]
        lo = -w1 if -w1 > -w2 else -w2
        if w3 - w_max > lo:
            lo = w3 - w_max
        if w4 - w_max > lo:
            lo = w4 - w_max
        hi = w_max - w1 if w_max - w1 < w_max - w2 else w_max - w2
        if w3 < hi:
            hi = w3
        if w4 < hi:
            hi = w4
        if lo >= hi:
            continue
        delta = lo + fracs[i] * (hi - lo)
        n1 = w1 + delta; n2 = w2 + delta
        n3 = w3 - delta; n4 = w4 - delta
        W[a, b] = n1; W[b, a] = n1
        W[c, d] = n2; W[d, c] = n2
        W[a, d] = n3; W[d, a] = n3
        W[c, b] = n4; W[b, c] = n4
    return W


def verify_strengths(
    W_obs: np.ndarray, W_surr: np.ndarray, tol: float = 1e-4
) -> bool:
    """True if max per-node strength delta < `tol`."""
    s_obs = W_obs.sum(axis=1)
    s_surr = W_surr.sum(axis=1)
    return bool(np.max(np.abs(s_obs - s_surr)) < tol)


def surrogate_cache_path(
    pat: str, band: str, phase: str,
    n_surr: int, swap_factor: int, seed: int,
    fc_method: str = "imcoh_abs",
) -> Path:
    """Canonical cache filename for one (patient, band, phase) ensemble."""
    return (SURROGATE_LRG_CACHE / pat
            / f"{band}_{phase}_R{n_surr}_swap{swap_factor}"
              f"_seed{seed}_{fc_method}.npz")


def _laplacian_eig(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return np.linalg.eigh(L)


def load_or_compute_surrogate_eigs(
    pat: str, band: str, phase: str,
    W: np.ndarray,
    n_surr: int, swap_factor: int, seed: int,
    rng: np.random.Generator,
    fc_method: str = "imcoh_abs",
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(eigvals[R, N], eigvecs[R, N, N])`` for the surrogate ensemble.

    If the canonical cache file exists with matching shape, load it.
    Otherwise generate `n_surr` matched-strength surrogates of `W` via
    `strength_preserving_shuffle` (with `n_swaps = swap_factor · N(N−1)/2`),
    eigendecompose each, save to cache, and return.

    Surrogates that fail the strength-preservation tolerance check
    (`verify_strengths` with tol=1e-4) get NaN entries — downstream code
    should mask these.

    Notes
    -----
    The seed in the filename pins the rng state implicitly: re-running with
    the same seed but consuming the rng differently (e.g. iterating bands
    in a different order) will produce a different surrogate ensemble.
    Cache hits therefore require both the seed and the per-cell rng-draw
    history to match — in practice, run the script with the same seed and
    iteration order to reuse the cache.
    """
    path = surrogate_cache_path(pat, band, phase, n_surr, swap_factor,
                                  seed, fc_method)
    N = W.shape[0]
    if path.exists():
        with np.load(path) as data:
            evals = data["eigvals"]
            evecs = data["eigvecs"]
        if (evals.shape == (n_surr, N)
                and evecs.shape == (n_surr, N, N)):
            return evals.astype(np.float64), evecs.astype(np.float64)
        if verbose:
            print(f"[matched_strength] cache shape mismatch at {path}; recomputing")

    n_swaps = swap_factor * (N * (N - 1)) // 2
    evals = np.empty((n_surr, N), dtype=np.float64)
    evecs = np.empty((n_surr, N, N), dtype=np.float64)
    ok_count = 0
    for r in range(n_surr):
        W_s = strength_preserving_shuffle(W, n_swaps, rng)
        if not verify_strengths(W, W_s, tol=1e-4):
            evals[r, :] = np.nan
            evecs[r, :, :] = np.nan
            continue
        ev, ec = _laplacian_eig(W_s)
        evals[r, :] = ev
        evecs[r, :, :] = ec
        ok_count += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, eigvals=evals, eigvecs=evecs)
    if verbose:
        print(f"[matched_strength] cached {ok_count}/{n_surr} surrogates → {path.name}")
    return evals, evecs
