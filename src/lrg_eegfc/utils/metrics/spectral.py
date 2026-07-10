"""Spectral / Grassmann subspace primitives on the LRG eigenstructure.

Promoted to the library on 2026-05-08 to back

- ``audit_37_e1_grassmann_triangle.py`` (private ``chordal``);
- ``audit_46_grassmann_principal_angles.py`` (private ``principal_angles``,
  ``chordal_from_angles``);
- ``audit_61_epi_grassmann_compute.py`` (Direction E — epi Grassmann
  embedding scope at
  ``.agents/guides/task-persistence-investigation/2026-05-08_epi-grassmann-embedding.md``).

References
----------
Björck & Golub (1973) — numerical principal-angle algorithm via SVD.
Edelman, Arias & Smith (1998) — chordal vs geodesic Grassmann metrics.
"""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

__all__ = [
    "principal_angles",
    "chordal_distance",
    "chordal_distances_for_k_grid",
    "chordal_from_angles",
    "grassmann_to_coord_subspace",
    "chordal_full_vs_resect",
    "participation_number",
    "node_set_mode_mass",
    "subspace_displacement",
    "laplacian_eig",
    "topk_basis",
]


def laplacian_eig(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Combinatorial-Laplacian eigendecomposition ``L = diag(W·1) − W``.

    Villegas-canonical combinatorial Laplacian (never ``L_rw`` / ``L_sym``;
    see ``feedback_combinatorial_laplacian_only``). Returns ``(eigvals,
    eigvecs)`` from ``numpy.linalg.eigh`` — eigenvalues ascending, so the
    trivial constant zero-mode is column 0 of ``eigvecs`` (drop it with
    :func:`topk_basis`). ``W`` must be symmetric, non-negative, zero-diagonal.

    Promoted 2026-07-09 (second caller: ``audit_165_grassmann_inference_arc``)
    from the private copies in ``audit_66_grassmann_matched_strength_surrogate``
    and ``matched_strength._laplacian_eig``.
    """
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return np.linalg.eigh(L)


def topk_basis(eigvecs: np.ndarray, k_max: int) -> np.ndarray:
    """Slowest ``k_max`` non-trivial eigenmodes: ``eigvecs[:, 1:k_max+1]``.

    Drops the trivial zero-mode at column 0 (``eigh`` ascending order) and
    returns a contiguous ``(N, k_max)`` orthonormal block — the leading
    invariant subspace used by the chordal Grassmann distance.

    Promoted 2026-07-09 (second caller: ``audit_165_grassmann_inference_arc``)
    from the private copy in ``audit_66_grassmann_matched_strength_surrogate``.
    """
    return np.ascontiguousarray(eigvecs[:, 1:k_max + 1])


def participation_number(V: np.ndarray) -> np.ndarray:
    """Per-mode participation number (inverse participation ratio reciprocal).

    For each column ``v_k`` of an eigenvector matrix ``V`` (``N × M``, columns
    unit-norm), returns

        PR_k = (Σ_i v_k(i)²)² / Σ_i v_k(i)⁴

    the standard solid-state *participation number* (Edwards & Thouless 1972).
    Range ``[1, N]``: ``PR_k = 1`` ⇔ all mass on one node (maximally
    localized); ``PR_k = N`` ⇔ uniform spread (maximally extended). For
    unit-norm columns ``Σ v_k(i)² = 1`` so ``PR_k = 1 / Σ_i v_k(i)⁴``; the
    explicit numerator keeps the definition correct for non-normalized inputs.

    The raw inverse participation ratio ``IPR_k = Σ_i v_k(i)⁴`` is
    ``1 / PR_k`` for unit-norm columns; ``PR`` is reported because it is in
    interpretable node units (small = localized) and is comparable across
    patients with different ``N`` only after normalization by ``N`` if needed.

    Parameters
    ----------
    V : (N, M) ndarray
        Eigenvector matrix; each column an eigenmode. The trivial constant
        zero-mode (if present) should be dropped by the caller before
        interpreting localization.

    Returns
    -------
    (M,) ndarray
        ``PR_k`` per column; ``nan`` where a column has zero ℓ⁴ norm.
    """
    if V.ndim != 2:
        raise ValueError(f"V must be 2-D (N, M); got shape {V.shape}")
    v2 = V * V
    s2 = v2.sum(axis=0)
    s4 = (v2 * v2).sum(axis=0)
    return np.where(s4 > 0, (s2 * s2) / s4, np.nan)


def node_set_mode_mass(V: np.ndarray, node_idx: np.ndarray) -> np.ndarray:
    """Per-mode squared mass on a node subset: ``m^S_k = Σ_{i∈S} v_k(i)²``.

    For unit-norm columns ``m^S_k ∈ [0, 1]`` is the fraction of mode ``k``'s
    squared amplitude on the nodes ``node_idx``. Sign-invariant by
    construction (squared amplitude); invariant to within-degenerate-block
    eigenvector rotation only when summed over the whole degenerate block.

    Parameters
    ----------
    V : (N, M) ndarray
        Eigenvector matrix (columns = modes).
    node_idx : array-like of int or bool
        Row indices (or boolean row mask) of the node subset ``S``.

    Returns
    -------
    (M,) ndarray
        ``m^S_k`` per mode.
    """
    if V.ndim != 2:
        raise ValueError(f"V must be 2-D (N, M); got shape {V.shape}")
    idx = np.asarray(node_idx)
    block = V[idx, :]                       # |S| × M
    return (block * block).sum(axis=0)


def subspace_displacement(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    """Per-mode cross-phase displacement of each column of ``V_b`` from the
    span of ``V_a`` (both ``N × M``, orthonormal columns in a common ``R^N``).

    For each mode ``k`` (column of ``V_b``):

        displacement_k = 1 − max_j |⟨v_b^k, v_a^j⟩|

    i.e. one minus the best single-mode overlap with the other phase's
    eigenbasis. Range ``[0, 1]``: ``0`` ⇔ the mode is unchanged (a perfect
    match exists in ``V_a``); large ⇔ the mode reorganized across phases.
    Sign-invariant (absolute overlap). This is a per-mode reduction of the
    cross-phase rotation that the leading-subspace chordal Grassmann distance
    aggregates over a fixed top-k block; here it is resolved mode-by-mode so
    it can be correlated against per-mode localization.

    Parameters
    ----------
    V_a, V_b : (N, M) ndarray
        Eigenvector matrices for the two phases (e.g. rest_pre and rest_post),
        nontrivial modes only, in a common ambient ``R^N``.

    Returns
    -------
    (M,) ndarray
        Per-mode displacement of ``V_b`` columns.
    """
    if V_a.ndim != 2 or V_b.ndim != 2:
        raise ValueError("V_a and V_b must be 2-D (N, M)")
    if V_a.shape[0] != V_b.shape[0]:
        raise ValueError(
            f"ambient mismatch: V_a rows {V_a.shape[0]} != V_b rows {V_b.shape[0]}"
        )
    overlap = np.abs(V_b.T @ V_a)           # M_b × M_a
    return 1.0 - overlap.max(axis=1)


def principal_angles(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    """Principal angles (radians) between the column spans of ``V_a`` and ``V_b``.

    Both inputs must live in a common ambient space ``R^D`` and have
    orthonormal columns (the Laplacian-eigenvector blocks satisfy this by
    construction). Internally: ``sigma = svd(V_a.T @ V_b)``,
    ``theta = arccos(clip(sigma, 0, 1))``.

    Returns
    -------
    np.ndarray
        ``min(p, q)``-long array of angles in ``[0, pi/2]`` ordered
        small → large (i.e. ``sigma`` ordered large → small).
    """
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return np.arccos(sigma)


def chordal_distance(V_a: np.ndarray, V_b: np.ndarray) -> float:
    """Chordal Grassmann distance ``sqrt(min(p, q) − Σ σ_i²)``.

    Inputs ``V_a`` (``D × p``) and ``V_b`` (``D × q``) must have orthonormal
    columns in a common ambient ``R^D``. Range: ``[0, sqrt(min(p, q))]``;
    ``0`` means one subspace contains the other; ``sqrt(min(p, q))`` means
    they are mutually orthogonal.
    """
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    pq = min(V_a.shape[1], V_b.shape[1])
    return float(np.sqrt(max(pq - float((sigma ** 2).sum()), 0.0)))


def chordal_distances_for_k_grid(
    eigvecs_a: np.ndarray, eigvecs_b: np.ndarray, k_grid: Sequence[int]
) -> np.ndarray:
    """Chordal Grassmann distance ``d(k)`` for every ``k`` in ``k_grid`` at once.

    Both inputs are full ``(N, N)`` eigenvector matrices in ``eigh`` ascending
    order (column 0 = trivial zero-mode). The leading-``k`` non-trivial subspace
    is ``eigvecs[:, 1:k+1]`` and

        d(k) = sqrt( k − ‖ U_a^{(k)ᵀ} U_b^{(k)} ‖_F² ),   domain ``[0, sqrt(k)]``.

    Computed once via the cumulative-Frobenius trick: form the
    ``(N−1)×(N−1)`` cross-Gram ``M = A^T B`` (``A = eigvecs_a[:, 1:]``), then
    ``‖U_a^{(k)ᵀ} U_b^{(k)}‖_F²`` is the 2-D cumulative sum of ``M∘M`` read at
    ``[k-1, k-1]`` — so all ``len(k_grid)`` cutoffs come from a single matmul.
    Mathematically identical to :func:`chordal_distance` on ``k``-sliced blocks
    (``‖M‖_F² = Σ σ_i²``); ``nan`` where ``k`` exceeds ``N−1``.

    Promoted 2026-07-09 (second caller: ``audit_165_grassmann_inference_arc``)
    from the private ``chordal_distances_for_k_grid`` in
    ``audit_70_grassmann_cluster_extent``.
    """
    A = eigvecs_a[:, 1:]  # (N, N-1)
    B = eigvecs_b[:, 1:]  # (N, N-1)
    M = A.T @ B           # (N-1, N-1)
    M2 = M * M
    csum = np.cumsum(np.cumsum(M2, axis=0), axis=1)
    out = np.full(len(k_grid), np.nan)
    n_max = csum.shape[0]
    for i, k in enumerate(k_grid):
        if k > n_max:
            continue
        frob2 = csum[k - 1, k - 1]
        out[i] = np.sqrt(max(0.0, k - frob2))
    return out


def chordal_from_angles(theta: np.ndarray) -> float:
    """Chordal distance recomposed from principal angles: ``sqrt(Σ sin²(θ_i))``."""
    return float(np.sqrt(float((np.sin(theta) ** 2).sum())))


def grassmann_to_coord_subspace(
    V_k: np.ndarray, idx: Sequence[int]
) -> Tuple[float, float]:
    """Chordal distance from ``span(V_k)`` to the coordinate subspace
    ``S = span{e_i : i ∈ idx}`` in the ambient ``R^N``.

    Mathematically equivalent to checking how much of the top-k eigenspace
    sits on the rows ``idx``. Let ``B = V_k[idx, :] ∈ R^{|idx| × k}``. The
    cross-Gram of ``V_k`` with the canonical basis of ``S`` is exactly
    ``B^T``, so ``Σ σ_i² = ‖B‖_F²``. Returns

    - ``eps_align`` = ``sqrt(min(k, |idx|) − ‖B‖_F²)``
        — chordal Grassmann distance, range ``[0, sqrt(min(k, |idx|))]``.
    - ``f_idx`` = ``‖B‖_F² / k``
        — average fraction of top-k mode mass on rows ``idx``,
        range ``[0, 1]``.

    Both numbers are reparametrisations of the same scalar; ``eps_align`` is
    the Grassmann form, ``f_idx`` is the interpretable mass-fraction form.
    """
    if V_k.ndim != 2:
        raise ValueError(f"V_k must be 2-D; got shape {V_k.shape}")
    idx_arr = np.asarray(idx, dtype=int)
    if idx_arr.size == 0:
        raise ValueError("idx must be non-empty")
    k = V_k.shape[1]
    block = V_k[idx_arr, :]
    fro2 = float(np.sum(block * block))
    m = min(k, idx_arr.size)
    eps_align = float(np.sqrt(max(m - fro2, 0.0)))
    f_idx = fro2 / k
    return eps_align, f_idx


def chordal_full_vs_resect(
    V_k_full: np.ndarray,
    V_k_resect: np.ndarray,
    retained_idx: Sequence[int],
    rank_tol: float = 1e-10,
) -> Tuple[float, int]:
    """Chordal distance between the full-graph eigenspace restricted to
    ``retained_idx`` rows (orthonormalised via QR) and the eigenspace
    computed from scratch on the resected graph.

    Parameters
    ----------
    V_k_full : ``N × k``
        Top-k non-trivial eigenvectors of the full FC graph's Laplacian.
        Orthonormal columns, ambient ``R^N``.
    V_k_resect : ``M × k``
        Top-k non-trivial eigenvectors of the resected graph's Laplacian.
        Orthonormal columns, ambient ``R^M`` where ``M = |retained_idx|``.
    retained_idx : sequence of int
        Row indices into ``V_k_full`` corresponding to nodes kept in the
        resected graph.
    rank_tol : float
        QR pivot tolerance for detecting rank-drop in the restricted block.

    Returns
    -------
    (delta, rank_block) : (float, int)
        - ``delta`` = ``d_chord(Q_full, V_k_resect)`` after orthonormalising
          the restricted block via QR. Range ``[0, sqrt(min(rank_block, k))]``.
          Small ⇒ resected modes match the non-epi rows of the full modes;
          large ⇒ resection induces a different mode set.
        - ``rank_block`` = rank of ``V_k_full[retained_idx, :]`` (≤ k).
          A rank-drop is itself a finding (the full top-k eigenspace was
          not entirely supported on the retained rows).
    """
    retained = np.asarray(retained_idx, dtype=int)
    if V_k_full.ndim != 2 or V_k_resect.ndim != 2:
        raise ValueError("Eigenvector blocks must be 2-D")
    if V_k_resect.shape[0] != retained.size:
        raise ValueError(
            f"V_k_resect rows ({V_k_resect.shape[0]}) must equal "
            f"|retained_idx| ({retained.size})"
        )
    if V_k_full.shape[1] != V_k_resect.shape[1]:
        raise ValueError(
            f"k mismatch: V_k_full has {V_k_full.shape[1]} cols, "
            f"V_k_resect has {V_k_resect.shape[1]}"
        )

    block = V_k_full[retained, :]                       # M × k
    Q, R_qr = np.linalg.qr(block, mode="reduced")       # Q: M × k, R: k × k
    diag_abs = np.abs(np.diag(R_qr))
    rank_block = int((diag_abs > rank_tol * (diag_abs.max() if diag_abs.size else 1.0)).sum())
    if rank_block < block.shape[1]:
        # truncate Q to the well-conditioned columns
        Q = Q[:, :rank_block]
    delta = chordal_distance(Q, V_k_resect)
    return delta, rank_block
