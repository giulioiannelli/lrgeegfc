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
    "chordal_from_angles",
    "grassmann_to_coord_subspace",
    "chordal_full_vs_resect",
]


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
