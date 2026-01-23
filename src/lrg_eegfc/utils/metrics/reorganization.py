"""Reorganization metrics computed from cached LRG results."""

from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

from lrgsglib.utils.basic.linalg import (
    ultrametric_matrix_distance,
    ultrametric_scaled_distance,
    ultrametric_rank_correlation,
    ultrametric_quantile_rmse,
    ultrametric_distance_permutation_robust,
    tree_robinson_foulds_distance,
    tree_cophenetic_correlation,
    tree_baker_gamma,
    tree_fowlkes_mallows_index,
)

from .comparison import compute_phase_distance_matrix

__all__ = [
    "build_metric_specs",
    "compute_metric_matrix",
    "compute_cluster_labels",
    "cluster_swap_coefficient",
    "compute_cluster_swap_matrix",
    "compute_ari_matrix",
]


def _ensure_square(matrix: np.ndarray) -> np.ndarray:
    if matrix.ndim == 1:
        return squareform(matrix)
    return matrix


def _safe_same_shape(matrix_a: np.ndarray, matrix_b: np.ndarray) -> bool:
    return (
        matrix_a is not None
        and matrix_b is not None
        and hasattr(matrix_a, "shape")
        and hasattr(matrix_b, "shape")
        and matrix_a.shape == matrix_b.shape
    )


def build_metric_specs(distance_metric: str = "euclidean") -> Dict[str, Dict[str, Callable]]:
    """Return metric specs keyed by short name.

    Each spec contains:
      - label: human-readable name
      - fn: callable (U1, U2, Z1, Z2) -> float
    """

    def _matrix_distance(U1: np.ndarray, U2: np.ndarray, *_):
        if not _safe_same_shape(U1, U2):
            return np.nan
        return ultrametric_matrix_distance(U1, U2, metric=distance_metric)

    def _scaled_distance(U1: np.ndarray, U2: np.ndarray, *_):
        if not _safe_same_shape(U1, U2):
            return np.nan
        return ultrametric_scaled_distance(
            U1, U2, metric=distance_metric, scale="log", normalize=True
        )

    def _rank_distance(U1: np.ndarray, U2: np.ndarray, *_):
        if not _safe_same_shape(U1, U2):
            return np.nan
        corr = ultrametric_rank_correlation(U1, U2, method="spearman")
        return 1.0 - corr if corr is not None else np.nan

    def _quantile_rmse(U1: np.ndarray, U2: np.ndarray, *_):
        if not _safe_same_shape(U1, U2):
            return np.nan
        return ultrametric_quantile_rmse(U1, U2)

    def _perm_robust(_, __, Z1: np.ndarray, Z2: np.ndarray):
        if Z1 is None or Z2 is None:
            return np.nan
        return ultrametric_distance_permutation_robust(Z1, Z2, metric=distance_metric)

    def _rf_distance(_, __, Z1: np.ndarray, Z2: np.ndarray):
        if Z1 is None or Z2 is None:
            return np.nan
        return tree_robinson_foulds_distance(Z1, Z2, normalized=True)

    def _cophenetic_corr(_, __, Z1: np.ndarray, Z2: np.ndarray):
        if Z1 is None or Z2 is None:
            return np.nan
        return tree_cophenetic_correlation(Z1, Z2)

    def _baker_gamma(_, __, Z1: np.ndarray, Z2: np.ndarray):
        if Z1 is None or Z2 is None:
            return np.nan
        return tree_baker_gamma(Z1, Z2)

    def _fowlkes_mallows(_, __, Z1: np.ndarray, Z2: np.ndarray):
        if Z1 is None or Z2 is None:
            return np.nan
        return tree_fowlkes_mallows_index(Z1, Z2)

    return {
        "matrix_distance": {
            "label": "Ultrametric matrix distance",
            "fn": _matrix_distance,
        },
        "scaled_distance": {
            "label": "Ultrametric scaled distance (log)",
            "fn": _scaled_distance,
        },
        "rank_distance": {
            "label": "Ultrametric rank distance (1 - Spearman)",
            "fn": _rank_distance,
        },
        "quantile_rmse": {
            "label": "Ultrametric quantile RMSE (log)",
            "fn": _quantile_rmse,
        },
        "permutation_robust": {
            "label": "Permutation-robust ultrametric distance",
            "fn": _perm_robust,
        },
        "tree_robinson_foulds": {
            "label": "Robinson-Foulds distance",
            "fn": _rf_distance,
        },
        "tree_cophenetic_corr": {
            "label": "Cophenetic correlation",
            "fn": _cophenetic_corr,
        },
        "tree_baker_gamma": {
            "label": "Baker gamma",
            "fn": _baker_gamma,
        },
        "tree_fowlkes_mallows": {
            "label": "Fowlkes-Mallows index",
            "fn": _fowlkes_mallows,
        },
    }


def compute_metric_matrix(
    phases: Iterable[str],
    results_by_phase: Dict[str, object],
    metric_fn: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], float],
) -> np.ndarray:
    """Compute phase-by-phase distance matrix for a single metric."""

    phases = list(phases)

    def _compute(pi: str, pj: str) -> float:
        res_i = results_by_phase.get(pi)
        res_j = results_by_phase.get(pj)
        if res_i is None or res_j is None:
            return np.nan
        U1 = _ensure_square(res_i.ultrametric_matrix)
        U2 = _ensure_square(res_j.ultrametric_matrix)
        return metric_fn(U1, U2, res_i.linkage_matrix, res_j.linkage_matrix)

    return compute_phase_distance_matrix(phases, _compute)


def compute_cluster_labels(linkage: np.ndarray, threshold: float) -> np.ndarray:
    """Compute flat cluster labels from a linkage matrix and cut threshold."""
    return fcluster(linkage, threshold, criterion="distance")


def _align_cluster_labels(labels_ref: np.ndarray, labels_target: np.ndarray) -> np.ndarray:
    """Align target labels to reference labels via maximum overlap."""
    ref = np.asarray(labels_ref)
    tgt = np.asarray(labels_target)
    ref_ids = np.unique(ref)
    tgt_ids = np.unique(tgt)
    cost = np.zeros((len(ref_ids), len(tgt_ids)), dtype=int)

    for i, rid in enumerate(ref_ids):
        for j, tid in enumerate(tgt_ids):
            cost[i, j] = np.sum((ref == rid) & (tgt == tid))

    row_ind, col_ind = linear_sum_assignment(cost.max() - cost)
    mapping = {tgt_ids[j]: ref_ids[i] for i, j in zip(row_ind, col_ind)}
    return np.array([mapping.get(label, label) for label in tgt])


def cluster_swap_coefficient(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    """Return swap coefficient (1.0 = identical, 0.0 = all nodes swap)."""
    aligned = _align_cluster_labels(labels_a, labels_b)
    return float(np.mean(np.asarray(labels_a) == aligned))


def compute_cluster_swap_matrix(
    phases: Iterable[str],
    labels_by_phase: Dict[str, np.ndarray],
) -> np.ndarray:
    phases = list(phases)
    matrix = np.zeros((len(phases), len(phases)), dtype=float)
    for i, pi in enumerate(phases):
        for j, pj in enumerate(phases):
            if i == j:
                matrix[i, j] = 1.0
            else:
                matrix[i, j] = cluster_swap_coefficient(labels_by_phase[pi], labels_by_phase[pj])
    return matrix


def compute_ari_matrix(
    phases: Iterable[str],
    labels_by_phase: Dict[str, np.ndarray],
) -> np.ndarray:
    phases = list(phases)
    matrix = np.zeros((len(phases), len(phases)), dtype=float)
    for i, pi in enumerate(phases):
        for j, pj in enumerate(phases):
            if i == j:
                matrix[i, j] = 1.0
            else:
                matrix[i, j] = float(adjusted_rand_score(labels_by_phase[pi], labels_by_phase[pj]))
    return matrix
