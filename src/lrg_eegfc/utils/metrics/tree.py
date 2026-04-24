"""Tree-level utilities for LRG dendrograms.

Operates on scipy ``linkage`` matrices ``Z`` (shape ``(n_leaves-1, 4)``) with
the internal-node representation from ``scipy.cluster.hierarchy.to_tree``.

Used by the cluster-birth-retention test (`h2_cluster_birth_retention`) and
the h-parametrized multiscale landscape (`h2_partition_multiscale_h`). See
`.agents/guides/02_methods/H2_METRICS.md` for the scientific role.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
from scipy.cluster.hierarchy import fcluster, to_tree


__all__ = [
    "dmax_from_Z",
    "tree_internal_nodes",
    "jaccard_leafsets",
    "h_log_grid",
    "fcluster_at_h_rel",
]


def dmax_from_Z(Z: np.ndarray) -> float:
    """Top merge height of a scipy linkage matrix — normalization reference."""
    return float(Z[-1, 2])


def tree_internal_nodes(Z: np.ndarray) -> list[dict]:
    """Enumerate every internal (non-leaf) node of ``Z``.

    Returns one dict per internal node with keys ``h`` (merge height, float),
    ``size`` (number of leaves in the subtree), ``leaves`` (frozenset of leaf
    ids), and ``h_rel`` (``h / dmax``).

    There are ``n_leaves - 1`` internal nodes.
    """
    dmax = dmax_from_Z(Z)
    root = to_tree(Z, rd=False)
    out: list[dict] = []

    def walk(node) -> frozenset[int]:
        if node.is_leaf():
            return frozenset([node.id])
        left = walk(node.get_left())
        right = walk(node.get_right())
        leaves = left | right
        out.append({
            "h": float(node.dist),
            "h_rel": float(node.dist) / dmax if dmax > 0 else 0.0,
            "size": len(leaves),
            "leaves": leaves,
        })
        return leaves

    walk(root)
    return out


def jaccard_leafsets(A: frozenset[int], B: frozenset[int]) -> float:
    """Symmetric Jaccard |A∩B| / |A∪B|. Returns 0 if both sets are empty."""
    if not A and not B:
        return 0.0
    inter = len(A & B)
    union = len(A | B)
    return inter / union if union > 0 else 0.0


def h_log_grid(dmaxes: Iterable[float], n: int = 60,
               h_min_floor: float = 1e-3) -> np.ndarray:
    """Geometric grid of relative heights ``h_rel ∈ [h_min, 1]``.

    ``h_min = max(h_min_floor, min(dmax) / max(dmax))`` — clamps the grid so
    the smallest scale is representable across all trees fed in. ``dmaxes``
    is only used to derive ``h_min``; the returned grid is in units of
    ``h_rel`` (fractional merge height), not absolute distance.
    """
    dm = np.array(list(dmaxes), dtype=float)
    if dm.size == 0 or not np.any(dm > 0):
        return np.geomspace(h_min_floor, 1.0, num=n)
    dm_pos = dm[dm > 0]
    data_floor = float(dm_pos.min() / dm_pos.max())
    h_min = max(h_min_floor, data_floor)
    if h_min >= 1.0:
        h_min = h_min_floor
    return np.geomspace(h_min, 1.0, num=n)


def fcluster_at_h_rel(Z: np.ndarray, h_rel: float) -> np.ndarray:
    """Flat partition at fractional cophenetic-distance threshold.

    ``h_rel ∈ (0, 1]`` is multiplied by ``dmax(Z)`` and fed to scipy's
    ``fcluster(..., criterion='distance')``. Equivalent to "cut this
    dendrogram at ``h_rel``·(its own tree height)".
    """
    dmax = dmax_from_Z(Z)
    return fcluster(Z, t=h_rel * dmax, criterion="distance")
