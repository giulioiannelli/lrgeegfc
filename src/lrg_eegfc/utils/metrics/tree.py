"""Tree-level utilities for LRG dendrograms.

Operates on scipy ``linkage`` matrices ``Z`` (shape ``(n_leaves-1, 4)``) with
the internal-node representation from ``scipy.cluster.hierarchy.to_tree``.

Used by the cluster-birth-retention test (`h2_cluster_birth_retention`) and
the h-parametrized multiscale landscape (`h2_partition_multiscale_h`). See
`.agents/guides/02_methods/H2_METRICS.md` for the scientific role.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
from scipy.cluster.hierarchy import cophenet, fcluster, linkage, to_tree
from scipy.spatial.distance import squareform


__all__ = [
    "dmax_from_Z",
    "tree_internal_nodes",
    "jaccard_leafsets",
    "h_log_grid",
    "fcluster_at_h_rel",
    "simpson_neff",
    "cluster_size_stats",
    "partition_vi_on_subset",
    "induced_linkage",
]


def induced_linkage(Z: np.ndarray, leaf_indices: Sequence[int],
                    method: str = "average") -> np.ndarray:
    """Linkage matrix of the dendrogram induced on a leaf subset.

    Given a scipy linkage matrix ``Z`` over leaves ``0..n-1`` and a subset
    ``leaf_indices ⊆ {0..n-1}`` of size ``k ≥ 2``, returns the linkage
    matrix of the induced subtree on those leaves — i.e. the dendrogram
    restricted to ``leaf_indices`` after pruning non-included leaves and
    contracting degree-2 internal nodes.

    Implementation: restricts the cophenetic distance matrix of ``Z`` to
    the chosen leaves and re-runs ``scipy.cluster.hierarchy.linkage`` on
    the condensed sub-matrix. For UPGMA (``method='average'``), the
    cophenetic-restriction approach is exact: every MRCA height in the
    induced tree equals the original MRCA height, because cophenetic
    distance ``d(i,j)`` IS the merge height of the MRCA of ``(i, j)``.

    Parameters
    ----------
    Z : np.ndarray, shape (n-1, 4)
        Linkage matrix over leaves ``0..n-1``.
    leaf_indices : sequence of int
        Leaves to keep, indices into ``0..n-1``. Order is preserved
        in the output: leaf ``i`` of the induced linkage corresponds
        to ``leaf_indices[i]`` of the original tree.
    method : str, default ``"average"``
        Linkage method. Must match the method used to build ``Z`` for
        exact MRCA-height preservation.

    Returns
    -------
    Z_sub : np.ndarray, shape (k-1, 4)
        Linkage matrix of the induced subtree.
    """
    leaves = np.asarray(leaf_indices, dtype=int)
    if leaves.size < 2:
        raise ValueError(f"induced_linkage requires |leaf_indices| ≥ 2; got {leaves.size}")
    coph = squareform(cophenet(Z))
    sub = coph[np.ix_(leaves, leaves)]
    sub_condensed = squareform(sub, checks=False)
    return linkage(sub_condensed, method=method)


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


def simpson_neff(labels: np.ndarray) -> float:
    """Simpson effective number of clusters: ``1 / Σᵢ pᵢ²``.

    For a flat partition (1-D label vector), ``pᵢ`` is the fraction of
    elements in cluster ``i``. ``n_eff = N`` when every element is alone;
    ``n_eff = 1`` when one cluster swallows everything. Used as a k-axis
    diagnostic — collapses near k=N (singleton domination) or near k=O(1)
    (giant domination).
    """
    labels = np.asarray(labels)
    if labels.size == 0:
        return 0.0
    _, counts = np.unique(labels, return_counts=True)
    p = counts / counts.sum()
    s = float((p * p).sum())
    return 1.0 / s if s > 0 else 0.0


def partition_vi_on_subset(labels1: np.ndarray, labels2: np.ndarray,
                              mask: np.ndarray) -> dict:
    """VI(P₁|M, P₂|M) — restrict two flat partitions to the leaves in ``mask``.

    Used to test whether a partition-level reorganization claim (whole-brain
    Δ_VI silent / strong) holds within an anatomically or functionally
    selected subset of leaves. Composes with any partition source — pass
    `fcluster(Z, k, 'maxclust')` for a per-k restriction, or
    `fcluster_at_h_rel(Z, h)` for an h-parametrized one.

    Parameters
    ----------
    labels1, labels2
        Cluster labels for the same N items. The full-tree partitions, not
        already restricted.
    mask
        Boolean (or integer index) array selecting the leaves to keep.

    Returns
    -------
    dict with keys
        ``vi`` — Variation of Information on the restricted labels (nats).
        ``n_subset`` — number of leaves in the subset (= ``mask.sum()``).
        ``k1_subset``, ``k2_subset`` — distinct labels seen in the subset.
        ``n_eff_1``, ``n_eff_2`` — Simpson effective number of clusters
            in the restricted partitions (giant-component diagnostic).
        ``max_cluster_fraction_1``, ``max_cluster_fraction_2`` — largest
            restricted-cluster size / ``n_subset`` (≈1 means the subset
            has all-collapsed-to-one-cluster, VI is trivially 0).
    """
    from .vi import compute_vi  # local import to avoid module-level cycles

    labels1 = np.asarray(labels1)
    labels2 = np.asarray(labels2)
    mask = np.asarray(mask)
    if mask.dtype == bool:
        sel = mask
    else:
        sel = np.zeros(labels1.size, dtype=bool)
        sel[mask] = True
    n_subset = int(sel.sum())
    if n_subset == 0:
        return {"vi": float("nan"), "n_subset": 0,
                 "k1_subset": 0, "k2_subset": 0,
                 "n_eff_1": 0.0, "n_eff_2": 0.0,
                 "max_cluster_fraction_1": float("nan"),
                 "max_cluster_fraction_2": float("nan")}
    sub1 = labels1[sel]
    sub2 = labels2[sel]
    return {
        "vi": float(compute_vi(sub1, sub2)),
        "n_subset": n_subset,
        "k1_subset": int(np.unique(sub1).size),
        "k2_subset": int(np.unique(sub2).size),
        "n_eff_1": simpson_neff(sub1),
        "n_eff_2": simpson_neff(sub2),
        "max_cluster_fraction_1": float(
            np.unique(sub1, return_counts=True)[1].max() / n_subset),
        "max_cluster_fraction_2": float(
            np.unique(sub2, return_counts=True)[1].max() / n_subset),
    }


def cluster_size_stats(labels: np.ndarray) -> dict:
    """Cluster-size summary for a flat partition.

    Returns ``{n, k, mean_size, max_size, min_size, singleton_fraction,
    max_cluster_fraction, n_eff}``. ``singleton_fraction`` = fraction of
    clusters of size 1, ``max_cluster_fraction`` = ``max|C| / N``.
    """
    labels = np.asarray(labels)
    n = int(labels.size)
    if n == 0:
        return {"n": 0, "k": 0, "mean_size": 0.0, "max_size": 0,
                "min_size": 0, "singleton_fraction": 0.0,
                "max_cluster_fraction": 0.0, "n_eff": 0.0}
    _, counts = np.unique(labels, return_counts=True)
    k = int(counts.size)
    return {
        "n": n,
        "k": k,
        "mean_size": float(counts.mean()),
        "max_size": int(counts.max()),
        "min_size": int(counts.min()),
        "singleton_fraction": float((counts == 1).sum()) / k,
        "max_cluster_fraction": float(counts.max()) / n,
        "n_eff": simpson_neff(labels),
    }
