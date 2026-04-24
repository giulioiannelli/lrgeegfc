"""Scalar tree-distance metrics for rooted dendrograms on aligned leaves.

Implements three scalar measures missing from the 14-metric MSC-era sweep
(see ``data/reports/imcoh_vi/stage0a_metric_catalog.md``). All three combine
branch-length AND topology information into one number per tree pair,
which none of the 14 did.

- **Kendall-Colijn** (Kendall & Colijn 2016, MBE 33(10):2735-2743) — per
  leaf-pair (m, M) vectors: ``m(i,j)`` = edges from root to MRCA (pure
  topology), ``M(i,j)`` = summed branch lengths along that path (pure
  heights). λ-blended L2 distance.
- **Matching Cluster** (Bogdanowicz & Giaro 2013, IJAMCS 23(3):669-684) —
  bipartite Hungarian matching between internal clusters of the two trees;
  cost = leaf-set symmetric difference (optionally plus height difference).
- **Weighted Robinson-Foulds** — bipartition-weighted edit cost; standard
  branch-length baseline from the phylogenetics literature. Implemented
  directly from scipy linkages; optional dendropy cross-check.

All functions take scipy ``linkage`` matrices ``Z`` (shape ``(n-1, 4)``)
with leaf IDs 0..n-1 and cluster IDs n..2n-2. Trees are assumed rooted
(UPGMA). Leaves are aligned across trees (same channel identities).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.cluster.hierarchy import to_tree
from scipy.optimize import linear_sum_assignment


__all__ = [
    "kc_vectors",
    "kc_distance",
    "matching_cluster_distance",
    "weighted_rf_distance",
]


# ---------------------------------------------------------------------------
# Kendall-Colijn (Kendall & Colijn 2016, Mol Biol Evol)
# ---------------------------------------------------------------------------

def _build_parent_height_depth(Z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(parent, height, depth)`` arrays of length ``2n-1``.

    - ``parent[c]`` = cluster ID of the merge that absorbs cluster ``c``,
      or ``-1`` if ``c`` is the root.
    - ``height[c]`` = merge height of cluster ``c`` (0 for leaves, ``Z[k,2]``
      for internal cluster ``n+k``).
    - ``depth[c]`` = number of edges from the root to node ``c``
      (root: 0; root's children: 1; etc.).
    """
    n = Z.shape[0] + 1
    total = 2 * n - 1
    parent = np.full(total, -1, dtype=np.int64)
    height = np.zeros(total, dtype=float)
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        cid = n + k
        parent[a] = cid
        parent[b] = cid
        height[cid] = float(Z[k, 2])
    # Root has id 2n-2 (last merge)
    depth = np.zeros(total, dtype=np.int64)
    # Process nodes in decreasing id order: the root (2n-2) is processed
    # first with depth 0; every child sees parent already depth-assigned.
    for c in range(total - 1, -1, -1):
        p = parent[c]
        if p >= 0:
            depth[c] = depth[p] + 1
    return parent, height, depth


def _ancestor_chain(leaf: int, parent: np.ndarray) -> list[int]:
    """Ancestor chain from root DOWN to the leaf itself (inclusive)."""
    chain = []
    node = leaf
    while node != -1:
        chain.append(node)
        node = int(parent[node])
    chain.reverse()
    return chain


def kc_vectors(Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the Kendall-Colijn ``(m, M)`` pair-vectors for a UPGMA tree.

    For each leaf pair ``(i, j)`` with ``i < j``, in scipy ``squareform``
    ordering, the function returns:

    - ``m[k]`` = number of edges from the root to the MRCA of ``(i, j)``
      (depth of the MRCA in edges, ``0`` at the root).
    - ``M[k]`` = merge height of that MRCA (in units of the ultrametric).

    Both arrays have length ``n*(n-1)/2`` with the standard scipy upper-
    triangle-row-major iteration over leaf pairs.
    """
    n = Z.shape[0] + 1
    parent, height, depth = _build_parent_height_depth(Z)
    chains = [_ancestor_chain(i, parent) for i in range(n)]
    chain_len = np.array([len(c) for c in chains], dtype=np.int64)

    npairs = n * (n - 1) // 2
    m = np.zeros(npairs, dtype=np.int64)
    M = np.zeros(npairs, dtype=float)
    idx = 0
    for i in range(n):
        ci = chains[i]
        for j in range(i + 1, n):
            cj = chains[j]
            # Walk from root until chains diverge — last common = MRCA.
            limit = min(chain_len[i], chain_len[j])
            mrca = ci[0]  # root always shared
            for k in range(limit):
                if ci[k] == cj[k]:
                    mrca = ci[k]
                else:
                    break
            m[idx] = depth[mrca]
            M[idx] = height[mrca]
            idx += 1
    return m, M


def kc_distance(Z1: np.ndarray, Z2: np.ndarray, lam: float = 0.5,
                 normalize: bool = True) -> float:
    """Kendall-Colijn L2 distance between two trees on the same leaves.

    ``d_KC(T1, T2; λ) = || v_λ(T1) − v_λ(T2) ||_2``, where
    ``v_λ = (1-λ) * m + λ * M``. ``λ = 0`` is pure topology;
    ``λ = 1`` is pure heights; ``λ = 0.5`` is a balanced blend.

    If ``normalize`` is True, the m-vector is divided by its max and the
    M-vector by its max **pooled across the two trees** before blending,
    so the two components contribute on comparable scales. This prevents
    dominance of whichever vector has the larger magnitude by default
    (KC without normalization is scale-sensitive; Kendall & Colijn 2016
    recommend reporting the curve d(λ) rather than a single d at λ=0.5).
    """
    m1, M1 = kc_vectors(Z1)
    m2, M2 = kc_vectors(Z2)
    if normalize:
        mmax = max(m1.max(), m2.max(), 1)
        Mmax = max(M1.max(), M2.max(), 1e-12)
        m1n, m2n = m1 / mmax, m2 / mmax
        M1n, M2n = M1 / Mmax, M2 / Mmax
    else:
        m1n, m2n = m1.astype(float), m2.astype(float)
        M1n, M2n = M1, M2
    v1 = (1 - lam) * m1n + lam * M1n
    v2 = (1 - lam) * m2n + lam * M2n
    return float(np.linalg.norm(v1 - v2))


# ---------------------------------------------------------------------------
# Matching Cluster distance (Bogdanowicz & Giaro 2013)
# ---------------------------------------------------------------------------

def _internal_clusters(Z: np.ndarray) -> list[tuple[frozenset[int], float]]:
    """Enumerate (leafset, merge_height) for every internal node of Z.

    Returns ``n-1`` entries, one per internal node. Order does not matter
    for MC — the Hungarian matching handles reordering.
    """
    root = to_tree(Z, rd=False)
    out: list[tuple[frozenset[int], float]] = []

    def walk(node) -> frozenset[int]:
        if node.is_leaf():
            return frozenset([int(node.id)])
        left = walk(node.get_left())
        right = walk(node.get_right())
        leaves = left | right
        out.append((leaves, float(node.dist)))
        return leaves

    walk(root)
    return out


def matching_cluster_distance(
    Z1: np.ndarray, Z2: np.ndarray, *,
    weighted: bool = False, height_scale: Optional[float] = None,
) -> float:
    """Bogdanowicz-Giaro 2013 Matching Cluster distance, rooted variant.

    Construction: every internal node defines a cluster (= leafset under
    the node). Build the bipartite graph between T1's and T2's internal
    clusters; edge cost is the symmetric-difference size
    ``|A Δ B|``. Hungarian min-weight perfect matching gives the MC
    distance = sum of edge costs over matched pairs.

    If ``weighted`` is True, the cost becomes
    ``|A Δ B| + height_scale * |h_A − h_B|`` where ``height_scale`` defaults
    to ``n / max_height`` so the two components contribute on comparable
    scales. This makes the metric topology + heights simultaneously. The
    Bogdanowicz-Giaro 2013 Theorem 4 shows weighted MC remains a metric.

    Trees must be rooted binary on the same leaf set. For UPGMA the number
    of internal nodes is ``n-1`` for both, so the bipartite graph is
    balanced.
    """
    c1 = _internal_clusters(Z1)
    c2 = _internal_clusters(Z2)
    n1, n2 = len(c1), len(c2)
    if n1 != n2:
        raise ValueError(f"MC requires same-size leaf sets; got {n1} vs {n2} internal nodes")
    n_leaves = Z1.shape[0] + 1

    cost = np.zeros((n1, n2), dtype=float)
    if weighted and height_scale is None:
        hmax = max(max(h for _, h in c1), max(h for _, h in c2), 1e-12)
        height_scale = n_leaves / hmax
    for i, (l1, h1) in enumerate(c1):
        for j, (l2, h2) in enumerate(c2):
            sym = len(l1 ^ l2)
            if weighted:
                cost[i, j] = sym + height_scale * abs(h1 - h2)
            else:
                cost[i, j] = float(sym)

    row, col = linear_sum_assignment(cost)
    return float(cost[row, col].sum())


# ---------------------------------------------------------------------------
# Weighted Robinson-Foulds distance (branch-length weighted)
# ---------------------------------------------------------------------------

def _bipartitions_with_weights(Z: np.ndarray) -> dict[frozenset[int], float]:
    """Map each non-root internal bipartition (leafset) -> branch length.

    Branch length = ``height(parent) - height(node)``. The root's bipartition
    is the full leaf set; it is excluded because every same-leaf tree has
    the same root bipartition and it would contribute zero to any
    bipartition-based comparison.
    """
    n = Z.shape[0] + 1
    parent, height, _ = _build_parent_height_depth(Z)
    root = 2 * n - 2
    out: dict[frozenset[int], float] = {}
    # Collect leafset per internal node
    leaves_under: dict[int, frozenset[int]] = {i: frozenset([i]) for i in range(n)}
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        cid = n + k
        leaves_under[cid] = leaves_under[a] | leaves_under[b]
    # Branch length above each internal node (except root)
    for k in range(n - 1):
        cid = n + k
        if cid == root:
            continue
        p = int(parent[cid])
        bl = float(height[p] - height[cid])
        out[leaves_under[cid]] = bl
    return out


def weighted_rf_distance(Z1: np.ndarray, Z2: np.ndarray) -> float:
    """Branch-length-weighted Robinson-Foulds distance.

    For each internal bipartition (leafset) ``c``, let ``w_T(c)`` be the
    branch length above the corresponding node in tree ``T`` (or ``0`` if
    ``c`` is absent in ``T``). Distance:

    ``d_wRF(T1, T2) = sum_c |w_T1(c) - w_T2(c)|``

    A shared bipartition with different heights contributes ``|h-diff|``;
    a bipartition present in only one tree contributes its own branch
    length. Simultaneously scores topology and heights in one scalar.

    Notoriously noisy between close trees (one leaf swap can flip many
    bipartitions at once) — kept as the standard phylogenetic baseline
    paired with KC and MC.
    """
    b1 = _bipartitions_with_weights(Z1)
    b2 = _bipartitions_with_weights(Z2)
    keys = set(b1) | set(b2)
    total = 0.0
    for c in keys:
        w1 = b1.get(c, 0.0)
        w2 = b2.get(c, 0.0)
        total += abs(w1 - w2)
    return total
