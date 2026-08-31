"""Classical low-level graph descriptors for weighted, non-negative FC adjacencies.

These are the *usual* network-science scalars — node strength, weighted
clustering, centralities, closeness — computed directly in numpy/scipy (no
networkx in the hot path, except the optional Brandes betweenness which imports
it lazily). They exist so a "pairwise-blind ladder" audit can ask whether any of
them, run through the *same* cross-phase ρ_sym estimator and matched-strength
null as the LRG cophenetic trace, reproduces the multiscale discrimination
(band-selectivity, inference-β-specificity, →OFC). See the scope report
``.agents/guides/task-persistence-investigation/2026-07-09_pairwise-descriptor-ladder.md``.

Every descriptor is a pure map ``W ↦ v`` from a symmetric, non-negative,
zero-diagonal adjacency ``W ∈ ℝ^{N×N}`` to a vector on a common index set:

- **pair-level** (length ``N(N−1)/2``): :func:`raw_edges` (the FC edge weights
  themselves — the bottom rung of the ladder).
- **node-level** (length ``N``): :func:`node_strength`,
  :func:`weighted_clustering_onnela`, :func:`eigenvector_centrality`,
  :func:`pagerank`, :func:`closeness_centrality`, :func:`betweenness_centrality`.

Names are general graph concepts (never manuscript-local), reusable by any future
figure or audit. ``node_strength`` unifies the ``W.sum(axis=1)`` one-liner that
was previously inlined across the matched-strength surrogate, disparity-filter,
and figure scripts (library-first rule).

On a *dense* FC graph (every ``|ImCoh| > 0``) each node's degree is ``N−1`` and
these reduce to cheap linear algebra: strength ``O(N²)``; Onnela clustering one
``O(N³)`` matmul; eigenvector centrality one ``eigh``; PageRank a power iteration;
closeness one all-pairs shortest path. Betweenness (Brandes) is the only
super-cubic term and is provided but not used by default.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.csgraph import shortest_path


def raw_edges(W: np.ndarray) -> np.ndarray:
    """Upper-triangular edge-weight vector (condensed), length ``N(N−1)/2``.

    The lowest-order representation: the pairwise FC values themselves, with no
    graph-theoretic transform. This is the ladder's bottom rung — a cross-phase
    ρ_sym on ``raw_edges`` is the raw-FC trace.
    """
    iu = np.triu_indices(W.shape[0], k=1)
    return np.asarray(W, float)[iu]


def node_strength(W: np.ndarray) -> np.ndarray:
    """Weighted degree ``s_i = Σ_j W_ij`` (length ``N``).

    The quantity the matched-strength null preserves *exactly*; hence a
    cross-phase trace built on ``node_strength`` is reproduced by the surrogate
    by construction and cannot clear the null. Unifies the inline
    ``W.sum(axis=1)`` scattered across the codebase.
    """
    return np.asarray(W, float).sum(axis=1)


def weighted_clustering_onnela(W: np.ndarray, w_max: float = 1.0) -> np.ndarray:
    """Onnela weighted clustering coefficient (length ``N``).

    ``C_i = diag(Â³)_i / (k_i (k_i − 1))`` with ``Â = (W / w_max)^{1/3}`` and
    ``k_i`` the (unweighted) degree = number of non-zero incident edges. This is
    the ``nx.clustering(G, weight='weight')`` definition computed in closed form
    (``diag(Â³)_i = Σ_k (Â²)_{ik} Â_{ik}`` for symmetric ``Â``), one matmul.

    On a fully-connected FC graph ``k_i = N−1`` for all nodes; the descriptor is
    then dominated by the weight geometry and empirically tracks strength at
    ``r ≈ 0.96`` — the canonical "strength-slaving" showcase.
    """
    W = np.asarray(W, float)
    A = np.cbrt(np.clip(W, 0.0, None) / w_max)
    np.fill_diagonal(A, 0.0)
    # diag(A^3) for symmetric A, without forming the full cube
    diag_a3 = ((A @ A) * A).sum(axis=1)
    k = (W > 0).sum(axis=1).astype(float)
    denom = k * (k - 1.0)
    out = np.zeros_like(diag_a3)
    nz = denom > 0
    out[nz] = diag_a3[nz] / denom[nz]
    return out


def eigenvector_centrality(W: np.ndarray) -> np.ndarray:
    """Eigenvector centrality = leading eigenvector of ``W`` (length ``N``).

    For a symmetric non-negative ``W`` the Perron eigenvector (largest
    eigenvalue) is non-negative; returned as its absolute value, L2-normalized.
    """
    W = np.asarray(W, float)
    _, vecs = np.linalg.eigh(W)
    v = np.abs(vecs[:, -1])
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def pagerank(W: np.ndarray, alpha: float = 0.85,
             tol: float = 1e-12, max_iter: int = 500) -> np.ndarray:
    """Weighted PageRank stationary distribution (length ``N``).

    Row-stochastic transition ``P_ij = W_ij / s_i`` (dangling rows → uniform),
    ``v = (1−α)/N · 1 + α Pᵀ v`` by power iteration. ``alpha`` is the damping.
    """
    W = np.asarray(W, float)
    N = W.shape[0]
    s = W.sum(axis=1)
    P = np.zeros_like(W)
    nz = s > 0
    P[nz] = W[nz] / s[nz, None]
    P[~nz] = 1.0 / N
    Pt = P.T.copy()
    v = np.full(N, 1.0 / N)
    teleport = (1.0 - alpha) / N
    for _ in range(max_iter):
        v_new = teleport + alpha * (Pt @ v)
        v_new /= v_new.sum()
        if np.abs(v_new - v).sum() < tol:
            v = v_new
            break
        v = v_new
    return v


def closeness_centrality(W: np.ndarray) -> np.ndarray:
    """Weighted closeness on the ``1/W`` distance graph (length ``N``).

    Edge distance ``d_ij = 1/W_ij`` (∞ where ``W_ij = 0``); all-pairs shortest
    paths (Dijkstra); ``closeness_i = (n_reachable) / Σ_j d_ij`` over reachable
    ``j ≠ i``. On a dense FC graph the direct edge is usually the geodesic.
    """
    W = np.asarray(W, float)
    N = W.shape[0]
    with np.errstate(divide="ignore"):
        D = 1.0 / W
    D[W <= 0] = np.inf
    np.fill_diagonal(D, 0.0)
    sp = shortest_path(D, method="D", directed=False)
    out = np.zeros(N)
    idx = np.arange(N)
    for i in range(N):
        d = sp[i]
        reach = np.isfinite(d) & (idx != i)
        tot = d[reach].sum()
        if tot > 0:
            out[i] = reach.sum() / tot
    return out


def betweenness_centrality(W: np.ndarray) -> np.ndarray:
    """Weighted betweenness on the ``1/W`` distance graph (length ``N``).

    Brandes algorithm via networkx (lazy import). Provided for completeness;
    on a fully-connected weighted graph most geodesics are direct edges, so
    betweenness is near-degenerate (most nodes ≈ 0) and expensive — the ladder
    audit gates its inclusion on a runtime check.
    """
    import networkx as nx  # lazy: keep the library import-light

    W = np.asarray(W, float)
    N = W.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))
    iu = np.triu_indices(N, k=1)
    for i, j, w in zip(iu[0].tolist(), iu[1].tolist(), W[iu].tolist()):
        if w > 0:
            G.add_edge(i, j, distance=1.0 / w)
    bc = nx.betweenness_centrality(G, weight="distance", normalized=True)
    return np.array([bc[i] for i in range(N)], dtype=float)


def gini_coefficient(x: np.ndarray) -> float:
    """Gini coefficient of a non-negative vector, in ``[0, 1]``.

    ``0`` = every entry equal, ``1`` = all the mass on one entry. Computed from
    the sorted vector as ``G = (2 Σ i x_(i)) / (n Σ x) − (n+1)/n``. A
    scale-invariant concentration measure: unlike the coefficient of variation
    it is bounded and insensitive to the units of ``x``, which matters when
    comparing weight fields produced by different transforms (e.g. ``⟨|Im C|⟩_f``
    vs ``⟨(Im C)²⟩_f``).

    Negative entries are not meaningful for a concentration measure and are
    clipped at ``0``.
    """
    v = np.sort(np.maximum(np.asarray(x, float).ravel(), 0.0))
    n = v.size
    tot = v.sum()
    if n == 0 or tot <= 0:
        return float("nan")
    idx = np.arange(1, n + 1)
    return float((2.0 * (idx * v).sum()) / (n * tot) - (n + 1.0) / n)


def weight_heterogeneity(W: np.ndarray) -> dict:
    """Scale-free descriptors of how unevenly a graph's edge weight is spread.

    On a **dense** weighted graph the heat-kernel propagator is degenerate
    (single specific-heat peak) precisely when the weight field is strongly
    heterogeneous, so these are the quantities that say how far a substrate sits
    from the collapsed regime -- independent of any cross-phase hypothesis.

    Returns ``gini``, ``cv`` (coefficient of variation), ``p99_over_p50``,
    ``max_over_p50`` and ``participation_ratio`` -- the last being
    ``(Σw)² / (Σw² · n_edges) ∈ (0, 1]``, the fraction of edges that
    effectively carry the weight (``1`` = perfectly uniform).
    """
    iu = np.triu_indices(W.shape[0], k=1)
    w = np.maximum(np.asarray(W, float)[iu], 0.0)
    mu = w.mean()
    med = float(np.median(w))
    return dict(
        gini=gini_coefficient(w),
        cv=float(w.std() / mu) if mu > 0 else float("nan"),
        p99_over_p50=float(np.percentile(w, 99) / med) if med > 0 else float("nan"),
        max_over_p50=float(w.max() / med) if med > 0 else float("nan"),
        participation_ratio=float(w.sum() ** 2 / ((w ** 2).sum() * w.size))
        if (w ** 2).sum() > 0 else float("nan"),
    )


# Registry: name -> (callable, level). Level 'pair' = length N(N-1)/2,
# 'node' = length N. Used by the pairwise-descriptor-ladder audit.
DESCRIPTORS = {
    "raw_fc": (raw_edges, "pair"),
    "strength": (node_strength, "node"),
    "clustering_onnela": (weighted_clustering_onnela, "node"),
    "eigcent": (eigenvector_centrality, "node"),
    "pagerank": (pagerank, "node"),
    "closeness": (closeness_centrality, "node"),
    "betweenness": (betweenness_centrality, "node"),
}
