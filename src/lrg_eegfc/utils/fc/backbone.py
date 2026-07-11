"""Parameter-free connected backbones for dense weighted graphs.

The LRG diffusion propagator is degenerate on a **fully-connected** weighted
graph (every pair has a non-zero edge): at ``tau = 1/lambda_max`` the cophenetic
tree is a near-monotone map of the raw FC, so ``D = 1/K`` and ``D = 1/A`` give
the same hierarchy (audit_174). LRG is built for **sparse / topological**
graphs, where the propagator integrates multi-step reachability that a single
adjacency scale cannot resolve (Villegas 2023/2025; guide section 6).

This module extracts a sparse **backbone** of a dense FC matrix so the LRG
machinery operates in the regime it is designed for. Every backbone here is:

- **parameter-free** (or swept over a pre-committed density grid for
  robustness) -- no significance level or magic threshold to tune;
- **spanning** -- every node is retained, so a per-phase backbone never drops a
  contact and cross-phase per-pair alignment (``rho_sym``) is preserved for
  free (no giant-component gymnastics);
- **connected** -- a single component, so graph geodesic / diffusion distances
  are finite and the UPGMA tree is well-defined.

Backbones
---------
- :func:`maximum_spanning_tree` -- the minimal connected backbone, ``N-1``
  edges. Zero parameters, guaranteed connected and spanning.
- :func:`tmfg_backbone` -- Triangulated Maximally Filtered Graph (Massara, Di
  Matteo & Aste 2016). Chordal, ``3N-6`` edges, contains 3- and 4-cliques so it
  carries genuine mesoscale structure, yet is still parameter-free,
  deterministic, connected and spanning. The canonical parameter-free choice.
- :func:`mst_union_top_fraction` -- ``MST`` union the strongest ``frac`` of all
  edges. Guarantees connectivity at any target density; the density grid used
  for robustness sweeps (``frac = 1.0`` recovers the fully-connected graph).

References
----------
Massara, G. P., Di Matteo, T. & Aste, T. "Network filtering for big data:
Triangulated Maximally Filtered Graph." *J. Complex Netw.* **5**, 161 (2016).

Tumminello, M., Aste, T., Di Matteo, T. & Mantegna, R. N. "A tool for
filtering information in complex systems." *PNAS* **102**, 10421 (2005).
"""
from __future__ import annotations

import numpy as np
import numba
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import (
    connected_components,
    minimum_spanning_tree,
    shortest_path,
)

__all__ = [
    "maximum_spanning_tree",
    "tmfg_backbone",
    "mst_union_top_fraction",
    "percolation_backbone",
    "backbone_density",
    "giant_component",
    "geodesic_distance",
]


def _clean(W: NDArray) -> NDArray:
    """Non-negative, symmetric, zero-diagonal copy of ``W``."""
    A = np.asarray(W, dtype=float).copy()
    A = np.maximum(A, 0.0)
    np.fill_diagonal(A, 0.0)
    return 0.5 * (A + A.T)


def backbone_density(A: NDArray) -> float:
    """Edge density = kept edges / ``N(N-1)/2`` (upper triangle)."""
    N = A.shape[0]
    if N < 2:
        return 0.0
    r, c = np.triu_indices(N, k=1)
    return float((A[r, c] > 0).sum()) / (N * (N - 1) / 2.0)


def maximum_spanning_tree(W: NDArray) -> NDArray:
    """Maximum spanning tree backbone (``N-1`` edges), weights preserved.

    Computed as the minimum spanning tree of ``-A`` (scipy). Connected and
    spanning whenever the input graph is connected (a fully-connected FC matrix
    always is). Zero parameters.
    """
    A = _clean(W)
    N = A.shape[0]
    mst = minimum_spanning_tree(csr_matrix(-A)).toarray()
    keep = mst != 0.0
    keep = keep | keep.T
    B = np.zeros((N, N))
    B[keep] = A[keep]
    np.fill_diagonal(B, 0.0)
    return B


@numba.njit(cache=True)
def _tmfg_edges(W):
    """Greedy TMFG (gain-vector implementation). Returns a boolean keep-mask.

    Start from the 4 highest-strength vertices (a tetrahedron), then repeatedly
    insert the remaining vertex/face pair of maximum edge-weight gain, splitting
    the chosen triangular face into three. O(N^2) amortised via a per-vertex
    best-face gain cache. Deterministic; ties broken by lowest index.
    """
    N = W.shape[0]
    keep = np.zeros((N, N), np.bool_)
    if N <= 4:
        for i in range(N):
            for j in range(i + 1, N):
                keep[i, j] = True
                keep[j, i] = True
        return keep

    # initial tetrahedron: 4 vertices of highest total strength
    strength = np.zeros(N)
    for i in range(N):
        s = 0.0
        for j in range(N):
            s += W[i, j]
        strength[i] = s
    order = np.argsort(strength)[::-1]
    c0 = order[0]
    c1 = order[1]
    c2 = order[2]
    c3 = order[3]

    in_graph = np.zeros(N, np.bool_)
    for v in (c0, c1, c2, c3):
        in_graph[v] = True
    for a in (c0, c1, c2, c3):
        for b in (c0, c1, c2, c3):
            if a < b:
                keep[a, b] = True
                keep[b, a] = True

    cap = 2 * N  # tetra has 4 faces; each insertion nets +2 faces => <= 2N-4
    faces = np.empty((cap, 3), np.int64)
    tetra = ((c0, c1, c2), (c0, c1, c3), (c0, c2, c3), (c1, c2, c3))
    nf = 0
    for t in tetra:
        faces[nf, 0] = t[0]
        faces[nf, 1] = t[1]
        faces[nf, 2] = t[2]
        nf += 1

    gain = np.full(N, -1.0)
    bestf = np.full(N, -1, np.int64)
    for v in range(N):
        if in_graph[v]:
            continue
        bg = -1.0
        bf = -1
        for f in range(nf):
            g = W[v, faces[f, 0]] + W[v, faces[f, 1]] + W[v, faces[f, 2]]
            if g > bg:
                bg = g
                bf = f
        gain[v] = bg
        bestf[v] = bf

    n_added = 4
    while n_added < N:
        vbest = -1
        gbest = -1e18
        for v in range(N):
            if in_graph[v]:
                continue
            if gain[v] > gbest:
                gbest = gain[v]
                vbest = v
        v = vbest
        f = bestf[v]
        x = faces[f, 0]
        y = faces[f, 1]
        z = faces[f, 2]
        keep[v, x] = True; keep[x, v] = True
        keep[v, y] = True; keep[y, v] = True
        keep[v, z] = True; keep[z, v] = True
        in_graph[v] = True
        n_added += 1

        # split face f = (x,y,z) into (v,x,y),(v,x,z),(v,y,z)
        faces[f, 0] = v; faces[f, 1] = x; faces[f, 2] = y
        faces[nf, 0] = v; faces[nf, 1] = x; faces[nf, 2] = z
        f1 = nf; nf += 1
        faces[nf, 0] = v; faces[nf, 1] = y; faces[nf, 2] = z
        f2 = nf; nf += 1

        for u in range(N):
            if in_graph[u]:
                continue
            if bestf[u] == f:
                # its best face was consumed -> recompute over all faces
                bg = -1.0
                bf = -1
                for ff in range(nf):
                    g = W[u, faces[ff, 0]] + W[u, faces[ff, 1]] + W[u, faces[ff, 2]]
                    if g > bg:
                        bg = g
                        bf = ff
                gain[u] = bg
                bestf[u] = bf
            else:
                for ff in (f, f1, f2):
                    g = W[u, faces[ff, 0]] + W[u, faces[ff, 1]] + W[u, faces[ff, 2]]
                    if g > gain[u]:
                        gain[u] = g
                        bestf[u] = ff
    return keep


def tmfg_backbone(W: NDArray) -> NDArray:
    """Triangulated Maximally Filtered Graph backbone, weights preserved.

    Parameter-free, deterministic, connected, spanning; ``3N-6`` edges for
    ``N > 3``. The canonical parameter-free multiscale backbone (Massara 2016).
    """
    A = _clean(W)
    keep = _tmfg_edges(A)
    B = np.zeros_like(A)
    B[keep] = A[keep]
    np.fill_diagonal(B, 0.0)
    return B


def percolation_backbone(W: NDArray) -> tuple[NDArray, float]:
    """Parameter-free percolation backbone: keep every edge above the connectivity
    bottleneck.

    Raising a global weight threshold ``theta`` drops edges; the giant component
    stays P_inf = 1 (all N nodes in one component) until ``theta`` passes the
    *bottleneck* -- the weakest edge that connectivity requires. That bottleneck is
    exactly ``min`` over the **maximum spanning tree**'s edges (any threshold at or
    below it keeps the graph connected; any threshold above it disconnects it).

    Keeping ALL edges ``>= theta*`` (not just the tree) yields the *sparsest
    globally-connected* graph, which is **cycle-rich** (every edge stronger than the
    bottleneck survives), not a tree. No free parameter -- the density is *set by the
    data* (the surviving-edge fraction) and reported, not chosen.

    Returns ``(backbone, theta_star)``. Density = :func:`backbone_density`.
    """
    A = _clean(W)
    mst = maximum_spanning_tree(A)
    pos = mst[mst > 0]
    theta = float(pos.min()) if pos.size else 0.0
    B = np.where(A >= theta, A, 0.0)
    np.fill_diagonal(B, 0.0)
    return B, theta


def mst_union_top_fraction(W: NDArray, frac: float) -> NDArray:
    """``MST`` union the strongest ``frac`` of upper-triangle edges.

    Connected at every ``frac`` (the MST spans all nodes); density increases
    monotonically with ``frac``. ``frac >= 1`` returns the full (cleaned)
    graph, recovering the fully-connected anchor. Weights preserved.
    """
    A = _clean(W)
    N = A.shape[0]
    if frac >= 1.0:
        return A
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    keep = np.zeros((N, N), bool)
    if w.size:
        n_keep = max(1, int(np.ceil(frac * w.size)))
        thr = np.partition(w, -n_keep)[-n_keep]
        sel = w >= thr
        keep[r[sel], c[sel]] = True
        keep[c[sel], r[sel]] = True
    mst = maximum_spanning_tree(A) > 0.0
    keep = keep | mst
    B = np.zeros((N, N))
    B[keep] = A[keep]
    np.fill_diagonal(B, 0.0)
    return B


def giant_component(A: NDArray) -> tuple[NDArray, NDArray]:
    """Largest connected component. Returns ``(sub_A, node_index)``.

    Only needed for non-spanning filters (disparity / hard threshold); the
    spanning backbones above never disconnect. ``node_index`` maps rows of
    ``sub_A`` back to the original node ids.
    """
    A = _clean(A)
    n_comp, labels = connected_components(csr_matrix(A > 0), directed=False)
    if n_comp == 1:
        return A, np.arange(A.shape[0])
    sizes = np.bincount(labels)
    big = int(np.argmax(sizes))
    idx = np.where(labels == big)[0]
    return A[np.ix_(idx, idx)], idx


def geodesic_distance(A: NDArray, mode: str = "inverse") -> NDArray:
    """Shortest-path (geodesic) distance on a weighted backbone.

    Non-diffusion topological baseline for the propagator ablation: edge length
    is ``1/w`` (``mode='inverse'``) or ``-log(w/w_max)`` (``mode='logratio'``),
    then Dijkstra all-pairs shortest paths. Finite for a connected backbone.
    """
    A = _clean(A)
    N = A.shape[0]
    L = np.zeros((N, N))
    pos = A > 0
    if mode == "inverse":
        L[pos] = 1.0 / A[pos]
    elif mode == "logratio":
        wmax = A[pos].max() if pos.any() else 1.0
        L[pos] = -np.log(A[pos] / (wmax * (1.0 + 1e-12)))
    else:
        raise ValueError(f"unknown mode {mode!r}")
    D = shortest_path(csr_matrix(L), method="D", directed=False)
    D = np.asarray(D, float)
    np.fill_diagonal(D, 0.0)
    return D
