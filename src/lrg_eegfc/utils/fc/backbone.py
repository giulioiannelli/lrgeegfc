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
    "pmfg_backbone",
    "disparity_backbone",
    "mst_union_top_fraction",
    "percolation_backbone",
    "percolation_sweep",
    "select_backbone",
    "backbone_density",
    "giant_component",
    "geodesic_distance",
]


def select_backbone(
    W: NDArray, kind: str, frac: float = 0.20, disparity_alpha: float = 0.20
) -> NDArray:
    """Dispatch to a named sparsifier -- one entry point for all pipelines.

    ``kind`` (matches the ``SA_BACKBONE`` env used across ``scripts/01_compute/
    sparsified_arc``):

    - ``"dense"``      -- the cleaned fully-connected graph (no sparsification);
    - ``"mst020"`` / ``"mst"`` -- :func:`mst_union_top_fraction` at the ``frac``
      argument (the callers pass ``frac`` explicitly; the name only selects the
      family);
    - ``"perc"``       -- :func:`percolation_backbone`;
    - ``"tmfg"``       -- :func:`tmfg_backbone` (fast chordal planar filter);
    - ``"pmfg"``       -- :func:`pmfg_backbone` (exact planar filter, slow);
    - ``"disparity"``  -- :func:`disparity_backbone` at ``disparity_alpha``
      (spanning via MST-union).

    Every returned backbone is connected + spanning, so cross-phase per-pair
    alignment holds without giant-component bookkeeping.
    """
    k = kind.lower()
    if k == "dense":
        return _clean(W)
    if k.startswith("mst"):                       # name selects family; frac is explicit
        return mst_union_top_fraction(W, frac)
    if k == "perc":
        return percolation_backbone(W)[0]
    if k == "tmfg":
        return tmfg_backbone(W)
    if k == "pmfg":
        return pmfg_backbone(W)
    if k == "disparity":
        return disparity_backbone(W, alpha=disparity_alpha, ensure_connected=True)
    raise ValueError(f"unknown backbone kind {kind!r}")


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


def pmfg_backbone(W: NDArray) -> NDArray:
    """Planar Maximally Filtered Graph backbone (Tumminello 2005), weights preserved.

    The *exact* planar-maximal filter, of which :func:`tmfg_backbone` is a fast
    chordal (greedy triangulation) approximation. Insert edges in **decreasing
    weight order**, keeping an edge iff the graph stays planar; stop at the
    maximal-planar budget ``3N-6``. The result is parameter-free, connected,
    spanning (it provably **contains the maximum spanning tree** -- the top-weight
    edges are inserted first and a forest is always planar) and cycle-rich (a
    maximal planar graph is a triangulation).

    Unlike TMFG, PMFG is **not** constrained to be chordal, so it is the
    principled parent method; TMFG trades that generality for an O(N^2) numba
    computation. This exact version uses incremental planarity testing
    (``networkx.check_planarity``, Left-Right algorithm) and costs O(N^3) --
    seconds per graph at N ~ 120 -- so it is intended for **observed graphs**,
    not inside a per-surrogate null (use TMFG there and verify PMFG ~ TMFG on the
    observed readout).

    Reference: Tumminello, Aste, Di Matteo & Mantegna, *PNAS* **102**, 10421 (2005).
    """
    import networkx as nx

    A = _clean(W)
    N = A.shape[0]
    if N <= 4:                                  # already planar (K4 is planar)
        B = A.copy()
        np.fill_diagonal(B, 0.0)
        return B
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    order = np.argsort(-w)                        # decreasing weight
    G = nx.Graph()
    G.add_nodes_from(range(N))
    emax = 3 * N - 6
    added = 0
    for k in order:
        if w[k] <= 0.0:                           # never add zero/negative edges
            break
        i, j = int(r[k]), int(c[k])
        G.add_edge(i, j)
        planar, _ = nx.check_planarity(G, counterexample=False)
        if not planar:
            G.remove_edge(i, j)
        else:
            added += 1
            if added >= emax:                    # maximal planar: no further edge fits
                break
    B = np.zeros((N, N))
    for i, j in G.edges():
        B[i, j] = A[i, j]
        B[j, i] = A[i, j]
    np.fill_diagonal(B, 0.0)
    return B


def disparity_backbone(
    W: NDArray, alpha: float = 0.05, ensure_connected: bool = True
) -> NDArray:
    """Disparity-filter (multiscale) backbone (Serrano 2009), weights preserved.

    The one principled sparsifier that drops the planarity prior of PMFG/TMFG:
    it is purely **statistical and local**. For each node ``i`` of degree
    ``k_i``, normalise its incident weights ``p_ij = w_ij / s_i`` and keep edge
    ``(i, j)`` iff it is significant against the null "the node's strength is
    distributed uniformly at random over its edges" from **either** endpoint::

        (1 - p_ij) ** (k_i - 1) < alpha        (Serrano et al. eq. for the p-value)

    This preserves edges that carry a disproportionate share of a node's
    strength at *every scale* of strength -- hence "multiscale backbone" -- with
    no global threshold. Its costs (relative to PMFG/TMFG): it carries the
    significance level ``alpha`` (swept, not tuned, in our use) and it can
    **fragment** the graph. To keep the backbone spanning -- required for
    cross-phase per-pair alignment (``rho_sym``) -- we optionally union it with
    the maximum spanning tree (``ensure_connected=True``); the MST adds at most
    ``N-1`` edges and never removes a disparity-significant one.

    O(N^2) pure-numpy: fast enough to recompute inside a per-surrogate null.

    Reference: Serrano, Boguna & Vespignani, "Extracting the multiscale backbone
    of complex weighted networks." *PNAS* **106**, 6483 (2009).
    """
    A = _clean(W)
    N = A.shape[0]
    if N < 2:
        return A
    strength = A.sum(1)
    deg = (A > 0).sum(1)                          # k_i (dense FC: N-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = A / strength[:, None]                 # p_ij from i's perspective (row-normalised)
        p = np.where(np.isfinite(p), p, 0.0)
        # p-value that edge (i,j) is compatible with the null, from i's view
        pval = np.where(A > 0, (1.0 - p) ** (deg[:, None] - 1), 1.0)
    sig = (A > 0) & (pval < alpha)               # significant from i's (row) view
    keep = sig | sig.T                           # keep if significant from EITHER endpoint
    if ensure_connected:
        keep = keep | (maximum_spanning_tree(A) > 0.0)
    np.fill_diagonal(keep, False)
    B = np.zeros((N, N))
    B[keep] = A[keep]
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


def percolation_sweep(W: NDArray, n_thresholds: int = 200) -> dict:
    """Percolation curve of a dense weighted graph vs a global weight threshold.

    Raises a global threshold ``theta`` over the observed edge-weight range and,
    at each level, keeps edges with ``weight >= theta`` and measures how the giant
    (largest connected) component erodes:

    - ``p_inf`` -- node fraction in the LCC (``= 1`` while the graph spans all N
      nodes, drops once ``theta`` passes the connectivity **bottleneck**);
    - ``e_inf`` -- fraction of the graph's edges that lie **inside** the LCC;
    - ``edge_frac`` -- surviving-edge fraction (kept / total), the density.

    Thresholds are placed at the empirical edge-weight quantiles (plus ``0`` and
    ``w_max``) so the curve is resolved where edges actually drop -- most edges
    fall in the first few percent of a linear grid, which would hide the giant
    component's collapse. ``theta_star`` is the largest threshold with
    ``p_inf == 1`` (the bottleneck), identical to :func:`percolation_backbone`.

    Returns a dict with arrays ``theta, p_inf, e_inf, edge_frac`` (each length
    ``n_kept <= n_thresholds+2``), scalars ``theta_star, edge_frac_star,
    p_inf_star`` (``== 1`` by construction), ``n_nodes, n_edges``, and
    ``w_max`` for normalising the threshold axis to ``theta/w_max in [0,1]``.
    """
    A = _clean(W)
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    w = w[w > 0]
    total_edges = float(w.size)
    if total_edges == 0:
        return dict(theta=np.zeros(1), p_inf=np.zeros(1), e_inf=np.zeros(1),
                    edge_frac=np.zeros(1), theta_star=0.0, edge_frac_star=0.0,
                    p_inf_star=0.0, n_nodes=int(N), n_edges=0, w_max=0.0)
    w_max = float(w.max())
    # thresholds at empirical quantiles (dense near the weak-edge mass) + endpoints
    qs = np.linspace(0.0, 1.0, n_thresholds)
    theta = np.unique(np.concatenate(([0.0], np.quantile(w, qs), [w_max])))

    p_inf = np.empty(theta.size)
    e_inf = np.empty(theta.size)
    edge_frac = np.empty(theta.size)
    for i, th in enumerate(theta):
        keep = A >= th
        np.fill_diagonal(keep, False)
        kept_edges = float(keep[r, c].sum())
        edge_frac[i] = kept_edges / total_edges
        n_comp, labels = connected_components(csr_matrix(keep), directed=False)
        if n_comp == N:                       # fully disconnected
            p_inf[i] = 1.0 / N
            e_inf[i] = 0.0
            continue
        sizes = np.bincount(labels)
        big = int(np.argmax(sizes))
        mask = labels == big
        p_inf[i] = float(mask.sum()) / N
        sub = keep[np.ix_(mask, mask)]
        e_inf[i] = float(np.triu(sub, 1).sum()) / total_edges

    connected = p_inf >= 1.0 - 1e-12
    if connected.any():
        j = int(np.where(connected)[0].max())
        theta_star = float(theta[j])
        edge_frac_star = float(edge_frac[j])
        p_inf_star = float(p_inf[j])
    else:                                     # never fully connected (should not happen for FC)
        theta_star, edge_frac_star, p_inf_star = 0.0, 1.0, float(p_inf[0])
    return dict(theta=theta, p_inf=p_inf, e_inf=e_inf, edge_frac=edge_frac,
                theta_star=theta_star, edge_frac_star=edge_frac_star,
                p_inf_star=p_inf_star, n_nodes=int(N), n_edges=int(total_edges),
                w_max=w_max)


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
