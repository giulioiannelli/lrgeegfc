"""Network-layout algorithms for FC adjacency matrices.

Promoted from ``scripts/01_compute/figures_for_notes/_shared.py`` on
2026-05-29 (per the library-naming meta-rule locked 2026-05-28: names
reflect general graph/network/statistics concepts, never manuscript-local
scope tokens).

This module is the canonical home for two-dimensional layout helpers that
operate on an FC adjacency matrix ``A`` and (optionally) the cached LRG
result for that ``(patient, phase, band, fc_method)`` cell. The dispatch
``compute_network_layout`` selects among:

- ``"mds_ultrametric"`` (default): classical MDS / PCoA on the LRG
  cophenetic distance matrix.
- ``"mds_log_ultrametric"``: log-compressed MDS for heavy-tailed D.
- ``"mds_lrg_continuous"``: classical MDS on the continuous-LRG
  propagator distance ``1/rho_ij`` at ``tau = 1 / lambda_max``.
- ``"kk_ultrametric"``: Kamada-Kawai using D_ij directly.
- ``"spring_lrg_distance"``: spring layout with edge weight ``~1/D_ij``.
- ``"sbm_nested"``: graph-tool nested SBM hierarchical layout.
- ``"spring"`` (fallback): spring on raw FC weights.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import NDArray

from lrg_eegfc.config.const import (
    LAYOUT_METHOD,
    LAYOUT_SEED_DEFAULT,
    SPRING_ITER_DEFAULT,
    spring_k_for,
)


__all__ = [
    "compute_network_layout",
    "compute_percolation_threshold",
]


def _seeded_initial_pos(N: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    theta = 2 * np.pi * np.arange(N) / N + rng.uniform(-0.2, 0.2, N)
    pos0 = np.column_stack([np.cos(theta), np.sin(theta)])
    return {i: pos0[i] for i in range(N)}


def compute_network_layout(
    A: NDArray,
    fc_method: str,
    patient: Optional[str] = None,
    phase: Optional[str] = None,
    band: Optional[str] = None,
    seed: int = LAYOUT_SEED_DEFAULT,
    spring_iter: int = SPRING_ITER_DEFAULT,
    spring_k: Optional[float] = 0.07,
    layout_method: Optional[str] = None,
) -> NDArray:
    """Return an (N, 2) layout matrix for an FC adjacency.

    Method chosen by :data:`LAYOUT_METHOD` per fc_method:
      - ``"mds_ultrametric"`` (default): classical MDS / PCoA on the LRG
        ultrametric distance matrix D_ij. Eigendecomposition of the
        double-centered ``-D²/2`` matrix; the first two eigenvectors
        scaled by ``√λ`` are the 2D coordinates. Deterministic, no
        iterative optimisation, distances preserved in least-squares.
      - ``"spring_lrg_distance"``: spring with edge weight ``~1/D_ij``.
      - ``"kk_ultrametric"``: Kamada-Kawai using D_ij directly.
      - ``"spring"``: spring on raw FC weights.
      - Falls back to spring on raw FC if the LRG cache is missing.
    """
    import networkx as nx
    from scipy.spatial.distance import squareform
    from lrg_eegfc.workflow.lrg import load_lrg_result

    N = A.shape[0]
    method = (layout_method if layout_method is not None
              else LAYOUT_METHOD.get(fc_method, "mds_ultrametric"))

    needs_lrg = method in (
        "kk_ultrametric", "spring_lrg_distance",
        "mds_ultrametric", "mds_log_ultrametric",
        "mds_lrg_continuous",
    )
    if needs_lrg and all(x is not None for x in (patient, phase, band)):
        lrg = load_lrg_result(patient, phase, band, fc_method)
        if lrg is not None:
            um = lrg.ultrametric_matrix
            U = squareform(um) if um.ndim == 1 else um
            M = min(N, U.shape[0])

            if method == "mds_lrg_continuous":
                from lrgsglib.utils.lrg.infocomm import lapl_dists
                from scipy.spatial.distance import squareform as _sf
                Adense = np.abs(A.copy())
                np.fill_diagonal(Adense, 0.0)
                if Adense.shape[0] != M:
                    Adense = Adense[:M, :M]
                L = np.diag(Adense.sum(axis=1)) - Adense

                # tau = 1 / lambda_max — inverse fastest-mode timescale.
                # DO NOT substitute lrg.optimal_threshold here: that is a
                # dendrogram cut HEIGHT in the cophenetic ultrametric,
                # NOT a diffusion time.
                eigvals_L = np.linalg.eigvalsh(L)
                lam_max = float(np.max(eigvals_L))
                if lam_max <= 0:
                    raise ValueError(
                        "Cannot pick tau for layout: Laplacian lambda_max <= 0; "
                        "graph likely empty or disconnected."
                    )
                tau_star = 1.0 / lam_max

                Dcond = lapl_dists(L, tau=tau_star)
                cap = float(np.percentile(Dcond, 99))
                Dcond = np.minimum(Dcond, cap)
                Dwork = _sf(Dcond)
                D2 = Dwork.astype(np.float64) ** 2
                n = M
                J = np.eye(n) - np.full((n, n), 1.0 / n)
                B = -0.5 * J @ D2 @ J
                B = 0.5 * (B + B.T)
                eigvals, eigvecs = np.linalg.eigh(B)
                idx = np.argsort(eigvals)[::-1][:2]
                lam = np.maximum(eigvals[idx], 0.0)
                coords = eigvecs[:, idx] * np.sqrt(lam)
                pos = {i: coords[i] for i in range(M)}

            elif method in ("mds_ultrametric", "mds_log_ultrametric"):
                if method == "mds_log_ultrametric":
                    Dwork = np.log1p(U.astype(np.float64))
                else:
                    Dwork = U.astype(np.float64)
                D2 = Dwork ** 2
                n = M
                J = np.eye(n) - np.full((n, n), 1.0 / n)
                B = -0.5 * J @ D2 @ J
                B = 0.5 * (B + B.T)
                eigvals, eigvecs = np.linalg.eigh(B)
                idx = np.argsort(eigvals)[::-1][:2]
                lam = np.maximum(eigvals[idx], 0.0)
                coords = eigvecs[:, idx] * np.sqrt(lam)
                pos = {i: coords[i] for i in range(M)}

            elif method == "kk_ultrametric":
                G = nx.complete_graph(M)
                dist = {i: {j: float(U[i, j]) for j in range(M) if j != i}
                        for i in range(M)}
                pos = nx.kamada_kawai_layout(
                    G, dist=dist, pos=_seeded_initial_pos(M, seed),
                )
            else:  # spring_lrg_distance
                with np.errstate(divide="ignore"):
                    W = np.where(U > 0, 1.0 / U, 0.0)
                np.fill_diagonal(W, 0.0)
                G = nx.from_numpy_array(W)
                k_val = spring_k if spring_k is not None else spring_k_for(
                    fc_method, band or "beta", M,
                )
                pos = nx.spring_layout(
                    G, k=k_val, iterations=spring_iter, seed=seed,
                    weight="weight",
                    pos=_seeded_initial_pos(M, seed),
                )
            arr = np.zeros((N, 2))
            for i in range(M):
                arr[i] = pos[i]
            return arr

    if method == "sbm_nested":
        data = _compute_sbm_data(A, seed=seed)
        return data["pos"]

    G = nx.from_numpy_array(A)
    k_val = spring_k if spring_k is not None else spring_k_for(
        fc_method, band or "beta", N,
    )
    pos = nx.spring_layout(
        G, k=k_val, iterations=spring_iter, seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


def _compute_sbm_data(A: NDArray, seed: int = LAYOUT_SEED_DEFAULT,
                      beta: float = 0.8) -> dict:
    """Nested-SBM inference + hierarchy tree + bezier control points.

    Returns a dict with keys:
        ``pos``   : (N, 2) vertex positions (leaves of the hierarchy tree).
        ``g``     : graph-tool Graph.
        ``weight``: edge property map with |FC| weights.
        ``state`` : fitted NestedBlockState.
        ``t``     : hierarchy tree (graph).
        ``tpos``  : tree-vertex positions (radial).
        ``cts``   : edge control points in the LOCAL edge frame.
    """
    import graph_tool.all as gt
    gt.seed_rng(seed)
    N = A.shape[0]

    g = gt.Graph(directed=False)
    g.add_vertex(N)
    weight = g.new_edge_property("double")
    ri, ci = np.triu_indices(N, k=1)
    vals = A[ri, ci]
    mask = vals > 0
    elist = np.column_stack([ri[mask], ci[mask], vals[mask]])
    g.add_edge_list(elist, eprops=[weight])

    state = gt.minimize_nested_blockmodel_dl(
        g, state_args=dict(recs=[weight], rec_types=["real-exponential"])
    )

    t, _tb, _tpos_base = gt.get_hierarchy_tree(state)
    root = next(v for v in t.vertices() if v.in_degree() == 0)
    tpos = gt.radial_tree_layout(t, root)
    cts = gt.get_hierarchy_control_points(g, t, tpos, beta=beta)

    pos = np.array(
        [[tpos[t.vertex(int(v))][0], tpos[t.vertex(int(v))][1]]
         for v in g.vertices()]
    )
    return {"pos": pos, "g": g, "weight": weight, "state": state,
            "t": t, "tpos": tpos, "cts": cts}


def compute_percolation_threshold(A: NDArray) -> float:
    """Minimum edge weight needed to keep the (undirected) graph connected.

    Defined as the minimum weight in the *maximum* spanning tree of |A|.
    Keeping edges with weight ``>= threshold`` yields the minimally
    connected sub-graph (the classical bond-percolation threshold from
    above).
    """
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.sparse import csr_matrix

    A_pos = np.where(A > 0, A, 0.0)
    np.fill_diagonal(A_pos, 0.0)
    mst_neg = minimum_spanning_tree(csr_matrix(-A_pos)).toarray()
    mst = np.where(mst_neg < 0, -mst_neg, 0.0)
    mst_w = mst[mst > 0]
    return float(mst_w.min()) if len(mst_w) else 0.0


# ---------------------------------------------------------------------------
# Network-template layouts (promoted from network_templates.py 2026-05-29,
# Phase 4-B split 7/7). Pure ``(A, ...) -> (N, 2)`` functions: spring, KK,
# spectral, Laplacian-PCA, circular-by-shaft, community-grouped,
# backbone-guided, LRG-seeded KK, SFDP / ARF / LRG-SFDP (graph-tool).
# Composed by ``compute_layout`` via ``LAYOUT_REGISTRY``.
# ---------------------------------------------------------------------------

# Imports required by these layouts that aren't already in this module
from typing import Callable, Dict, Sequence, Tuple
import networkx as nx
from lrg_eegfc.utils.fc.msc.sparsify import disparity_filter
from lrg_eegfc.workflow.lrg import load_lrg_result as _load_lrg_result_for_layouts  # noqa: F401


# ---------------------------------------------------------------------------
# Layouts — return (N, 2) np.ndarray.  Pure: no drawing, no side
# effects on inputs.
# ---------------------------------------------------------------------------

def _build_graph_with_distance(A: NDArray) -> nx.Graph:
    G = nx.from_numpy_array(np.asarray(A))
    for _, _, d in G.edges(data=True):
        w = d.get("weight", 1.0)
        d["distance"] = 1.0 / (w + 1e-6)
    return G


def layout_spring(
    A: NDArray,
    *,
    k: Optional[float] = None,
    k_base: float = 5.0,
    iterations: int = 600,
    seed: int = 42,
) -> NDArray:
    """Fruchterman–Reingold (NetworkX ``spring_layout``).

    Pass ``k`` (raw spring constant, fed to ``nx.spring_layout``
    verbatim) for direct control.  When ``k`` is ``None`` (default) the
    canonical dense-graph tuning ``k_base / √N`` is used instead.
    Stochastic; set ``seed`` for reproducibility.
    """
    A = np.asarray(A)
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    k_eff = k if k is not None else k_base / np.sqrt(N)
    pos = nx.spring_layout(
        G, k=k_eff, iterations=iterations,
        seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


def layout_kk(A: NDArray) -> NDArray:
    """Kamada–Kawai with ``distance = 1 / weight``.  Deterministic."""
    A = np.asarray(A)
    N = A.shape[0]
    G = _build_graph_with_distance(A)
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


def layout_spectral(A: NDArray) -> NDArray:
    """Project onto Laplacian eigenvectors λ₂, λ₃.  Deterministic.

    Falls flat for nearly-fully-connected uniform-weight networks
    (no clear bipartition); use :func:`layout_laplacian_pca` for a
    richer projection in that regime.
    """
    A = np.asarray(A)
    D = np.diag(A.sum(axis=1))
    L = D - A
    _, eigvecs = np.linalg.eigh(L)
    coords = eigvecs[:, 1:3]
    coords = coords - coords.min(axis=0)
    denom = coords.max(axis=0) + 1e-12
    return coords / denom


def layout_laplacian_pca(A: NDArray, n_components: int = 2) -> NDArray:
    """First five non-trivial Laplacian eigenvectors → 2-D PCA.
    Deterministic.  Captures multiscale structure invisible to the
    pure (λ₂, λ₃) projection.
    """
    A = np.asarray(A)
    D = np.diag(A.sum(axis=1))
    L = D - A
    _, eigvecs = np.linalg.eigh(L)
    emb = eigvecs[:, 1:6]
    emb_centered = emb - emb.mean(axis=0)
    U, S, _ = np.linalg.svd(emb_centered, full_matrices=False)
    coords = U[:, :n_components] * S[:n_components]
    coords = coords - coords.min(axis=0)
    denom = coords.max(axis=0) + 1e-12
    return coords / denom


def layout_circular_by_shaft(probes: Sequence[str]) -> NDArray:
    """Nodes on a unit circle, sorted by electrode probe.

    Deterministic baseline.  Use to rule out probe-trivial modules:
    if the same module structure appears here, it isn't real geometry.
    """
    N = len(probes)
    uniq = sorted(set(probes))
    probe_to_idx = {p: i for i, p in enumerate(uniq)}
    order = np.argsort([probe_to_idx[p] * 100000 + i for i, p in enumerate(probes)])
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    coords = np.column_stack([np.cos(theta), np.sin(theta)])
    inverse = np.empty(N, dtype=int)
    inverse[order] = np.arange(N)
    return coords[inverse]


def layout_community_grouped(
    A: NDArray,
    *,
    seed: int = 42,
    R_out: float = 3.0,
) -> NDArray:
    """Louvain communities on an outer ring; sub-spring within each.

    Reads communities from ``A`` directly (Louvain on the weighted
    graph).  Sub-spring is stochastic; pin ``seed`` for reproducibility.
    """
    A = np.asarray(A)
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    try:
        import networkx.algorithms.community as nx_comm
        communities = list(
            nx_comm.louvain_communities(
                G, weight="weight", seed=seed, resolution=1.0,
            )
        )
    except Exception:
        communities = [set(range(N))]

    K = len(communities)
    coords = np.zeros((N, 2))
    for ci, nodes in enumerate(communities):
        theta_c = 2 * np.pi * ci / max(K, 1)
        center = np.array([R_out * np.cos(theta_c), R_out * np.sin(theta_c)])
        nodes = sorted(nodes)
        nn = len(nodes)
        if nn == 1:
            coords[nodes[0]] = center
            continue
        sub = A[np.ix_(nodes, nodes)]
        sub_G = nx.from_numpy_array(sub)
        sub_pos = nx.spring_layout(
            sub_G, k=1.5 / np.sqrt(nn),
            iterations=200, seed=seed, weight="weight",
        )
        local = np.array([sub_pos[i] for i in range(nn)])
        span = np.abs(local).max() + 1e-9
        local = local / span
        for idx, node in enumerate(nodes):
            coords[node] = center + local[idx]
    return coords


def layout_backbone_guided(
    A: NDArray,
    *,
    alpha: float = 0.05,
    iterations: int = 400,
    seed: int = 42,
) -> NDArray:
    """Lay out the disparity-filter backbone; full graph drawn on top.

    The layout sees only the statistically-significant subset (Serrano
    et al. 2009) so modules emerge; the caller still draws every
    edge of ``A``.
    """
    B = disparity_filter(np.abs(A), alpha=alpha)
    N = A.shape[0]
    G_back = nx.from_numpy_array(B)
    pos = nx.spring_layout(
        G_back, k=2.0 / np.sqrt(N),
        iterations=iterations, seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


def layout_lrg_kk(
    A: NDArray,
    ultrametric_matrix: NDArray,
) -> NDArray:
    """Kamada–Kawai on the LRG ultrametric distance ``T_ρ = 1/ρ̂_ij``.

    The LRG ultrametric is a single mathematical object — it doesn't
    depend on any community-label cut.  We feed it directly to KK as
    the pairwise distance the layout should preserve.

    Parameters
    ----------
    A
        FC matrix (only used for shape; layout uses ultrametric).
    ultrametric_matrix
        LRG ultrametric distances in scipy *condensed* form (length
        N·(N−1)/2), as cached on :class:`LRGResult.ultrametric_matrix`.

    Notes
    -----
    Modules emerge naturally because the LRG ultrametric clusters
    nodes that share a high communication probability — there is no
    need for an artificial inter-cluster distance multiplier.
    """
    from scipy.spatial.distance import squareform
    A = np.asarray(A)
    N = A.shape[0]
    expected = N * (N - 1) // 2
    if ultrametric_matrix.shape[0] != expected:
        raise ValueError(
            f"ultrametric_matrix has length {ultrametric_matrix.shape[0]}, "
            f"expected {expected} for N={N}"
        )
    D = squareform(np.asarray(ultrametric_matrix))  # (N, N)
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    d_flat = D[r, c]
    for i, j, d in zip(r, c, d_flat):
        G.add_edge(int(i), int(j), distance=float(d))
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


def layout_sfdp(
    A: NDArray,
    *,
    gamma: float = 0.1,
    mu: float = 0.0,
    C: float = 0.2,
    p: float = 2.0,
    max_iter: int = 0,
) -> NDArray:
    """Scalable Force-Directed Placement (graph-tool).  Multilevel,
    fast on dense graphs (≥ 50× faster than NetworkX on N=115).
    """
    import graph_tool.all as gt
    g, ew = matrix_to_gt(np.abs(A))
    pos = gt.sfdp_layout(
        g, eweight=ew, gamma=gamma, mu=mu, C=C, p=p, max_iter=max_iter,
    )
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


def layout_arf(A: NDArray, *, d: float = 0.5, a: float = 10.0) -> NDArray:
    """Attractive-Repulsive Force layout (Geipel 2007, graph-tool).
    Robust on very dense graphs.  Deterministic post-init.
    """
    import graph_tool.all as gt
    g, ew = matrix_to_gt(np.abs(A))
    pos = gt.arf_layout(g, weight=ew, d=d, a=a)
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


def layout_lrg_sfdp(
    A: NDArray,
    community_labels: Sequence[int],
    *,
    gamma: float = 0.1,
    mu: float = 0.0,
    C: float = 0.2,
    p: float = 2.0,
    max_iter: int = 0,
) -> NDArray:
    """SFDP with ``groups=`` from an LRG dendrogram cut.  graph-tool's
    SFDP gives same-group vertices an extra attractive force; ``gamma``
    controls module-separation strength (0.03 → weak, 0.6 → strong).
    """
    import graph_tool.all as gt
    g, ew = matrix_to_gt(np.abs(A))
    grp = g.new_vertex_property("int")
    grp.a = np.asarray(community_labels, dtype=int)
    pos = gt.sfdp_layout(
        g, eweight=ew, groups=grp,
        gamma=gamma, mu=mu, C=C, p=p, max_iter=max_iter,
    )
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


# Layout registry — string name → callable.  ``compute_layout`` is the
# entry point templates use; it dispatches on the name and threads
# extra kwargs (``probes`` for circular_by_shaft, ``community_labels``
# for LRG-seeded layouts).

LAYOUT_REGISTRY: Dict[str, Callable] = {
    "spring": layout_spring,
    "kk": layout_kk,
    "spectral": layout_spectral,
    "laplacian_pca": layout_laplacian_pca,
    "circular_by_shaft": layout_circular_by_shaft,
    "community_grouped": layout_community_grouped,
    "backbone_guided": layout_backbone_guided,
    "lrg_kk": layout_lrg_kk,
    "sfdp": layout_sfdp,
    "arf": layout_arf,
    "lrg_sfdp": layout_lrg_sfdp,
}

# Default gallery for ``mosaic_layout_compare``.  Includes one of every
# family.  LRG-seeded layouts only render when ``community_labels`` is
# supplied; absent labels they're skipped silently.
DEFAULT_GALLERY: Tuple[str, ...] = (
    "spring",
    "kk",
    "spectral",
    "laplacian_pca",
    "community_grouped",
    "backbone_guided",
    "circular_by_shaft",
    "sfdp",
    "arf",
    "lrg_kk",
    "lrg_sfdp",
)


def is_lrg_layout(name: str) -> bool:
    return name.startswith("lrg_")


def compute_layout(
    name: str,
    A: NDArray,
    *,
    probes: Optional[Sequence[str]] = None,
    community_labels: Optional[Sequence[int]] = None,
    ultrametric_matrix: Optional[NDArray] = None,
    **kwargs,
) -> NDArray:
    """Dispatch layout by name; thread the right extra args.

    - ``circular_by_shaft`` consumes ``probes``
    - ``lrg_kk`` consumes ``ultrametric_matrix`` (LRG ultrametric ``T_ρ``,
      condensed) — does NOT need community_labels
    - ``lrg_sfdp`` consumes ``community_labels`` (cut of the LRG dendrogram)
    - Other layouts take only ``A`` plus layout-specific kwargs.
    """
    if name not in LAYOUT_REGISTRY:
        raise ValueError(f"Unknown layout {name!r}; choose from {list(LAYOUT_REGISTRY)}")
    fn = LAYOUT_REGISTRY[name]
    if name == "circular_by_shaft":
        if probes is None:
            raise ValueError("circular_by_shaft requires probes=")
        return fn(probes, **kwargs)
    if name == "lrg_kk":
        if ultrametric_matrix is None:
            raise ValueError(
                "lrg_kk requires ultrametric_matrix= (condensed T_ρ from LRGResult)"
            )
        return fn(A, ultrametric_matrix, **kwargs)
    if name == "lrg_sfdp":
        if community_labels is None:
            raise ValueError("lrg_sfdp requires community_labels=")
        return fn(A, community_labels, **kwargs)
    return fn(A, **kwargs)
