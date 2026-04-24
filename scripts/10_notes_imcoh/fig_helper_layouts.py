#!/usr/bin/env python3
"""Gallery of network-layout algorithms for dense weighted ImCoh networks.

The spring layout is stochastic and parameter-sensitive — it also needs
per-case tuning for dense uniform networks like ImCoh. This helper produces
a gallery of alternative layouts, most of them deterministic and parameter-
free, so we can pick one that reveals structure without exploration.

Layouts included:
  - spring (Fruchterman-Reingold) with 3 k values
  - kamada_kawai      — weighted shortest-path based, deterministic
  - spectral           — Laplacian eigenvectors (2nd & 3rd smallest)
  - circular by shaft  — deterministic baseline
  - community-grouped  — Louvain communities placed on a ring, nodes
                          within each community on an inner ring
  - PCA on Laplacian   — project first 2 non-trivial Laplacian eigenvectors
  - backbone-guided    — lay out the disparity-filter backbone; draw all edges

Run:
  python scripts/10_notes_imcoh/fig_helper_layouts.py [-v]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import networkx as nx

from lrg_eegfc.utils.metrics.hypothesis import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _shaft_colors(probes):
    """Color each node by its electrode shaft (tab20)."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _draw_gamma_edges(ax, pos_arr, mat, probes, gamma: float = 6.0,
                      width_range=(0.15, 4.0), alpha_range=(0.03, 0.9)):
    N = mat.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = mat[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    t = (w_act / w_act.max()) ** gamma
    widths = width_range[0] + t * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])
    order = np.argsort(w_act)
    segments, colors, linew = [], [], []
    for oi in order:
        idx = active[oi]
        i, j = r[idx], c[idx]
        segments.append([pos_arr[i], pos_arr[j]])
        linew.append(widths[oi])
        if probes[i] == probes[j]:
            colors.append((0.8, 0.2, 0.2, alphas[oi]))
        else:
            colors.append((0.35, 0.35, 0.35, alphas[oi]))
    lc = LineCollection(segments, colors=colors, linewidths=linew,
                        zorder=1, rasterized=True)
    ax.add_collection(lc)


def _finalize_axes(ax, pos_arr, title: str):
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.axis("off")
    margin = 0.08
    xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
    ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    ax.set_xlim(xmin - margin * span, xmax + margin * span)
    ax.set_ylim(ymin - margin * span, ymax + margin * span)


# ---------------------------------------------------------------------------
# Layout implementations (all return an (N, 2) array)
# ---------------------------------------------------------------------------

def layout_spring(G: nx.Graph, k_base: float, seed: int = 42) -> np.ndarray:
    """Spring layout from an existing graph (no re-conversion)."""
    N = G.number_of_nodes()
    pos = nx.spring_layout(G, k=k_base / np.sqrt(N), iterations=600,
                           seed=seed, weight="weight")
    return np.array([pos[i] for i in range(N)])


def layout_kamada_kawai(G: nx.Graph) -> np.ndarray:
    """Kamada-Kawai via the built-in edge weight (no O(N²) distance dict).

    Deterministic. NetworkX's KK uses a numpy-based FW shortest path
    internally, which is lighter on memory than building an explicit
    distance dict.
    """
    N = G.number_of_nodes()
    # Invert similarity to distance via the 'weight' attribute in-place
    # by using a wrapper weight function: d = 1 / w.
    # NetworkX KK accepts a custom weight attribute name.
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


def layout_spectral(A: np.ndarray) -> np.ndarray:
    """Project onto the 2nd and 3rd smallest Laplacian eigenvectors."""
    N = A.shape[0]
    D = np.diag(A.sum(axis=1))
    L = D - A
    eigvals, eigvecs = np.linalg.eigh(L)
    # Skip the zero eigenvalue (constant vector)
    coords = eigvecs[:, 1:3]
    # Normalize to [0, 1] per axis for consistent display
    coords = coords - coords.min(axis=0)
    denom = coords.max(axis=0) + 1e-12
    coords = coords / denom
    return coords


def layout_circular_by_shaft(probes: list[str]) -> np.ndarray:
    """Nodes sorted by shaft on a circle. Deterministic baseline."""
    N = len(probes)
    uniq = sorted(set(probes))
    probe_to_idx = {p: i for i, p in enumerate(uniq)}
    order = np.argsort([probe_to_idx[p] * 10000 + i for i, p in enumerate(probes)])
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    coords = np.column_stack([np.cos(theta), np.sin(theta)])
    # Place node i at the position corresponding to its rank in order
    inverse = np.empty(N, dtype=int)
    inverse[order] = np.arange(N)
    return coords[inverse]


def layout_community_grouped(G: nx.Graph, A: np.ndarray,
                              seed: int = 42) -> np.ndarray:
    """Louvain communities on an outer ring, sub-spring within each.

    Uses the shared graph G to avoid re-conversion. Keeps the per-community
    induced subgraph small to limit memory.
    """
    N = A.shape[0]
    try:
        import networkx.algorithms.community as nx_comm
        communities = list(nx_comm.louvain_communities(G, weight="weight",
                                                       seed=seed, resolution=1.0))
    except Exception:
        communities = [set(range(N))]

    K = len(communities)
    coords = np.zeros((N, 2))
    R_out = 3.0
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
        sub_pos = nx.spring_layout(sub_G, k=1.5 / np.sqrt(nn),
                                   iterations=200, seed=seed,
                                   weight="weight")
        local = np.array([sub_pos[i] for i in range(nn)])
        del sub_G, sub
        span = np.abs(local).max() + 1e-9
        local = local / span
        for idx, node in enumerate(nodes):
            coords[node] = center + local[idx]
    return coords


def layout_laplacian_pca(A: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Low-dim projection from the first few non-trivial Laplacian eigenvectors.
    Very similar to spectral but uses more components + PCA for the final 2D.
    """
    N = A.shape[0]
    D = np.diag(A.sum(axis=1))
    L = D - A
    eigvals, eigvecs = np.linalg.eigh(L)
    emb = eigvecs[:, 1:6]  # skip trivial, take next 5
    # Further reduce to 2D via PCA
    emb_centered = emb - emb.mean(axis=0)
    U, S, Vt = np.linalg.svd(emb_centered, full_matrices=False)
    coords = U[:, :n_components] * S[:n_components]
    # Normalize
    coords = coords - coords.min(axis=0)
    denom = coords.max(axis=0) + 1e-12
    return coords / denom


def layout_backbone_guided(A: np.ndarray, alpha: float = 0.05,
                            seed: int = 42) -> np.ndarray:
    """Use the disparity-filter (Serrano et al.) statistically-significant
    backbone to drive the layout. Not a top-X% threshold — it's a null-model
    backbone that keeps every locally-significant edge regardless of weight.
    The drawn network uses the FULL graph; only the layout is driven by the
    backbone so modules can emerge.
    """
    from lrg_eegfc.utils.fc.msc.sparsify import disparity_filter
    B = disparity_filter(A.copy(), alpha=alpha)
    N = A.shape[0]
    G_back = nx.from_numpy_array(B)
    pos = nx.spring_layout(G_back, k=2.0 / np.sqrt(N), iterations=400,
                           seed=seed, weight="weight")
    coords = np.array([pos[i] for i in range(N)])
    del G_back, B
    return coords


# ---------------------------------------------------------------------------
# Gallery
# ---------------------------------------------------------------------------

def _build_graph_with_distance(A: np.ndarray) -> nx.Graph:
    """Build one weighted graph with both 'weight' and 'distance' edge attrs.
    Done once per case to avoid repeated numpy→nx conversions.
    """
    G = nx.from_numpy_array(A)
    # Attach distance = 1/weight for algorithms that want a shortest-path metric
    for u, v, data in G.edges(data=True):
        data["distance"] = 1.0 / (data["weight"] + 1e-6)
    return G


def run_gallery(
    patient: str, band: str, phase: str, fc_method: str,
    output_dir: Path, gamma: float = 6.0,
    width_range: tuple[float, float] = (0.15, 4.0),
    verbose: bool = False,
):
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        print(f"  SKIP: no {fc_method} for {patient} {band} {phase}")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    G = _build_graph_with_distance(A)

    # Build layout closures with access to the shared G
    layout_specs = [
        ("Kamada-Kawai",          lambda: layout_kamada_kawai(G)),
        ("Spectral (λ₂, λ₃)",     lambda: layout_spectral(A)),
        ("Laplacian PCA (5→2)",   lambda: layout_laplacian_pca(A)),
        ("Community-grouped",     lambda: layout_community_grouped(G, A)),
        ("Backbone-guided",       lambda: layout_backbone_guided(A)),
        ("Circular by shaft",     lambda: layout_circular_by_shaft(probes)),
        ("Spring k_base=5",       lambda: layout_spring(G, 5.0)),
        ("Spring k_base=15",      lambda: layout_spring(G, 15.0)),
    ]

    n = len(layout_specs)
    cols = 4
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 4.0 * rows))
    axes = axes.flatten()

    for idx, (name, fn) in enumerate(layout_specs):
        ax = axes[idx]
        try:
            pos = fn()
        except Exception as exc:
            ax.text(0.5, 0.5, f"{name}\nFAILED:\n{exc}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="red")
            ax.axis("off")
            if verbose:
                print(f"    {name}: {exc}")
            continue
        _draw_gamma_edges(ax, pos, A, probes,
                          gamma=gamma, width_range=width_range)
        ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=20,
                   edgecolors="white", linewidths=0.3, zorder=5)
        _finalize_axes(ax, pos, name)
        if verbose:
            print(f"    {name}: ok")

    for idx in range(n, rows * cols):
        axes[idx].axis("off")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Layout gallery — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"full weighted graph, γ={gamma}, nodes colored by electrode shaft",
        fontsize=13, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_fig(
        fig,
        output_dir / f"layout_gallery_{fc_method}_{patient}_{band}_{phase}",
    )

    # Free graph/matrix before next case
    del G, A, mat
    import gc; gc.collect()


def main():
    parser = argparse.ArgumentParser(
        description="Layout-algorithm gallery for dense ImCoh networks.",
    )
    parser.add_argument("--output-dir", type=Path,
                        default=SECTION2_ROOT / "fig_helper" / "layouts")
    parser.add_argument("--gamma", type=float, default=6.0)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Sample variety: 3 patients × 3 bands × (rest + task)
    cases = [
        ("Pat_05", "beta",  "rest_pre"),
        ("Pat_05", "alpha", "task_learn"),
        ("Pat_02", "theta", "rest_post"),
        ("Pat_02", "alpha", "task_learn"),
        ("Pat_08", "beta",  "task_test"),
        ("Pat_08", "low_gamma", "rest_pre"),
    ]
    for patient, band, phase in cases:
        print(f"--- Gallery: {patient} {band} {phase} ---")
        run_gallery(
            patient, band, phase, fc_method="imcoh_abs",
            output_dir=args.output_dir, gamma=args.gamma,
            verbose=args.verbose,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
