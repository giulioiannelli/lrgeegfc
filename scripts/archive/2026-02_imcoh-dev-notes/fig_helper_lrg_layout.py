#!/usr/bin/env python3
"""LRG-seeded Kamada-Kawai layout prototype for dense ImCoh networks.

Principle: use the ImCoh LRG dendrogram cut at n communities to inflate
inter-community distances, so KK is mathematically forced to pull modules
apart. Intra-community distances stay at 1/w_ij; inter-community distances
are multiplied by a separation factor.

Run:
  python scripts/10_notes_imcoh/fig_helper_lrg_layout.py [-v]
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
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.metrics.hypothesis import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result


def _shaft_colors(probes):
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _comm_colors(labels):
    K = int(labels.max())
    cmap = plt.get_cmap("tab20", max(K, 3))
    return [cmap(int(l) - 1) for l in labels]


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


def lrg_seeded_layout(
    A: np.ndarray,
    community_labels: np.ndarray,
    separation: float = 5.0,
) -> np.ndarray:
    """Kamada-Kawai on a distance matrix where inter-community distances are
    multiplied by `separation`. Forces modular layout without ring placement.
    """
    N = A.shape[0]
    # Base distance: inverse of similarity
    D = 1.0 / (A + 1e-6)
    np.fill_diagonal(D, 0.0)

    # Inter-community mask
    cl = np.asarray(community_labels)
    inter = cl[:, None] != cl[None, :]
    D[inter] *= separation

    # Build graph with both weight and modified distance attributes
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    for i, j in zip(r, c):
        if A[i, j] > 0:
            G.add_edge(int(i), int(j),
                       weight=float(A[i, j]),
                       distance=float(D[i, j]))
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


def prototype(
    patient: str = "Pat_05", band: str = "beta", phase: str = "rest_pre",
    fc_method: str = "imcoh_abs",
    n_communities_list: tuple = (5, 10, 15, 20),
    separations: tuple = (1.5, 3.0, 5.0, 10.0),
    output_dir: Path = SECTION2_ROOT / "fig_helper",
    verbose: bool = False,
):
    """Grid: rows = n_communities (LRG cut level), cols = separation factor.
    Top row uses shaft coloring; rest use community coloring so you can see
    whether modules are genuinely separated or just probe-driven.
    """
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        print(f"  SKIP: no {fc_method} matrix")
        return
    lrg = load_lrg_result(patient, phase, band, fc_method)
    if lrg is None:
        print(f"  SKIP: no {fc_method} LRG")
        return

    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    shaft_colors = _shaft_colors(probes)
    N = mat.shape[0]
    if lrg.n_nodes != N:
        print(f"  WARN: N_lrg={lrg.n_nodes} != N_mat={N}; using min")
        N = min(N, lrg.n_nodes)

    A = np.abs(mat[:N, :N].copy())
    np.fill_diagonal(A, 0.0)
    shaft_colors = shaft_colors[:N]
    probes = probes[:N]

    rows = len(n_communities_list)
    cols = len(separations)

    fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 4.0 * rows))
    if rows == 1:
        axes = axes.reshape(1, -1)

    for i, nc in enumerate(n_communities_list):
        comm = fcluster(lrg.linkage_matrix, t=nc, criterion="maxclust")[:N]
        comm_cols = _comm_colors(comm)
        for j, sep in enumerate(separations):
            ax = axes[i, j]
            pos = lrg_seeded_layout(A, comm, separation=sep)
            _draw_gamma_edges(ax, pos, A, probes, gamma=6.0)
            ax.scatter(pos[:, 0], pos[:, 1], c=comm_cols, s=20,
                       edgecolors="white", linewidths=0.3, zorder=5)
            ax.axis("off")
            if i == 0:
                ax.set_title(f"sep×{sep}", fontsize=11, fontweight="bold")
            if j == 0:
                ax.text(-0.05, 0.5, f"n_comm={nc}", rotation=90,
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=11, fontweight="bold")
            margin = 0.08
            xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
            ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
            span = max(xmax - xmin, ymax - ymin, 1e-9)
            ax.set_xlim(xmin - margin * span, xmax + margin * span)
            ax.set_ylim(ymin - margin * span, ymax + margin * span)
            if verbose:
                print(f"    n_comm={nc} sep={sep}: done")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"LRG-seeded Kamada-Kawai — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"rows = LRG cut level, cols = inter-community separation factor; "
        f"nodes colored by LRG community, edges: red=same-probe, grey=cross",
        fontsize=13, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_fig(fig, output_dir / f"lrg_layout_{fc_method}_{patient}_{band}_{phase}")

    # Second figure: the chosen cut n=10 with shaft coloring for comparison
    nc = 10
    comm = fcluster(lrg.linkage_matrix, t=nc, criterion="maxclust")[:N]
    fig2, axes2 = plt.subplots(1, len(separations),
                               figsize=(4.0 * len(separations), 4.5))
    for j, sep in enumerate(separations):
        ax = axes2[j]
        pos = lrg_seeded_layout(A, comm, separation=sep)
        _draw_gamma_edges(ax, pos, A, probes, gamma=6.0)
        ax.scatter(pos[:, 0], pos[:, 1], c=shaft_colors, s=20,
                   edgecolors="white", linewidths=0.3, zorder=5)
        ax.set_title(f"sep×{sep}", fontsize=11, fontweight="bold")
        ax.axis("off")
        margin = 0.08
        xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
        ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
        span = max(xmax - xmin, ymax - ymin, 1e-9)
        ax.set_xlim(xmin - margin * span, xmax + margin * span)
        ax.set_ylim(ymin - margin * span, ymax + margin * span)
    fig2.suptitle(
        f"LRG-seeded KK, n_comm=10, colored by ELECTRODE SHAFT — "
        f"{patient} {band_tex} {phase}  (is the layout probe-independent?)",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig2.tight_layout()
    save_fig(fig2, output_dir / f"lrg_layout_shaft_{fc_method}_{patient}_{band}_{phase}")


def main():
    parser = argparse.ArgumentParser(description="LRG-seeded KK layout prototype.")
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--fc-method", default="imcoh_abs")
    parser.add_argument("--output-dir", type=Path,
                        default=SECTION2_ROOT / "fig_helper")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prototype(
        args.patient, args.band, args.phase, args.fc_method,
        output_dir=args.output_dir, verbose=args.verbose,
    )
    print("\nDone. Inspect the two figures and pick (n_comm, sep).")


if __name__ == "__main__":
    main()
