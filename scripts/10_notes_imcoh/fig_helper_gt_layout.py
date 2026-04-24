#!/usr/bin/env python3
"""Graph-tool based layouts for dense ImCoh networks.

graph-tool's SFDP and ARF layouts are production-quality force-directed
algorithms that handle dense weighted graphs far better than NetworkX.
SFDP additionally accepts a ``groups`` argument that forces same-group
nodes to attract more strongly — exactly what we need to reveal LRG
community structure.

Run:
  python scripts/10_notes_imcoh/fig_helper_gt_layout.py [-v]
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
from scipy.cluster.hierarchy import fcluster

import graph_tool.all as gt

from _shared import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _shaft_colors(probes):
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _comm_colors(labels):
    K = max(int(np.max(labels)), 2)
    cmap = plt.get_cmap("tab20", K)
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


def _gt_graph_from_matrix(A: np.ndarray):
    """Build a graph-tool Graph with a float edge weight property.
    Only stores upper triangle as undirected edges."""
    N = A.shape[0]
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    ew = g.new_edge_property("double")
    r, c = np.triu_indices(N, k=1)
    mask = A[r, c] > 0
    if mask.any():
        edge_list = np.column_stack([r[mask], c[mask]])
        g.add_edge_list(edge_list)
        ew.a = A[r[mask], c[mask]]
    return g, ew


def _positions_to_array(pos_prop, N: int) -> np.ndarray:
    return np.array([list(pos_prop[v]) for v in range(N)], dtype=float)


# ---------------------------------------------------------------------------
# Layouts
# ---------------------------------------------------------------------------

def layout_sfdp(A: np.ndarray, groups: np.ndarray = None,
                gamma: float = 0.3, mu: float = 0.0,
                C: float = 0.2, p: float = 2.0) -> np.ndarray:
    """SFDP force-directed layout.
    - groups: per-vertex int community label → same-group attraction.
    - gamma: controls the ratio of repulsion vs attraction (smaller = tighter).
    - C: scale of the repulsive force.
    - p: exponent of the attractive force (1 = softer, 2 = stiffer).
    Weighted via edge weights (stronger edges pull harder).
    """
    g, ew = _gt_graph_from_matrix(A)
    grp = None
    if groups is not None:
        grp = g.new_vertex_property("int")
        grp.a = np.asarray(groups, dtype=int)
    pos = gt.sfdp_layout(
        g, eweight=ew, groups=grp,
        gamma=gamma, mu=mu, C=C, p=p,
        max_iter=0,
    )
    return _positions_to_array(pos, A.shape[0])


def layout_arf(A: np.ndarray, d: float = 0.5, a: float = 10.0) -> np.ndarray:
    """Attractive-repulsive force layout (Geipel 2007).
    Robust on dense weighted graphs. Deterministic given the initial positions.
    """
    g, ew = _gt_graph_from_matrix(A)
    pos = gt.arf_layout(g, weight=ew, d=d, a=a, max_iter=0)
    return _positions_to_array(pos, A.shape[0])


def layout_fruchterman(A: np.ndarray, a: float = 1.0, r: float = 1.0,
                        scale: float = 1.0) -> np.ndarray:
    """Fruchterman-Reingold in graph-tool (numerically solid version)."""
    g, ew = _gt_graph_from_matrix(A)
    pos = gt.fruchterman_reingold_layout(g, weight=ew, a=a, r=r,
                                         scale=scale, n_iter=200)
    return _positions_to_array(pos, A.shape[0])


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_prototype(
    patient: str, band: str, phase: str, fc_method: str,
    output_dir: Path, verbose: bool = False,
):
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        print(f"  SKIP: no {fc_method} for {patient}")
        return
    lrg = load_lrg_result(patient, phase, band, fc_method)
    if lrg is None:
        print(f"  SKIP: no LRG")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    N = min(mat.shape[0], lrg.n_nodes)
    A = np.abs(mat[:N, :N].copy())
    np.fill_diagonal(A, 0)
    probes = probes[:N]
    shaft_cols = _shaft_colors(probes)

    # LRG communities at n=10 (tunable)
    n_comm = 10
    comm = fcluster(lrg.linkage_matrix, t=n_comm, criterion="maxclust")[:N]
    comm_cols = _comm_colors(comm)

    # ---- Layout variants to compare ----
    variants = [
        ("SFDP no groups",       lambda: layout_sfdp(A, groups=None)),
        ("SFDP + LRG groups",    lambda: layout_sfdp(A, groups=comm)),
        ("SFDP tight (γ=0.1)",   lambda: layout_sfdp(A, groups=comm, gamma=0.1)),
        ("SFDP stiff p=4",       lambda: layout_sfdp(A, groups=comm, p=4.0)),
        ("ARF layout",           lambda: layout_arf(A)),
        ("ARF (d=2, a=20)",      lambda: layout_arf(A, d=2.0, a=20.0)),
        ("Fruchterman-Reingold", lambda: layout_fruchterman(A)),
        ("FR weak a=0.3",        lambda: layout_fruchterman(A, a=0.3, r=2.0)),
    ]

    cols = 4
    rows = (len(variants) + cols - 1) // cols

    # ---- figure 1: nodes colored by LRG community ----
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 4.2 * rows))
    axes = axes.flatten()
    for idx, (name, fn) in enumerate(variants):
        ax = axes[idx]
        try:
            pos = fn()
        except Exception as exc:
            ax.text(0.5, 0.5, f"{name}\nFAILED: {exc}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="red")
            ax.axis("off")
            if verbose:
                print(f"    {name}: {exc}")
            continue
        _draw_gamma_edges(ax, pos, A, probes, gamma=6.0)
        ax.scatter(pos[:, 0], pos[:, 1], c=comm_cols, s=22,
                   edgecolors="white", linewidths=0.4, zorder=5)
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.axis("off")
        margin = 0.08
        xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
        ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
        span = max(xmax - xmin, ymax - ymin, 1e-9)
        ax.set_xlim(xmin - margin * span, xmax + margin * span)
        ax.set_ylim(ymin - margin * span, ymax + margin * span)
        if verbose:
            print(f"    {name}: ok")
    for i in range(len(variants), rows * cols):
        axes[i].axis("off")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"graph-tool layout gallery — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"nodes colored by LRG community (n={n_comm}); edges: red=same-probe, grey=cross",
        fontsize=13, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_fig(fig, output_dir / f"gt_layout_comm_{fc_method}_{patient}_{band}_{phase}")

    # ---- figure 2: same layouts but nodes colored by electrode shaft ----
    fig2, axes2 = plt.subplots(rows, cols, figsize=(4.2 * cols, 4.2 * rows))
    axes2 = axes2.flatten()
    for idx, (name, fn) in enumerate(variants):
        ax = axes2[idx]
        try:
            pos = fn()
        except Exception:
            ax.axis("off")
            continue
        _draw_gamma_edges(ax, pos, A, probes, gamma=6.0)
        ax.scatter(pos[:, 0], pos[:, 1], c=shaft_cols, s=22,
                   edgecolors="white", linewidths=0.4, zorder=5)
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.axis("off")
        margin = 0.08
        xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
        ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
        span = max(xmax - xmin, ymax - ymin, 1e-9)
        ax.set_xlim(xmin - margin * span, xmax + margin * span)
        ax.set_ylim(ymin - margin * span, ymax + margin * span)
    for i in range(len(variants), rows * cols):
        axes2[i].axis("off")

    fig2.suptitle(
        f"Same layouts, nodes colored by ELECTRODE SHAFT — {patient} {band_tex} {phase}\n"
        f"(shaft coloring shows whether modules are shaft-independent)",
        fontsize=13, fontweight="bold", y=1.00,
    )
    fig2.tight_layout()
    save_fig(fig2, output_dir / f"gt_layout_shaft_{fc_method}_{patient}_{band}_{phase}")


def main():
    parser = argparse.ArgumentParser()
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

    run_prototype(
        args.patient, args.band, args.phase, args.fc_method,
        args.output_dir, args.verbose,
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
