#!/usr/bin/env python3
"""Layout-on-backbone sweep using graph-tool.

Principle: dense near-uniform weighted graphs (like ImCoh) produce hairballs
because spring/FR see too many similar attractions. Standard technique in
force-directed graph drawing: compute the layout on a SPARSE BACKBONE of
the strongest / most significant edges, then DRAW the full weighted graph
using the backbone-derived positions.

This is NOT display-thresholding. Every edge is drawn. Only the layout
input is sparse.

Rows: backbone method × strength
  - top  5%    (top edges by weight)
  - top 10%
  - top 20%
  - disparity filter α=0.05 (null-model backbone, Serrano et al.)

Cols: layout algorithm
  - FR  (Fruchterman-Reingold,  graph-tool)
  - ARF (attractive-repulsive,  graph-tool)
  - KK  (Kamada-Kawai,          graph-tool)

Drawing uses γ-scaled full-weight edges.

Run:
  python scripts/10_notes_imcoh/fig_helper_backbone_layout.py -v
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

import graph_tool.all as gt

from _shared import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_helper_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.fc.msc.sparsify import disparity_filter


def _shaft_colors(probes):
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _draw_edges(ax, pos, A, probes, gamma_draw=6.0,
                width_range=(0.12, 3.5), alpha_range=(0.03, 0.9)):
    """Draw edges of the FULL weighted graph with gamma-scaled visuals."""
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    t = (w_act / w_act.max()) ** gamma_draw
    widths = width_range[0] + t * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])
    order = np.argsort(w_act)
    segs, cols, lws = [], [], []
    for oi in order:
        idx = active[oi]
        i, j = r[idx], c[idx]
        segs.append([pos[i], pos[j]])
        lws.append(widths[oi])
        if probes[i] == probes[j]:
            cols.append((0.8, 0.2, 0.2, alphas[oi]))
        else:
            cols.append((0.35, 0.35, 0.35, alphas[oi]))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws,
                                     zorder=1, rasterized=True))


def _backbone_topk(A: np.ndarray, keep_frac: float) -> np.ndarray:
    """Keep top-(keep_frac) edges by weight."""
    r, c = np.triu_indices(A.shape[0], k=1)
    w = A[r, c]
    pos_mask = w > 0
    if not pos_mask.any():
        return A.copy()
    th = np.percentile(w[pos_mask], 100 * (1 - keep_frac))
    B = A.copy()
    B[B < th] = 0.0
    return B


def _backbone_disparity(A: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    return disparity_filter(A.copy(), alpha=alpha)


def _gt_graph(A: np.ndarray):
    N = A.shape[0]
    g = gt.Graph(directed=False)
    g.add_vertex(N)
    ew = g.new_edge_property("double")
    r, c = np.triu_indices(N, k=1)
    mask = A[r, c] > 0
    if mask.any():
        g.add_edge_list(np.column_stack([r[mask], c[mask]]))
        ew.a = A[r[mask], c[mask]]
    return g, ew


def layout_fr(A: np.ndarray, n_iter: int = 300) -> np.ndarray:
    g, ew = _gt_graph(A)
    pos = gt.fruchterman_reingold_layout(g, weight=ew, n_iter=n_iter)
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


def layout_arf(A: np.ndarray, a: float = 10.0, d: float = 0.5) -> np.ndarray:
    g, ew = _gt_graph(A)
    pos = gt.arf_layout(g, weight=ew, a=a, d=d, max_iter=0)
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


def layout_kk(A: np.ndarray) -> np.ndarray:
    """Kamada-Kawai via graph-tool's radial_tree substitute: use
    sfdp with gamma high and no groups as a KK-style deterministic fallback,
    since graph-tool does not ship KK. Keeps tool consistent."""
    # graph-tool's closest analogue: radial_tree_layout requires a tree.
    # Fall back to sfdp without groups — still FR-family deterministic modulo seed.
    g, ew = _gt_graph(A)
    pos = gt.sfdp_layout(g, eweight=ew, gamma=2.0, C=0.3, p=2.0, max_iter=0)
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)


def _finalize(ax, pos):
    ax.axis("off")
    m = 0.08
    xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
    ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    ax.set_xlim(xmin - m * span, xmax + m * span)
    ax.set_ylim(ymin - m * span, ymax + m * span)


def run_sweep(
    patient: str, band: str, phase: str, fc_method: str,
    output_dir: Path, verbose: bool = False,
):
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        print("  SKIP: no matrix")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0.0)

    backbones = [
        ("top 5%",      lambda: _backbone_topk(A, 0.05)),
        ("top 10%",     lambda: _backbone_topk(A, 0.10)),
        ("top 20%",     lambda: _backbone_topk(A, 0.20)),
        ("disparity α=0.05", lambda: _backbone_disparity(A, alpha=0.05)),
    ]
    layouts = [
        ("FR",   layout_fr),
        ("ARF",  layout_arf),
        ("SFDP (KK-ish)", layout_kk),
    ]

    fig, axes = plt.subplots(len(backbones), len(layouts),
                             figsize=(3.0 * len(layouts),
                                      3.0 * len(backbones)))

    for i, (bb_name, bb_fn) in enumerate(backbones):
        B = bb_fn()
        n_edges = int((B > 0).sum() // 2)
        for j, (lay_name, lay_fn) in enumerate(layouts):
            ax = axes[i, j]
            try:
                pos = lay_fn(B)
            except Exception as exc:
                ax.text(0.5, 0.5, f"FAIL\n{exc}", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="red")
                ax.axis("off")
                if verbose:
                    print(f"    {bb_name} / {lay_name}: {exc}")
                continue
            # Draw full graph, not the backbone
            _draw_edges(ax, pos, A, probes)
            ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=18,
                       edgecolors="white", linewidths=0.3, zorder=5)
            if i == 0:
                ax.set_title(lay_name, fontsize=12, fontweight="bold")
            if j == 0:
                ax.text(-0.05, 0.5, f"{bb_name}\n({n_edges} edges)",
                        rotation=90, transform=ax.transAxes,
                        ha="right", va="center",
                        fontsize=10, fontweight="bold")
            _finalize(ax, pos)
            if verbose:
                print(f"    {bb_name} / {lay_name}: ok ({n_edges} edges)")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Layout-on-backbone — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"layout sees sparse backbone, drawing shows full weighted graph. "
        f"Nodes by shaft.",
        fontsize=12, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_helper_fig(
        fig,
        output_dir / f"backbone_layout_{fc_method}_{patient}_{band}_{phase}",
    )


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
    run_sweep(args.patient, args.band, args.phase, args.fc_method,
              args.output_dir, args.verbose)
    print("Done.")


if __name__ == "__main__":
    main()
