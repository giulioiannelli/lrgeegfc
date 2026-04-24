#!/usr/bin/env python3
"""LRG-seeded SFDP layout prototype using graph-tool.

Idea: SFDP is a multi-level force-directed algorithm that accepts a per-
vertex ``groups`` property. Vertices sharing a group attract each other
more strongly, which forces modules apart on the canvas. We feed the LRG
dendrogram cut (at several n_communities) as groups.

Graph-tool is used for both the graph object and the layout — faster and
better suited to dense weighted graphs than NetworkX.

Run:
  python scripts/10_notes_imcoh/fig_helper_lrg_sfdp.py [-v]
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

from lrg_eegfc.utils.metrics.hypothesis import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_helper_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result


def _shaft_colors(probes):
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _comm_colors(labels):
    K = max(int(np.max(labels)), 2)
    cmap = plt.get_cmap("tab20", K)
    return [cmap(int(l) - 1) for l in labels]


def _draw_edges(ax, pos, A, probes, gamma_draw=6.0,
                width_range=(0.12, 3.5), alpha_range=(0.03, 0.9)):
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


def sfdp_lrg(A, groups=None, gamma=0.1, mu=0.0, C=0.2, p=2.0):
    g, ew = _gt_graph(A)
    grp_prop = None
    if groups is not None:
        grp_prop = g.new_vertex_property("int")
        grp_prop.a = np.asarray(groups, dtype=int)
    pos = gt.sfdp_layout(g, eweight=ew, groups=grp_prop,
                         gamma=gamma, mu=mu, C=C, p=p, max_iter=0)
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
    lrg = load_lrg_result(patient, phase, band, fc_method)
    if lrg is None:
        print("  SKIP: no LRG")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    N = min(mat.shape[0], lrg.n_nodes)
    A = np.abs(mat[:N, :N].copy())
    np.fill_diagonal(A, 0.0)
    probes = probes[:N]

    # Rows: n_communities from LRG cut; Cols: gamma (SFDP module-separation knob)
    n_comm_list = [5, 10, 15, 25]
    gamma_vals = [0.03, 0.08, 0.15, 0.3, 0.6]

    # Two figures: one colored by LRG community, one by shaft
    for color_by, color_fn_name in [("comm", "comm"), ("shaft", "shaft")]:
        fig, axes = plt.subplots(len(n_comm_list), len(gamma_vals),
                                 figsize=(2.6 * len(gamma_vals),
                                          2.6 * len(n_comm_list)))
        for i, nc in enumerate(n_comm_list):
            comm = fcluster(lrg.linkage_matrix, t=nc,
                            criterion="maxclust")[:N]
            for j, gamma in enumerate(gamma_vals):
                ax = axes[i, j]
                try:
                    pos = sfdp_lrg(A, groups=comm, gamma=gamma)
                except Exception as exc:
                    ax.text(0.5, 0.5, f"FAIL\n{exc}", ha="center", va="center",
                            transform=ax.transAxes, fontsize=8, color="red")
                    ax.axis("off")
                    continue
                _draw_edges(ax, pos, A, probes)
                if color_by == "comm":
                    node_colors = _comm_colors(comm)
                else:
                    node_colors = _shaft_colors(probes)
                ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=15,
                           edgecolors="white", linewidths=0.3, zorder=5)
                if i == 0:
                    ax.set_title(f"γ={gamma}", fontsize=11, fontweight="bold")
                if j == 0:
                    ax.text(-0.04, 0.5, f"n_comm={nc}", rotation=90,
                            transform=ax.transAxes, ha="right", va="center",
                            fontsize=10, fontweight="bold")
                _finalize(ax, pos)
                if verbose:
                    print(f"    n_comm={nc} γ={gamma} ({color_by}): ok")

        band_tex = BRAIN_BAND_TEX_DICT[band]
        fig.suptitle(
            f"LRG-seeded SFDP — {fc_method.upper()} {patient} {band_tex} {phase}\n"
            f"rows=LRG cut, cols=γ (SFDP module-separation). "
            f"Nodes colored by {'LRG community' if color_by == 'comm' else 'electrode shaft'}.",
            fontsize=12, fontweight="bold", y=1.00,
        )
        fig.tight_layout()
        save_helper_fig(
            fig,
            output_dir
            / f"lrg_sfdp_{color_by}_{fc_method}_{patient}_{band}_{phase}",
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
