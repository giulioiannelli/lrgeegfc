#!/usr/bin/env python3
"""Layout diagnostic: do strong edges connect nearby nodes?

For each candidate layout recipe we draw the SAME node positions twice:

  left column : only the top 5% edges (strong-link test)
  right column: all edges with γ-scaled styling

A good layout places the endpoints of the top-5% edges close together on
the canvas. A bad one shows them crossing the canvas (the current failure
mode of spring k=0.01 — too little repulsion, nodes frozen at seed).

Candidates tested:
  - KK on rank(A)^α, α ∈ {1, 2, 4}
  - spring k ∈ {0.2, 0.5, 1.0} on rank(A)^1, 2000 iter
  - spring k ∈ {0.2, 0.5} on rank(A)^2, 2000 iter

Run:
  python scripts/10_notes_imcoh/fig_helper_layout_diagnostic.py -v
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

from _shared import (
    BRAIN_BAND_TEX_DICT,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_helper_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix


def _shaft_colors(probes):
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _top_frac(A: np.ndarray, frac: float) -> np.ndarray:
    """Keep only the top `frac` fraction of upper-triangle edges."""
    r, c = np.triu_indices(A.shape[0], k=1)
    w = A[r, c]
    pos = w > 0
    if not pos.any():
        return A.copy()
    th = np.percentile(w[pos], 100 * (1 - frac))
    B = A.copy()
    B[B < th] = 0.0
    return B


def _boost_top(A: np.ndarray, frac: float, factor: float) -> np.ndarray:
    """Multiply edges above the (1-frac) percentile by `factor`. All weak
    edges are preserved so every node stays connected → no periphery ring.
    """
    r, c = np.triu_indices(A.shape[0], k=1)
    w = A[r, c]
    pos = w > 0
    if not pos.any():
        return A.copy()
    th = np.percentile(w[pos], 100 * (1 - frac))
    B = A.copy()
    mask = B >= th
    B = np.where(mask, B * factor, B)
    return B


def rank_transform(A, alpha):
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = w > 0
    new = np.zeros_like(w)
    if active.any():
        rk = np.argsort(np.argsort(w[active])).astype(float)
        new[active] = ((rk + 1) / len(rk)) ** alpha
    B = np.zeros_like(A)
    B[r, c] = new
    return B + B.T


def layout_spring(A, k, iters=2000, seed=42):
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    pos = nx.spring_layout(G, k=k, iterations=iters, seed=seed,
                           weight="weight")
    return np.array([pos[i] for i in range(N)])


def layout_kk(A):
    """KK with distance = 1 / weight (upper triangle)."""
    N = A.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    m = A[r, c] > 0
    for i, j in zip(r[m], c[m]):
        w = A[i, j]
        G.add_edge(int(i), int(j), weight=float(w),
                   distance=1.0 / (w + 1e-12))
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


def _draw_top_edges_only(ax, pos, A, probes, keep_frac=0.05,
                         lw=1.4, alpha=0.85):
    """Only the top keep_frac of edges."""
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    if not (w > 0).any():
        return
    th = np.percentile(w[w > 0], 100 * (1 - keep_frac))
    idx = np.where(w >= th)[0]
    segs, cols = [], []
    for k in idx:
        i, j = r[k], c[k]
        segs.append([pos[i], pos[j]])
        if probes[i] == probes[j]:
            cols.append((0.8, 0.2, 0.2, alpha))
        else:
            cols.append((0.10, 0.35, 0.75, alpha))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw,
                                     zorder=1, rasterized=True))


def _draw_all_edges(ax, pos, A, probes, gamma=1.5,
                    wr=(0.2, 4.0), ar=(0.08, 0.9)):
    """Rank-based width/alpha scaling so the visible gradient does NOT depend
    on whether the weight distribution is heavy-tailed or near-uniform.

    t_k = rank(w_k) / n_edges  (uniform in (0, 1] by construction)
    width = wr[0] + (t**γ) * (wr[1] - wr[0])
    alpha = ar[0] + (t**γ) * (ar[1] - ar[0])
    """
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    # Rank-based t, uniform in (0, 1]
    ranks = np.argsort(np.argsort(w_act)).astype(float)
    t = ((ranks + 1) / len(ranks)) ** gamma
    widths = wr[0] + t * (wr[1] - wr[0])
    alphas = ar[0] + t * (ar[1] - ar[0])
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


def _finalize(ax, pos):
    ax.axis("off")
    m = 0.08
    xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
    ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    ax.set_xlim(xmin - m * span, xmax + m * span)
    ax.set_ylim(ymin - m * span, ymax + m * span)


def run_sweep(patient, band, phase, fc_method, output_dir, verbose=False):
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    # Full graph layout with top-5% edges boosted by a factor M, so strong
    # edges dominate placement but weak edges still attract isolated nodes.
    recipes = []
    for M in [5, 20, 100, 500]:
        Aboost = _boost_top(A, frac=0.05, factor=M)
        recipes.append((f"KK on A(top5%·{M})",
                        lambda Ab=Aboost: layout_kk(Ab)))
        recipes.append((f"spring k=0.5 on A(top5%·{M})",
                        lambda Ab=Aboost: layout_spring(Ab, 0.5)))

    # Column layout: [placement test: top-5% only] + [γ sweep for drawing]
    gammas = [1.0, 1.5, 2.0, 3.0]
    ncols = 1 + len(gammas)
    fig, axes = plt.subplots(len(recipes), ncols,
                             figsize=(2.8 * ncols, 2.8 * len(recipes)))

    for i, (name, fn) in enumerate(recipes):
        try:
            pos = fn()
        except Exception as exc:
            for c_idx in range(ncols):
                axes[i, c_idx].text(0.5, 0.5, f"FAIL\n{exc}",
                                    ha="center", va="center",
                                    transform=axes[i, c_idx].transAxes,
                                    fontsize=8, color="red")
                axes[i, c_idx].axis("off")
            continue

        # Col 0: top 5% only (placement test)
        ax0 = axes[i, 0]
        _draw_top_edges_only(ax0, pos, A, probes, keep_frac=0.05)
        ax0.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=16,
                    edgecolors="white", linewidths=0.3, zorder=5)
        ax0.set_title(f"{name}\ntop 5% only", fontsize=9, fontweight="bold")
        _finalize(ax0, pos)

        # Col 1..: all edges with rank-based γ sweep
        for j, g in enumerate(gammas):
            ax = axes[i, 1 + j]
            _draw_all_edges(ax, pos, A, probes, gamma=g,
                            wr=(0.15, 5.0), ar=(0.06, 0.95))
            ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=16,
                       edgecolors="white", linewidths=0.3, zorder=5)
            ax.set_title(f"γ={g}", fontsize=10, fontweight="bold")
            _finalize(ax, pos)

        if verbose:
            print(f"    {name}: ok")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Layout diagnostic — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"left = top 5% edges only (strong-link placement test); "
        f"right = all edges. Nodes by shaft.",
        fontsize=12, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_helper_fig(
        fig,
        output_dir / f"layout_diag_{fc_method}_{patient}_{band}_{phase}",
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
