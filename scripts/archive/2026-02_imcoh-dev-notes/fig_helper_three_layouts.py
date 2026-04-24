#!/usr/bin/env python3
"""Three-panel layout comparison — FULL VECTOR PDF output.

Layouts are computed with networkx (spring, KK, KK-on-ultrametric).
Drawing is done with matplotlib but with the edge styling QUANTIZED into
a small number of bins, so the PDF stores O(n_bins) graphic primitives
instead of O(n_edges). That keeps the PDF compact AND fully vectorial
(no rasterized images).

Run:
  python scripts/10_notes_imcoh/fig_helper_three_layouts.py \
      --gamma 2.0 --spring-k 0.01 --spring-iter 1000
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

from scipy.spatial.distance import squareform

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


def _probe_color_map(probes):
    """Return {probe_id: rgba} using the same tab20 palette as the nodes."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return {p: cmap(i) for i, p in enumerate(uniq)}


# ---------------------------------------------------------------------------
# Binned vector edge drawing
# ---------------------------------------------------------------------------

def draw_edges_binned(ax, pos, A, gamma=2.0,
                      wr=(0.2, 1.5), ar=(0.08, 0.9),
                      n_bins: int = 6,
                      probes=None, highlight_same_probe: bool = False,
                      probe_color_map: dict | None = None,
                      base_color=(0.0, 0.0, 0.0)):
    """Rank-based width/alpha. Each edge k gets:
          t_k   = ( (rank_k + 1) / n_edges ) ** gamma
          width = wmin + t_k * (wmax - wmin)
          alpha = amin + t_k * (amax - amin)
    If *highlight_same_probe* and *probes* is given, same-probe edges are
    coloured red (keeping their rank-based width and alpha); cross-probe
    edges stay grey.
    """
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    ranks = np.argsort(np.argsort(w_act)).astype(float)
    t = ((ranks + 1) / len(ranks)) ** gamma
    widths = wr[0] + t * (wr[1] - wr[0])
    alphas = ar[0] + t * (ar[1] - ar[0])

    ii, jj = r[active], c[active]
    order = np.argsort(t)  # strong edges painted last
    segs = np.stack([pos[ii[order]], pos[jj[order]]], axis=1)
    widths_o = widths[order]
    alphas_o = alphas[order]

    # Base colour per edge: cross-probe = `base_color` (default black);
    # same-probe = that shaft's node colour.
    grey = np.array(base_color, dtype=float)
    if highlight_same_probe and probes is not None:
        same = np.array(
            [probes[ii[k]] == probes[jj[k]] for k in order], dtype=bool,
        )
    else:
        same = np.zeros(len(order), dtype=bool)

    # Compute RGB per edge
    rgb = np.tile(grey, (len(order), 1))
    if highlight_same_probe and probes is not None and probe_color_map is not None:
        for idx in np.where(same)[0]:
            k = order[idx]
            probe_id = probes[ii[k]]
            col = probe_color_map.get(probe_id)
            if col is not None:
                rgb[idx] = col[:3]

    # Same-probe edges keep rank-based width but use full alpha; cross-probe
    # edges fade with rank.
    alpha_out = np.where(same, 1.0, alphas_o)
    colors = np.concatenate([rgb, alpha_out[:, None]], axis=1)

    ax.add_collection(LineCollection(
        segs,
        colors=colors,
        linewidths=widths_o,
        zorder=1,
        rasterized=False,
    ))


# ---------------------------------------------------------------------------
# Layouts
# ---------------------------------------------------------------------------

def _seeded_initial_pos(N: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    theta = 2 * np.pi * np.arange(N) / N + rng.uniform(-0.2, 0.2, N)
    pos0 = np.column_stack([np.cos(theta), np.sin(theta)])
    return {i: pos0[i] for i in range(N)}


def layout_spring(A, k=None, iterations=1000, seed=42):
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    pos = nx.spring_layout(G, k=k, iterations=iterations, seed=seed,
                           weight="weight")
    return np.array([pos[i] for i in range(N)])


def layout_kk_invweight(A, seed=42):
    N = A.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    m = A[r, c] > 0
    for i, j in zip(r[m], c[m]):
        w = A[i, j]
        G.add_edge(int(i), int(j), distance=1.0 / (w + 1e-12))
    pos = nx.kamada_kawai_layout(G, weight="distance",
                                 pos=_seeded_initial_pos(N, seed))
    return np.array([pos[i] for i in range(N)])


def layout_kk_ultrametric(U, seed=42):
    N = U.shape[0]
    G = nx.complete_graph(N)
    dist = {i: {j: float(U[i, j]) for j in range(N) if j != i}
            for i in range(N)}
    pos = nx.kamada_kawai_layout(G, dist=dist,
                                 pos=_seeded_initial_pos(N, seed))
    return np.array([pos[i] for i in range(N)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--fc-method", default="imcoh_abs")
    parser.add_argument("--gamma", type=float, default=2.0,
                        help="rank-based edge scaling exponent")
    parser.add_argument("--spring-k", type=float, default=None)
    parser.add_argument("--spring-iter", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--wmin", type=float, default=0.2)
    parser.add_argument("--wmax", type=float, default=1.5)
    parser.add_argument("--amin", type=float, default=0.08)
    parser.add_argument("--amax", type=float, default=0.9)
    parser.add_argument("--n-bins", type=int, default=6,
                        help="number of (width, alpha) bins for edges")
    parser.add_argument("--highlight-same-probe", action="store_true",
                        help="Color same-probe edges red (keeping rank-based width/alpha)")
    parser.add_argument("--node-size", type=float, default=22.0,
                        help="Marker size for nodes (matplotlib s=...)")
    parser.add_argument("--edge-color", type=str, default="black",
                        help="Color for cross-probe (regular) edges. "
                             "Any matplotlib color name or #RRGGBB.")
    parser.add_argument("--output-dir", type=Path,
                        default=SECTION2_ROOT / "fig_helper")
    args = parser.parse_args()

    apply_pub_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mat = load_fc_matrix(args.patient, args.phase, args.band, args.fc_method)
    if mat is None:
        raise SystemExit("no matrix")
    lrg = load_lrg_result(args.patient, args.phase, args.band, args.fc_method)
    if lrg is None:
        raise SystemExit("no LRG")
    ch = load_channel_labels(args.patient)
    probes = extract_probe_labels(ch)
    N = min(mat.shape[0], lrg.n_nodes)
    A = np.abs(mat[:N, :N].copy())
    np.fill_diagonal(A, 0)
    probes = probes[:N]
    node_colors = _shaft_colors(probes)
    probe_cmap = _probe_color_map(probes)

    um = lrg.ultrametric_matrix
    U = squareform(um) if um.ndim == 1 else um
    U = U[:N, :N]

    print("  KK (1/w) ...")
    pos_kk = layout_kk_invweight(A, seed=args.seed)
    print(f"  spring (k={args.spring_k}, iter={args.spring_iter}) ...")
    pos_sp = layout_spring(A, k=args.spring_k, iterations=args.spring_iter,
                           seed=args.seed)
    print("  KK on LRG ultrametric ...")
    pos_lrg = layout_kk_ultrametric(U, seed=args.seed)

    band_tex = BRAIN_BAND_TEX_DICT[args.band]
    spring_k_str = (f"{args.spring_k:g}" if args.spring_k is not None
                    else "1/√N")
    labels = [
        f"KK (dist=1/w, seed={args.seed})",
        f"Spring (k={spring_k_str}, it={args.spring_iter}, seed={args.seed})",
        f"KK on LRG ultram. (seed={args.seed})",
    ]
    positions = [pos_kk, pos_sp, pos_lrg]

    panel_w, panel_h = 3.2, 3.4
    fig, axes = plt.subplots(1, 3,
                             figsize=(panel_w * 3, panel_h),
                             constrained_layout=True)

    wr, ar = (args.wmin, args.wmax), (args.amin, args.amax)
    for ax, pos, lab in zip(axes, positions, labels):
        import matplotlib.colors as _mcolors
        base_rgb = _mcolors.to_rgb(args.edge_color)
        draw_edges_binned(ax, pos, A, gamma=args.gamma,
                          wr=wr, ar=ar, n_bins=args.n_bins,
                          probes=probes,
                          highlight_same_probe=args.highlight_same_probe,
                          probe_color_map=probe_cmap,
                          base_color=base_rgb)
        ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=args.node_size,
                   edgecolors="white", linewidths=0.4, zorder=5)
        ax.set_title(f"{lab}    rank-γ={args.gamma}",
                     fontsize=9, fontweight="bold")
        ax.axis("off")
        m = 0.08
        xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
        ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
        span = max(xmax - xmin, ymax - ymin, 1e-9)
        ax.set_xlim(xmin - m * span, xmax + m * span)
        ax.set_ylim(ymin - m * span, ymax + m * span)

    out_stem = (args.output_dir
                / f"three_layouts_{args.fc_method}_{args.patient}_"
                  f"{args.band}_{args.phase}")
    save_fig(fig, out_stem, dpi=150, fmt="pdf")
    print("Done.")


if __name__ == "__main__":
    main()
