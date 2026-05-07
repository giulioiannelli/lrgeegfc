#!/usr/bin/env python3
"""diag_edge_suppression — pick a scaling that reveals structure.

Step 1: histogram + quantile annotations of the upper-triangle edge
weight distribution.

Step 2: same KK layout drawn with five suppression rules, so we can
choose visually:

    A. linear:        t = w/w_max
    B. power γ=4:     t = (w/w_max)^4
    C. power γ=8:     t = (w/w_max)^8
    D. Q75 threshold: t = max(0, w-Q75) / (w_max-Q75)
    E. Q90 threshold: t = max(0, w-Q90) / (w_max-Q90)

Width and alpha both scale with t.  Edges with t=0 are not drawn at all
(both width and alpha effectively vanish).

Run:
    conda activate lapbrain
    python scripts/01_compute/diagnostics/diag_edge_suppression.py
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from matplotlib.collections import LineCollection

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import (
    _load_inputs,
    draw_nodes,
    shaft_colors,
    _shaft_color_map,
    CROSS_PROBE_RGB,
)


def _kk_layout(A: np.ndarray) -> np.ndarray:
    G = nx.from_numpy_array(A)
    for _, _, d in G.edges(data=True):
        d["distance"] = 1.0 / (d.get("weight", 1.0) + 1e-6)
    pos_dict = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos_dict[i] for i in range(A.shape[0])])


def _draw_with_t(ax, pos, A, probes, t_vec, *,
                 width_range=(0.15, 4.0), alpha_range=(0.0, 0.95)):
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    keep = (w > 0) & (t_vec > 0.0)
    if not keep.any():
        return
    rr = r[keep]
    cc = c[keep]
    t_kept = t_vec[keep]
    # Draw strong edges last so they sit on top.
    order = np.argsort(t_kept)
    rr = rr[order]
    cc = cc[order]
    t_kept = t_kept[order]

    widths = width_range[0] + t_kept * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t_kept * (alpha_range[1] - alpha_range[0])

    shaft_map = _shaft_color_map(probes)
    segments = []
    colors = []
    linew = []
    for k in range(len(t_kept)):
        i, j = int(rr[k]), int(cc[k])
        segments.append([pos[i], pos[j]])
        linew.append(float(widths[k]))
        if probes[i] == probes[j]:
            rgba = shaft_map[probes[i]]
            colors.append((rgba[0], rgba[1], rgba[2], float(alphas[k])))
        else:
            colors.append((*CROSS_PROBE_RGB, float(alphas[k])))
    ax.add_collection(LineCollection(segments, colors=colors,
                                     linewidths=linew, zorder=1))


def _finalize(ax, pos, *, title: str):
    ax.axis("off")
    ax.set_title(title, fontsize=8)
    if pos.size:
        m = 0.08
        xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
        ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
        span = max(xmax - xmin, ymax - ymin, 1e-9)
        ax.set_xlim(xmin - m * span, xmax + m * span)
        ax.set_ylim(ymin - m * span, ymax + m * span)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--fc-method", default="imcoh_abs")
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "_diagnostic"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    np.fill_diagonal(A, 0.0)
    N = A.shape[0]

    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    pos_w = w[w > 0]
    w_max = pos_w.max()

    # Quantiles of the positive-weight upper-triangle edges.
    q50, q75, q90, q95, q99 = np.quantile(pos_w, [0.5, 0.75, 0.9, 0.95, 0.99])
    print(f"N = {N},  edges = {pos_w.size}")
    print(f"  w_min = {pos_w.min():.4f}  w_max = {w_max:.4f}")
    print(f"  Q50={q50:.4f}  Q75={q75:.4f}  Q90={q90:.4f}  "
          f"Q95={q95:.4f}  Q99={q99:.4f}")
    print(f"  mean = {pos_w.mean():.4f}  std = {pos_w.std():.4f}  "
          f"CV = {pos_w.std()/pos_w.mean():.3f}")

    pos = _kk_layout(A)

    # Define the 5 suppression rules as functions on full w array.
    def rule_linear(w):
        return w / w_max

    def rule_gamma(g):
        return lambda w: (w / w_max) ** g

    def rule_threshold(thr):
        return lambda w: np.where(
            w > thr,
            (w - thr) / (w_max - thr + 1e-12),
            0.0,
        )

    rules = [
        ("A. linear  t=w/w_max",            rule_linear),
        ("B. power γ=4",                     rule_gamma(4.0)),
        ("C. power γ=8",                     rule_gamma(8.0)),
        ("D. threshold Q75, linear above",   rule_threshold(q75)),
        ("E. threshold Q90, linear above",   rule_threshold(q90)),
    ]

    fig = plt.figure(figsize=(15, 7.5))
    gs = fig.add_gridspec(2, 3)

    # Top-left: histogram with quantile lines.
    ax_hist = fig.add_subplot(gs[0, 0])
    ax_hist.hist(pos_w, bins=80, color="steelblue", edgecolor="white", linewidth=0.3)
    for q, lbl in [(q50, "Q50"), (q75, "Q75"), (q90, "Q90"),
                   (q95, "Q95"), (q99, "Q99")]:
        ax_hist.axvline(q, color="crimson", linestyle="--", linewidth=0.8)
        ax_hist.text(q, ax_hist.get_ylim()[1] * 0.95, lbl,
                     rotation=90, fontsize=6, color="crimson",
                     ha="right", va="top")
    ax_hist.set_xlabel(r"$|w_{ij}|$", fontsize=9)
    ax_hist.set_ylabel("count", fontsize=9)
    ax_hist.set_title(
        f"weight distribution    n_edges={pos_w.size},  "
        f"CV={pos_w.std()/pos_w.mean():.2f}",
        fontsize=8,
    )
    ax_hist.tick_params(labelsize=7)

    # Five suppression rule panels — order: top row B/C, bottom row A/D/E.
    panel_axes = [
        fig.add_subplot(gs[0, 1]),  # rule A (linear current)
        fig.add_subplot(gs[0, 2]),  # rule B (γ=4)
        fig.add_subplot(gs[1, 0]),  # rule C (γ=8)
        fig.add_subplot(gs[1, 1]),  # rule D (Q75)
        fig.add_subplot(gs[1, 2]),  # rule E (Q90)
    ]
    for ax, (label, fn) in zip(panel_axes, rules):
        t = fn(w)
        n_drawn = int((t > 0.0).sum())
        ttl = f"{label}\n{n_drawn} edges drawn"
        _draw_with_t(ax, pos, A, probes, t)
        draw_nodes(ax, pos, c=shaft_colors(probes), s=10)
        _finalize(ax, pos, title=ttl)
        print(f"  {label}: {n_drawn} edges drawn")

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    head = (
        f"{args.patient}  ·  {band_tex}  ·  {args.phase}  ·  "
        f"{args.fc_method}    KK layout, edge-color: shaft"
    )
    fig.text(0.5, 0.99, head, ha="center", va="top",
             fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out = (
        Path(args.out_dir)
        / f"edge_suppression_{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {out}")
    size_kb = out.stat().st_size / 1024
    print(f"size: {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
