#!/usr/bin/env python3
"""diag_spring_k_sweep — sweep spring layout k over orders of magnitude.

One graph, raw weights (no rank-transform), spring layout with k
varied across orders of magnitude.  Each panel labelled with the k
value used and the resulting pairwise-distance CV (proxy for
"non-uniform spread").

Edges below ``--min-alpha`` are dropped from rendering only (the
layout still sees the full weighted graph) so the vector PDF stays
light-weight.

Run:
    conda activate lapbrain
    python scripts/01_compute/diagnostics/diag_spring_k_sweep.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import (
    _load_inputs,
    draw_gamma_edges,
    draw_nodes,
    shaft_colors,
)


def _spring(A: np.ndarray, *, k: float, iterations: int, seed: int = 42) -> np.ndarray:
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    pos = nx.spring_layout(
        G, k=k, iterations=iterations, seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


def _pairwise_distance_cv(pos: np.ndarray) -> float:
    diffs = pos[:, None, :] - pos[None, :, :]
    d = np.sqrt((diffs ** 2).sum(axis=-1))
    r, c = np.triu_indices(pos.shape[0], k=1)
    flat = d[r, c]
    if flat.mean() <= 0:
        return 0.0
    return float(flat.std() / flat.mean())


def _draw_panel(ax, pos, A, probes, *, title: str, min_alpha: float):
    # alpha_range floor = 0.0 so the weakest edges genuinely scale to
    # zero opacity; min_alpha then prunes them from the LineCollection.
    draw_gamma_edges(
        ax, pos, A,
        coloring="shaft", probes=probes,
        gamma=1.0, width_range=(0.15, 4.0), alpha_range=(0.0, 0.9),
        min_alpha=min_alpha,
    )
    draw_nodes(ax, pos, c=shaft_colors(probes), s=12)
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
        "--ks", default="0.001,0.01,0.1,1.0,10.0",
        help="Comma-separated k values (NetworkX absolute, not k/√N). "
             "Default spans 4 orders of magnitude.",
    )
    parser.add_argument("--iterations", type=int, default=600)
    parser.add_argument(
        "--min-alpha", type=float, default=0.05,
        help="Drop edges below this rendered alpha — keeps file size "
             "light without changing the layout.  Default 0.05.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "_diagnostic"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    np.fill_diagonal(A, 0.0)
    N = A.shape[0]

    ks = [float(s) for s in args.ks.split(",")]
    n = len(ks)

    fig, axes = plt.subplots(1, n, figsize=(2.4 * n, 2.6))
    if n == 1:
        axes = [axes]

    for k, ax in zip(ks, axes):
        pos = _spring(A, k=k, iterations=args.iterations)
        cv_d = _pairwise_distance_cv(pos)
        ttl = f"k = {k:g}\nCV(d) = {cv_d:.2f}"
        _draw_panel(ax, pos, A, probes, title=ttl, min_alpha=args.min_alpha)
        print(f"  k={k}: CV(d)={cv_d:.3f}")

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    head = (
        f"{args.patient}  ·  {band_tex}  ·  {args.phase}  ·  "
        f"{args.fc_method}    spring, iter={args.iterations}, "
        f"min-α(render)={args.min_alpha:g},  N={N}"
    )
    fig.text(0.5, 1.02, head, ha="center", va="bottom",
             fontsize=9, fontweight="bold")
    fig.tight_layout()

    out = (
        Path(args.out_dir)
        / f"spring_k_sweep_{args.patient}_{args.band}_{args.phase}_"
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
