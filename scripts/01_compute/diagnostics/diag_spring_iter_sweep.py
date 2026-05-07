#!/usr/bin/env python3
"""diag_spring_iter_sweep — sweep spring iterations at fixed k.

Fixes a single k value and varies the iteration count to see
whether the layout converges to richer structure given more time.

Run:
    conda activate lapbrain
    python scripts/01_compute/diagnostics/diag_spring_iter_sweep.py \
        --k 0.01 --iterations 100,500,2000,10000,50000
"""
from __future__ import annotations

import argparse
import time
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
    parser.add_argument("--k", type=float, default=0.01,
                        help="Fixed spring k.  Default 0.01.")
    parser.add_argument(
        "--iterations", default="100,500,2000,10000,50000",
        help="Comma-separated iteration counts to sweep.",
    )
    parser.add_argument(
        "--min-alpha", type=float, default=0.10,
        help="Drop edges below this rendered alpha for file weight.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "_diagnostic"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    np.fill_diagonal(A, 0.0)
    N = A.shape[0]

    iters = [int(s) for s in args.iterations.split(",")]
    n = len(iters)
    fig, axes = plt.subplots(1, n, figsize=(2.4 * n, 2.6))
    if n == 1:
        axes = [axes]

    for it, ax in zip(iters, axes):
        t0 = time.perf_counter()
        pos = _spring(A, k=args.k, iterations=it)
        elapsed = time.perf_counter() - t0
        cv_d = _pairwise_distance_cv(pos)
        ttl = f"iter = {it}\nCV(d) = {cv_d:.2f}   ({elapsed:.1f}s)"
        _draw_panel(ax, pos, A, probes, title=ttl, min_alpha=args.min_alpha)
        print(f"  iter={it:>6d}: CV(d)={cv_d:.3f}, {elapsed:.1f}s")

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    head = (
        f"{args.patient}  ·  {band_tex}  ·  {args.phase}  ·  "
        f"{args.fc_method}    spring, k={args.k:g}, "
        f"min-α(render)={args.min_alpha:g},  N={N}"
    )
    fig.text(0.5, 1.02, head, ha="center", va="bottom",
             fontsize=9, fontweight="bold")
    fig.tight_layout()

    out = (
        Path(args.out_dir)
        / f"spring_iter_sweep_{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}_k{args.k:g}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {out}")
    size_kb = out.stat().st_size / 1024
    print(f"size: {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
