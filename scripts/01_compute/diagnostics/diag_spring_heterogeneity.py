#!/usr/bin/env python3
"""diag_spring_heterogeneity — algorithmic k+α tuning for spring layout.

Spring layout produces a uniform bubble on a fully-connected
near-uniform-weight graph because pairwise equilibrium distances all
collapse to ~k/√w_ij when w_ij is roughly constant.  This script
demonstrates the algorithmic fix:

    1. Diagnose weight heterogeneity via CV(w) = std(w)/mean(w).
    2. If CV(w) < 1, rank-transform the weights:
       w' = (rank(w) / E)^α.  Sweep α to amplify heterogeneity.
    3. Run spring on the transformed weights; draw with the original
       weights so the visualization stays comparable.
    4. Pick (α, k) by reading off the resulting node-distance CV — a
       proxy for "how non-uniform is the layout".

Output: one PDF, 2 rows × 5 cols, sweeping α ∈ {1, 2, 4, 8, 16} at
two k regimes.  Each panel annotated with (CV(w'), CV(d)).

Run:
    conda activate lapbrain
    python scripts/01_compute/diagnostics/diag_spring_heterogeneity.py
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
    rank_transform,
    shaft_colors,
)


def _weight_cv(A: np.ndarray) -> float:
    """Coefficient of variation of upper-triangle weights."""
    r, c = np.triu_indices(A.shape[0], k=1)
    w = A[r, c]
    w = w[w > 0]
    if w.size == 0 or w.mean() <= 0:
        return 0.0
    return float(w.std() / w.mean())


def _pairwise_distance_cv(pos: np.ndarray) -> float:
    """CV of pairwise Euclidean distances between layout positions.
    A bubble has CV(d) ≈ 0.2-0.3; a heterogeneous layout has CV(d) > 0.5.
    """
    N = pos.shape[0]
    diffs = pos[:, None, :] - pos[None, :, :]
    d = np.sqrt((diffs ** 2).sum(axis=-1))
    r, c = np.triu_indices(N, k=1)
    flat = d[r, c]
    if flat.mean() <= 0:
        return 0.0
    return float(flat.std() / flat.mean())


def _spring(A: np.ndarray, *, k_base: float, iterations: int, seed: int = 42) -> np.ndarray:
    N = A.shape[0]
    G = nx.from_numpy_array(A)
    pos = nx.spring_layout(
        G, k=k_base / np.sqrt(N), iterations=iterations,
        seed=seed, weight="weight",
    )
    return np.array([pos[i] for i in range(N)])


def _draw_panel(ax, pos, A_orig, probes, *, title: str):
    draw_gamma_edges(
        ax, pos, A_orig,
        coloring="shaft", probes=probes,
        gamma=1.0, width_range=(0.15, 4.0), alpha_range=(0.03, 0.9),
    )
    draw_nodes(ax, pos, c=shaft_colors(probes), s=14)
    ax.axis("off")
    ax.set_title(title, fontsize=8, fontweight="normal")
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
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument(
        "--alphas", default="1,2,4,8,16",
        help="Comma-separated rank-transform α values.",
    )
    parser.add_argument(
        "--k-bases", default="1.0,5.0",
        help="Comma-separated k_base values; spring k = k_base/√N.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "_diagnostic"),
    )
    args = parser.parse_args()

    A, probes = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    np.fill_diagonal(A, 0.0)
    cv_raw = _weight_cv(A)
    print(f"raw CV(w) = {cv_raw:.3f}")

    alphas = [float(s) for s in args.alphas.split(",")]
    k_bases = [float(s) for s in args.k_bases.split(",")]

    n_rows = len(k_bases)
    n_cols = len(alphas)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.8 * n_cols, 2.8 * n_rows),
    )
    if n_rows == 1:
        axes = np.atleast_2d(axes)

    for i, k_base in enumerate(k_bases):
        for j, alpha in enumerate(alphas):
            ax = axes[i, j]
            if alpha == 1.0:
                A_layout = A.copy()
            else:
                A_layout = rank_transform(A, alpha=alpha)
            cv_w = _weight_cv(A_layout)

            pos = _spring(
                A_layout,
                k_base=k_base, iterations=args.iterations,
            )
            cv_d = _pairwise_distance_cv(pos)

            label_alpha = "raw" if alpha == 1.0 else f"α={alpha:g}"
            ttl = (
                f"{label_alpha}, k={k_base:g}/√N\n"
                f"CV(w')={cv_w:.2f}  CV(d)={cv_d:.2f}"
            )
            _draw_panel(ax, pos, A, probes, title=ttl)
            print(f"  {label_alpha} k_base={k_base}: CV(w')={cv_w:.3f} CV(d)={cv_d:.3f}")

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    head = (
        f"{args.patient}  ·  {band_tex}  ·  {args.phase}  ·  {args.fc_method}    "
        f"raw CV(w)={cv_raw:.3f},  iter={args.iterations}"
    )
    fig.text(
        0.5, 0.99, head,
        ha="center", va="top", fontsize=10, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out = (
        Path(args.out_dir)
        / f"spring_heterogeneity_{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
