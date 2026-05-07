#!/usr/bin/env python3
"""Rank-transform weight sweep for spring / KK on ImCoh.

The hypothesis: raw ImCoh has a narrow-range weight distribution, so
spring/KK can't separate modules. Rank-transforming the weights maps
them to exactly uniform over [0, 1] (by rank), then raising to a power α
redistributes mass toward the top — mimicking MSC's heavy-tailed shape
by construction.

Rows: spring k={0.01, 0.03, 0.1} and kamada-kawai.
Cols: α ∈ {1, 2, 4, 8, 16}.

Drawing uses the ORIGINAL weights so the visualization stays comparable.

Run:
  python scripts/10_notes_imcoh/fig_helper_rank_transform.py
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


def _draw_edges(ax, pos, A, probes, gamma_draw=6.0,
                width_range=(0.12, 3.5), alpha_range=(0.03, 0.9)):
    """Draw edges using ORIGINAL weights with gamma-scaled width/alpha."""
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


def rank_transform(A: np.ndarray, alpha: float) -> np.ndarray:
    """Replace upper-triangle weights with (rank / n_edges) ** alpha.
    Lower triangle mirrored. Diagonal zero. Preserves edge ordering,
    remaps distribution to uniform ranks on [0, 1], then amplifies top.
    """
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = w > 0
    if not active.any():
        return A.copy()
    # Rank among positive edges only
    w_pos = w[active]
    ranks = np.argsort(np.argsort(w_pos)).astype(float)
    ranks_norm = (ranks + 1) / len(ranks)  # in (0, 1]
    new_w = np.zeros_like(w)
    new_w[active] = ranks_norm ** alpha
    # Write back (upper triangle)
    B = np.zeros_like(A)
    B[r, c] = new_w
    B = B + B.T
    return B


def _spring_layout(A_rank, k, iterations=2000, seed=42):
    N = A_rank.shape[0]
    G = nx.from_numpy_array(A_rank)
    pos = nx.spring_layout(G, k=k, iterations=iterations, seed=seed,
                           weight="weight")
    return np.array([pos[i] for i in range(N)])


def _kk_layout(A_rank):
    N = A_rank.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    m = A_rank[r, c] > 0
    for i, j in zip(r[m], c[m]):
        G.add_edge(int(i), int(j),
                   weight=float(A_rank[i, j]),
                   distance=1.0 / (A_rank[i, j] + 1e-12))
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


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
        print("  SKIP")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    alphas = [1, 2, 4, 8, 16]
    layout_rows = [
        ("KK",          lambda Ar: _kk_layout(Ar)),
        ("spring k=0.01", lambda Ar: _spring_layout(Ar, 0.01)),
        ("spring k=0.03", lambda Ar: _spring_layout(Ar, 0.03)),
        ("spring k=0.1",  lambda Ar: _spring_layout(Ar, 0.1)),
    ]

    fig, axes = plt.subplots(len(layout_rows), len(alphas),
                             figsize=(2.6 * len(alphas),
                                      2.6 * len(layout_rows)))

    for i, (row_name, row_fn) in enumerate(layout_rows):
        for j, alpha in enumerate(alphas):
            ax = axes[i, j]
            A_rank = rank_transform(A, alpha)
            try:
                pos = row_fn(A_rank)
            except Exception as exc:
                ax.text(0.5, 0.5, f"FAIL\n{exc}", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="red")
                ax.axis("off")
                if verbose:
                    print(f"    {row_name} α={alpha}: {exc}")
                continue
            _draw_edges(ax, pos, A, probes)
            ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=14,
                       edgecolors="white", linewidths=0.3, zorder=5)
            if i == 0:
                ax.set_title(f"α={alpha}", fontsize=11, fontweight="bold")
            if j == 0:
                ax.text(-0.04, 0.5, row_name, rotation=90,
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=10, fontweight="bold")
            _finalize(ax, pos)
            if verbose:
                print(f"    {row_name} α={alpha}: ok")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Rank-transform sweep — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"layout input = (rank/E)^α; edges drawn with original weights",
        fontsize=12, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_helper_fig(
        fig,
        output_dir / f"rank_sweep_{fc_method}_{patient}_{band}_{phase}",
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
