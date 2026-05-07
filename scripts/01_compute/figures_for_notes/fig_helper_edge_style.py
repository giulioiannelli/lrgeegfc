#!/usr/bin/env python3
"""Edge-drawing style sweep with a LOCKED placement.

The placement (rank-transform α=1, spring k=0.01, many iterations) was
visually the most promising ImCoh layout so far. This script locks that
placement and sweeps only how edges are drawn: γ of the width/alpha curve,
the width range, the alpha range.

Layout input: the rank-transformed matrix A' = rank(A)/E (α=1).
Drawing input: the ORIGINAL weights A.

Grid:
  rows = width range ∈ {(0.1, 3), (0.2, 5), (0.3, 8)}
  cols = γ_draw   ∈ {4, 6, 8, 12, 20, 40, 80}
  alpha range fixed at (0.02, 0.9) for row 0 / (0.01, 0.95) row 1 / (0.005, 1.0) row 2

Run:
  python scripts/10_notes_imcoh/fig_helper_edge_style.py -v
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


def rank_transform(A: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = w > 0
    new_w = np.zeros_like(w)
    if active.any():
        ranks = np.argsort(np.argsort(w[active])).astype(float)
        new_w[active] = ((ranks + 1) / len(ranks)) ** alpha
    B = np.zeros_like(A)
    B[r, c] = new_w
    return B + B.T


def locked_layout(A: np.ndarray, k: float = 0.01, iterations: int = 2000):
    """Reproduce the exact (k=0.01, α=1, 2000 iter) placement that was visually
    promising in the earlier rank_transform sweep."""
    Ar = rank_transform(A, alpha=1.0)
    G = nx.from_numpy_array(Ar)
    pos = nx.spring_layout(G, k=k, iterations=iterations, seed=42,
                           weight="weight")
    return np.array([pos[i] for i in range(A.shape[0])])


def _draw_edges(ax, pos, A, probes, gamma, width_range, alpha_range):
    N = A.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = A[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    t = (w_act / w_act.max()) ** gamma
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
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    # LOCK the layout (once)
    print("  computing locked layout (rank α=1, spring k=0.01, 2000 iter)...")
    pos = locked_layout(A, k=0.01, iterations=2000)

    gammas = [4, 6, 8, 12, 20, 40, 80]
    width_ranges = [(0.10, 3.0), (0.20, 5.0), (0.30, 8.0)]
    alpha_ranges = [(0.02, 0.70), (0.01, 0.90), (0.005, 1.00)]

    fig, axes = plt.subplots(len(width_ranges), len(gammas),
                             figsize=(2.6 * len(gammas),
                                      2.6 * len(width_ranges)))

    for i, (wr, ar) in enumerate(zip(width_ranges, alpha_ranges)):
        for j, g in enumerate(gammas):
            ax = axes[i, j]
            _draw_edges(ax, pos, A, probes, gamma=g,
                        width_range=wr, alpha_range=ar)
            ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=14,
                       edgecolors="white", linewidths=0.3, zorder=5)
            if i == 0:
                ax.set_title(f"γ={g}", fontsize=11, fontweight="bold")
            if j == 0:
                ax.text(-0.05, 0.5,
                        f"w {wr[0]}-{wr[1]}\nα {ar[0]}-{ar[1]}",
                        rotation=90, transform=ax.transAxes,
                        ha="right", va="center", fontsize=9,
                        fontweight="bold")
            _finalize(ax, pos)
            if verbose:
                print(f"    wr={wr} γ={g}: ok")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Edge-drawing sweep (locked layout) — {fc_method.upper()} "
        f"{patient} {band_tex} {phase}\n"
        f"layout fixed: spring k=0.01, rank-α=1, 2000 iter. "
        f"Rows=(width,alpha) range; cols=γ.",
        fontsize=12, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_helper_fig(
        fig,
        output_dir / f"edge_style_{fc_method}_{patient}_{band}_{phase}",
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
