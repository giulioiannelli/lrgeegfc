#!/usr/bin/env python3
"""Weight-transform sweep for spring / kamada-kawai on ImCoh networks.

Hypothesis: ImCoh weights are near-uniform, so spring & KK can't separate
modules because all attractions are ~equal. Raising the weights to a power
α amplifies strong edges and suppresses weak ones, effectively sharpening
the attraction landscape.

We sweep α ∈ {1, 2, 4, 8, 16, 32} against:
  - spring (FR) with 3 representative k values
  - kamada-kawai (deterministic) with d = 1 / w^α as distance

Only ONE case is produced: Pat_05 β rest_pre ImCoh.
Edges are drawn the same way in every panel so we compare layouts.

Run:
  python scripts/10_notes_imcoh/fig_helper_weight_transform.py [-v]
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

from lrg_eegfc.utils.metrics.hypothesis import (
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
    """Draw edges with gamma-scaled width/alpha on the ORIGINAL weights
    (so the visualization doesn't change between panels)."""
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


def _finalize(ax, pos, title):
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.axis("off")
    m = 0.08
    xmin, xmax = pos[:, 0].min(), pos[:, 0].max()
    ymin, ymax = pos[:, 1].min(), pos[:, 1].max()
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    ax.set_xlim(xmin - m * span, xmax + m * span)
    ax.set_ylim(ymin - m * span, ymax + m * span)


def _spring_layout(A_powered, k, iterations=800, seed=42):
    """Spring layout with an ABSOLUTE k value (not k_base/sqrt(N))."""
    N = A_powered.shape[0]
    G = nx.from_numpy_array(A_powered)
    pos = nx.spring_layout(G, k=k, iterations=iterations,
                           seed=seed, weight="weight")
    return np.array([pos[i] for i in range(N)])


def _kk_layout(A_powered):
    """Kamada-Kawai with distance = 1 / w (powered weights → powered distances)."""
    N = A_powered.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))
    r, c = np.triu_indices(N, k=1)
    m = A_powered[r, c] > 0
    if m.any():
        for i, j in zip(r[m], c[m]):
            G.add_edge(int(i), int(j),
                       weight=float(A_powered[i, j]),
                       distance=1.0 / (A_powered[i, j] + 1e-12))
    pos = nx.kamada_kawai_layout(G, weight="distance")
    return np.array([pos[i] for i in range(N)])


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
    np.fill_diagonal(A, 0)
    # Normalize so pow() behaves comparably across bands/patients
    A_norm = A / (A.max() + 1e-30)

    # Sweep rows: spring k (absolute), cols: iterations count.
    # Fix α=2 (moderate weight amplification).
    alpha = 2
    Ap = A_norm ** alpha

    k_values = [0.005, 0.01, 0.02, 0.04, 0.07, 0.1]
    iters_values = [500, 2000, 5000, 10000]

    fig, axes = plt.subplots(len(k_values), len(iters_values),
                             figsize=(2.6 * len(iters_values),
                                      2.6 * len(k_values)))

    for i, k in enumerate(k_values):
        for j, it in enumerate(iters_values):
            ax = axes[i, j]
            try:
                pos = _spring_layout(Ap, k=k, iterations=it)
            except Exception as exc:
                ax.text(0.5, 0.5, f"FAIL\n{exc}", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="red")
                ax.axis("off")
                if verbose:
                    print(f"    k={k} iter={it}: {exc}")
                continue
            _draw_edges(ax, pos, A, probes)
            ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=16,
                       edgecolors="white", linewidths=0.3, zorder=5)
            if i == 0:
                ax.set_title(f"iter={it}", fontsize=12, fontweight="bold")
            if j == 0:
                ax.text(-0.05, 0.5, f"k={k}", rotation=90,
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=11, fontweight="bold")
            _finalize(ax, pos, "")
            if verbose:
                print(f"    k={k} iter={it}: ok")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Spring k × iterations sweep — {fc_method.upper()} {patient} {band_tex} {phase}  "
        f"(w^{alpha} weights; nodes by shaft; red = same-probe)",
        fontsize=13, fontweight="bold", y=1.00,
    )
    fig.tight_layout()
    save_helper_fig(
        fig,
        output_dir / f"weight_transform_{fc_method}_{patient}_{band}_{phase}",
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
