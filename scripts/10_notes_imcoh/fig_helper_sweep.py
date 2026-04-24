#!/usr/bin/env python3
"""Helper sweep figures for tuning network layout parameters.

Two sweep grids, both for one representative (patient, band, phase, fc_method):

  1. Spring-k x edge-pct sweep  — find the k range where modular structure emerges
  2. Edge-width gamma sweep     — find the exponent that makes strong edges stand out

Run:
  python scripts/10_notes_imcoh/fig_helper_sweep.py [-v]
  python scripts/10_notes_imcoh/fig_helper_sweep.py --patient Pat_05 --band beta \\
      --phase rest_pre --fc-method imcoh
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
    apply_pub_style, save_fig, SECTION2_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix


def _shaft_colors(probes):
    """Color each node by its electrode shaft (A, B, G, ...). Each shaft gets
    a distinct tab20 color so probe geometry is visually readable."""
    uniq = sorted(set(probes))
    cmap = plt.get_cmap("tab20", len(uniq))
    return [cmap(uniq.index(p)) for p in probes]


def _draw_gamma_edges(ax, pos_arr, mat, probes, gamma: float,
                      width_range=(0.15, 5.0), alpha_range=(0.05, 0.95)):
    """Draw all edges with gamma-scaled width and alpha.

    Uses the FULL weighted graph (no top-X% thresholding).
    t = (w / w_max) ** gamma, then linearly mapped to width and alpha.
    gamma > 1 accentuates the top edges and makes weak ones thin/faint.
    """
    N = mat.shape[0]
    r, c = np.triu_indices(N, k=1)
    w = mat[r, c]
    active = np.where(w > 0)[0]
    if active.size == 0:
        return
    w_act = w[active]
    t = (w_act / w_act.max()) ** gamma
    widths = width_range[0] + t * (width_range[1] - width_range[0])
    alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])

    order = np.argsort(w_act)
    segments, colors, linew = [], [], []
    for oi in order:
        idx = active[oi]
        i, j = r[idx], c[idx]
        segments.append([pos_arr[i], pos_arr[j]])
        linew.append(widths[oi])
        if probes[i] == probes[j]:
            colors.append((0.8, 0.2, 0.2, alphas[oi]))
        else:
            colors.append((0.35, 0.35, 0.35, alphas[oi]))
    lc = LineCollection(segments, colors=colors, linewidths=linew,
                        zorder=1, rasterized=True)
    ax.add_collection(lc)


# ---------------------------------------------------------------------------
# Sweep 1: spring-k on the FULL weighted graph, various gamma for drawing
# ---------------------------------------------------------------------------

def sweep_spring_k(
    patient: str, band: str, phase: str, fc_method: str,
    output_dir: Path, verbose: bool = False,
):
    """Layout uses the full weighted graph (all edges).
    Rows = draw gamma (how aggressively to emphasize strong edges).
    Cols = k_base (spring repulsion).
    """
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        print(f"  SKIP: no {fc_method} data for {patient} {phase} {band}")
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    N = len(ch)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    k_bases = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    gammas = [1.0, 3.0, 6.0, 10.0]

    G = nx.from_numpy_array(A)

    fig, axes = plt.subplots(
        len(gammas), len(k_bases),
        figsize=(2.8 * len(k_bases), 2.8 * len(gammas)),
    )

    for j, kb in enumerate(k_bases):
        k_val = kb / np.sqrt(N)
        pos = nx.spring_layout(G, k=k_val, iterations=800, seed=42,
                               weight="weight")
        pos_arr = np.array([pos[k] for k in range(N)])
        for i, g in enumerate(gammas):
            ax = axes[i, j]
            _draw_gamma_edges(ax, pos_arr, A, probes, gamma=g)
            ax.scatter(pos_arr[:, 0], pos_arr[:, 1], c=node_colors, s=16,
                       edgecolors="white", linewidths=0.3, zorder=5)
            ax.axis("off")

            if i == 0:
                ax.set_title(f"k_base={kb}\nk={k_val:.3f}", fontsize=10)
            if j == 0:
                ax.text(-0.05, 0.5, f"γ={g}", rotation=90,
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=11, fontweight="bold")

            margin = 0.08
            xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
            ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
            span = max(xmax - xmin, ymax - ymin)
            ax.set_xlim(xmin - margin * span, xmax + margin * span)
            ax.set_ylim(ymin - margin * span, ymax + margin * span)

        if verbose:
            print(f"    k_base={kb}: layout done")

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Spring-k × edge-gamma sweep — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"full weighted graph (no edge dropping). "
        f"Nodes colored by electrode shaft.",
        fontsize=13, y=1.01,
    )
    fig.tight_layout()
    save_fig(
        fig,
        output_dir / f"sweep_k_{fc_method}_{patient}_{band}_{phase}",
    )


# ---------------------------------------------------------------------------
# Sweep 2: edge-width gamma
# ---------------------------------------------------------------------------

def sweep_edge_gamma(
    patient: str, band: str, phase: str, fc_method: str,
    k_base: float,
    output_dir: Path, verbose: bool = False,
):
    """Gamma sweep on the FULL weighted graph (no thresholding).
    Rows = width range, cols = gamma exponent.
    Higher gamma makes weak edges thin/faint and strong edges stand out.
    """
    mat = load_fc_matrix(patient, phase, band, fc_method)
    if mat is None:
        return
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    node_colors = _shaft_colors(probes)
    N = len(ch)
    A = np.abs(mat.copy())
    np.fill_diagonal(A, 0)

    G = nx.from_numpy_array(A)
    k_val = k_base / np.sqrt(N)
    pos = nx.spring_layout(G, k=k_val, iterations=800, seed=42, weight="weight")
    pos_arr = np.array([pos[k] for k in range(N)])

    gammas = [1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 20.0]
    width_ranges = [(0.15, 4.0), (0.2, 6.0), (0.3, 9.0)]

    fig, axes = plt.subplots(
        len(width_ranges), len(gammas),
        figsize=(2.8 * len(gammas), 2.8 * len(width_ranges)),
    )

    for i, wr in enumerate(width_ranges):
        for j, g in enumerate(gammas):
            ax = axes[i, j]
            _draw_gamma_edges(ax, pos_arr, A, probes, gamma=g,
                              width_range=wr, alpha_range=(0.03, 0.95))
            ax.scatter(pos_arr[:, 0], pos_arr[:, 1], c=node_colors, s=16,
                       edgecolors="white", linewidths=0.3, zorder=5)
            ax.axis("off")

            if i == 0:
                ax.set_title(f"γ={g}", fontsize=11)
            if j == 0:
                ax.text(-0.05, 0.5, f"width\n{wr[0]}-{wr[1]}", rotation=90,
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=10, fontweight="bold")

            margin = 0.08
            xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
            ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
            span = max(xmax - xmin, ymax - ymin)
            ax.set_xlim(xmin - margin * span, xmax + margin * span)
            ax.set_ylim(ymin - margin * span, ymax + margin * span)

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Edge-width γ sweep — {fc_method.upper()} {patient} {band_tex} {phase}\n"
        f"full weighted graph (no thresholding), k_base={k_base}. "
        f"widths = w_min + (w/w_max)^γ · (w_max - w_min). Nodes colored by shaft.",
        fontsize=13, y=1.01,
    )
    fig.tight_layout()
    save_fig(
        fig,
        output_dir / f"sweep_gamma_{fc_method}_{patient}_{band}_{phase}",
    )


def main():
    parser = argparse.ArgumentParser(description="Layout parameter sweep figures.")
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--fc-method", default="imcoh_abs")
    parser.add_argument("--k-base", type=float, default=10.0,
                        help="k_base for the gamma sweep (choose after the k sweep)")
    parser.add_argument("--output-dir", type=Path,
                        default=SECTION2_ROOT / "fig_helper")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("--- Sweep 1: spring-k × edge-pct ---")
    sweep_spring_k(
        args.patient, args.band, args.phase, args.fc_method,
        args.output_dir, args.verbose,
    )

    print(f"\n--- Sweep 2: edge-width γ (k_base={args.k_base}) ---")
    sweep_edge_gamma(
        args.patient, args.band, args.phase, args.fc_method,
        args.k_base,
        args.output_dir, args.verbose,
    )

    print("\nDone. Inspect figures and pick parameters that reveal modular structure.")


if __name__ == "__main__":
    main()
