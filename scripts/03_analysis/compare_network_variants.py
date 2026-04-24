#!/usr/bin/env python3
"""Compare correlation vs dense/validated MSC networks for one case."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from lrg_eegfc.workflow.corr import load_corr_matrix
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.paths import CORR_CACHE, MSC_CACHE, FIGURES_ROOT


def _threshold_matrix(matrix: np.ndarray, quantile: float) -> Tuple[np.ndarray, float]:
    if quantile <= 0:
        return matrix.copy(), 0.0
    triu = matrix[np.triu_indices_from(matrix, k=1)]
    triu = triu[np.isfinite(triu)]
    if triu.size == 0:
        return matrix.copy(), 0.0
    threshold = float(np.quantile(triu, quantile))
    masked = matrix.copy()
    masked[masked < threshold] = 0.0
    return masked, threshold


def _draw_network(ax, matrix: np.ndarray, title: str, quantile: float) -> None:
    matrix = matrix.copy()
    np.fill_diagonal(matrix, 0)
    masked, threshold = _threshold_matrix(matrix, quantile)
    graph = nx.from_numpy_array(masked)
    if graph.number_of_edges() == 0:
        ax.text(0.5, 0.5, "No edges", ha="center", va="center")
        ax.set_title(title)
        ax.axis("off")
        return

    try:
        pos = nx.spring_layout(graph, seed=42, k=0.2, iterations=50)
    except Exception:
        pos = nx.circular_layout(graph)

    widths = [graph[u][v]["weight"] for u, v in graph.edges()]
    nx.draw(
        graph,
        pos=pos,
        ax=ax,
        node_size=60,
        node_color="lightblue",
        edge_color="gray",
        width=widths,
        with_labels=False,
    )
    ax.set_title(f"{title}\nq={quantile:.2f}, thr={threshold:.3f}", fontsize=9)
    ax.axis("off")


def _imshow(ax, matrix: np.ndarray, title: str) -> None:
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1, interpolation="none")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare correlation vs MSC network variants."
    )
    parser.add_argument("--patient", required=True)
    parser.add_argument("--phase", required=True)
    parser.add_argument("--band", required=True)
    parser.add_argument("--nperseg", type=int, default=1024)
    parser.add_argument("--n-surrogates", type=int, default=200)
    parser.add_argument(
        "--edge-quantile",
        type=float,
        default=0.9,
        help="Quantile threshold for network visualization (default: 0.9)",
    )
    parser.add_argument("--cache-root-corr", type=Path, default=CORR_CACHE)
    parser.add_argument("--cache-root-msc", type=Path, default=MSC_CACHE)
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output file path (default: data/figures/comparison/<patient>/<band>_<phase>_network_variants.png)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    corr_matrix = load_corr_matrix(
        args.patient,
        args.phase,
        args.band,
        cache_root=args.cache_root_corr,
        filter_type="abs",
        zero_diagonal=True,
    )
    if corr_matrix is None:
        print("Missing correlation cache.")
        return 1

    dense_msc = load_msc_matrix(
        args.patient,
        args.phase,
        args.band,
        cache_root=args.cache_root_msc,
        sparsify="none",
        n_surrogates=0,
        nperseg=args.nperseg,
    )
    if dense_msc is None:
        print("Missing dense MSC cache.")
        return 1

    validated_msc = load_msc_matrix(
        args.patient,
        args.phase,
        args.band,
        cache_root=args.cache_root_msc,
        sparsify="soft",
        n_surrogates=args.n_surrogates,
        nperseg=args.nperseg,
    )
    if validated_msc is None:
        print("Missing validated MSC cache.")
        return 1

    corr_matrix = np.clip(corr_matrix, 0, None)
    dense_msc = np.clip(dense_msc, 0, None)
    validated_msc = np.clip(validated_msc, 0, None)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    _imshow(axes[0, 0], corr_matrix, "Correlation (abs)")
    _imshow(axes[0, 1], dense_msc, "MSC dense")
    _imshow(axes[0, 2], validated_msc, "MSC validated")

    _draw_network(axes[1, 0], corr_matrix, "Correlation network", args.edge_quantile)
    _draw_network(axes[1, 1], dense_msc, "MSC dense network", args.edge_quantile)
    _draw_network(axes[1, 2], validated_msc, "MSC validated network", args.edge_quantile)

    fig.suptitle(
        f"{args.patient} {args.band} {args.phase} - Network variants",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if args.output_path is None:
        output_dir = FIGURES_ROOT / "comparison" / args.patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{args.band}_{args.phase}_network_variants.png"
    else:
        output_path = args.output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if args.verbose:
        print(f"Saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
