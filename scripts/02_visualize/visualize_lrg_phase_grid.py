#!/usr/bin/env python3
"""Plot a phase grid of LRG dendrograms for one patient/band."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, LRG_CACHE
from lrg_eegfc.workflow.lrg import load_lrg_result


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot dendrograms across phases for a patient/band."
    )
    parser.add_argument("--patient", required=True)
    parser.add_argument("--band", required=True)
    parser.add_argument("--fc-method", required=True, choices=["msc", "corr"])
    parser.add_argument("--phases", nargs="+", default=None)
    parser.add_argument("--cache-root", type=Path, default=LRG_CACHE)
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output file path (default: data/figures/lrg/<patient>/<band>_<fc>_phase_dendrograms.png)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    phases = _resolve_list(args.phases, list(PHASE_LABELS))
    n_phases = len(phases)
    ncols = 2
    nrows = int(np.ceil(n_phases / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows))
    axes = np.array(axes).reshape(-1)

    for idx, phase in enumerate(phases):
        ax = axes[idx]
        result = load_lrg_result(
            args.patient,
            phase,
            args.band,
            args.fc_method,
            cache_root=args.cache_root,
        )
        if result is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center")
            ax.set_title(f"{phase} (missing)")
            ax.axis("off")
            continue

        linkage = result.linkage_matrix
        optimal_threshold = result.optimal_threshold
        dendrogram(
            linkage,
            ax=ax,
            color_threshold=optimal_threshold,
            no_labels=(result.n_nodes > 40),
        )
        ax.axhline(optimal_threshold, color="red", linestyle="--", linewidth=1)
        ax.set_title(f"{phase}", fontsize=10)
        ax.set_xlabel("Cluster distance", fontsize=8)
        ax.set_ylabel("Nodes", fontsize=8)

    for idx in range(n_phases, len(axes)):
        axes[idx].axis("off")

    fig.suptitle(
        f"{args.patient} {args.band} ({args.fc_method}) - Phase dendrograms",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if args.output_path is None:
        output_dir = FIGURES_ROOT / "lrg" / args.patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{args.band}_{args.fc_method}_phase_dendrograms.png"
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
