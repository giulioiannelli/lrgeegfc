#!/usr/bin/env python3
"""Create a phase x band grid of MSC matrices for a single patient."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, MSC_CACHE
from lrg_eegfc.workflow.msc import load_msc_matrix


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot a phase x band grid of cached MSC matrices."
    )
    parser.add_argument("--patient", required=True, help="Patient ID (e.g., Pat_02)")
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--phases", nargs="+", default=None)
    parser.add_argument(
        "--sparsify",
        choices=["none", "soft", "fdr", "disparity", "hybrid", "ecm"],
        default="none",
        help="MSC sparsify mode (default: none)",
    )
    parser.add_argument("--n-surrogates", type=int, default=0)
    parser.add_argument("--nperseg", type=int, default=1024)
    parser.add_argument("--cache-root", type=Path, default=MSC_CACHE)
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output file path (default: data/figures/msc/<patient>/phase_band_grid.png)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    bands = _resolve_list(args.bands, list(BRAIN_BANDS.keys()))
    phases = _resolve_list(args.phases, list(PHASE_LABELS))

    n_rows = len(phases)
    n_cols = len(bands)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.2 * n_cols, 3.2 * n_rows),
        squeeze=False,
    )

    missing = 0
    for i, phase in enumerate(phases):
        for j, band in enumerate(bands):
            ax = axes[i][j]
            matrix = load_msc_matrix(
                args.patient,
                phase,
                band,
                cache_root=args.cache_root,
                sparsify=args.sparsify,
                n_surrogates=args.n_surrogates,
                nperseg=args.nperseg,
            )
            if matrix is None:
                missing += 1
                ax.set_facecolor("#f0f0f0")
                ax.text(
                    0.5,
                    0.5,
                    "missing",
                    ha="center",
                    va="center",
                    fontsize=10,
                )
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                ax.imshow(
                    matrix,
                    cmap="viridis",
                    vmin=0,
                    vmax=1,
                    interpolation="none",
                    aspect="equal",
                )
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(band, fontsize=10)
                if j == 0:
                    ax.set_ylabel(phase, fontsize=10)

    fig.suptitle(
        f"{args.patient} MSC grid (sparsify={args.sparsify}, nperseg={args.nperseg})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if args.output_path is None:
        output_dir = FIGURES_ROOT / "msc" / args.patient
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "phase_band_grid.png"
    else:
        output_path = args.output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if args.verbose:
        print(f"Saved: {output_path}")
        if missing:
            print(f"Missing matrices: {missing}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
