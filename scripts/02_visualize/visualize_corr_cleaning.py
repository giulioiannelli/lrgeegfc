#!/usr/bin/env python3
"""Visualize Marchenko-Pastur cleaning outputs from cached results."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import CORR_CACHE, FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.workflow.cleaning import load_cleaned_corr_matrix
from lrg_eegfc.workflow.corr import load_corr_matrix


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def _plot_cleaning(
    corr_matrix: np.ndarray,
    cleaned_matrix: np.ndarray,
    metadata: dict,
    title: str,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    im0 = axes[0].imshow(corr_matrix, cmap="viridis", vmin=0, vmax=1, interpolation="none")
    axes[0].set_title("Original |corr|", fontsize=10)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(cleaned_matrix, cmap="viridis", vmin=0, vmax=1, interpolation="none")
    axes[1].set_title("Cleaned (MP + threshold)", fontsize=10)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    eigvals = np.asarray(metadata.get("eigenvalues", []))
    lambda_min = metadata.get("lambda_min")
    lambda_max = metadata.get("lambda_max")
    axes[2].hist(eigvals, bins=40, color="#4c72b0", alpha=0.8)
    if lambda_min is not None:
        axes[2].axvline(lambda_min, color="black", linestyle="--", label="lambda_min")
    if lambda_max is not None:
        axes[2].axvline(lambda_max, color="red", linestyle="--", label="lambda_max")
    axes[2].set_title("Eigenvalue spectrum", fontsize=10)
    axes[2].set_xlabel("Eigenvalue")
    axes[2].set_ylabel("Count")
    axes[2].legend(fontsize=8)

    threshold = metadata.get("threshold")
    if threshold is not None:
        fig.suptitle(f"{title} (threshold={threshold:.4f})", fontsize=12)
    else:
        fig.suptitle(title, fontsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot correlation cleaning diagnostics from cached results."
    )
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--phases", nargs="+", default=None)
    parser.add_argument("--cache-root", type=Path, default=CORR_CACHE)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=FIGURES_ROOT / "cleaning",
    )
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    from lrg_eegfc.utils.io import list_patients

    patients = _resolve_list(args.patients, list_patients(SEEG_DATAPATH))
    bands = _resolve_list(args.bands, list(BRAIN_BANDS.keys()))
    phases = _resolve_list(args.phases, list(PHASE_LABELS))

    for patient in patients:
        for band in bands:
            for phase in phases:
                output_path = (
                    args.output_dir
                    / patient
                    / f"{band}_{phase}_corr_cleaning.png"
                )
                if args.skip_existing and output_path.exists():
                    if args.verbose:
                        print(f"• Skipping (exists): {output_path}")
                    continue

                cleaned, metadata = load_cleaned_corr_matrix(
                    patient, phase, band, cache_root=args.cache_root, load_metadata=True
                )
                if cleaned is None or metadata is None:
                    if args.verbose:
                        print(f"[skip] missing cleaned cache for {patient} {band} {phase}")
                    continue

                corr_matrix = load_corr_matrix(
                    patient,
                    phase,
                    band,
                    cache_root=args.cache_root,
                    filter_type="abs",
                    zero_diagonal=True,
                )
                if corr_matrix is None:
                    if args.verbose:
                        print(f"[skip] missing corr cache for {patient} {band} {phase}")
                    continue

                _plot_cleaning(
                    corr_matrix,
                    cleaned,
                    metadata,
                    f"{patient} {band} {phase}",
                    output_path,
                )
                if args.verbose:
                    print(f"✓ Saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
