#!/usr/bin/env python3
"""Plot the same band/phase MSC matrices across all patients."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.utils.io import list_patients
from lrg_eegfc.workflow.msc import load_msc_matrix


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def _grid_dims(n_items: int, max_cols: int = 4) -> tuple[int, int]:
    cols = min(max_cols, max(1, int(math.ceil(math.sqrt(n_items)))))
    rows = int(math.ceil(n_items / cols))
    return rows, cols


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot MSC matrices across all patients for a band/phase."
    )
    parser.add_argument("--band", required=True, choices=list(BRAIN_BANDS.keys()))
    parser.add_argument("--phase", required=True, choices=list(PHASE_LABELS))
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument(
        "--sparsify",
        choices=["none", "soft"],
        default="none",
        help="MSC sparsify mode (default: none)",
    )
    parser.add_argument("--n-surrogates", type=int, default=0)
    parser.add_argument("--nperseg", type=int, default=1024)
    parser.add_argument("--cache-root", type=Path, default=Path("data/msc_cache"))
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output file path (default: data/figures/msc/all_patients_<band>_<phase>.png)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    patients = _resolve_list(
        args.patients, list_patients(Path("data/stereoeeg_patients"))
    )
    if not patients:
        print("No patients found.")
        return 1

    rows, cols = _grid_dims(len(patients))
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(3.2 * cols, 3.2 * rows),
        squeeze=False,
    )

    missing = 0
    for idx, patient in enumerate(patients):
        r, c = divmod(idx, cols)
        ax = axes[r][c]
        matrix = load_msc_matrix(
            patient,
            args.phase,
            args.band,
            cache_root=args.cache_root,
            sparsify=args.sparsify,
            n_surrogates=args.n_surrogates,
            nperseg=args.nperseg,
        )
        if matrix is None:
            missing += 1
            ax.set_facecolor("#f0f0f0")
            ax.text(0.5, 0.5, "missing", ha="center", va="center", fontsize=10)
        else:
            ax.imshow(
                matrix,
                cmap="viridis",
                vmin=0,
                vmax=1,
                interpolation="none",
                aspect="equal",
            )
        ax.set_title(patient, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

    for idx in range(len(patients), rows * cols):
        r, c = divmod(idx, cols)
        axes[r][c].axis("off")

    fig.suptitle(
        f"MSC {args.band} {args.phase} (sparsify={args.sparsify})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if args.output_path is None:
        output_dir = Path("data/figures/msc")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"all_patients_{args.band}_{args.phase}.png"
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
