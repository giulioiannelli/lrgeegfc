#!/usr/bin/env python3
"""Visualize cached time-window FC matrices as a grid."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault("JOBLIB_MULTIPROCESSING", "0")

import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.time_windows import build_window_run_id, get_window_cache_dir


def _load_windows_csv(path: Path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Visualize cached time-window FC matrices as a grid."
    )
    parser.add_argument("--patient", required=True, help="Patient ID")
    parser.add_argument("--phase", required=True, help="Phase name")
    parser.add_argument("--band", required=True, help="Band name")
    parser.add_argument(
        "--fc-method",
        choices=["corr", "msc"],
        default="corr",
        help="FC method for cached windows",
    )

    parser.add_argument("--window-sec", type=float, required=True, help="Window length in seconds")
    parser.add_argument("--overlap", type=float, default=0.25, help="Window overlap fraction")
    parser.add_argument(
        "--filter-time",
        type=int,
        default=None,
        help="Filter-time tag used in cache naming (dev runs)",
    )

    parser.add_argument(
        "--filter-type",
        choices=["abs", "pos", "neg", "none"],
        default="abs",
        help="Correlation filter type (corr only)",
    )
    parser.add_argument(
        "--keep-diagonal",
        action="store_true",
        help="Keep diagonal in correlation matrices (corr only)",
    )
    parser.add_argument("--filter-order", type=int, default=4, help="Bandpass filter order")

    parser.add_argument("--sparsify", choices=["none", "soft", "fdr", "disparity", "hybrid", "ecm"], default="none")
    parser.add_argument("--n-surrogates", type=int, default=0)
    parser.add_argument("--nperseg", type=int, default=1024)
    parser.add_argument("--noverlap", type=int, default=None)

    parser.add_argument(
        "--max-windows",
        type=int,
        default=9,
        help="Max windows to plot (default: 9)",
    )
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--cmap", type=str, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for figures",
    )

    args = parser.parse_args()

    zero_diagonal = not args.keep_diagonal

    run_id = build_window_run_id(
        fc_method=args.fc_method,
        window_sec=args.window_sec,
        overlap=args.overlap,
        filter_time=args.filter_time,
        filter_type=args.filter_type,
        zero_diagonal=zero_diagonal,
        filter_order=args.filter_order,
        sparsify=args.sparsify,
        n_surrogates=args.n_surrogates,
        nperseg=args.nperseg,
        noverlap=args.noverlap,
    )

    run_dir = get_window_cache_dir(
        patient=args.patient,
        phase=args.phase,
        band=args.band,
        fc_method=args.fc_method,
        window_sec=args.window_sec,
        overlap=args.overlap,
        filter_time=args.filter_time,
        filter_type=args.filter_type,
        zero_diagonal=zero_diagonal,
        filter_order=args.filter_order,
        sparsify=args.sparsify,
        n_surrogates=args.n_surrogates,
        nperseg=args.nperseg,
        noverlap=args.noverlap,
        ensure_dir=False,
    )

    if not run_dir.exists():
        raise SystemExit(f"Run directory not found: {run_dir}")

    windows_csv = run_dir / "windows.csv"
    if not windows_csv.exists():
        raise SystemExit(f"Missing windows.csv in {run_dir}")

    rows = _load_windows_csv(windows_csv)
    if not rows:
        raise SystemExit(f"No windows listed in {windows_csv}")

    rows = rows[: args.max_windows]

    if args.fc_method == "msc":
        vmin = 0.0 if args.vmin is None else args.vmin
        vmax = 1.0 if args.vmax is None else args.vmax
        cmap = args.cmap or "viridis"
    else:
        if args.filter_type == "none":
            vmin = -1.0 if args.vmin is None else args.vmin
        else:
            vmin = 0.0 if args.vmin is None else args.vmin
        vmax = 1.0 if args.vmax is None else args.vmax
        cmap = args.cmap or "coolwarm"

    n_plots = len(rows)
    n_cols = math.ceil(math.sqrt(n_plots))
    n_rows = math.ceil(n_plots / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.flatten()

    for ax in axes[n_plots:]:
        ax.axis("off")

    for ax, row in zip(axes, rows):
        idx = int(row["window_index"])
        start = float(row["start_sec"])
        stop = float(row["stop_sec"])
        matrix_path = run_dir / f"win-{idx:04d}.npy"
        if not matrix_path.exists():
            ax.set_title(f"win {idx} missing")
            ax.axis("off")
            continue

        matrix = np.load(matrix_path)
        im = ax.imshow(matrix, vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_title(f"win {idx} ({start:.1f}-{stop:.1f}s)")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        f"{args.fc_method.upper()} {args.patient} {args.phase} {args.band} ({run_id})",
        fontsize=12,
    )
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    cbar = fig.colorbar(im, ax=axes[:n_plots], shrink=0.7)
    cbar.ax.tick_params(labelsize=8)

    output_dir = args.output_dir or FIGURES_ROOT / "time_windows" / args.patient / args.phase / args.band
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{args.fc_method}_windows_{run_id}_grid.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
