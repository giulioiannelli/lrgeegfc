#!/usr/bin/env python3
"""Visualize cached reorganization metric matrices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.utils.metrics.reorganization import build_metric_specs


def _parse_metrics_file(path: Path) -> tuple[str, str, str] | None:
    if path.suffix != ".npz":
        return None
    patient = path.parent.name
    if not patient.startswith("Pat_"):
        return None
    stem = path.stem
    if not stem.endswith("_metrics"):
        return None
    prefix = stem[: -len("_metrics")]
    parts = prefix.split("_")
    if len(parts) < 2:
        return None
    fc_method = parts[-1]
    band = "_".join(parts[:-1])
    return patient, band, fc_method


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def _plot_metrics(
    phases: List[str],
    matrices: Dict[str, np.ndarray],
    labels: Dict[str, str],
    output_path: Path,
    title: str,
) -> None:
    metric_keys = list(matrices.keys())
    n_metrics = len(metric_keys)
    if n_metrics == 0:
        raise ValueError("No metric matrices to plot.")

    ncols = 3
    nrows = int(np.ceil(n_metrics / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = np.array(axes).reshape(-1)

    for ax in axes[n_metrics:]:
        ax.axis("off")

    for idx, key in enumerate(metric_keys):
        mat = matrices[key]
        ax = axes[idx]
        finite_vals = mat[np.isfinite(mat) & ~np.eye(len(phases), dtype=bool)]
        if finite_vals.size:
            vmin, vmax = float(finite_vals.min()), float(finite_vals.max())
        else:
            vmin, vmax = 0.0, 1.0
        im = ax.imshow(mat, cmap="YlOrRd", vmin=vmin, vmax=vmax, interpolation="none")
        ax.set_title(labels.get(key, key), fontsize=11)
        ax.set_xticks(range(len(phases)))
        ax.set_yticks(range(len(phases)))
        ax.set_xticklabels(phases, rotation=45, ha="right", fontsize=9)
        ax.set_yticklabels(phases, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot reorganization metric matrices from cached outputs."
    )
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--fc-method", choices=["msc", "corr"], default="msc")
    parser.add_argument("--metrics", nargs="+", default=None)
    parser.add_argument(
        "--results-root", type=Path, default=Path("results/reorganization")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/figures/reorganization")
    )
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    metric_specs = build_metric_specs()
    default_metrics = list(metric_specs.keys()) + ["cluster_swap", "cluster_ari"]
    metric_keys = args.metrics or default_metrics

    metric_labels = {key: spec["label"] for key, spec in metric_specs.items()}
    metric_labels.update(
        {
            "cluster_swap": "Cluster swap coefficient",
            "cluster_ari": "Adjusted Rand index",
        }
    )

    patients_filter = set(args.patients) if args.patients else None
    bands = _resolve_list(args.bands, list(BRAIN_BANDS.keys()))

    for path in sorted(args.results_root.glob("Pat_*/*.npz")):
        parsed = _parse_metrics_file(path)
        if parsed is None:
            continue
        patient, band, fc_method = parsed
        if fc_method != args.fc_method:
            continue
        if patients_filter and patient not in patients_filter:
            continue
        if band not in bands:
            continue

        output_path = args.output_dir / patient / f"{band}_{fc_method}_metrics_grid.png"
        if args.skip_existing and output_path.exists():
            if args.verbose:
                print(f"• Skipping (exists): {output_path}")
            continue

        with np.load(path) as data:
            phases = data.get("phases")
            if phases is None:
                if args.verbose:
                    print(f"[skip] missing phases in {path.name}")
                continue
            phases = [str(p) for p in phases.tolist()]

            matrices: Dict[str, np.ndarray] = {}
            for key in metric_keys:
                if key in data:
                    matrix = np.asarray(data[key])
                    if matrix.ndim == 2:
                        matrices[key] = matrix

        if not matrices:
            if args.verbose:
                print(f"[skip] no metrics to plot for {path.name}")
            continue

        title = f"{patient} {band} ({fc_method}) reorganization metrics"
        _plot_metrics(phases, matrices, metric_labels, output_path, title)

        if args.verbose:
            print(f"✓ Saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
