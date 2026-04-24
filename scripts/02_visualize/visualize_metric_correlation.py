#!/usr/bin/env python3
"""Summarize cross-metric agreement for reorganization metrics.

Examples
--------
python scripts/py/visualize_metric_correlation.py --fc-method msc --verbose
python scripts/py/visualize_metric_correlation.py --fc-method corr --method pearson --verbose
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.utils.metrics.reorganization import build_metric_specs


def _parse_result_path(path: Path) -> Optional[Tuple[str, str, str]]:
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
    if not parts:
        return None
    fc_method = parts[-1]
    if fc_method not in {"msc", "corr"}:
        return None
    band = "_".join(parts[:-1])
    if not band:
        return None
    return patient, band, fc_method


def _phase_pairs(phases: Iterable[str]) -> List[Tuple[int, int]]:
    phase_list = list(phases)
    pairs = []
    for i in range(len(phase_list)):
        for j in range(i + 1, len(phase_list)):
            pairs.append((i, j))
    return pairs


def _pairwise_correlation(
    values_a: np.ndarray,
    values_b: np.ndarray,
    method: str,
    min_samples: int,
) -> Tuple[float, int]:
    mask = np.isfinite(values_a) & np.isfinite(values_b)
    count = int(mask.sum())
    if count < min_samples:
        return np.nan, count

    data_a = values_a[mask]
    data_b = values_b[mask]
    if np.nanstd(data_a) == 0 or np.nanstd(data_b) == 0:
        return np.nan, count

    if method == "pearson":
        corr = float(np.corrcoef(data_a, data_b)[0, 1])
        return corr, count

    corr = spearmanr(data_a, data_b).correlation
    return (float(corr) if np.isfinite(corr) else np.nan), count


def _collect_metric_values(
    paths: List[Path],
    metric_keys: List[str],
    verbose: bool,
) -> Dict[str, np.ndarray]:
    values: Dict[str, List[float]] = {key: [] for key in metric_keys}
    sample_count = 0
    for path in paths:
        parsed = _parse_result_path(path)
        if parsed is None:
            continue
        if verbose:
            print(f"• Reading {path}")
        with np.load(path) as data:
            if "phases" not in data:
                if verbose:
                    print(f"  [skip] missing phases in {path.name}")
                continue
            phases = data["phases"].tolist()
            pairs = _phase_pairs(phases)
            for i, j in pairs:
                sample_count += 1
                for key in metric_keys:
                    if key not in data:
                        values[key].append(np.nan)
                        continue
                    matrix = np.asarray(data[key])
                    if matrix.ndim != 2:
                        values[key].append(np.nan)
                        continue
                    if i >= matrix.shape[0] or j >= matrix.shape[1]:
                        values[key].append(np.nan)
                        continue
                    values[key].append(float(matrix[i, j]))
    if verbose:
        print(f"Collected {sample_count} phase-pair samples.")
    if sample_count == 0:
        return {}
    return {key: np.asarray(vals, dtype=float) for key, vals in values.items()}


def _plot_correlation_matrix(
    matrix: np.ndarray,
    labels: List[str],
    output_path: Path,
    title: str,
) -> None:
    n = len(labels)
    fig_size = max(6.0, 0.6 * n)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    cmap = plt.cm.get_cmap("coolwarm").copy()
    cmap.set_bad(color="#d9d9d9")
    masked = np.ma.masked_invalid(matrix)

    im = ax.imshow(masked, cmap=cmap, vmin=-1, vmax=1, interpolation="none")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            if not np.isfinite(val):
                continue
            color = "white" if abs(val) > 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)

    ax.set_title(title, fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Correlation")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute cross-metric correlation summary for reorganization metrics."
    )
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--fc-method", choices=["msc", "corr"], default="msc")
    parser.add_argument("--method", choices=["spearman", "pearson"], default="spearman")
    parser.add_argument("--min-samples", type=int, default=3)
    parser.add_argument("--results-root", type=Path, default=Path("results/reorganization"))
    parser.add_argument("--summary-dir", type=Path, default=Path("results/reorganization/summary"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=FIGURES_ROOT / "reorganization",
    )
    parser.add_argument("--metrics", nargs="+", default=None)
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

    bands = args.bands or list(BRAIN_BANDS.keys())
    patients_filter = set(args.patients) if args.patients else None
    bands_filter = set(bands)

    if args.verbose:
        print("=== Metric correlation summary ===")
        print(f"  fc method: {args.fc_method}")
        print(f"  method: {args.method}")
        print(f"  results root: {args.results_root}")

    result_paths: List[Path] = []
    for path in args.results_root.glob("Pat_*/*.npz"):
        parsed = _parse_result_path(path)
        if parsed is None:
            continue
        patient, band, fc_method = parsed
        if fc_method != args.fc_method:
            continue
        if patients_filter and patient not in patients_filter:
            continue
        if band not in bands_filter:
            continue
        result_paths.append(path)

    if not result_paths:
        print("No metric files found for the requested filters.")
        return 1

    values_by_metric = _collect_metric_values(result_paths, metric_keys, args.verbose)
    if not values_by_metric:
        print("No metric values collected; check input files.")
        return 1

    if args.metrics is None:
        present_keys = [
            key for key in metric_keys if np.isfinite(values_by_metric[key]).any()
        ]
        missing_keys = [key for key in metric_keys if key not in present_keys]
        if args.verbose and missing_keys:
            print(f"[warn] No data for metrics: {', '.join(missing_keys)}")
        if not present_keys:
            print("No metrics with finite values found; check inputs.")
            return 1
        metric_keys = present_keys
        values_by_metric = {key: values_by_metric[key] for key in metric_keys}

    metric_values = [values_by_metric[key] for key in metric_keys]
    n_metrics = len(metric_keys)
    corr_matrix = np.full((n_metrics, n_metrics), np.nan, dtype=float)
    count_matrix = np.zeros((n_metrics, n_metrics), dtype=int)

    for i in range(n_metrics):
        for j in range(i, n_metrics):
            if i == j:
                count = int(np.isfinite(metric_values[i]).sum())
                count_matrix[i, j] = count
                corr_matrix[i, j] = 1.0 if count >= args.min_samples else np.nan
                continue
            corr, count = _pairwise_correlation(
                metric_values[i],
                metric_values[j],
                args.method,
                args.min_samples,
            )
            corr_matrix[i, j] = corr
            corr_matrix[j, i] = corr
            count_matrix[i, j] = count
            count_matrix[j, i] = count

    labels = [metric_labels.get(key, key) for key in metric_keys]
    corr_df = pd.DataFrame(corr_matrix, index=labels, columns=labels)
    count_df = pd.DataFrame(count_matrix, index=labels, columns=labels)

    summary_prefix = f"metric_correlation_{args.fc_method}_{args.method}"
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    corr_df.to_csv(args.summary_dir / f"{summary_prefix}.csv")
    count_df.to_csv(args.summary_dir / f"{summary_prefix}_counts.csv")

    figure_path = args.output_dir / f"{summary_prefix}.png"
    title = f"Metric Correlation ({args.fc_method}, {args.method})"
    _plot_correlation_matrix(corr_matrix, labels, figure_path, title)

    if args.verbose:
        print(f"Saved: {args.summary_dir / f'{summary_prefix}.csv'}")
        print(f"Saved: {figure_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
