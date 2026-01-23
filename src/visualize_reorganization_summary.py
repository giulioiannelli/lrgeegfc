#!/usr/bin/env python3
"""Band-dependent summary of reorganization metrics."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def _summarize(values: List[float]) -> Dict[str, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {"n": 0, "mean": np.nan, "sem": np.nan}
    mean = float(np.mean(arr))
    sem = float(np.std(arr, ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else 0.0
    return {"n": int(arr.size), "mean": mean, "sem": sem}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize reorganization metrics by band."
    )
    parser.add_argument("--fc-method", choices=["msc", "corr"], default="msc")
    parser.add_argument("--metric", default="cluster_swap")
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument(
        "--results-root", type=Path, default=Path("results/reorganization")
    )
    parser.add_argument(
        "--summary-dir", type=Path, default=Path("results/reorganization/summary")
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/figures/summary")
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    metric_specs = build_metric_specs()
    metric_labels = {key: spec["label"] for key, spec in metric_specs.items()}
    metric_labels.update(
        {
            "cluster_swap": "Cluster swap coefficient",
            "cluster_ari": "Adjusted Rand index",
        }
    )

    bands = _resolve_list(args.bands, list(BRAIN_BANDS.keys()))
    patients_filter = set(args.patients) if args.patients else None

    values_by_band: Dict[str, List[float]] = {band: [] for band in bands}

    for path in sorted(args.results_root.glob("Pat_*/*.npz")):
        parsed = _parse_metrics_file(path)
        if parsed is None:
            continue
        patient, band, fc_method = parsed
        if fc_method != args.fc_method:
            continue
        if band not in values_by_band:
            continue
        if patients_filter and patient not in patients_filter:
            continue
        with np.load(path) as data:
            if args.metric not in data:
                continue
            matrix = np.asarray(data[args.metric])
            if matrix.ndim != 2:
                continue
            n = matrix.shape[0]
            triu_vals = matrix[np.triu_indices(n, k=1)]
            values_by_band[band].extend(triu_vals.tolist())

    summary_rows: List[Dict[str, object]] = []
    for band in bands:
        stats = _summarize(values_by_band[band])
        summary_rows.append(
            {
                "band": band,
                "metric": args.metric,
                "metric_label": metric_labels.get(args.metric, args.metric),
                **stats,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary_dir / f"band_summary_{args.metric}_{args.fc_method}.csv"
    summary_df.to_csv(summary_path, index=False)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(
        summary_df["band"],
        summary_df["mean"],
        yerr=summary_df["sem"],
        color="#4c72b0",
        alpha=0.8,
        capsize=4,
    )
    ax.set_ylabel(metric_labels.get(args.metric, args.metric))
    ax.set_title(f"Band summary ({args.fc_method})")
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / f"band_summary_{args.metric}_{args.fc_method}.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if args.verbose:
        print(f"Saved: {summary_path}")
        print(f"Saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
