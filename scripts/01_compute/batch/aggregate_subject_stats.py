#!/usr/bin/env python3
"""Aggregate reorganization metrics across subjects."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.config.paths import FIGURES_ROOT
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


def _summarize_values(values: np.ndarray) -> Dict[str, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"n": 0, "mean": np.nan, "sem": np.nan, "median": np.nan, "iqr": np.nan}
    mean = float(np.mean(finite))
    std = float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0
    sem = float(std / np.sqrt(finite.size)) if finite.size > 1 else 0.0
    median = float(np.median(finite))
    q75, q25 = np.percentile(finite, [75, 25])
    iqr = float(q75 - q25)
    return {"n": int(finite.size), "mean": mean, "sem": sem, "median": median, "iqr": iqr}


def _plot_mean_matrix(
    matrix: np.ndarray,
    phases: List[str],
    title: str,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    masked = np.ma.masked_invalid(matrix)
    finite_vals = matrix[np.isfinite(matrix) & ~np.eye(len(phases), dtype=bool)]
    if finite_vals.size:
        vmin, vmax = float(finite_vals.min()), float(finite_vals.max())
    else:
        vmin, vmax = 0.0, 1.0
    im = ax.imshow(masked, cmap="YlOrRd", vmin=vmin, vmax=vmax, interpolation="none")
    ax.set_xticks(range(len(phases)))
    ax.set_yticks(range(len(phases)))
    ax.set_xticklabels(phases, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(phases, fontsize=9)
    ax.set_title(title, fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate reorganization metrics across patients."
    )
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--fc-method", choices=["msc", "corr"], default="msc")
    parser.add_argument("--metrics", nargs="+", default=None)
    parser.add_argument(
        "--results-root", type=Path, default=Path("results/reorganization")
    )
    parser.add_argument(
        "--summary-dir", type=Path, default=Path("results/reorganization/summary")
    )
    parser.add_argument(
        "--figures-dir", type=Path, default=FIGURES_ROOT / "summary"
    )
    parser.add_argument("--no-figures", action="store_true")
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

    matrices_by_band: Dict[str, Dict[str, List[np.ndarray]]] = {
        band: {key: [] for key in metric_keys} for band in bands
    }
    phases_by_band: Dict[str, List[str]] = {}

    for path in sorted(args.results_root.glob("Pat_*/*.npz")):
        parsed = _parse_metrics_file(path)
        if parsed is None:
            continue
        patient, band, fc_method = parsed
        if fc_method != args.fc_method:
            continue
        if band not in matrices_by_band:
            continue
        if patients_filter and patient not in patients_filter:
            continue
        with np.load(path) as data:
            phases = data.get("phases")
            if phases is None:
                continue
            phases_list = [str(p) for p in phases.tolist()]
            if band not in phases_by_band:
                phases_by_band[band] = phases_list
            elif phases_by_band[band] != phases_list:
                if args.verbose:
                    print(f"[skip] phase mismatch in {path.name}")
                continue

            for key in metric_keys:
                if key not in data:
                    continue
                matrix = np.asarray(data[key])
                if matrix.ndim != 2:
                    continue
                matrices_by_band[band][key].append(matrix)

    summary_rows: List[Dict[str, object]] = []

    for band in bands:
        phases = phases_by_band.get(band)
        if not phases:
            continue
        for key in metric_keys:
            mats = matrices_by_band[band].get(key, [])
            if not mats:
                continue
            stack = np.stack(mats, axis=0)
            mean_matrix = np.nanmean(stack, axis=0)
            sem_matrix = np.nanstd(stack, axis=0, ddof=1) / np.sqrt(stack.shape[0])

            for i, phase_i in enumerate(phases):
                for j, phase_j in enumerate(phases):
                    if i >= j:
                        continue
                    values = stack[:, i, j]
                    stats = _summarize_values(values)
                    summary_rows.append(
                        {
                            "band": band,
                            "metric": key,
                            "metric_label": metric_labels.get(key, key),
                            "phase_i": phase_i,
                            "phase_j": phase_j,
                            **stats,
                        }
                    )

            band_dir = args.summary_dir / band
            band_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(mean_matrix, index=phases, columns=phases).to_csv(
                band_dir / f"{key}_{args.fc_method}_mean.csv"
            )
            pd.DataFrame(sem_matrix, index=phases, columns=phases).to_csv(
                band_dir / f"{key}_{args.fc_method}_sem.csv"
            )

            if not args.no_figures:
                fig_path = (
                    args.figures_dir
                    / band
                    / f"{key}_{args.fc_method}_mean.png"
                )
                _plot_mean_matrix(
                    mean_matrix,
                    phases,
                    f"{band} {metric_labels.get(key, key)} (mean)",
                    fig_path,
                )

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        args.summary_dir.mkdir(parents=True, exist_ok=True)
        summary_path = args.summary_dir / f"{args.fc_method}_metrics_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        if args.verbose:
            print(f"Saved summary: {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
