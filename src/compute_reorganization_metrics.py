#!/usr/bin/env python3
"""Compute reorganization metrics from cached LRG results."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.utils.io import list_patients
from lrg_eegfc.utils.metrics.reorganization import (
    build_metric_specs,
    compute_metric_matrix,
    compute_cluster_labels,
    compute_cluster_swap_matrix,
    compute_ari_matrix,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


def _load_phase_results(
    patient: str,
    band: str,
    phases: List[str],
    fc_method: str,
    cache_root: Path,
) -> Dict[str, object]:
    results: Dict[str, object] = {}
    for phase in phases:
        result = load_lrg_result(patient, phase, band, fc_method, cache_root)
        if result is None:
            return {}
        results[phase] = result
    return results


def _save_npz(output_path: Path, payload: Dict[str, np.ndarray], phases: List[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, phases=np.array(phases), **payload)


def _save_csvs(
    output_dir: Path,
    prefix: str,
    payload: Dict[str, np.ndarray],
    phases: List[str],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, matrix in payload.items():
        if matrix.ndim != 2:
            continue
        df = pd.DataFrame(matrix, index=phases, columns=phases)
        df.to_csv(output_dir / f"{prefix}_{key}.csv")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute reorganization metrics from cached LRG results."
    )
    parser.add_argument("--patients", nargs="+", default=None, help="Patient IDs (default: auto-detect)")
    parser.add_argument("--bands", nargs="+", default=None, help="Band names (default: all)")
    parser.add_argument("--phases", nargs="+", default=None, help="Phase names (default: all)")
    parser.add_argument("--fc-method", choices=["msc", "corr"], default="msc")
    parser.add_argument("--cache-root", type=Path, default=Path("data/lrg_cache"))
    parser.add_argument("--output-root", type=Path, default=Path("results/reorganization"))
    parser.add_argument("--distance-metric", default="euclidean")
    parser.add_argument("--skip-csv", action="store_true")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    patients = args.patients or list_patients(Path("data/stereoeeg_patients"))
    bands = args.bands or list(BRAIN_BANDS.keys())
    phases = args.phases or list(PHASE_LABELS)

    metric_specs = build_metric_specs(distance_metric=args.distance_metric)

    for patient in patients:
        if args.verbose:
            print(f"\nPatient: {patient}")
        for band in bands:
            results_by_phase = _load_phase_results(
                patient, band, phases, args.fc_method, args.cache_root
            )
            if not results_by_phase:
                if args.verbose:
                    print(f"  [skip] Missing LRG cache for {band}")
                continue

            payload: Dict[str, np.ndarray] = {}

            for key, spec in metric_specs.items():
                matrix = compute_metric_matrix(phases, results_by_phase, spec["fn"])
                payload[key] = matrix

            labels_by_phase = {
                phase: compute_cluster_labels(
                    results_by_phase[phase].linkage_matrix,
                    results_by_phase[phase].optimal_threshold,
                )
                for phase in phases
            }

            payload["cluster_swap"] = compute_cluster_swap_matrix(phases, labels_by_phase)
            payload["cluster_ari"] = compute_ari_matrix(phases, labels_by_phase)

            output_path = (
                args.output_root / patient / f"{band}_{args.fc_method}_metrics.npz"
            )
            _save_npz(output_path, payload, phases)
            if not args.skip_csv:
                _save_csvs(
                    args.output_root / patient,
                    f"{band}_{args.fc_method}",
                    payload,
                    phases,
                )
            if args.verbose:
                print(f"  saved: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
