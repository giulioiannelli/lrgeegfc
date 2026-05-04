#!/usr/bin/env python
"""Functional tree distance (δ_S, δ_P) curves on the primitive 4-patient set.

Scope:
    .agents/guides/task-persistence-investigation/2026-04-28_functional-tree-distance.md

For each (patient ∈ {Pat_02, Pat_03, Pat_06, Pat_13}, band, phase-pair
(A, B)):
- Compute δ_S(τ; A, B) — Spearman rank-distance on `D = 1/ρ`.
- Compute δ_P(τ; A, B) — Pearson value-distance on `ρ`.
- Snapshot scalars at `τ = τ_min = max(1/λ_max^A, 1/λ_max^B)`.
- Log-τ-normalised integrals over `I = [τ_min, τ_max]`,
  `τ_max = min(τ*_A, τ*_B)`.

Cache layout (one .npz per (patient, band, phase-pair)):
    data/cache/functional_tree_distance/<patient>/<band>_<phaseA>-<phaseB>.npz

Summary CSV (one row per evaluated cell):
    data/cache/functional_tree_distance/summary_n4_imcoh_abs.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import ALL_PHASE_PAIRS, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.metrics import compute_functional_tree_distance

PRIMITIVE_PATIENTS: Tuple[str, ...] = ("Pat_02", "Pat_03", "Pat_06", "Pat_13")
DEFAULT_FC_METHOD = "imcoh_abs"
OUT_ROOT = CACHE_ROOT / "functional_tree_distance"


def _save_npz(out_path: Path, result) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_path,
        patient=result.patient,
        band=result.band,
        phase_A=result.phase_A,
        phase_B=result.phase_B,
        fc_method=result.fc_method,
        tau_grid=result.tau_grid,
        delta_S=result.delta_S,
        delta_P=result.delta_P,
        delta_hat_S=result.delta_hat_S,
        delta_hat_P=result.delta_hat_P,
        Delta_S=result.Delta_S,
        Delta_P=result.Delta_P,
        tau_min=result.tau_min,
        tau_max=result.tau_max,
        tau_star_A=result.tau_star_A,
        tau_star_B=result.tau_star_B,
        lambda_max_A=result.lambda_max_A,
        lambda_max_B=result.lambda_max_B,
        n_nodes=result.n_nodes,
        interval_compatible=result.interval_compatible,
        n_tau=result.n_tau,
        savgol_window=result.savgol_window,
        common_giant_phases=np.array(result.extra["common_giant_phases"]),
        interior_A=result.extra["interior_A"],
        interior_B=result.extra["interior_B"],
    )


def _row_from_result(result) -> dict:
    return {
        "patient": result.patient,
        "band": result.band,
        "phase_A": result.phase_A,
        "phase_B": result.phase_B,
        "fc_method": result.fc_method,
        "n_nodes": result.n_nodes,
        "lambda_max_A": result.lambda_max_A,
        "lambda_max_B": result.lambda_max_B,
        "tau_star_A": result.tau_star_A,
        "tau_star_B": result.tau_star_B,
        "tau_min": result.tau_min,
        "tau_max": result.tau_max,
        "interval_compatible": int(result.interval_compatible),
        "delta_hat_S": result.delta_hat_S,
        "delta_hat_P": result.delta_hat_P,
        "Delta_S": result.Delta_S,
        "Delta_P": result.Delta_P,
        "interior_A": int(result.extra["interior_A"]),
        "interior_B": int(result.extra["interior_B"]),
        "n_tau": result.n_tau,
        "savgol_window": result.savgol_window,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--patients",
        nargs="+",
        default=list(PRIMITIVE_PATIENTS),
        help="Patient identifiers (default: primitive set Pat_02/03/06/13).",
    )
    parser.add_argument(
        "--bands",
        nargs="+",
        default=list(BRAIN_BANDS_NAMES),
        help="Bands (default: all in BRAIN_BANDS_NAMES).",
    )
    parser.add_argument(
        "--fc-method",
        default=DEFAULT_FC_METHOD,
        help="FC method (default: imcoh_abs).",
    )
    parser.add_argument(
        "--n-tau",
        type=int,
        default=100,
        help="Log-spaced τ-grid size in the meaningful interval (default: 100).",
    )
    parser.add_argument(
        "--savgol-window",
        type=int,
        default=5,
        help="Savitzky-Golay smoothing window for τ* (default: 5).",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=OUT_ROOT,
        help=f"Cache root (default: {OUT_ROOT}).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute and overwrite existing .npz files.",
    )
    parser.add_argument(
        "--common-giant-all-phases",
        action="store_true",
        help=(
            "Restrict to the giant-component intersection across all four "
            "phases (default: only the two phases in each pair)."
        ),
    )
    args = parser.parse_args()

    pairs: List[Tuple[str, str]] = list(ALL_PHASE_PAIRS)
    out_root: Path = args.out_root
    out_root.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    common_phases = (
        ("rest_pre", "task_learn", "task_test", "rest_post")
        if args.common_giant_all_phases
        else None
    )

    total = len(args.patients) * len(args.bands) * len(pairs)
    counter = 0
    for patient in args.patients:
        for band in args.bands:
            for phase_A, phase_B in pairs:
                counter += 1
                out_path = (
                    out_root / patient / f"{band}_{phase_A}-{phase_B}.npz"
                )
                if out_path.exists() and not args.overwrite:
                    print(
                        f"[{counter}/{total}] {patient}/{band}/"
                        f"{phase_A}-{phase_B}: cached"
                    )
                    data = np.load(out_path, allow_pickle=True)
                    rows.append(
                        {
                            "patient": str(data["patient"]),
                            "band": str(data["band"]),
                            "phase_A": str(data["phase_A"]),
                            "phase_B": str(data["phase_B"]),
                            "fc_method": str(data["fc_method"]),
                            "n_nodes": int(data["n_nodes"]),
                            "lambda_max_A": float(data["lambda_max_A"]),
                            "lambda_max_B": float(data["lambda_max_B"]),
                            "tau_star_A": float(data["tau_star_A"]),
                            "tau_star_B": float(data["tau_star_B"]),
                            "tau_min": float(data["tau_min"]),
                            "tau_max": float(data["tau_max"]),
                            "interval_compatible": int(
                                data["interval_compatible"]
                            ),
                            "delta_hat_S": float(data["delta_hat_S"]),
                            "delta_hat_P": float(data["delta_hat_P"]),
                            "Delta_S": float(data["Delta_S"]),
                            "Delta_P": float(data["Delta_P"]),
                            "interior_A": int(data["interior_A"]),
                            "interior_B": int(data["interior_B"]),
                            "n_tau": int(data["n_tau"]),
                            "savgol_window": int(data["savgol_window"]),
                        }
                    )
                    continue

                print(
                    f"[{counter}/{total}] {patient}/{band}/"
                    f"{phase_A}-{phase_B}: computing"
                )
                try:
                    result = compute_functional_tree_distance(
                        patient=patient,
                        band=band,
                        phase_A=phase_A,
                        phase_B=phase_B,
                        fc_method=args.fc_method,
                        n_tau=args.n_tau,
                        savgol_window=args.savgol_window,
                        common_giant_phases=common_phases,
                    )
                except Exception as exc:
                    print(f"    FAILED: {exc}")
                    rows.append(
                        {
                            "patient": patient,
                            "band": band,
                            "phase_A": phase_A,
                            "phase_B": phase_B,
                            "fc_method": args.fc_method,
                            "n_nodes": -1,
                            "lambda_max_A": np.nan,
                            "lambda_max_B": np.nan,
                            "tau_star_A": np.nan,
                            "tau_star_B": np.nan,
                            "tau_min": np.nan,
                            "tau_max": np.nan,
                            "interval_compatible": 0,
                            "delta_hat_S": np.nan,
                            "delta_hat_P": np.nan,
                            "Delta_S": np.nan,
                            "Delta_P": np.nan,
                            "interior_A": 0,
                            "interior_B": 0,
                            "n_tau": 0,
                            "savgol_window": args.savgol_window,
                            "error": str(exc),
                        }
                    )
                    continue

                _save_npz(out_path, result)
                rows.append(_row_from_result(result))

    summary = pd.DataFrame(rows)
    summary_path = (
        out_root
        / f"summary_n{len(args.patients)}_{args.fc_method}.csv"
    )
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote summary CSV: {summary_path}  ({len(summary)} rows)")


if __name__ == "__main__":
    main()
