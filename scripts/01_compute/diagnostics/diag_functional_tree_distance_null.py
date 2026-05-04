#!/usr/bin/env python3
"""Within-baseline split-half null for the functional tree distance.

For each (patient, band) in the primitive set:

1. Load the full-data FTD result to get the anchor τ-grid `[τ_min, τ_max]`.
2. Load each phase's timeseries, split at T//2, compute |ImCoh| FC for each
   half (existing `compute_imcoh_abs_halves` helper, halved nperseg).
3. Restrict each half-FC to the common giant-component electrodes already
   used by the full-data computation (read from the cached `.npz`).
4. Run the FTD machinery on the half-Laplacians at the *anchor* τ-grid:
   - **Within-phase Δ:** `Δ(half_A, half_B)` per phase  → 4 noise-floor values.
   - **Between-phase Δ at half level:** `Δ(phase_X half_a, phase_Y half_a)`
     for both `a ∈ {A, B}` and all 6 phase-pairs → 12 between-phase values.
5. Save per-cell .npz + summary CSV.

Compares against the observed full-data Δ values (already cached) so the
question becomes: is the cross-phase Δ above the within-phase noise floor?
For the contrast `C = Δ(test, post) − Δ(pre, post)`, also report the
half-level contrast under the within-phase scale.
"""
from __future__ import annotations

import argparse
import gc
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import networkx as nx
import pandas as pd

import sys
HYP_DIR = Path(__file__).resolve().parents[1] / "hypothesis_tests"
if str(HYP_DIR) not in sys.path:
    sys.path.insert(0, str(HYP_DIR))

from lrg_eegfc.config.const import (
    ALL_PHASE_PAIRS,
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    FS_OVERRIDES,
    PHASE_LABELS,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.functional_tree_distance import (
    delta_integrals_from_adjacencies,
)
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore

PRIMITIVE_PATIENTS: Tuple[str, ...] = ("Pat_02", "Pat_03", "Pat_06", "Pat_13")
DEFAULT_FC_METHOD = "imcoh_abs"
FTD_ROOT = CACHE_ROOT / "functional_tree_distance"
NULL_ROOT = CACHE_ROOT / "functional_tree_distance_null"


def _load_anchor(patient: str, band: str, phase_A: str, phase_B: str) -> Optional[dict]:
    """Pull the full-data anchor τ-grid + node indices for a phase-pair cell."""
    path = FTD_ROOT / patient / f"{band}_{phase_A}-{phase_B}.npz"
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    return {
        "tau_grid": np.asarray(data["tau_grid"]),
        "n_nodes": int(data["n_nodes"]),
        "Delta_S_obs": float(data["Delta_S"]),
        "Delta_P_obs": float(data["Delta_P"]),
        "tau_min": float(data["tau_min"]),
        "tau_max": float(data["tau_max"]),
    }


def _common_giant_indices_from_phases(
    fc_per_phase: Dict[str, np.ndarray],
) -> np.ndarray:
    sets = []
    for adj in fc_per_phase.values():
        G = nx.from_numpy_array(adj)
        components = list(nx.connected_components(G))
        if not components:
            sets.append(set())
            continue
        giant = max(components, key=len)
        sets.append(set(giant))
    if not sets:
        return np.array([], dtype=int)
    common = set.intersection(*sets)
    return np.array(sorted(common), dtype=int)


def _load_X(patient: str, phase: str) -> Optional[Tuple[np.ndarray, float]]:
    try:
        X = load_timeseries(patient, phase, SEEG_DATAPATH)
    except (FileNotFoundError, OSError):
        return None
    if X is None:
        return None
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(patient, 2048.0)
    return X, fs


def _compute_halves_all_phases(
    patient: str,
    bands: List[str],
    verbose: bool = False,
) -> Optional[Dict[str, Dict[Tuple[str, str], np.ndarray]]]:
    """Return {phase: {(band, tag): FC}}; tag ∈ {A, B}."""
    out: Dict[str, Dict[Tuple[str, str], np.ndarray]] = {}
    band_dict = {b: BRAIN_BANDS[b] for b in bands}
    for phase in PHASE_LABELS:
        loaded = _load_X(patient, phase)
        if loaded is None:
            if verbose:
                print(f"  {patient}/{phase}: timeseries missing, skip")
            return None
        X, fs = loaded
        nperseg_half = max(256, nperseg_for_fs(fs) // 2)
        if verbose:
            print(
                f"  {patient}/{phase}: X.shape={X.shape}  "
                f"nperseg_half={nperseg_half}"
            )
        out[phase] = compute_imcoh_abs_halves(X, fs, nperseg_half, band_dict)
        del X
        gc.collect()
    return out


def _restrict(adj: np.ndarray, idx: np.ndarray) -> np.ndarray:
    return adj[np.ix_(idx, idx)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--patients", nargs="+", default=list(PRIMITIVE_PATIENTS))
    parser.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS_NAMES))
    parser.add_argument(
        "--out-root",
        type=Path,
        default=NULL_ROOT,
        help=f"Cache root (default: {NULL_ROOT}).",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []

    for patient in args.patients:
        print(f"=== {patient}: computing halves for all 4 phases ===")
        halves_by_phase = _compute_halves_all_phases(
            patient, args.bands, verbose=args.verbose
        )
        if halves_by_phase is None:
            print(f"  {patient}: missing timeseries; skip")
            continue

        for band in args.bands:
            try:
                fc_per_phase: Dict[str, np.ndarray] = {}
                for phase in PHASE_LABELS:
                    fc_A = halves_by_phase[phase].get((band, "A"))
                    fc_B = halves_by_phase[phase].get((band, "B"))
                    if fc_A is None or fc_B is None:
                        raise KeyError(f"{patient}/{phase}/{band}: half FC missing")
                    fc_per_phase[(phase, "A")] = fc_A
                    fc_per_phase[(phase, "B")] = fc_B
            except KeyError as exc:
                print(f"  {patient}/{band}: {exc}; skip")
                continue

            common_idx = _common_giant_indices_from_phases(fc_per_phase)
            if common_idx.size < 4:
                print(
                    f"  {patient}/{band}: common giant component too small "
                    f"({common_idx.size}); skip"
                )
                continue

            fc_restricted: Dict[Tuple[str, str], np.ndarray] = {
                k: _restrict(v, common_idx) for k, v in fc_per_phase.items()
            }

            anchor = _load_anchor(patient, band, "rest_pre", "task_test")
            if anchor is None:
                print(
                    f"  {patient}/{band}: full-data anchor missing; skip"
                )
                continue
            tau_grid = anchor["tau_grid"]

            within_phase: Dict[str, dict] = {}
            for phase in PHASE_LABELS:
                res = delta_integrals_from_adjacencies(
                    fc_restricted[(phase, "A")],
                    fc_restricted[(phase, "B")],
                    tau_grid,
                )
                within_phase[phase] = {
                    "Delta_S": res["Delta_S"],
                    "Delta_P": res["Delta_P"],
                    "delta_S": res["delta_S"],
                    "delta_P": res["delta_P"],
                }

            between_phase: Dict[Tuple[str, str, str], dict] = {}
            for phase_A, phase_B in ALL_PHASE_PAIRS:
                for tag in ("A", "B"):
                    res = delta_integrals_from_adjacencies(
                        fc_restricted[(phase_A, tag)],
                        fc_restricted[(phase_B, tag)],
                        tau_grid,
                    )
                    between_phase[(phase_A, phase_B, tag)] = {
                        "Delta_S": res["Delta_S"],
                        "Delta_P": res["Delta_P"],
                    }

            out_path = args.out_root / patient / f"{band}_null.npz"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                out_path,
                patient=patient,
                band=band,
                fc_method=DEFAULT_FC_METHOD,
                tau_grid=tau_grid,
                tau_min=anchor["tau_min"],
                tau_max=anchor["tau_max"],
                n_nodes=int(common_idx.size),
                within_phase_Delta_S=np.array(
                    [within_phase[p]["Delta_S"] for p in PHASE_LABELS]
                ),
                within_phase_Delta_P=np.array(
                    [within_phase[p]["Delta_P"] for p in PHASE_LABELS]
                ),
                within_phase_labels=np.array(list(PHASE_LABELS)),
                between_phase_keys=np.array(
                    [f"{a}|{b}|{t}" for (a, b, t) in between_phase.keys()]
                ),
                between_phase_Delta_S=np.array(
                    [v["Delta_S"] for v in between_phase.values()]
                ),
                between_phase_Delta_P=np.array(
                    [v["Delta_P"] for v in between_phase.values()]
                ),
            )

            within_S_med = float(
                np.median([within_phase[p]["Delta_S"] for p in PHASE_LABELS])
            )
            within_P_med = float(
                np.median([within_phase[p]["Delta_P"] for p in PHASE_LABELS])
            )

            anchor_pp = _load_anchor(patient, band, "rest_pre", "rest_post")
            anchor_tp = _load_anchor(patient, band, "task_test", "rest_post")

            row = {
                "patient": patient,
                "band": band,
                "n_nodes": int(common_idx.size),
                "within_S_median": within_S_med,
                "within_P_median": within_P_med,
                "obs_Delta_S_pre_post": (
                    anchor_pp["Delta_S_obs"] if anchor_pp else np.nan
                ),
                "obs_Delta_P_pre_post": (
                    anchor_pp["Delta_P_obs"] if anchor_pp else np.nan
                ),
                "obs_Delta_S_test_post": (
                    anchor_tp["Delta_S_obs"] if anchor_tp else np.nan
                ),
                "obs_Delta_P_test_post": (
                    anchor_tp["Delta_P_obs"] if anchor_tp else np.nan
                ),
            }
            row["obs_C_S"] = (
                row["obs_Delta_S_test_post"] - row["obs_Delta_S_pre_post"]
            )
            row["obs_C_P"] = (
                row["obs_Delta_P_test_post"] - row["obs_Delta_P_pre_post"]
            )
            for phase in PHASE_LABELS:
                row[f"within_S_{phase}"] = within_phase[phase]["Delta_S"]
                row[f"within_P_{phase}"] = within_phase[phase]["Delta_P"]

            tt_pp_AA = (
                between_phase[("task_test", "rest_post", "A")]["Delta_S"]
                - between_phase[("rest_pre", "rest_post", "A")]["Delta_S"]
            )
            tt_pp_BB = (
                between_phase[("task_test", "rest_post", "B")]["Delta_S"]
                - between_phase[("rest_pre", "rest_post", "B")]["Delta_S"]
            )
            row["half_C_S_A"] = tt_pp_AA
            row["half_C_S_B"] = tt_pp_BB

            tt_pp_AA_p = (
                between_phase[("task_test", "rest_post", "A")]["Delta_P"]
                - between_phase[("rest_pre", "rest_post", "A")]["Delta_P"]
            )
            tt_pp_BB_p = (
                between_phase[("task_test", "rest_post", "B")]["Delta_P"]
                - between_phase[("rest_pre", "rest_post", "B")]["Delta_P"]
            )
            row["half_C_P_A"] = tt_pp_AA_p
            row["half_C_P_B"] = tt_pp_BB_p

            rows.append(row)
            print(
                f"  {patient}/{band}: "
                f"within_S_med={within_S_med:.3f}  C_S_obs={row['obs_C_S']:+.3f}  "
                f"C_S_halves=({tt_pp_AA:+.3f}, {tt_pp_BB:+.3f})"
            )

    summary = pd.DataFrame(rows)
    summary_path = (
        args.out_root
        / f"summary_n{len(args.patients)}_{DEFAULT_FC_METHOD}.csv"
    )
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote summary CSV: {summary_path}  ({len(summary)} rows)")


if __name__ == "__main__":
    main()
