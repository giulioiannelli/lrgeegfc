#!/usr/bin/env python3
"""Audit Step 23 — per-measure double-check protocol (Phase 0).

For each surviving measure (verdict ``current`` in
``measure_audit_n10_imcoh_abs.csv``), hand-recompute one cell from raw
LRG inputs and compare to the cached row to within float-tolerance. The
recomputation does NOT call the producer script; it uses only library
helpers (``compute_vi``, ``tree_internal_nodes``, ``jaccard_leafsets``,
``scipy.stats.spearmanr``) to derive the predicate from scratch.

Smoke cell: ``Pat_02 / alpha`` (cohort-typical). For sparse-by-design
measures (e.g. trace_modules J=0.9), uses Pat_02 / delta (the only cell
known to have a non-null T-regime row).

Outputs:
    data/audit/measure_audit/double_check_n10_imcoh_abs.csv

Each row: (measure_id, cell, hand_value, cached_value, abs_diff,
rel_diff, tolerance, passed, notes).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.utils.metrics.tree import (
    dmax_from_Z, fcluster_at_h_rel, jaccard_leafsets, tree_internal_nodes,
)
from lrg_eegfc.utils.metrics.vi import compute_vi
from lrg_eegfc.workflow.lrg import load_lrg_result


HALVES_CACHE = ROOT / "data/cache/imcoh_lrg_halves"


def _load_halves(patient: str, phase: str, band: str, half: str) -> tuple[np.ndarray, np.ndarray]:
    """Load split-half LRG result. half ∈ {'A', 'B'}, phase ∈ {'rest_pre', 'rest_post'}."""
    path = HALVES_CACHE / patient / f"{band}_{phase}_{half}_lrg_imcoh-abs.npz"
    if not path.exists():
        raise FileNotFoundError(f"halves cache miss: {path}")
    z = np.load(path, allow_pickle=False)
    D = np.asarray(z["ultrametric_matrix"])
    if D.ndim == 1:
        D = squareform(D)
    Z = np.asarray(z["linkage_matrix"])
    return D, Z


# ---------------------------------------------------------------------------
# Smoke cell + per-measure double-check protocols
# ---------------------------------------------------------------------------

DEFAULT_PATIENT = "Pat_02"
DEFAULT_BAND = "alpha"
TOL_ABS = 1e-6
TOL_REL = 1e-4


def _load_D_Z(patient: str, phase: str, band: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (square ultrametric matrix D, scipy linkage Z)."""
    r = load_lrg_result(patient, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    if r is None:
        raise FileNotFoundError(f"LRG cache miss for {patient}/{band}/{phase}")
    D = np.asarray(r.ultrametric_matrix)
    if D.ndim == 1:
        D = squareform(D)
    Z = np.asarray(r.linkage_matrix)
    return D, Z


def _upper_tri(M: np.ndarray) -> np.ndarray:
    n = M.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    return M[iu, ju]


# --- L1 H2c shared --------------------------------------------------------

def doublecheck_h2c(patient: str, band: str) -> tuple[float, float, str]:
    """Hand-recompute Spearman ρ on (Δ_task, Δ_rest) shared baseline.

    Returns (hand, cached, cell_label).
    """
    D_pre, _ = _load_D_Z(patient, "rest_pre", band)
    D_test, _ = _load_D_Z(patient, "task_test", band)
    D_post, _ = _load_D_Z(patient, "rest_post", band)
    delta_task = _upper_tri(D_test - D_pre)
    delta_rest = _upper_tri(D_post - D_pre)
    rho_hand, _ = spearmanr(delta_task, delta_rest)

    cached = pd.read_csv(ROOT / "data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv")
    row = cached[(cached.patient == patient) & (cached.band == band)]
    rho_cached = float(row["rho_task"].iloc[0])
    return float(rho_hand), rho_cached, f"{patient}/{band}"


# --- L5(k) partition multiscale -------------------------------------------

def doublecheck_partition_multiscale(patient: str, band: str, k: int = 28) -> tuple[float, float, str]:
    """Hand-recompute Δ_VI(k) = VI(c_pre, c_post) − VI(c_test, c_post).

    Returns (hand, cached, cell_label).
    """
    _, Z_pre = _load_D_Z(patient, "rest_pre", band)
    _, Z_test = _load_D_Z(patient, "task_test", band)
    _, Z_post = _load_D_Z(patient, "rest_post", band)
    c_pre = fcluster(Z_pre, k, criterion="maxclust")
    c_test = fcluster(Z_test, k, criterion="maxclust")
    c_post = fcluster(Z_post, k, criterion="maxclust")
    d_VI_hand = compute_vi(c_pre, c_post) - compute_vi(c_test, c_post)

    cached = pd.read_csv(ROOT / "data/reports/imcoh_vi/h2_partition_multiscale_raw.csv")
    row = cached[(cached.patient == patient) & (cached.band == band) & (cached.k == k)]
    if row.empty:
        return float(d_VI_hand), float("nan"), f"{patient}/{band}/k={k} (no cache row)"
    d_VI_cached = float(row["d_VI"].iloc[0])
    return float(d_VI_hand), d_VI_cached, f"{patient}/{band}/k={k}"


# --- L5(k) drift floor ----------------------------------------------------

def doublecheck_dvi_split_baseline(patient: str, band: str, k: int = 28) -> tuple[float, float, str]:
    """Hand-recompute d_VI_full at (patient, band, k=28).

    The split-baseline CSV's `d_VI_full` column should equal the headline
    Δ_VI(k) (same predicate, just paired against the half-data variants
    in the same file).
    """
    _, Z_pre = _load_D_Z(patient, "rest_pre", band)
    _, Z_test = _load_D_Z(patient, "task_test", band)
    _, Z_post = _load_D_Z(patient, "rest_post", band)
    c_pre = fcluster(Z_pre, k, criterion="maxclust")
    c_test = fcluster(Z_test, k, criterion="maxclust")
    c_post = fcluster(Z_post, k, criterion="maxclust")
    d_VI_hand = compute_vi(c_pre, c_post) - compute_vi(c_test, c_post)

    cached = pd.read_csv(ROOT / "data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv")
    row = cached[(cached.patient == patient) & (cached.band == band) & (cached.k == k)]
    if row.empty:
        return float(d_VI_hand), float("nan"), f"{patient}/{band}/k={k} (no cache row)"
    d_VI_cached = float(row["d_VI_full"].iloc[0])
    return float(d_VI_hand), d_VI_cached, f"{patient}/{band}/k={k}"


# --- L4 trace-modules J=0.9 -----------------------------------------------

def doublecheck_trace_modules(patient: str, band: str, k_lo: int = 23, k_hi: int = 31, j_min: float = 0.9) -> tuple[float, float, str]:
    """Hand-recompute the count of T-regime subtrees at (patient, band, k_lo..k_hi).

    Mirrors the producer's scope from `audit_15_trace_modules.py:DEFAULT_K_BINS`:
    enumerates candidate subtrees from T_test ∪ T_learn via fcluster at k cuts
    in [k_lo, k_hi] (not all internal nodes), then classifies as T iff:
      max-Jaccard(V_task, T_pre) < j_min  AND  max-Jaccard(V_task, T_post) >= j_min.

    Hand value = number of T candidates. Cached value = number of rows in
    trace_subtrees_n10_imcoh_abs.csv at this (patient, band) cell.
    """
    _, Z_pre = _load_D_Z(patient, "rest_pre", band)
    _, Z_test = _load_D_Z(patient, "task_test", band)
    _, Z_post = _load_D_Z(patient, "rest_post", band)
    _, Z_learn = _load_D_Z(patient, "task_learn", band)

    nodes_pre = tree_internal_nodes(Z_pre)
    nodes_post = tree_internal_nodes(Z_post)
    leaves_pre = [n["leaves"] for n in nodes_pre if len(n["leaves"]) >= 3]
    leaves_post = [n["leaves"] for n in nodes_post if len(n["leaves"]) >= 3]

    def _fcluster_leafsets(Z: np.ndarray, k: int) -> list[frozenset[int]]:
        labels = fcluster(Z, k, criterion="maxclust")
        out: dict[int, set[int]] = {}
        for leaf_idx, lab in enumerate(labels):
            out.setdefault(int(lab), set()).add(int(leaf_idx))
        return [frozenset(s) for s in out.values() if len(s) >= 3]

    n_T = 0
    for Z in (Z_test, Z_learn):
        for k_cut in range(k_lo, k_hi + 1):
            for S in _fcluster_leafsets(Z, k_cut):
                jp_max = max((jaccard_leafsets(S, T) for T in leaves_pre), default=0.0)
                jq_max = max((jaccard_leafsets(S, T) for T in leaves_post), default=0.0)
                if jp_max < j_min and jq_max >= j_min:
                    n_T += 1

    cached = pd.read_csv(ROOT / "data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv")
    row = cached[(cached.patient == patient) & (cached.band == band)]
    n_T_cached = len(row)
    return float(n_T), float(n_T_cached), f"{patient}/{band} k=[{k_lo},{k_hi}]"


# --- L1 continuous-trace per_cell (identity check vs h2c_shared) ---------

def doublecheck_continuous_trace_per_cell(patient: str, band: str) -> tuple[float, float, str]:
    """Hand-recompute Spearman ρ and compare to per_cell_summary.csv `rho` column.

    `per_cell_summary.csv` `rho` is the SAME predicate as `h2c_ultrametric_drift_raw.csv`
    `rho_task` — Spearman on (Δ_task, Δ_rest) over upper-triangle, shared baseline.
    """
    D_pre, _ = _load_D_Z(patient, "rest_pre", band)
    D_test, _ = _load_D_Z(patient, "task_test", band)
    D_post, _ = _load_D_Z(patient, "rest_post", band)
    rho_hand, _ = spearmanr(_upper_tri(D_test - D_pre), _upper_tri(D_post - D_pre))

    cached = pd.read_csv(ROOT / "data/reports/imcoh_continuous_trace/per_cell_summary.csv")
    row = cached[(cached.patient == patient) & (cached.band == band)]
    rho_cached = float(row["rho"].iloc[0])
    return float(rho_hand), rho_cached, f"{patient}/{band}"


# --- L1 continuous-trace split-baseline ----------------------------------

def doublecheck_continuous_trace_split_baseline(patient: str, band: str) -> tuple[float, float, str]:
    """Hand-recompute Spearman ρ with disjoint half-baselines (rpre_A vs rpre_B).

    Predicate: Δ_task = D_test − D_rpre_A; Δ_rest = D_rpost − D_rpre_B; ρ = Spearman(...)
    """
    D_test, _ = _load_D_Z(patient, "task_test", band)
    D_post, _ = _load_D_Z(patient, "rest_post", band)
    D_pre_A, _ = _load_halves(patient, "rest_pre", band, "A")
    D_pre_B, _ = _load_halves(patient, "rest_pre", band, "B")
    rho_hand, _ = spearmanr(
        _upper_tri(D_test - D_pre_A), _upper_tri(D_post - D_pre_B),
    )

    cached = pd.read_csv(ROOT / "data/reports/imcoh_continuous_trace/per_cell_summary_split.csv")
    row = cached[(cached.patient == patient) & (cached.band == band)]
    rho_cached = float(row["rho"].iloc[0])
    return float(rho_hand), rho_cached, f"{patient}/{band}"


# --- L1 H2e split-half drift floor + within-rpre reliability -------------

def doublecheck_h2e_rho_null_drift(patient: str, band: str) -> tuple[float, float, str]:
    """Hand-recompute ρ_null_drift = Spearman(D_post_B − D_post_A, D_pre_B − D_pre_A)."""
    D_pre_A, _ = _load_halves(patient, "rest_pre", band, "A")
    D_pre_B, _ = _load_halves(patient, "rest_pre", band, "B")
    D_post_A, _ = _load_halves(patient, "rest_post", band, "A")
    D_post_B, _ = _load_halves(patient, "rest_post", band, "B")
    rho_hand, _ = spearmanr(
        _upper_tri(D_post_B - D_post_A), _upper_tri(D_pre_B - D_pre_A),
    )

    cached = pd.read_csv(ROOT / "data/reports/imcoh_vi/h2e_split_half_rho_raw.csv")
    row = cached[(cached.patient == patient) & (cached.band == band)]
    rho_cached = float(row["rho_null_drift"].iloc[0])
    return float(rho_hand), rho_cached, f"{patient}/{band} ρ_null_drift"


# --- L5(h_rel) partition Δ_VI ---------------------------------------------

def doublecheck_dvi_hrel(patient: str, band: str, h_rel_target: float = 0.59) -> tuple[float, float, str]:
    """Hand-recompute Δ_VI(h_rel) via fcluster_at_h_rel.

    Uses the exact h_rel of the closest cached cell (linspace grid) so that
    the hand-recompute and the cached row evaluate the same h_rel.
    """
    cached = pd.read_csv(ROOT / "data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv")
    sub = cached[
        (cached.patient == patient) & (cached.band == band) & (cached.grid == "linspace")
    ]
    if sub.empty:
        return float("nan"), float("nan"), f"{patient}/{band} (no linspace cache rows)"
    closest_idx = (sub["h_rel"] - h_rel_target).abs().idxmin()
    row = sub.loc[closest_idx]
    h_rel_cached = float(row["h_rel"])  # use the cached cell's exact h_rel
    d_VI_cached = float(row["d_VI"])

    _, Z_pre = _load_D_Z(patient, "rest_pre", band)
    _, Z_test = _load_D_Z(patient, "task_test", band)
    _, Z_post = _load_D_Z(patient, "rest_post", band)
    c_pre = fcluster_at_h_rel(Z_pre, h_rel_cached)
    c_test = fcluster_at_h_rel(Z_test, h_rel_cached)
    c_post = fcluster_at_h_rel(Z_post, h_rel_cached)
    d_VI_hand = compute_vi(c_pre, c_post) - compute_vi(c_test, c_post)

    return float(d_VI_hand), d_VI_cached, f"{patient}/{band}/h_rel={h_rel_cached:.4f}"


# --- L5_kcut_heights: h(k) = Z[-k, 2] -------------------------------------

def doublecheck_kcut_heights(patient: str, band: str, phase: str = "rest_pre", k: int = 28) -> tuple[float, float, str]:
    """Hand-recompute h(k) = Z[-k, 2] (kth-largest merge height)."""
    _, Z = _load_D_Z(patient, phase, band)
    h_k_hand = float(Z[-k, 2])

    cached = pd.read_csv(ROOT / "data/audit/dvi_split_baseline/kcut_heights_n10_imcoh_abs.csv")
    row = cached[
        (cached.patient == patient) & (cached.band == band)
        & (cached.phase == phase) & (cached.k == k)
    ]
    if row.empty:
        return h_k_hand, float("nan"), f"{patient}/{band}/{phase}/k={k} (no cache row)"
    h_k_cached = float(row["h_k"].iloc[0])
    return h_k_hand, h_k_cached, f"{patient}/{band}/{phase}/k={k} h(k)"


# --- L5_dmin_dmax_inventory: dmin = Z[0, 2], dmax = Z[-1, 2] -------------

def doublecheck_dmin(patient: str, band: str, phase: str = "rest_pre") -> tuple[float, float, str]:
    """Hand-recompute dmin = first merge height = Z[0, 2]."""
    _, Z = _load_D_Z(patient, phase, band)
    dmin_hand = float(Z[0, 2])

    cached = pd.read_csv(ROOT / "data/audit/dvi_split_baseline/dmin_dmax_inventory.csv")
    row = cached[
        (cached.patient == patient) & (cached.band == band) & (cached.phase == phase)
    ]
    if row.empty:
        return dmin_hand, float("nan"), f"{patient}/{band}/{phase} (no cache row)"
    dmin_cached = float(row["dmin"].iloc[0])
    return dmin_hand, dmin_cached, f"{patient}/{band}/{phase} dmin"


# --- cohort_metadata: Pat_03 sampling rate is 1024 Hz --------------------

def doublecheck_cohort_metadata_pat03() -> tuple[float, float, str]:
    """Pat_03 sampling rate must be 1024 Hz per the always-flag rule."""
    cached = pd.read_csv(ROOT / "data/audit/cohort_metadata.csv")
    row = cached[cached.patient_id == "Pat_03"]
    fs_cached = float(row["sampling_rate_Hz"].iloc[0])
    return 1024.0, fs_cached, "Pat_03 sampling_rate_Hz"


# --- L7_mspc: structural row-count check -----------------------------------

def doublecheck_mspc_row_count(patient: str, band: str) -> tuple[float, float, str]:
    """Structural check: rows in MSPC = n_leaves × |k_grid| for this (patient, band)."""
    _, Z = _load_D_Z(patient, "rest_pre", band)
    n_leaves_hand = int(Z.shape[0]) + 1   # n-1 merges → n leaves

    cached = pd.read_csv(ROOT / "data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv")
    sub = cached[(cached.patient == patient) & (cached.band == band)]
    n_leaves_in_cache = sub["leaf_id"].nunique()
    return float(n_leaves_hand), float(n_leaves_in_cache), f"{patient}/{band} n_leaves"


# --- L6_cnp: structural row-count check -----------------------------------

def doublecheck_cnp_row_count(patient: str, band: str) -> tuple[float, float, str]:
    """Structural check: CNP has one row per leaf for this (patient, band)."""
    _, Z = _load_D_Z(patient, "rest_pre", band)
    n_leaves_hand = int(Z.shape[0]) + 1

    cached = pd.read_csv(ROOT / "data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv")
    sub = cached[(cached.patient == patient) & (cached.band == band)]
    n_leaves_in_cache = len(sub)
    return float(n_leaves_hand), float(n_leaves_in_cache), f"{patient}/{band} n_leaves"


# --- L4_aux_cohesion_cbr: structural row-count check ---------------------

def doublecheck_cohesion_cbr_row_count(patient: str, band: str) -> tuple[float, float, str]:
    """Structural check: Cohesion-CBR has one row per leaf for this (patient, band)."""
    _, Z = _load_D_Z(patient, "rest_pre", band)
    n_leaves_hand = int(Z.shape[0]) + 1

    cached = pd.read_csv(ROOT / "data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv")
    sub = cached[(cached.patient == patient) & (cached.band == band)]
    n_leaves_in_cache = len(sub)
    return float(n_leaves_hand), float(n_leaves_in_cache), f"{patient}/{band} n_leaves"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

CHECKS = [
    # Original 4
    ("L1_h2c_shared", doublecheck_h2c, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L5k_partition_multiscale", doublecheck_partition_multiscale, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND, k=28)),
    ("L5k_drift_floor", doublecheck_dvi_split_baseline, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND, k=28)),
    # trace-modules: delta k=23-31 (canonical bin from audit_15 DEFAULT_K_BINS)
    ("L4_trace_modules_J090", doublecheck_trace_modules, dict(patient=DEFAULT_PATIENT, band="delta", k_lo=23, k_hi=31)),
    # v2 extensions (Phase 0 0d follow-up)
    ("L1_continuous_trace_per_cell", doublecheck_continuous_trace_per_cell, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L1_continuous_trace_split_baseline", doublecheck_continuous_trace_split_baseline, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L1_h2e_rho_null_drift", doublecheck_h2e_rho_null_drift, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L5_hrel_partition", doublecheck_dvi_hrel, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND, h_rel_target=0.59)),
    ("L5_kcut_heights", doublecheck_kcut_heights, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND, phase="rest_pre", k=28)),
    ("L5_dmin_inventory", doublecheck_dmin, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND, phase="rest_pre")),
    ("cohort_metadata_pat03_fs", doublecheck_cohort_metadata_pat03, {}),
    ("L7_mspc_row_count", doublecheck_mspc_row_count, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L6_cnp_row_count", doublecheck_cnp_row_count, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
    ("L4_aux_cohesion_cbr_row_count", doublecheck_cohesion_cbr_row_count, dict(patient=DEFAULT_PATIENT, band=DEFAULT_BAND)),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/audit/measure_audit")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "double_check_n10_imcoh_abs.csv"

    rows = []
    for measure_id, fn, kwargs in CHECKS:
        try:
            hand, cached, label = fn(**kwargs)
            abs_diff = abs(hand - cached) if not np.isnan(cached) else float("nan")
            rel_diff = abs_diff / max(abs(cached), TOL_ABS) if not np.isnan(cached) else float("nan")
            passed = (
                not np.isnan(cached)
                and abs_diff <= max(TOL_ABS, TOL_REL * abs(cached))
            )
            note = "OK" if passed else (
                "no cached row" if np.isnan(cached)
                else f"diff above tol (|Δ|={abs_diff:.3g})"
            )
        except Exception as e:
            hand = cached = abs_diff = rel_diff = float("nan")
            label = ""
            passed = False
            note = f"exception: {type(e).__name__}: {e}"
        rows.append({
            "measure_id": measure_id,
            "cell": label,
            "hand_value": hand,
            "cached_value": cached,
            "abs_diff": abs_diff,
            "rel_diff": rel_diff,
            "tolerance": max(TOL_ABS, TOL_REL * max(abs(cached) if not np.isnan(cached) else 0, TOL_ABS)),
            "passed": passed,
            "notes": note,
        })
        if args.verbose:
            mark = "✓" if passed else "✗"
            print(f"  [{mark}] {measure_id:<32s}  hand={hand:.6g}  cached={cached:.6g}  diff={abs_diff:.3g}  ({note})")

    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    n_pass = int(df["passed"].sum())
    print(f"[audit_23] {n_pass}/{len(df)} double-checks passed")
    print(f"[audit_23] wrote: {csv_path.relative_to(ROOT)}")
    return 0 if n_pass == len(df) else 1


if __name__ == "__main__":
    sys.exit(main())
