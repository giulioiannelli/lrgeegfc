#!/usr/bin/env python3
"""Audit Step 24 — cohort-coverage matrix v1 (rebuild plan §10 step 3).

Builds the rebuild's central deliverable: three CSVs that materialise
the per-cell verdict, the per-(band, rung) aggregate, and the per-band
triangulation, plus the headline heatmap PDF.

    cohort_coverage_matrix_n10_imcoh_abs.csv  per-cell rows
    band_verdict_n10_imcoh_abs.csv            per (band, rung [, scale_value])
    triangulation_n10_imcoh_abs.csv           per band T(b)
    cohort_coverage_matrix_n10_imcoh_abs.pdf  band x rung heatmap

Implements the data-flow defined in
    .agents/guides/task-persistence-investigation/2026-04-29_cohort-coverage-matrix.md
and the predicates defined in
    .agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md

L2 (entropy curve / specific heat) is intentionally NOT a rung — for our
fully-connected weight-heterogeneous outlier case S(τ) and C(τ) are
non-informative (continuous spectrum, no discrete C(τ) peaks, see
.agents/guides/02_methods/lrg-framework-guide.md §6 + memory entry
lrg_outlier_case_fully_connected.md).

v1 caveat: cells without rung-matched nulls (L3, L4_null, L5_hrel_null,
L6_null, L7_null) emit pi="uncontrolled". v2 (rebuild plan §10 step 16)
replaces hatched cells after §7.2-7.6 controls land.

Naming: audit_24 to avoid collision with the existing
audit_18_taustar_four_phase.py.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, bh_fdr


# =============================================================================
# Decision-rules constants
# Pre-registered in 2026-04-29_decision-rules.md; immutable in v1.
# =============================================================================

COHORT_THRESHOLD = 0.8        # frac_pos >= 0.8  (>= 8/10)
NEGATIVE_THRESHOLD = 0.2      # frac_pos <= 0.2  -> negative
WILCOXON_Q_MAX = 0.05         # BH-FDR q-cutoff per rung, m=6 bands
RIDGE_LENGTH_MIN = 3          # contiguous frac_pos >= 0.8 cells along scale axis
ELIGIBILITY_THRESHOLD = 0.5   # cell ineligible if (n_eligible / |P|) < 0.5

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
RUNG_ORDER = ["L1", "L3", "L4", "L4_aux", "L5_k", "L5_hrel", "L6", "L7", "E1_subspace"]
DISSENTER_MAP = {"Pat_02": "x", "Pat_03": "+"}

# Eigenvector-direct E-series scale grid (pivot plan §3.4).
# Sparse log-ish: {2, 3, 5, 8, 13, 21}.
E1_K_GRID = [2, 3, 5, 8, 13, 21]

# Triangulation short codes for figure annotation
T_SHORT_CODES = {
    "eigenvector-confirmed":         "EC",
    "headline-triangulated":         "HT",
    "headline-triangulated-pending": "HT*",
    "partition-resolution-locked":   "PRL",
    "continuous-only":               "CO",
    "anatomy-suspect":               "AS",
    "ergodic":                       "ER",
    "uncontrolled-exploratory":      "UE",
    "inconsistent":                  "IN",
}
T_COLOURS = {
    "eigenvector-confirmed":         "#1a9850",
    "headline-triangulated":         "#2ca02c",
    "headline-triangulated-pending": "#7fcf78",
    "partition-resolution-locked":   "#d4a017",
    "continuous-only":               "#1f77b4",
    "anatomy-suspect":               "#ff7f0e",
    "ergodic":                       "#7f7f7f",
    "uncontrolled-exploratory":      "#9467bd",
    "inconsistent":                  "#d62728",
}


# =============================================================================
# Consumer paths — read contract from cohort-coverage-matrix scope
# =============================================================================

CONSUMERS: dict[str, dict] = {
    "L1": dict(
        measure="data/reports/imcoh_continuous_trace/per_cell_summary_split.csv",
        null="data/reports/imcoh_continuous_trace/controls_summary.csv",
        measure_audit_id="L1_continuous_trace_split_baseline",
        null_audit_id="L1_continuous_trace_controls",
    ),
    "L4": dict(
        measure="data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv",
        null=None,
        measure_audit_id="L4_trace_modules_J090",
        null_audit_id=None,
    ),
    "L4_aux": dict(
        measure="data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv",
        null=None,
        measure_audit_id="L4_aux_cohesion_cbr",
        null_audit_id=None,
    ),
    "L5_k": dict(
        measure="data/reports/imcoh_vi/h2_partition_multiscale_raw.csv",
        null="data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv",
        measure_audit_id="L5k_partition_multiscale",
        null_audit_id="L5k_drift_floor",
    ),
    "L5_hrel": dict(
        measure="data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv",
        null=None,
        measure_audit_id="L5_hrel_partition",
        null_audit_id=None,
    ),
    "L6": dict(
        measure="data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv",
        null=None,
        measure_audit_id="L6_cnp",
        null_audit_id=None,
    ),
    "L7": dict(
        measure="data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv",
        null=None,
        measure_audit_id="L7_mspc",
        null_audit_id=None,
    ),
    # E1 — first eigenvector-direct rung (pivot plan §3). Measure and null
    # live in the SAME CSV (delta_e1 vs delta_e1_null columns); we set
    # `null` to the same path as a sentinel for `rungs_with_null` and
    # bypass the per-cell file check inside populate_e1_subspace.
    "E1_subspace": dict(
        measure="data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv",
        null="data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv",
        measure_audit_id=None,   # New rung — Phase 0 audit row pending
        null_audit_id=None,
        skip_audit_check=True,
    ),
    # L3 has no producer in v1 -> emit uncontrolled skeleton rows.
    # L2 (entropy/specific heat) is permanently removed — non-informative
    # for our continuous-spectrum outlier case.
}

DMIN_DMAX_PATH          = "data/audit/dvi_split_baseline/dmin_dmax_inventory.csv"
COHORT_METADATA_PATH    = "data/audit/cohort_metadata.csv"
MEASURE_AUDIT_PATH      = "data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv"
AUDIT_RUN_METADATA_PATH = "data/audit/measure_audit/audit_run_metadata.json"

# Output paths
OUT_DIR     = "data/audit/cohort_coverage_matrix"
PER_CELL_CSV    = f"{OUT_DIR}/cohort_coverage_matrix_n10_imcoh_abs.csv"
BAND_VERDICT_CSV = f"{OUT_DIR}/band_verdict_n10_imcoh_abs.csv"
TRIANGULATION_CSV = f"{OUT_DIR}/triangulation_n10_imcoh_abs.csv"
RUN_METADATA    = f"{OUT_DIR}/audit_24_run_metadata.json"
FIGURE_PDF      = "data/outputs/figures/section6/cohort_coverage_matrix_n10_imcoh_abs.pdf"


# =============================================================================
# Preflight
# =============================================================================

def preflight(audit_df: pd.DataFrame, audit_ts: float) -> None:
    """Refuse to run if any consumed measure's audit verdict is not 'current'
    or any input CSV's mtime exceeds the Phase 0 audit timestamp."""
    needed_ids: set[str] = set()
    for rung, c in CONSUMERS.items():
        if c.get("skip_audit_check"):
            continue
        for key in ("measure_audit_id", "null_audit_id"):
            mid = c.get(key)
            if mid is not None:
                needed_ids.add(mid)
    audit_lookup = audit_df.set_index("measure_id")["audit_verdict"].to_dict()
    bad: list[str] = []
    for mid in needed_ids:
        v = audit_lookup.get(mid)
        if v != "current":
            bad.append(f"{mid}: verdict={v!r}")
    if bad:
        sys.exit("Preflight FAIL — measures with non-'current' audit verdict:\n  "
                 + "\n  ".join(bad)
                 + "\nRe-run audit_21_measure_correctness.py and fix.")
    # mtime gate
    stale: list[str] = []
    for rung, c in CONSUMERS.items():
        if c.get("skip_audit_check"):
            continue
        for path in (c["measure"], c["null"]):
            if path is None:
                continue
            p = ROOT / path
            if not p.exists():
                continue
            if p.stat().st_mtime > audit_ts + 1.0:
                stale.append(f"{path} (mtime > audit_ts)")
    if stale:
        sys.exit("Preflight FAIL — inputs newer than Phase 0 audit timestamp:\n  "
                 + "\n  ".join(stale)
                 + "\nRe-run audit_21_measure_correctness.py to refresh verdicts.")


# =============================================================================
# Cartesian skeleton
# =============================================================================

def k_grid_for_rung(rung: str, dmin_df: pd.DataFrame) -> list[int]:
    """Return the integer-k grid used by a scale-axis rung."""
    if rung == "L5_k":
        # h2_partition_multiscale ranges k=2..119 (variable per patient)
        return list(range(2, 120))
    if rung == "L7":
        # MSPC k axis matches L5_k
        return list(range(2, 120))
    return []


def hrel_grid(measure_path: Path) -> list[float]:
    """Return the unique h_rel values from the linspace grid of dvi_hrel."""
    df = pd.read_csv(measure_path)
    return sorted(df.loc[df["grid"] == "linspace", "h_rel"].unique().tolist())


def build_skeleton(patients: list[str], dmin_df: pd.DataFrame,
                   hrel_values: list[float]) -> pd.DataFrame:
    """Cartesian product (patient × band × rung × scale_value)."""
    rows: list[dict] = []
    for p in patients:
        for b in BAND_ORDER:
            for r in RUNG_ORDER:
                if r == "L5_k":
                    for k in k_grid_for_rung("L5_k", dmin_df):
                        rows.append(dict(patient=p, band=b, rung=r,
                                         scale_axis="k", scale_value=float(k)))
                elif r == "L5_hrel":
                    for h in hrel_values:
                        rows.append(dict(patient=p, band=b, rung=r,
                                         scale_axis="h_rel", scale_value=float(h)))
                elif r == "L7":
                    for k in k_grid_for_rung("L7", dmin_df):
                        rows.append(dict(patient=p, band=b, rung=r,
                                         scale_axis="k", scale_value=float(k)))
                elif r == "E1_subspace":
                    for k in E1_K_GRID:
                        rows.append(dict(patient=p, band=b, rung=r,
                                         scale_axis="k", scale_value=float(k)))
                else:
                    rows.append(dict(patient=p, band=b, rung=r,
                                     scale_axis="none", scale_value=np.nan))
    df = pd.DataFrame(rows)
    df["measure_value"] = np.nan
    df["null_value"] = np.nan
    df["passes_null"] = pd.NA
    df["passes_zero"] = pd.NA
    df["passes_count"] = pd.NA
    df["eligible"] = True
    df["pi"] = ""
    df["dissenter_marker"] = df["patient"].map(DISSENTER_MAP).fillna("none")
    df["producer_artefact"] = ""
    df["null_artefact"] = ""
    df["audit_verdict"] = ""
    df["notes"] = ""
    return df


# =============================================================================
# Per-rung populate
# =============================================================================

def populate_l1(skel: pd.DataFrame) -> pd.DataFrame:
    """L1 = ρ_split (Run A) vs ρ_null_drift (Run C)."""
    m = pd.read_csv(ROOT / CONSUMERS["L1"]["measure"])
    n = pd.read_csv(ROOT / CONSUMERS["L1"]["null"])
    assert len(m) == 60, f"L1 measure has {len(m)} rows, expected 60"
    assert len(n) == 60, f"L1 null has {len(n)} rows, expected 60"
    m_lookup = m.set_index(["patient", "band"])["rho"].to_dict()
    n_lookup = n.set_index(["patient", "band"])["rho_null_drift"].to_dict()
    mask = skel["rung"] == "L1"
    for idx in skel[mask].index:
        key = (skel.at[idx, "patient"], skel.at[idx, "band"])
        skel.at[idx, "measure_value"] = m_lookup.get(key, np.nan)
        skel.at[idx, "null_value"]    = n_lookup.get(key, np.nan)
        skel.at[idx, "producer_artefact"] = CONSUMERS["L1"]["measure"]
        skel.at[idx, "null_artefact"]     = CONSUMERS["L1"]["null"]
    return skel


def populate_l5_k(skel: pd.DataFrame, dmin_df: pd.DataFrame) -> pd.DataFrame:
    """L5_k = d_VI (test↔post − pre↔post) vs drift_dVI (half-baseline)."""
    m = pd.read_csv(ROOT / CONSUMERS["L5_k"]["measure"])
    n = pd.read_csv(ROOT / CONSUMERS["L5_k"]["null"])
    m_lookup = m.set_index(["patient", "band", "k"])["d_VI"].to_dict()
    n_lookup = n.set_index(["patient", "band", "k"])["drift_dVI"].to_dict()
    nleaf_lookup = (
        dmin_df[dmin_df["phase"] == "rest_post"]
        .set_index(["patient", "band"])["n_leaves"].to_dict()
    )
    mask = skel["rung"] == "L5_k"
    for idx in skel[mask].index:
        p = skel.at[idx, "patient"]; b = skel.at[idx, "band"]
        k = int(skel.at[idx, "scale_value"])
        skel.at[idx, "measure_value"] = m_lookup.get((p, b, k), np.nan)
        skel.at[idx, "null_value"]    = n_lookup.get((p, b, k), np.nan)
        skel.at[idx, "producer_artefact"] = CONSUMERS["L5_k"]["measure"]
        skel.at[idx, "null_artefact"]     = CONSUMERS["L5_k"]["null"]
        # Eligibility: k ≤ n_leaves(rest_post)
        nl = nleaf_lookup.get((p, b))
        if nl is None or k > nl:
            skel.at[idx, "eligible"] = False
            skel.at[idx, "notes"] = "k > n_leaves(rest_post)"
    return skel


def populate_l5_hrel(skel: pd.DataFrame, dmin_df: pd.DataFrame) -> pd.DataFrame:
    """L5_hrel = d_VI at fixed h_rel (linspace grid). Null missing in v1."""
    m = pd.read_csv(ROOT / CONSUMERS["L5_hrel"]["measure"])
    m = m[m["grid"] == "linspace"]
    m_lookup = m.set_index(["patient", "band", "h_rel"])["d_VI"].to_dict()
    hmin_lookup = (
        dmin_df[dmin_df["phase"] == "rest_post"]
        .set_index(["patient", "band"])["h_rel_min"].to_dict()
    )
    mask = skel["rung"] == "L5_hrel"
    for idx in skel[mask].index:
        p = skel.at[idx, "patient"]; b = skel.at[idx, "band"]
        h = float(skel.at[idx, "scale_value"])
        skel.at[idx, "measure_value"] = m_lookup.get((p, b, h), np.nan)
        skel.at[idx, "producer_artefact"] = CONSUMERS["L5_hrel"]["measure"]
        # Eligibility: h_rel ≥ h_rel_min(rest_post)
        hmin = hmin_lookup.get((p, b))
        if hmin is None or h < hmin:
            skel.at[idx, "eligible"] = False
            skel.at[idx, "notes"] = f"h_rel < h_rel_min={hmin}"
    return skel


def _count_T_per_pb(path: Path, regime_col: str, regime_value: str = "T") -> dict:
    """Count rows where regime_col == regime_value, grouped by (patient, band)."""
    df = pd.read_csv(path)
    df = df[df[regime_col] == regime_value]
    return df.groupby(["patient", "band"]).size().to_dict()


def populate_l4(skel: pd.DataFrame) -> pd.DataFrame:
    """L4 = count of T-modules in trace_subtrees CSV (sparse-by-design)."""
    counts = _count_T_per_pb(ROOT / CONSUMERS["L4"]["measure"], regime_col="regime")
    mask = skel["rung"] == "L4"
    for idx in skel[mask].index:
        key = (skel.at[idx, "patient"], skel.at[idx, "band"])
        skel.at[idx, "measure_value"] = float(counts.get(key, 0))
        skel.at[idx, "producer_artefact"] = CONSUMERS["L4"]["measure"]
    return skel


def populate_l4_aux(skel: pd.DataFrame) -> pd.DataFrame:
    """L4_aux = count of T-leaves in Cohesion-CBR leaf_assignment."""
    counts = _count_T_per_pb(ROOT / CONSUMERS["L4_aux"]["measure"], regime_col="dominant")
    mask = skel["rung"] == "L4_aux"
    for idx in skel[mask].index:
        key = (skel.at[idx, "patient"], skel.at[idx, "band"])
        skel.at[idx, "measure_value"] = float(counts.get(key, 0))
        skel.at[idx, "producer_artefact"] = CONSUMERS["L4_aux"]["measure"]
    return skel


def populate_l6(skel: pd.DataFrame) -> pd.DataFrame:
    """L6 = count of T-leaves in CNP leaf_assignment."""
    counts = _count_T_per_pb(ROOT / CONSUMERS["L6"]["measure"], regime_col="dominant")
    mask = skel["rung"] == "L6"
    for idx in skel[mask].index:
        key = (skel.at[idx, "patient"], skel.at[idx, "band"])
        skel.at[idx, "measure_value"] = float(counts.get(key, 0))
        skel.at[idx, "producer_artefact"] = CONSUMERS["L6"]["measure"]
    return skel


def populate_e1_subspace(skel: pd.DataFrame) -> pd.DataFrame:
    """E1 = Δ^E1_k via Grassmann chordal distance on top-k Laplacian
    eigenspaces. Per-cell positivity: ``Δ^E1 < 0 AND |Δ^E1| > |Δ^E1_null|``.
    The continuous-rungs branch in `assign_pi` evaluates `m > n` on the
    NEGATED magnitudes so the existing ``positive ⇔ m > n AND m > 0``
    contract holds: we map measure_value = −Δ^E1, null_value = |Δ^E1_null|.
    """
    path = ROOT / CONSUMERS["E1_subspace"]["measure"]
    if not path.exists():
        # CSV not yet produced — leave skeleton rows blank (will be flagged
        # ineligible by assign_pi).
        return skel
    df = pd.read_csv(path)
    delta_lookup     = df.set_index(["patient", "band", "k"])["delta_e1"].to_dict()
    delta_null_lookup = df.set_index(["patient", "band", "k"])["delta_e1_null"].to_dict()
    mask = skel["rung"] == "E1_subspace"
    for idx in skel[mask].index:
        p = skel.at[idx, "patient"]; b = skel.at[idx, "band"]
        k = int(skel.at[idx, "scale_value"])
        de = delta_lookup.get((p, b, k))
        dn = delta_null_lookup.get((p, b, k))
        if de is None or dn is None or pd.isna(de) or pd.isna(dn):
            skel.at[idx, "measure_value"] = np.nan
            skel.at[idx, "null_value"]    = np.nan
        else:
            # Trace direction is `Δ < 0`. Map to the continuous-rung
            # gate `m > n AND m > 0` by setting m = -Δ_E1 (positive when
            # trace), n = |Δ_E1_null| (always non-negative).
            skel.at[idx, "measure_value"] = float(-de)
            skel.at[idx, "null_value"]    = float(abs(dn))
        skel.at[idx, "producer_artefact"] = CONSUMERS["E1_subspace"]["measure"]
        skel.at[idx, "null_artefact"]     = CONSUMERS["E1_subspace"]["null"]
    return skel


def populate_l7(skel: pd.DataFrame, dmin_df: pd.DataFrame, chunksize: int = 100_000) -> pd.DataFrame:
    """L7 = per-(p, b, k) count of T-labelled leaves in MSPC, streamed."""
    counts: dict = {}
    for chunk in pd.read_csv(ROOT / CONSUMERS["L7"]["measure"], chunksize=chunksize):
        sub = chunk[chunk["dominant_at_k"] == "T"]
        for (p, b, k), n in sub.groupby(["patient", "band", "k"]).size().items():
            counts[(p, b, int(k))] = counts.get((p, b, int(k)), 0) + int(n)
    nleaf_lookup = (
        dmin_df[dmin_df["phase"] == "rest_post"]
        .set_index(["patient", "band"])["n_leaves"].to_dict()
    )
    mask = skel["rung"] == "L7"
    for idx in skel[mask].index:
        p = skel.at[idx, "patient"]; b = skel.at[idx, "band"]
        k = int(skel.at[idx, "scale_value"])
        skel.at[idx, "measure_value"] = float(counts.get((p, b, k), 0))
        skel.at[idx, "producer_artefact"] = CONSUMERS["L7"]["measure"]
        nl = nleaf_lookup.get((p, b))
        if nl is None or k > nl:
            skel.at[idx, "eligible"] = False
            skel.at[idx, "notes"] = "k > n_leaves(rest_post)"
    return skel


# =============================================================================
# Per-cell positivity (decision-rules §2.1)
# =============================================================================

def assign_pi(skel: pd.DataFrame) -> pd.DataFrame:
    """Decide per-cell pi ∈ {0, 1, ineligible, uncontrolled}.

    'uncontrolled' is reserved for rungs whose NULL ARTEFACT is entirely
    absent in v1 (L3, L4, L4_aux, L5_hrel, L6, L7). 'ineligible' is used
    for individual cells where the cell-level measure or null value is
    missing (e.g. k > n_leaves(half) in dvi_split_baseline)."""
    rungs_with_null = {r for r, c in CONSUMERS.items() if c.get("null") is not None}
    continuous_rungs = {"L1", "L5_k", "E1_subspace"}
    for idx in skel.index:
        if not skel.at[idx, "eligible"]:
            skel.at[idx, "pi"] = "ineligible"
            continue
        rung = skel.at[idx, "rung"]
        m = skel.at[idx, "measure_value"]
        n = skel.at[idx, "null_value"]

        if rung not in rungs_with_null:
            skel.at[idx, "pi"] = "uncontrolled"
            continue

        if pd.isna(m):
            skel.at[idx, "pi"] = "ineligible"
            existing = str(skel.at[idx, "notes"] or "")
            skel.at[idx, "notes"] = (existing + "; measure missing").lstrip("; ")
            continue
        if pd.isna(n):
            skel.at[idx, "pi"] = "ineligible"
            existing = str(skel.at[idx, "notes"] or "")
            skel.at[idx, "notes"] = (existing + "; null missing").lstrip("; ")
            continue

        if rung in continuous_rungs:
            passes_null = bool(m > n)
            passes_zero = bool(m > 0)
            skel.at[idx, "passes_null"] = passes_null
            skel.at[idx, "passes_zero"] = passes_zero
            skel.at[idx, "pi"] = "1" if (passes_null and passes_zero) else "0"
        else:
            skel.at[idx, "pi"] = "uncontrolled"
    return skel


# =============================================================================
# Aggregate to band_verdict (decision-rules §2.3)
# =============================================================================

def longest_contiguous_run(values: np.ndarray, threshold: float) -> int:
    """Longest contiguous run of (values >= threshold)."""
    sup = np.asarray(values) >= threshold
    if not sup.any():
        return 0
    runs = []
    cur = 0
    for s in sup:
        if s:
            cur += 1
        else:
            runs.append(cur); cur = 0
    runs.append(cur)
    return max(runs)


def aggregate_band_verdict(per_cell: pd.DataFrame) -> pd.DataFrame:
    """Per (band, rung [, scale_value]) aggregation; emits band-level row with
    scale_value=NaN that carries the V verdict (with ridge condition collapsed)."""
    rows: list[dict] = []
    for band in BAND_ORDER:
        for rung in RUNG_ORDER:
            sub = per_cell[(per_cell["band"] == band) & (per_cell["rung"] == rung)]
            if len(sub) == 0:
                continue
            scale_axis = sub["scale_axis"].iloc[0]

            if scale_axis == "none":
                # Single (band, rung) row
                rows.append(_aggregate_one(sub, band, rung, scale_axis, scale_value=np.nan,
                                          add_band_verdict=True))
            else:
                # Per-scale rows + band-level aggregate row
                per_scale_rows = []
                for s_val, ss in sub.groupby("scale_value"):
                    per_scale_rows.append(_aggregate_one(ss, band, rung, scale_axis,
                                                       scale_value=float(s_val),
                                                       add_band_verdict=False))
                # Band-level aggregate uses scale-axis ridge condition
                frac_pos_arr = np.array([r["frac_pos"] for r in per_scale_rows])
                ridge_len = longest_contiguous_run(frac_pos_arr, COHORT_THRESHOLD)
                # Band-level Wilcoxon: take min raw p across scales as the
                # band-level p (most generous to positive call). BH-FDR
                # correction across the 6 bands (within rung) is applied
                # later by _bh_correct_within_rung.
                wilcoxon_ps = [r["wilcoxon_p"] for r in per_scale_rows
                               if pd.notna(r["wilcoxon_p"])]
                min_p = float(min(wilcoxon_ps)) if wilcoxon_ps else np.nan
                pi_uncontrolled = any(r["any_uncontrolled"] for r in per_scale_rows)
                # Provisional V (will be tightened by BH q-cutoff in _bh_correct_within_rung)
                if pi_uncontrolled:
                    V = "uncontrolled"
                elif (any(r["frac_pos"] >= COHORT_THRESHOLD for r in per_scale_rows)
                      and ridge_len >= RIDGE_LENGTH_MIN
                      and pd.notna(min_p) and min_p <= WILCOXON_Q_MAX):
                    V = "positive"
                elif all(r["frac_pos"] <= NEGATIVE_THRESHOLD for r in per_scale_rows):
                    V = "negative"
                else:
                    V = "silent"
                n_eligible_max = max(r["n_eligible"] for r in per_scale_rows)
                n_total = per_scale_rows[0]["n_total_cells"]
                if n_eligible_max / n_total < ELIGIBILITY_THRESHOLD:
                    V = "ineligible"
                rows.extend(per_scale_rows)
                rows.append(dict(
                    band=band, rung=rung, scale_axis=scale_axis, scale_value=np.nan,
                    n_pos=int(sum(r["n_pos"] for r in per_scale_rows)),
                    n_eligible=int(n_eligible_max),
                    frac_pos=float(np.nanmax([r["frac_pos"] for r in per_scale_rows])),
                    frac_pos_robust=float(np.nanmax([r["frac_pos_robust"] for r in per_scale_rows])),
                    wilcoxon_z=np.nan, wilcoxon_p=min_p, wilcoxon_q=np.nan,
                    ridge_length=ridge_len, V=V, n_total_cells=n_total,
                    any_uncontrolled=bool(pi_uncontrolled),
                    notes=f"band-level aggregate (max-frac across {len(per_scale_rows)} scales; ridge_len={ridge_len})",
                ))
    band_df = pd.DataFrame(rows)
    return _bh_correct_within_rung(band_df)


def _aggregate_one(sub: pd.DataFrame, band: str, rung: str, scale_axis: str,
                   scale_value: float, add_band_verdict: bool) -> dict:
    """Aggregate one (band, rung [, scale]) cell across patients."""
    n_total = sub["patient"].nunique()
    pi_vals = sub["pi"].tolist()
    n_eligible = sum(1 for p in pi_vals if p != "ineligible")
    n_pos = sum(1 for p in pi_vals if p == "1")
    any_uncontrolled = any(p == "uncontrolled" for p in pi_vals)
    frac_pos = (n_pos / n_eligible) if n_eligible > 0 else np.nan

    # Robust (drop dissenters)
    sub_robust = sub[~sub["patient"].isin(DISSENTER_MAP.keys())]
    pi_r = sub_robust["pi"].tolist()
    n_e_r = sum(1 for p in pi_r if p != "ineligible")
    n_p_r = sum(1 for p in pi_r if p == "1")
    frac_pos_robust = (n_p_r / n_e_r) if n_e_r > 0 else np.nan

    # Wilcoxon (only if both Q_full and Q_null available)
    wilcoxon_z_val = np.nan; wilcoxon_p = np.nan
    diffs = (sub["measure_value"] - sub["null_value"]).dropna().values
    if len(diffs) >= 3:
        wilcoxon_z_val, wilcoxon_p = wilcoxon_z(diffs)

    # Per-cell verdict V (only added when add_band_verdict=True; for per-scale rows we leave V empty)
    if add_band_verdict:
        if any_uncontrolled:
            V = "uncontrolled"
        elif n_eligible / max(n_total, 1) < ELIGIBILITY_THRESHOLD:
            V = "ineligible"
        elif (frac_pos >= COHORT_THRESHOLD
              and pd.notna(wilcoxon_p) and wilcoxon_p <= WILCOXON_Q_MAX):
            # Wilcoxon q computed at rung level later; here we mark with raw p,
            # bh_correct_within_rung() will revisit V if q > q_max.
            V = "positive"
        elif (pd.notna(frac_pos) and frac_pos <= NEGATIVE_THRESHOLD):
            V = "negative"
        else:
            V = "silent"
    else:
        V = ""

    return dict(
        band=band, rung=rung, scale_axis=scale_axis, scale_value=scale_value,
        n_pos=int(n_pos), n_eligible=int(n_eligible),
        frac_pos=frac_pos, frac_pos_robust=frac_pos_robust,
        wilcoxon_z=wilcoxon_z_val, wilcoxon_p=wilcoxon_p, wilcoxon_q=np.nan,
        ridge_length=np.nan, V=V, n_total_cells=int(n_total),
        any_uncontrolled=bool(any_uncontrolled),
        notes="",
    )


def _bh_correct_within_rung(df: pd.DataFrame) -> pd.DataFrame:
    """Apply BH-FDR within rung (m=6 bands), update wilcoxon_q + V."""
    for rung in RUNG_ORDER:
        # Only operate on band-level rows (scale_value NaN)
        if rung in {"L5_k", "L5_hrel", "L7", "E1_subspace"}:
            band_rows = df[(df["rung"] == rung) & df["scale_value"].isna()]
        else:
            band_rows = df[df["rung"] == rung]
        if len(band_rows) == 0:
            continue
        ps = band_rows["wilcoxon_p"].fillna(1.0).tolist()
        if all(np.isnan(band_rows["wilcoxon_p"].values)) or all(p >= 1.0 for p in ps):
            continue
        qs = bh_fdr(ps)
        for q, idx in zip(qs, band_rows.index):
            df.at[idx, "wilcoxon_q"] = q
            # Recompute V if we used raw p as the gate
            if df.at[idx, "V"] == "positive" and (q > WILCOXON_Q_MAX):
                df.at[idx, "V"] = "silent"
    # Cleanup helper column
    if "any_uncontrolled" in df.columns:
        df = df.drop(columns="any_uncontrolled")
    return df


# =============================================================================
# Triangulation (decision-rules §2.4)
# =============================================================================

def triangulate(band_verdict: pd.DataFrame) -> pd.DataFrame:
    """Per-band T(b)."""
    # Pull band-level verdicts (scale_value NaN row per (band, rung))
    rows: list[dict] = []
    for band in BAND_ORDER:
        verdicts: dict[str, str] = {}
        for rung in RUNG_ORDER:
            mask = (band_verdict["band"] == band) & (band_verdict["rung"] == rung) & (band_verdict["scale_value"].isna())
            if rung not in {"L5_k", "L5_hrel", "L7", "E1_subspace"}:
                mask = (band_verdict["band"] == band) & (band_verdict["rung"] == rung)
            sub = band_verdict[mask]
            if len(sub) == 0:
                verdicts[rung] = "uncontrolled"
            else:
                verdicts[rung] = sub["V"].iloc[0]

        T, pending = triangulate_band(verdicts)

        # Dissenter sensitivity
        diss = "stable"
        # Compare strict vs robust frac_pos for L1, L5_k (the controlled rungs)
        strict_pos = sum(1 for r in ["L1", "L5_k"] if verdicts.get(r) == "positive")
        # Recompute robust verdicts inline by re-evaluating frac_pos_robust against threshold
        # (approximation: use frac_pos_robust for the same band-level rows)
        robust_pos = 0
        for r in ["L1", "L5_k"]:
            if r in {"L5_k"}:
                row = band_verdict[(band_verdict["band"] == band) & (band_verdict["rung"] == r) & (band_verdict["scale_value"].isna())]
            else:
                row = band_verdict[(band_verdict["band"] == band) & (band_verdict["rung"] == r)]
            if len(row) and pd.notna(row["frac_pos_robust"].iloc[0]) \
                and row["frac_pos_robust"].iloc[0] >= COHORT_THRESHOLD:
                robust_pos += 1
        if robust_pos > strict_pos:
            diss = "dissenter-suppressed"
        elif robust_pos < strict_pos:
            diss = "dissenter-driven"

        rows.append(dict(
            band=band,
            **{f"V_{r}": verdicts[r] for r in RUNG_ORDER},
            T=T,
            T_pending_controls=",".join(pending) if pending else "",
            dissenter_sensitivity=diss,
            narrative_pointer=".agents/reports/2026-04-29_task-trace-rebuild-verdict.md",
        ))
    return pd.DataFrame(rows)


def triangulate_band(verdicts: dict[str, str]) -> tuple[str, list[str]]:
    """Per decision-rules §2.4 + the §"Triangulation under partial control
    coverage" open question (returns *-pending suffix for L1 + L5_k positive
    when L3 / L7 controls are still missing in v1)."""
    pending = [r for r in RUNG_ORDER if verdicts.get(r) == "uncontrolled"]
    g = verdicts.get  # alias

    # 0. Eigenvector-confirmed (pivot plan §6) — V_L5_k AND V_E1 positive.
    if g("L5_k") == "positive" and g("E1_subspace") == "positive":
        return ("eigenvector-confirmed", pending)

    # 1. Strict headline-triangulated
    if (g("L1") == "positive" and g("L5_k") == "positive"
            and (g("L3") == "positive" or g("L7") == "positive")):
        return ("headline-triangulated", pending)

    # 2. Headline-triangulated with L3/L7 pending control (v1 default)
    if (g("L1") == "positive" and g("L5_k") == "positive"
            and g("L3") in {"uncontrolled", "ineligible"}
            and g("L7") in {"uncontrolled", "ineligible"}):
        return ("headline-triangulated-pending", pending)

    # 3. Partition-resolution-locked
    if (g("L5_k") == "positive" and g("L5_hrel") != "positive"
            and g("L1") != "positive"):
        return ("partition-resolution-locked", pending)

    # 4. Continuous-only / non-modular
    if (g("L1") == "positive" and g("L5_k") != "positive"
            and g("L4") != "positive" and g("L7") != "positive"):
        return ("continuous-only", pending)

    # 5. Anatomy-suspect
    fine = {"L5_k", "L5_hrel", "L7"}
    if any(g(r) == "positive" for r in fine) and any(g(r) == "ineligible" for r in fine):
        return ("anatomy-suspect", pending)

    # 6. Ergodic — interpreted as "all CONTROLLED key rungs are silent/negative
    # AND no rung anywhere is positive". Uncontrolled rungs in v1 don't block
    # this verdict (they carry no information until §7 controls land).
    key = ["L1", "L5_k", "L5_hrel", "L4", "L7"]
    controlled_in_key = [r for r in key if g(r) not in {"uncontrolled", "ineligible"}]
    if (controlled_in_key
            and all(g(r) in {"silent", "negative"} for r in controlled_in_key)
            and not any(g(r) == "positive" for r in RUNG_ORDER)):
        return ("ergodic", pending)

    # 7. Uncontrolled-exploratory: positive somewhere with pending controls
    if any(g(r) == "positive" for r in RUNG_ORDER) and pending:
        return ("uncontrolled-exploratory", pending)

    # 8. Inconsistent: none of the above (rung disagreement, narrative needed)
    return ("inconsistent", pending)


# =============================================================================
# Figure
# =============================================================================

def render_figure(band_verdict: pd.DataFrame, triang: pd.DataFrame, out_path: Path) -> None:
    """Band x rung heatmap with triangulation strip on the right."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, (ax, ax_t) = plt.subplots(1, 2, figsize=(12, 5.5),
                                    gridspec_kw={"width_ratios": [9, 1]},
                                    constrained_layout=True)

    # Build display matrix: max frac_pos per (band, rung) across scale (band-level row)
    nb, nr = len(BAND_ORDER), len(RUNG_ORDER)
    fill = np.full((nb, nr), np.nan)
    Vmat = np.empty((nb, nr), dtype=object); Vmat[:] = ""
    label = np.empty((nb, nr), dtype=object); label[:] = ""
    is_scale = np.zeros((nb, nr), dtype=bool)
    for ib, band in enumerate(BAND_ORDER):
        for ir, rung in enumerate(RUNG_ORDER):
            if rung in {"L5_k", "L5_hrel", "L7", "E1_subspace"}:
                row = band_verdict[(band_verdict["band"] == band)
                                   & (band_verdict["rung"] == rung)
                                   & (band_verdict["scale_value"].isna())]
                is_scale[ib, ir] = True
            else:
                row = band_verdict[(band_verdict["band"] == band)
                                   & (band_verdict["rung"] == rung)]
            if len(row) == 0:
                Vmat[ib, ir] = "uncontrolled"
                continue
            fp = row["frac_pos"].iloc[0]
            fill[ib, ir] = fp if pd.notna(fp) else 0.0
            Vmat[ib, ir] = row["V"].iloc[0]
            n_pos = int(row["n_pos"].iloc[0])
            n_elig = int(row["n_eligible"].iloc[0])
            if Vmat[ib, ir] == "uncontrolled":
                label[ib, ir] = ""
            elif is_scale[ib, ir]:
                ridge = row["ridge_length"].iloc[0]
                ridge_str = f", L={int(ridge)}" if pd.notna(ridge) else ""
                label[ib, ir] = f"max {fp:.2f}{ridge_str}"
            else:
                label[ib, ir] = f"{n_pos}/{n_elig}"

    im = ax.imshow(fill, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    im.set_rasterized(True)
    ax.set_xticks(range(nr)); ax.set_xticklabels(RUNG_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(nb))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BAND_ORDER])

    for ib in range(nb):
        for ir in range(nr):
            V = Vmat[ib, ir]
            text = label[ib, ir]
            if text:
                ax.text(ir, ib, text,
                        ha="center", va="center", fontsize=7,
                        color="white" if (pd.notna(fill[ib, ir]) and fill[ib, ir] < 0.5) else "black")
            # Border per V verdict
            if V == "positive":
                ax.add_patch(Rectangle((ir-0.5, ib-0.5), 1, 1, fill=False,
                                       edgecolor="black", linewidth=2.0))
            elif V == "negative":
                ax.add_patch(Rectangle((ir-0.5, ib-0.5), 1, 1, fill=False,
                                       edgecolor="red", linewidth=2.0,
                                       linestyle="dashed"))
            elif V == "silent":
                ax.add_patch(Rectangle((ir-0.5, ib-0.5), 1, 1, fill=False,
                                       edgecolor="grey", linewidth=0.5))
            # Hatch for uncontrolled
            if V == "uncontrolled":
                ax.add_patch(Rectangle((ir-0.5, ib-0.5), 1, 1,
                                       fill=False, hatch="////",
                                       edgecolor="dimgrey", linewidth=0.0))

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("frac_pos = n_pos / n_eligible", fontsize=8)
    ax.set_xlabel("rung", fontsize=9)
    ax.set_ylabel("band", fontsize=9)

    # Triangulation strip
    Tmat = np.array([T_COLOURS.get(t, "#cccccc") for t in triang["T"]])
    ax_t.set_xlim(0, 1); ax_t.set_ylim(-0.5, nb - 0.5); ax_t.invert_yaxis()
    for ib in range(nb):
        T = triang["T"].iloc[ib]
        ax_t.add_patch(Rectangle((0, ib - 0.5), 1, 1,
                                  facecolor=T_COLOURS.get(T, "#cccccc"),
                                  edgecolor="white", linewidth=1.0))
        ax_t.text(0.5, ib, T_SHORT_CODES.get(T, "?"), ha="center", va="center",
                  fontsize=10, fontweight="bold",
                  color="white" if T not in {"partition-resolution-locked"} else "black")
        pend = triang["T_pending_controls"].iloc[ib]
        if pend:
            n_pending = len(pend.split(","))
            ax_t.text(0.5, ib + 0.32, f"+{n_pending} uncontrolled",
                      ha="center", va="center", fontsize=5,
                      color="white" if T not in {"partition-resolution-locked"} else "black")
    ax_t.set_yticks([]); ax_t.set_xticks([])
    for spine in ax_t.spines.values():
        spine.set_visible(False)
    ax_t.set_xlabel("T(b)", fontsize=9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Driver
# =============================================================================

def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--cohort", default="n10", choices=["n10"])
    p.add_argument("--fc-method", default="imcoh_abs", choices=["imcoh_abs"])
    p.add_argument("--hrel-grid", default="linspace", choices=["linspace", "logspace"])
    p.add_argument("--no-figure", action="store_true",
                   help="skip rendering the headline PDF")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    log = print if args.verbose else (lambda *a, **k: None)
    t0 = time.time()

    log(f"[audit_24] root = {ROOT}")
    log(f"[audit_24] cohort = {args.cohort}, fc_method = {args.fc_method}")

    # Preflight
    audit_df = pd.read_csv(ROOT / MEASURE_AUDIT_PATH)
    audit_meta = json.loads((ROOT / AUDIT_RUN_METADATA_PATH).read_text())
    audit_ts = float(audit_meta["audit_time_unix"])
    log(f"[audit_24] preflight: audit_ts = {audit_meta['audit_time_iso']}, "
        f"verdicts = {audit_meta['verdict_counts']}")
    preflight(audit_df, audit_ts)

    # Load cohort + dmin/dmax
    cohort = pd.read_csv(ROOT / COHORT_METADATA_PATH)
    patients = sorted(cohort["patient_id"].tolist())
    dmin_df = pd.read_csv(ROOT / DMIN_DMAX_PATH)
    log(f"[audit_24] cohort: {len(patients)} patients")

    # Build skeleton
    hrel_values = hrel_grid(ROOT / CONSUMERS["L5_hrel"]["measure"])
    log(f"[audit_24] h_rel grid ({args.hrel_grid}): {len(hrel_values)} values")
    skel = build_skeleton(patients, dmin_df, hrel_values)
    log(f"[audit_24] skeleton: {len(skel)} cells")

    # Populate
    log("[audit_24] populating L1...");      skel = populate_l1(skel)
    log("[audit_24] populating L5_k...");    skel = populate_l5_k(skel, dmin_df)
    log("[audit_24] populating L5_hrel..."); skel = populate_l5_hrel(skel, dmin_df)
    log("[audit_24] populating L4...");      skel = populate_l4(skel)
    log("[audit_24] populating L4_aux..."); skel = populate_l4_aux(skel)
    log("[audit_24] populating L6...");      skel = populate_l6(skel)
    log("[audit_24] populating L7 (streamed)..."); skel = populate_l7(skel, dmin_df)
    log("[audit_24] populating E1_subspace...");   skel = populate_e1_subspace(skel)

    # audit verdict bookkeeping
    audit_lookup = audit_df.set_index("measure_id")["audit_verdict"].to_dict()
    for rung, c in CONSUMERS.items():
        m_id = c["measure_audit_id"]
        skel.loc[skel["rung"] == rung, "audit_verdict"] = audit_lookup.get(m_id, "")

    # Per-cell positivity
    log("[audit_24] assigning pi...")
    skel = assign_pi(skel)

    # Emit per-cell
    out_per_cell = ROOT / PER_CELL_CSV
    out_per_cell.parent.mkdir(parents=True, exist_ok=True)
    skel.to_csv(out_per_cell, index=False)
    log(f"[audit_24] wrote {out_per_cell} ({len(skel)} rows)")

    # Aggregate to band_verdict
    log("[audit_24] aggregating band_verdict...")
    band_verdict = aggregate_band_verdict(skel)
    band_verdict.to_csv(ROOT / BAND_VERDICT_CSV, index=False)
    log(f"[audit_24] wrote {ROOT / BAND_VERDICT_CSV} ({len(band_verdict)} rows)")

    # Triangulate
    log("[audit_24] triangulating...")
    triang = triangulate(band_verdict)
    triang.to_csv(ROOT / TRIANGULATION_CSV, index=False)
    log(f"[audit_24] wrote {ROOT / TRIANGULATION_CSV} ({len(triang)} rows)")
    if args.verbose:
        for _, row in triang.iterrows():
            print(f"    {row['band']:>11s}  T = {row['T']}"
                  + (f"   pending: {row['T_pending_controls']}"
                     if row['T_pending_controls'] else ""))

    # Run metadata
    git_hash = ""
    try:
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        pass
    meta = dict(
        audit_24_time_unix=time.time(),
        audit_24_time_iso=time.strftime("%Y-%m-%dT%H:%M:%S"),
        git_hash=git_hash,
        phase0_audit_time_iso=audit_meta["audit_time_iso"],
        cohort=args.cohort, fc_method=args.fc_method,
        n_patients=len(patients), n_bands=len(BAND_ORDER), n_rungs=len(RUNG_ORDER),
        n_per_cell_rows=len(skel),
        n_band_verdict_rows=len(band_verdict),
        n_triangulation_rows=len(triang),
        hrel_grid=args.hrel_grid, n_hrel_values=len(hrel_values),
        dissenter_map=DISSENTER_MAP,
        decision_rules_constants=dict(
            cohort_threshold=COHORT_THRESHOLD,
            negative_threshold=NEGATIVE_THRESHOLD,
            wilcoxon_q_max=WILCOXON_Q_MAX,
            ridge_length_min=RIDGE_LENGTH_MIN,
            eligibility_threshold=ELIGIBILITY_THRESHOLD,
        ),
        elapsed_sec=time.time() - t0,
    )
    (ROOT / RUN_METADATA).write_text(json.dumps(meta, indent=2))
    log(f"[audit_24] wrote {ROOT / RUN_METADATA}")

    # Figure
    if not args.no_figure:
        log("[audit_24] rendering figure...")
        render_figure(band_verdict, triang, ROOT / FIGURE_PDF)
        log(f"[audit_24] wrote {ROOT / FIGURE_PDF}")

    log(f"[audit_24] done in {time.time() - t0:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
