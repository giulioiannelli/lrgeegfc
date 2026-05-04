#!/usr/bin/env python3
"""Audit Step 21 — per-measure correctness audit (Phase 0 of task-trace rebuild).

For every measure in the rebuild plan §6 reuse map, run six φ-checks and
emit one row to ``data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv``:

    phi_coh   cohort actual vs cohort_metadata.csv (auto)
    phi_fc    producer reads from imcoh_lrg/ or imcoh_lrg_halves/ (auto)
    phi_help  producer does NOT redefine library helpers locally (auto, grep)
    phi_name  artefact filename matches CLAUDE.md canonical patterns (auto)
    phi_fid   producer predicate matches scope-report predicate (MANUAL flag)
    phi_smoke producer reproduces one cell from cache to tolerance (MANUAL — see audit_23)

Verdict ∈ {current, stale-numeric, broken, duplicate, archive} per
``.agents/guides/task-persistence-investigation/2026-04-29_measure-correctness-audit.md``.

Outputs:
    data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv
    data/audit/measure_audit/audit_run_metadata.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()


# ---------------------------------------------------------------------------
# Reuse map — every (measure, producer, scope, artefact) tracked by the audit
# ---------------------------------------------------------------------------

@dataclass
class Measure:
    measure_id: str
    producer: str           # path relative to ROOT, may be "" if scope-only
    scope_report: str       # path under .agents/, may be "" if no scope yet
    artefact: str           # path relative to ROOT
    rung: str               # L1, L2, ..., L7 per rebuild plan §2
    artefact_type: str      # csv | npz | parquet
    expected_fc_method: str = "imcoh_abs"
    cohort_claim: str = "n10"   # n9 | n10
    canonical_name_regex: str = ""
    cohort_sparse: bool = False     # if True, allow < n10 patients (sparse-by-design measures)
    skip_phi_fc: bool = False       # metadata producers don't need imcoh cache
    notes: str = ""


REUSE_MAP: list[Measure] = [
    Measure(
        measure_id="L1_h2c_shared",
        producer="scripts/01_compute/hypothesis_tests/h2c_ultrametric_drift.py",
        scope_report="",
        artefact="data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv",
        rung="L1",
        artefact_type="csv",
    ),
    Measure(
        measure_id="L1_continuous_trace_per_cell",
        producer="scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md",
        artefact="data/reports/imcoh_continuous_trace/per_cell_summary.csv",
        rung="L1",
        artefact_type="csv",
    ),
    Measure(
        measure_id="L1_continuous_trace_split_baseline",
        producer="scripts/01_compute/hypothesis_tests/continuous_trace_matrix.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md",
        artefact="data/reports/imcoh_continuous_trace/per_cell_summary_split.csv",
        rung="L1",
        artefact_type="csv",
    ),
    Measure(
        measure_id="L1_continuous_trace_controls",
        producer="scripts/01_compute/hypothesis_tests/continuous_trace_controls.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md",
        artefact="data/reports/imcoh_continuous_trace/controls_summary.csv",
        rung="L1",
        artefact_type="csv",
        skip_phi_fc=True,  # consumes upstream CSVs (per_cell_summary*, h2e_split_half), not LRG cache
    ),
    Measure(
        measure_id="L1_h2e_split_half",
        producer="scripts/01_compute/hypothesis_tests/h2e_split_half.py",
        scope_report="",
        artefact="data/reports/imcoh_vi/h2e_split_half_rho_raw.csv",
        rung="L1_null",
        artefact_type="csv",
        cohort_claim="n10",  # ledger says "n=9 stale" but file actually has n=10
    ),
    Measure(
        measure_id="L5k_partition_multiscale",
        producer="scripts/01_compute/hypothesis_tests/h2_partition_multiscale.py",
        scope_report="",
        artefact="data/reports/imcoh_vi/h2_partition_multiscale_raw.csv",
        rung="L5(k)",
        artefact_type="csv",
    ),
    Measure(
        measure_id="L5k_drift_floor",
        producer="scripts/01_compute/diagnostics/diag_dvi_split_baseline.py",
        scope_report="",
        artefact="data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv",
        rung="L5(k)_null",
        artefact_type="csv",
        canonical_name_regex=r"^dvi_split_baseline_n10_imcoh_abs\.csv$",
    ),
    Measure(
        measure_id="L5_hrel_partition",
        producer="scripts/01_compute/diagnostics/diag_dvi_hrel.py",
        scope_report="",
        artefact="data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv",
        rung="L5(h_rel)",
        artefact_type="csv",
        canonical_name_regex=r"^dvi_hrel_n10_imcoh_abs\.csv$",
    ),
    Measure(
        measure_id="L5_kcut_heights",
        producer="scripts/01_compute/diagnostics/diag_kcut_heights.py",
        scope_report="",
        artefact="data/audit/dvi_split_baseline/kcut_heights_n10_imcoh_abs.csv",
        rung="L5_aux",
        artefact_type="csv",
        canonical_name_regex=r"^kcut_heights_n10_imcoh_abs\.csv$",
    ),
    Measure(
        measure_id="L5_dmin_dmax_inventory",
        producer="scripts/01_compute/diagnostics/diag_dvi_hrel.py",
        scope_report="",
        artefact="data/audit/dvi_split_baseline/dmin_dmax_inventory.csv",
        rung="L5_aux",
        artefact_type="csv",
    ),
    Measure(
        measure_id="L4_trace_modules_J090",
        producer="scripts/01_compute/audit/audit_15_trace_modules.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md",
        artefact="data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv",
        rung="L4",
        artefact_type="csv",
        cohort_sparse=True,  # T-regime J=0.9 returns ≈0 modules cohort-wide by design
    ),
    Measure(
        measure_id="L4_aux_cohesion_cbr",
        producer="scripts/01_compute/audit/audit_12_cohesion_cbr.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md",
        artefact="data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv",
        rung="L4_aux",
        artefact_type="csv",
        cohort_claim="n10",  # CSV actually has n=10 (ledger note "n=9 active" is stale)
        skip_phi_fc=True,    # consumes precomputed LRG linkages, not raw imcoh cache
    ),
    Measure(
        measure_id="L6_cnp",
        producer="scripts/01_compute/audit/audit_13_cnp.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-25_cophenetic-neighbourhood.md",
        artefact="data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv",
        rung="L6",
        artefact_type="csv",
        skip_phi_fc=True,  # consumes precomputed LRG cophenetic distances
    ),
    Measure(
        measure_id="L7_mspc",
        producer="scripts/01_compute/audit/audit_14_mspc.py",
        scope_report=".agents/guides/task-persistence-investigation/2026-04-26_multiscale-partition-coherence.md",
        artefact="data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv",
        rung="L7",
        artefact_type="csv",
        skip_phi_fc=True,  # consumes precomputed LRG linkages
    ),
    Measure(
        measure_id="cohort_metadata",
        producer="scripts/01_compute/audit/audit_01_cohort_metadata.py",
        scope_report="",
        artefact="data/audit/cohort_metadata.csv",
        rung="metadata",
        artefact_type="csv",
        canonical_name_regex=r"^cohort_metadata\.csv$",
        skip_phi_fc=True,  # metadata script reads raw .fif files, not imcoh cache
    ),
]


# ---------------------------------------------------------------------------
# Library-helper symbols that must NOT be locally redefined in producers
# ---------------------------------------------------------------------------

LIBRARY_HELPERS = [
    # vi.py
    "compute_vi", "conditional_entropy",
    # hypothesis.py
    "wilcoxon_z", "rank_biserial", "boot_ci_mean", "bh_fdr", "cluster_stats",
    # tree.py
    "tree_internal_nodes", "dmax_from_Z", "jaccard_leafsets", "h_log_grid",
    "fcluster_at_h_rel", "simpson_neff", "cluster_size_stats",
    "partition_vi_on_subset",
    # tree_distance.py
    "kc_distance", "matching_cluster_distance", "weighted_rf_distance",
    # workflow.fc / lrg
    "load_fc_matrix", "load_lrg_result",
]

# Canonical filename pattern fallbacks
CANONICAL_NAME_PATTERNS = {
    "imcoh_lrg": re.compile(r"^[a-z_]+_(rest_pre|rest_post|task_learn|task_test)_lrg_imcoh-(abs|sq)\.npz$"),
    "imcoh_lrg_halves": re.compile(r"^[a-z_]+_(rest_pre|rest_post)_(A|B)_lrg_imcoh-abs\.npz$"),
    "imcoh_freqresolved": re.compile(r"^[a-z_]+_(rest_pre|rest_post|task_learn|task_test)_imcoh_freqresolved_nperseg-(\d+)\.npy$"),
}


# ---------------------------------------------------------------------------
# φ-checks
# ---------------------------------------------------------------------------

def _read_artefact_patients(artefact: Path) -> list[str]:
    """Return the unique patient list in an artefact (CSV with `patient` col)."""
    if not artefact.exists():
        return []
    if artefact.suffix == ".csv":
        try:
            df = pd.read_csv(artefact, usecols=["patient"])
            return sorted(df["patient"].dropna().unique().tolist())
        except Exception:
            try:
                df = pd.read_csv(artefact)
                for col in ("patient", "patient_id", "Patient", "subject"):
                    if col in df.columns:
                        return sorted(df[col].dropna().astype(str).unique().tolist())
            except Exception:
                pass
            return []
    return []


def check_phi_coh(measure: Measure, cohort_n10: list[str]) -> tuple[int, str]:
    """Cohort actual vs claimed."""
    artefact = ROOT / measure.artefact
    actual = _read_artefact_patients(artefact)
    if not actual:
        return 0, "could not extract patient column from artefact"
    actual_set = set(actual)
    n10_set = set(cohort_n10)
    extra = sorted(actual_set - n10_set)
    if extra:
        return 0, f"unexpected patients in artefact: {extra}"
    if measure.cohort_sparse:
        # Sparse-by-design measures (e.g. T-regime J=0.9 trace_modules)
        # are allowed to omit patients with zero rows.
        return 1, f"sparse-by-design; {len(actual)}/{len(cohort_n10)} patients with rows"
    if measure.cohort_claim == "n10":
        if actual_set == n10_set:
            return 1, f"n=10 ok ({len(actual)} patients)"
        missing = sorted(n10_set - actual_set)
        return 0, f"missing={missing}"
    elif measure.cohort_claim == "n9":
        n9_set = n10_set - {"Pat_14"}
        if actual_set == n9_set:
            return 1, "n=9 ok (Pat_14 honestly excluded)"
        return 0, f"n=9 claim but actual={sorted(actual_set)}"
    else:
        return 0, f"unknown cohort_claim={measure.cohort_claim}"


def check_phi_fc(measure: Measure) -> tuple[int, str]:
    """Producer reads from imcoh_lrg{,_halves}/ and respects Jensen order.

    Skipped (returns 1) for measures whose producer is a metadata script or
    a downstream consumer of already-processed CSVs — those producers
    transitively respect the FC method via their upstream dependencies.
    """
    if measure.skip_phi_fc:
        return 1, "skipped (metadata or upstream-CSV consumer)"
    if not measure.producer:
        return 1, "no producer (scope-only)"
    src = ROOT / measure.producer
    if not src.exists():
        return 0, "producer file missing"
    body = src.read_text()
    reads_imcoh = (
        "IMCOH_LRG_CACHE" in body
        or "imcoh_lrg" in body
        or "load_fc_matrix" in body
        or "load_lrg_result" in body
        or "imcoh_abs" in body
    )
    if not reads_imcoh:
        return 0, "producer does not reference imcoh cache"
    # Heuristic Jensen-violation flag: literal `np.abs(np.mean(...))` is wrong
    # order, but `np.mean(np.abs(...))` is correct. Only the explicit reverse
    # order is a hard fail; everything else is OK because the canonical
    # transform happens inside `workflow.fc.load_fc_matrix`.
    if re.search(r"np\.abs\s*\(\s*np\.mean\(", body):
        return 0, "Jensen-violating np.abs(np.mean(...)) detected"
    return 1, "imcoh cache referenced; no Jensen-violation pattern"


def check_phi_help(measure: Measure) -> tuple[int, list[str], str]:
    """Producer does NOT redefine library helpers locally."""
    if not measure.producer:
        return 1, [], "no producer (scope-only)"
    src = ROOT / measure.producer
    if not src.exists():
        return 0, [], "producer file missing"
    body = src.read_text()
    duplicated = []
    for sym in LIBRARY_HELPERS:
        # detect `def sym(...)` or `def _sym(...)` at module level
        if re.search(rf"^def\s+_?{re.escape(sym)}\s*\(", body, re.MULTILINE):
            duplicated.append(sym)
    if duplicated:
        return 0, duplicated, f"local re-implementation: {duplicated}"
    return 1, [], "no library-helper duplicates"


def check_phi_name(measure: Measure) -> tuple[int, str]:
    """Artefact filename matches canonical patterns."""
    if not measure.artefact:
        return 0, "no artefact path"
    name = Path(measure.artefact).name
    if measure.canonical_name_regex:
        if re.match(measure.canonical_name_regex, name):
            return 1, "matches canonical regex"
        return 0, f"does not match {measure.canonical_name_regex!r}"
    # Loose default: report-style CSVs are OK as long as they end in .csv/.npz/.parquet
    if Path(measure.artefact).suffix not in {".csv", ".npz", ".parquet"}:
        return 0, f"unsupported extension {Path(measure.artefact).suffix}"
    return 1, "extension OK; no specific regex required"


def check_phi_fid_manual(measure: Measure) -> tuple[int, str]:
    """φ_fid is human-judged in v1: emit a flag for manual review.

    Returns 1 if scope-report exists (auditor can compare); 0 if missing.
    """
    if not measure.scope_report:
        # Acceptable for measures whose scope is the predicate's docstring
        # (e.g. h2c_ultrametric_drift, h2_partition_multiscale, h2e_split_half
        # — these predate the task-persistence-investigation folder rule).
        return 1, "no scope report; producer docstring is the contract"
    scope = ROOT / measure.scope_report
    if not scope.exists():
        return 0, f"scope-report missing: {measure.scope_report}"
    return 1, "scope-report exists; MANUAL φ_fid review owed (audit_23)"


def check_phi_smoke_manual(measure: Measure) -> tuple[int, str]:
    """φ_smoke is run in audit_23 (per-measure double-check). Skip in v1."""
    return 1, "deferred to audit_23 double-check protocol"


def classify_verdict(
    phi_coh: int, phi_fc: int, phi_help: int, phi_name: int,
    phi_fid: int, phi_smoke: int, cohort_claim: str,
) -> str:
    if phi_help == 0:
        return "duplicate"
    if phi_fid == 0 or phi_smoke == 0 or phi_fc == 0:
        return "broken"
    if phi_coh == 0:
        if cohort_claim == "n9":
            return "stale-numeric"
        return "broken"
    if phi_name == 0:
        return "current"  # soft failure; rename in 0e
    return "current"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", default="n10", help="cohort tag (default n10)")
    parser.add_argument("--fc-method", default="imcoh_abs",
                        help="FC method to audit (default imcoh_abs)")
    parser.add_argument("--out-dir", default="data/audit/measure_audit",
                        help="output directory for CSV + metadata")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"measure_audit_{args.cohort}_{args.fc_method}.csv"
    meta_path = out_dir / "audit_run_metadata.json"

    cohort_meta_csv = ROOT / "data/audit/cohort_metadata.csv"
    cohort = pd.read_csv(cohort_meta_csv)
    cohort_n10 = sorted(cohort["patient_id"].astype(str).tolist())
    if args.verbose:
        print(f"[audit_21] cohort_n10 ({len(cohort_n10)}): {cohort_n10}")
        print(f"[audit_21] auditing {len(REUSE_MAP)} measures")

    rows = []
    for m in REUSE_MAP:
        ph_coh, n_coh = check_phi_coh(m, cohort_n10)
        ph_fc, n_fc = check_phi_fc(m)
        ph_help, dups, n_help = check_phi_help(m)
        ph_name, n_name = check_phi_name(m)
        ph_fid, n_fid = check_phi_fid_manual(m)
        ph_smoke, n_smoke = check_phi_smoke_manual(m)
        verdict = classify_verdict(
            ph_coh, ph_fc, ph_help, ph_name, ph_fid, ph_smoke, m.cohort_claim,
        )
        notes = "; ".join(filter(None, [
            f"phi_coh={n_coh}",
            f"phi_fc={n_fc}",
            f"phi_help={n_help}",
            f"phi_name={n_name}",
            f"phi_fid={n_fid}",
            f"phi_smoke={n_smoke}",
        ]))
        rows.append({
            "measure_id": m.measure_id,
            "rung": m.rung,
            "producer_script": m.producer,
            "scope_report": m.scope_report,
            "artefact_path": m.artefact,
            "fc_method_claimed": m.expected_fc_method,
            "cohort_claimed": m.cohort_claim,
            "phi_fid": ph_fid,
            "phi_coh": ph_coh,
            "phi_fc": ph_fc,
            "phi_help": ph_help,
            "phi_smoke": ph_smoke,
            "phi_name": ph_name,
            "helpers_duplicated": ",".join(dups),
            "audit_verdict": verdict,
            "notes": notes,
        })
        if args.verbose:
            mark = "✓" if verdict == "current" else "·"
            print(f"  [{mark}] {m.measure_id:<40s} {verdict}")

    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    if args.verbose:
        print(f"[audit_21] wrote {csv_path} ({len(df)} rows)")

    # metadata sidecar
    try:
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(ROOT), stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        git_hash = "unknown"
    meta_path.write_text(json.dumps({
        "audit_time_unix": time.time(),
        "audit_time_iso": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": git_hash,
        "python_version": sys.version,
        "n_measures_audited": len(rows),
        "verdict_counts": df["audit_verdict"].value_counts().to_dict(),
        "csv_path": str(csv_path.relative_to(ROOT)),
    }, indent=2))

    # terse summary
    counts = df["audit_verdict"].value_counts().to_dict()
    print(f"[audit_21] verdict counts: {counts}")
    print(f"[audit_21] wrote: {csv_path.relative_to(ROOT)}")
    print(f"[audit_21] meta : {meta_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
