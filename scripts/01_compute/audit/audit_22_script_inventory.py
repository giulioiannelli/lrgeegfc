#!/usr/bin/env python3
"""Audit Step 22 — script inventory + dead-duplicate decisions (Phase 0).

Walks ``scripts/01_compute/{hypothesis_tests, diagnostics, batch, audit,
figures_embedded}/`` and emits one row per script to
``data/audit/measure_audit/script_inventory_n10.csv``:

    path, mtime, lines, predicate_summary, produces_artefact,
    artefact_status, supersedes, archive_destination

Decision rule for ``archive_destination``:

  - **keep** if the script's artefact appears in the rebuild §6 reuse map
    (data/reports/imcoh_*/, data/audit/dvi_split_baseline/,
    data/audit/trace_modules/, data/audit/per_patient_hierarchy_*/, etc.)
    OR if the script is ``current`` per
    ``.agents/reports/2026-04-25_measure-ledger.csv``.
  - **archive** to ``scripts/archive/2026-04_pre-rebuild/<subdir>/`` if
    the script's artefact appears only in the ``superseded`` rows of the
    ledger OR if no artefact path is detected and no other current
    script depends on this one.

Outputs:
    data/audit/measure_audit/script_inventory_n10.csv
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()


SCRIPT_SUBDIRS = [
    "scripts/01_compute/hypothesis_tests",
    "scripts/01_compute/diagnostics",
    "scripts/01_compute/batch",
    "scripts/01_compute/audit",
    "scripts/01_compute/figures_embedded",
]

LEDGER_CSV = ".agents/reports/2026-04-25_measure-ledger.csv"

# Reuse map paths from rebuild plan §6 (anything under these dirs is current input)
REUSE_DIRS = {
    "data/reports/imcoh_vi",
    "data/reports/imcoh_continuous_trace",
    "data/audit/dvi_split_baseline",
    "data/audit/dvi_spatial_scale",         # diag_dvi_spatial_scale.py + fig_dvi_spatial_cohort.py
    "data/audit/four_phase_taumin",         # audit_17_taumin_four_phase.py
    "data/audit/four_phase_tau_star",       # audit_18_taustar_four_phase.py
    "data/audit/psi_tau_scan",              # diag_psi_tau_scan.py (closed but recent)
    "data/audit/residual_subspace",         # diag_residual_subspace.py (closed but recent)
    "data/audit/trace_modules",
    "data/audit/per_patient_hierarchy",     # audit_08 ancestor figures
    "data/audit/per_patient_hierarchy_cohesion",
    "data/audit/per_patient_hierarchy_cnp",
    "data/audit/per_patient_hierarchy_mspc",
    "data/audit/measure_audit",
    "data/audit/cohort_metadata.csv",
    "data/audit/ridge_diagnostics",
    "data/outputs/figures/section6",
}

# These artefact substrings flag MSC-era / pre-2026-04-15 work
ARCHIVE_HINTS_MSC = ["compute_msc", "rescaled_msc", "bipolar"]
ARCHIVE_HINTS_DEAD = ["task_trace_per_patient", "task_trace_cross_patient",
                      "task_trace_old5_only", "containment_cbr",
                      "consensus_subtree", "task_anchored_cbr"]
# MRL-era scripts from the superseded scalar-trace session (per
# 2026-04-25_mrl-vs-cbr-reconciliation.md)
ARCHIVE_HINTS_MRL = ["trace_backward", "trace_null", "trace_sweep",
                     "fig_trace_anchor", "fig_trace_asymmetry",
                     "fig_trace_forward_vs_backward", "fig_trace_multipatient",
                     "fig_trace_per_band"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DOCSTRING_RE = re.compile(r'"""(.*?)"""', re.DOTALL)
ARTEFACT_PATTERNS = [
    re.compile(r'\.to_csv\s*\(\s*([^)]+)\)'),
    re.compile(r'\.to_parquet\s*\(\s*([^)]+)\)'),
    re.compile(r'np\.savez(?:_compressed)?\s*\(\s*([^,)]+)'),
    re.compile(r'np\.save\s*\(\s*([^,)]+)'),
    re.compile(r'\.savefig\s*\(\s*([^,)]+)'),
    re.compile(r'plt\.savefig\s*\(\s*([^,)]+)'),
]
PATH_LITERAL_RE = re.compile(r'(data/[a-zA-Z0-9_\-/]+\.(csv|parquet|npz|npy|pdf|png|md|json))')


def _docstring(body: str) -> str:
    """Return first triple-quoted module docstring, single-line summary."""
    m = DOCSTRING_RE.search(body)
    if not m:
        return ""
    text = m.group(1).strip()
    # First non-blank line
    for line in text.split("\n"):
        line = line.strip()
        if line:
            return line[:200]
    return ""


def _detect_artefacts(body: str) -> list[str]:
    """Heuristic: extract output-artefact paths from .to_csv / np.savez / savefig calls.

    Falls back to scanning for `data/...` literals near the end of the file
    (within the last 1/3) as a coarse approximation when the call pattern
    uses Path() composition.
    """
    found = set()
    for pat in ARTEFACT_PATTERNS:
        for m in pat.finditer(body):
            arg = m.group(1).strip().strip('"\'')
            # Try to extract a literal path inside the arg
            sub = PATH_LITERAL_RE.search(arg)
            if sub:
                found.add(sub.group(1))
            elif arg.startswith("data/"):
                found.add(arg)
    # Fallback: scan the whole body for any data/ literal
    for m in PATH_LITERAL_RE.finditer(body):
        found.add(m.group(1))
    return sorted(found)


def _classify_artefact(art: str, ledger_current: set[str], ledger_superseded: set[str]) -> str:
    """Return 'reuse-map', 'ledger-current', 'ledger-superseded', or 'unknown'."""
    if art in ledger_current:
        return "ledger-current"
    if art in ledger_superseded:
        return "ledger-superseded"
    for d in REUSE_DIRS:
        if art.startswith(d):
            return "reuse-map"
    return "unknown"


def _decide_archive(
    script: Path, artefacts: list[str], statuses: list[str],
) -> tuple[str, str]:
    """Return (archive_destination, reason)."""
    name = script.name
    rel = str(script.relative_to(ROOT))

    # MSC-era hints
    if any(h in name.lower() for h in ARCHIVE_HINTS_MSC):
        return ("scripts/archive/2026-04_pre-rebuild/msc-era/", "msc-era hint in script name")

    # Known-dead hypothesis branches (CBR variants superseded by Cohesion-CBR)
    if any(h in name.lower() for h in ARCHIVE_HINTS_DEAD):
        return ("scripts/archive/2026-04_pre-rebuild/dead-branches/", "dead-branch hint in script name")

    # MRL-era scalar-trace scripts superseded by 2026-04-25 reconciliation
    if any(h in name.lower() for h in ARCHIVE_HINTS_MRL):
        return ("scripts/archive/2026-04_pre-rebuild/mrl-era/", "mrl-era hint in script name (superseded by 2026-04-25 reconciliation)")

    # If at least one artefact is in reuse-map or ledger-current → KEEP
    if any(s in {"reuse-map", "ledger-current"} for s in statuses):
        return ("keep", "produces current artefact")

    # If all artefacts are ledger-superseded → ARCHIVE
    if statuses and all(s == "ledger-superseded" for s in statuses):
        return ("scripts/archive/2026-04_pre-rebuild/superseded/", "all artefacts in ledger.superseded")

    # If no artefacts detected → likely orchestration / library; KEEP unless older than 2026-02
    if not artefacts:
        return ("keep", "no artefact detected (orchestration?)")

    # All unknowns → MANUAL
    if statuses and all(s == "unknown" for s in statuses):
        return ("manual-review", "artefacts not in ledger or reuse-map; needs manual triage")

    return ("keep", "mixed; defer to manual review")


def _git_first_commit(script: Path) -> Optional[str]:
    try:
        out = subprocess.check_output(
            ["git", "log", "--diff-filter=A", "--follow", "--format=%ai|%s", "--", str(script.relative_to(ROOT))],
            cwd=str(ROOT), stderr=subprocess.DEVNULL,
        ).decode()
        lines = out.strip().splitlines()
        return lines[-1] if lines else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/audit/measure_audit")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "script_inventory_n10.csv"

    # Load ledger artefact sets (engine=python tolerates the one badly-quoted
    # row at line 30 of the ledger; flag for fix in step 0e).
    ledger = pd.read_csv(ROOT / LEDGER_CSV, engine="python", on_bad_lines="warn")
    ledger_current = set(
        ledger.loc[ledger.status == "current", "output_artifact"].dropna().astype(str).tolist()
    )
    ledger_superseded = set(
        ledger.loc[ledger.status == "superseded", "output_artifact"].dropna().astype(str).tolist()
    )
    if args.verbose:
        print(f"[audit_22] ledger: {len(ledger_current)} current, {len(ledger_superseded)} superseded artefacts")

    rows = []
    for subdir_rel in SCRIPT_SUBDIRS:
        subdir = ROOT / subdir_rel
        if not subdir.exists():
            continue
        for script in sorted(subdir.glob("*.py")):
            body = script.read_text(errors="ignore")
            lines = body.count("\n") + 1
            mtime = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(script.stat().st_mtime))
            doc = _docstring(body)
            artefacts = _detect_artefacts(body)
            statuses = [_classify_artefact(a, ledger_current, ledger_superseded) for a in artefacts]
            archive, reason = _decide_archive(script, artefacts, statuses)
            git_intro = _git_first_commit(script) or ""

            rows.append({
                "path": str(script.relative_to(ROOT)),
                "subdir": subdir_rel.split("/")[-1],
                "mtime": mtime,
                "lines": lines,
                "predicate_summary": doc,
                "produces_artefacts": "; ".join(artefacts),
                "artefact_status_summary": "; ".join(statuses),
                "archive_destination": archive,
                "decision_reason": reason,
                "git_intro": git_intro,
            })
            if args.verbose:
                glyph = {"keep": "K", "manual-review": "?"}.get(archive, "A")
                print(f"  [{glyph}] {script.relative_to(ROOT)} -> {archive}")

    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    if args.verbose:
        print()
    counts = df["archive_destination"].value_counts().to_dict()
    print(f"[audit_22] inventoried {len(df)} scripts; archive decisions: {counts}")
    print(f"[audit_22] wrote: {csv_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
