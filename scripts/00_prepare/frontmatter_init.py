#!/usr/bin/env python3
"""Prepend YAML frontmatter to every `.md` under `.agents/`.

Idempotent: files that already start with `---\\n` are skipped. Era and
status are inferred from path + filename using the RULES table below.
When you add new `.agents/` files, either write frontmatter by hand
(see `.agents/guides/04_rules/frontmatter-schema.md`) or extend RULES
and re-run this.
"""
from __future__ import annotations

import re
from datetime import date
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parents[2]
AGENTS = ROOT / ".agents"
TODAY = "2026-04-24"


class Meta(NamedTuple):
    era: str
    status: str
    type: str


# Order matters — first matching prefix wins.
RULES: list[tuple[str, Meta]] = [
    # Root
    (".agents/README.md", Meta("CROSS_ERA", "current", "guide")),
    (".agents/START_HERE.md", Meta("COHORT_N9", "current", "guide")),
    # Indices
    (".agents/guides/INDEX.md", Meta("CROSS_ERA", "current", "guide")),
    (".agents/plans/INDEX.md", Meta("CROSS_ERA", "current", "guide")),
    # Guides
    (".agents/guides/01_project/", Meta("CROSS_ERA", "current", "guide")),
    (".agents/guides/02_methods/MSC_METHOD_GUIDE.md", Meta("MSC", "superseded", "guide")),
    (".agents/guides/02_methods/", Meta("IMCOH_ABS", "current", "guide")),
    (".agents/guides/03_implementation/", Meta("CROSS_ERA", "current", "guide")),
    # Plans — developed
    (".agents/plans/developed/", Meta("MSC", "dead", "plan")),
    # Plans — archive (legacy pre-2026 = MSC; 2026-01 = IMCOH_SQ)
    (".agents/plans/archive/2025-", Meta("MSC", "superseded", "plan")),
    (".agents/plans/archive/2026-01", Meta("IMCOH_SQ", "superseded", "plan")),
    # Plans — active (pre-reset Jan 2026 era, selected quirks first)
    (".agents/plans/active/2026-01-10_stage-00", Meta("IMCOH_SQ", "dead", "plan")),
    (".agents/plans/active/2026-01-10_stage-02_fc-msc", Meta("MSC", "superseded", "plan")),
    (".agents/plans/active/2026-01-10_stage-05_figures", Meta("MSC", "superseded", "plan")),
    (".agents/plans/active/2026-01-10_stage-05u", Meta("IMCOH_SQ", "dead", "plan")),
    (".agents/plans/active/2026-01-10_fc_params", Meta("MSC", "superseded", "plan")),
    (".agents/plans/active/2026-01-10_notebook_audit", Meta("IMCOH_SQ", "dead", "plan")),
    (".agents/plans/active/2026-01-10_data_inventory", Meta("IMCOH_SQ", "superseded", "plan")),
    (".agents/plans/active/2026-01-24_spatial", Meta("IMCOH_SQ", "dead", "plan")),
    (".agents/plans/active/", Meta("IMCOH_SQ", "superseded", "plan")),  # catch-all
    # Reports — specific first
    (".agents/reports/IMCOH_RESET_REPORT", Meta("IMCOH_ABS", "current", "report")),
    (".agents/reports/IMCOH_RECOVERY_PLAN", Meta("IMCOH_ABS", "current", "plan")),
    (".agents/reports/IMCOH_PROCESS_REPORT", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/IMCOH_VERIFICATION_RESULTS", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/IMCOH_PAT02_AND_CONTROLS", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/IMCOH_GAP_ANALYSIS", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/SECTION2_FIGURES_HANDOFF", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/WRITING_AGENT_BRIEFING", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/IMCOH_RESULTS_FOR_WRITING", Meta("IMCOH_ABS", "superseded", "report")),
    (".agents/reports/EPILEPTIC_IMCOH_FINAL", Meta("IMCOH_ABS", "current", "report")),
    (".agents/reports/H1_H4_VI_RESULTS_POST_RESET", Meta("COHORT_N9", "current", "report")),
    (".agents/reports/MULTISCALE_TASK_TRACE_FOR_WRITING", Meta("COHORT_N9", "current", "report")),
    (".agents/reports/PIPELINE_STATUS", Meta("CROSS_ERA", "current", "guide")),
    (".agents/reports/MODULAR_TRACE_INVESTIGATION", Meta("COHORT_N9", "dead", "report")),
    (".agents/reports/stage2_literature_menu", Meta("COHORT_N9", "dead", "report")),
    (".agents/reports/fig_section2_descriptions", Meta("IMCOH_SQ", "superseded", "report")),
    (".agents/reports/REPORT_FIGURES", Meta("MSC", "superseded", "report")),
    (".agents/reports/migration_plan_dryrun", Meta("COHORT_N9", "current", "plan")),
    (".agents/reports/SESSION_zesty-roaming-lecun", Meta("IMCOH_SQ", "dead", "report")),
    # Reports — catch-all (archive subfolder etc.)
    (".agents/reports/", Meta("IMCOH_SQ", "superseded", "report")),
]


def has_frontmatter(text: str) -> bool:
    return text.startswith("---\n")


def infer_meta(rel_path: str) -> Meta | None:
    for pattern, meta in RULES:
        if rel_path.startswith(pattern):
            return meta
    return None


def infer_name(path: Path) -> str:
    stem = path.stem
    stem = re.sub(r"^\d{4}-\d{2}-\d{2}_", "", stem)
    return stem.lower().replace("_", "-")


def infer_created(path: Path) -> str:
    m = re.match(r"^(\d{4}-\d{2}-\d{2})_", path.name)
    return m.group(1) if m else TODAY


def main() -> None:
    n_total = 0
    n_skipped = 0
    n_written = 0
    unmatched: list[str] = []

    for md in sorted(AGENTS.rglob("*.md")):
        if md.name == ".gitkeep":
            continue
        n_total += 1
        rel = str(md.relative_to(ROOT))
        text = md.read_text(encoding="utf-8")
        if has_frontmatter(text):
            n_skipped += 1
            continue
        meta = infer_meta(rel)
        if meta is None:
            unmatched.append(rel)
            continue
        name = infer_name(md)
        created = infer_created(md)
        fm = (
            "---\n"
            f"name: {name}\n"
            f"type: {meta.type}\n"
            f"era: {meta.era}\n"
            f"status: {meta.status}\n"
            f"created: {created}\n"
            f"updated: {TODAY}\n"
            "pointers: []\n"
            "---\n\n"
        )
        md.write_text(fm + text, encoding="utf-8")
        n_written += 1

    print(f"Scanned: {n_total} | had frontmatter: {n_skipped} | wrote: {n_written}")
    if unmatched:
        print(f"\nUnmatched ({len(unmatched)}):")
        for p in unmatched:
            print(f"  {p}")


if __name__ == "__main__":
    main()
