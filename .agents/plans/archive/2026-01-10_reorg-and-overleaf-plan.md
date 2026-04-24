---
name: reorg-and-overleaf-plan
type: plan
era: IMCOH_SQ
status: superseded
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# Plan (2026-01-10): Repo reorg + pipeline prep + Overleaf bundle

Note: superseded by the stage plans in `.agents/plans/active/`.

## Scope
Prepare the codebase and notebooks for a clean, reproducible pipeline by Tuesday
and an Overleaf-ready results bundle by Thursday. Development runs should use a
single-patient skeleton (default Pat_02), cache-first, and low surrogate counts
to avoid long runtimes while keeping everything patient-agnostic.

## Constraints
- Dev mode: single patient by default (Pat_02), but no hard-coded patient names.
- Cache-first everywhere.
- Use low surrogate counts during development.
- Long-running validation steps deferred to final run.
- Notebook archiving before deletion.

## Plan

1) Notebooks and index (done)
- Create numbered notebook folders under `ipynb/`.
- Move active notebooks to the new structure.
- Archive legacy notebooks under `ipynb/90_archive/.old/`.
- Create `ipynb/INDEX.md` with stable IDs and summaries.

2) Dev-mode runbook and pipeline defaults (in progress)
- Add a developer runbook under `.agents/`.
- Ensure pipelines support single-patient dev runs (added bands/phases/filter-time to FC scripts).
- Make timeouts explicit for long stages.
- Keep surrogate counts low in dev.

3) Source refactor with compatibility (pending)
- Create clear subpackages (io/fc/lrg/metrics/viz/pipelines).
- Keep legacy imports working via re-exports.
- Update docs to reference the new layout.

4) Overleaf-ready bundle (pending)
- Define final figures, tables, and metrics.
- Generate `outputs/overleaf/` with figures/tables/summaries.
- Validate full-run instructions for all patients.

## Current status
- Step 1 is complete.
- Step 2 is in progress (FC scripts updated; single-patient dev notebooks added).

## Notes
- Final full-run (all patients) is deferred until the pipeline is stable.
- Notebook deletions only after functionality is migrated.
