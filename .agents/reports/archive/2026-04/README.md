---
name: reports-archive-2026-04-readme
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - ../../2026-04-24_pipeline-status.md
  - ../../../era-map.md
  - ../../2026-04-24_post-mortem-scalar-session.md
---

# Archive — 2026-04

**Reports moved here during the April 2026 reorg. Two kinds live
together: SUPERSEDED (pre-reset |ImCoh|² outputs, invalidated
numerically) and DEAD (scalar-session artifacts, experimentally
complete but not cited). Era and status live in each file's
frontmatter; this bucket is "archived in April 2026", not "belongs to
April 2026".**

## Contents

### Superseded pre-reset reports (era: IMCOH_SQ)

- `2026-04-15_imcoh-process.md`
- `2026-04-15_imcoh-verification.md`
- `2026-04-15_imcoh-pat02-controls.md`
- `2026-04-15_imcoh-gap-analysis.md`
- `2026-04-15_section2-figures-handoff.md`
- `2026-04-15_writing-agent-briefing.md`
- `2026-04-15_fig-section2-descriptions.md`

All produced during the |ImCoh|² era before the ordering-bug reset on
2026-04-15. Their numerical claims are invalid; their methodology
descriptions survive only inasmuch as the `imcoh-guide.md` supersedes
them.

### Superseded writing handoff (era: IMCOH_ABS, status: superseded)

- `2026-04-24_imcoh-results-for-writing.md` — replaced by the
  canonical `2026-04-24_multiscale-task-trace.md` in `../..`.

### Dead scalar-session artifacts (era: COHORT_N9, status: dead)

- `2026-04-24_modular-trace-investigation.md` — handoff log for the
  failed Stage 0b → 4 gauntlet. Preserved for historical trace; its
  STOP directive stands.
- `2026-04-24_stage2-literature-menu.md` — tree-distance literature
  search output that fed Stage 3 (KC / MC / wRF). Methodology useful;
  the gate-testing framing is dead.

See `../../2026-04-24_post-mortem-scalar-session.md` for the canonical
lesson.

### Historical session tag + MSC-era figure registry

- `2026-04-13_session-reorganize-scripts.md` — session label from the
  scripts/data reorganization pass.
- `2026-02-26_report-figures.md` — MSC-era figure registry, superseded
  by `2026-04-24_multiscale-task-trace.md` + the per-section handoffs.

## Policy

- **Don't cite** files in this directory as current results.
- **Don't delete** — archives encode lineage. Use `git mv` or add a new
  bucket `archive/YYYY-MM/` for future rounds.
- **Do reference** them from post-mortems and era-map entries.
