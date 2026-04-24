---
name: era-map
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - .agents/reports/2026-04-24_pipeline-status.md
  - .agents/guides/02_methods/imcoh-guide.md
  - .agents/guides/02_methods/probe-bias-guide.md
  - .agents/guides/03_implementation/data-layout.md
---

# Era map

**Four eras. Each changed what we could claim. Each is a load-bearing
frontmatter tag and a filter for what's citable.**

## MSC — 2025-11 → 2026-04-10

- **Carrier:** magnitude-squared coherence (real part of cross-spectrum).
- **Invalidated by:** same-probe volume-conduction bias dominating
  community structure at coarse scales
  (`.agents/guides/02_methods/probe-bias-guide.md`).
- **Landmark docs:**
  - `.agents/guides/02_methods/msc-method-guide.md` (methodology)
  - `.agents/reports/archive/2026-04/` (pre-reset reports)
- **What survives:** the general LRG workflow, scripts scaffolding,
  cache layout, CLI.
- **What dies:** every numerical claim, every H2a result, every figure.

## IMCOH_SQ — 2026-04-10 → 2026-04-15

- **Carrier:** `|ImCoh|²` (squared imaginary coherence, band-averaged).
- **Invalidated by:** the 2026-04-15 reset. An ordering bug in
  `compute_lrg_analysis` silently produced wrong LRG output for this
  FC variant. Full diagnosis:
  `.agents/reports/2026-04-15_imcoh-reset.md`.
- **Lifetime:** 5 days.
- **What dies:** all quantitative statements.

## IMCOH_ABS — 2026-04-15 → present

- **Carrier:** `<|ImCoh|>_f` — frequency-resolved imaginary coherence,
  magnitude-averaged over the band (Jensen-safe: `mean(|·|)`, not
  `|mean(·)|`).
- **Trigger:** reset fix + corrected order-of-operations.
- **State:** current production FC. Post-reset LRG dendrograms are the
  numerical ground truth.
- **Landmark docs:**
  - `.agents/reports/2026-04-15_imcoh-reset.md` (diagnosis)
  - `.agents/guides/02_methods/imcoh-guide.md` (methodology)
  - `.agents/reports/2026-04-24_pipeline-status.md` (era index)

## COHORT_N9 — 2026-04-22 → present

- **Carrier:** IMCOH_ABS + 9-patient locked cohort. Members: Pat_02,
  03, 05, 06, 07, 08, 10, 13, 15. Excluded: Pat_14 (corrupt
  `task_test`).
- **Trigger:** data normalization + quality pass on 2026-04-22.
- **State:** the cohort any cohort-wide claim must reference.
- **Quirks:**
  - Pat_03 at 1024 Hz — negative control, flagged distinctly.
  - Pat_10 task rows [53, 54, 55] dropped at load.
- **Landmark docs:**
  - `.agents/guides/03_implementation/data-layout.md`
  - `.agents/reports/2026-04-22_migration-dryrun.md`

## CROSS_ERA

- **Meaning:** methodology, tooling, rules — not bound to an epoch.
- **Examples:** coding rules, naming conventions, frontmatter schema,
  renormalization style, never/always list.

## Sub-episode (not a separate era)

**April 2026 scalar-testing gauntlet** — ran under IMCOH_ABS +
COHORT_N9, produced nothing publishable. Documented in
`.agents/reports/2026-04-24_post-mortem-scalar-session.md`. Not a
separate era; a lesson inside COHORT_N9.
