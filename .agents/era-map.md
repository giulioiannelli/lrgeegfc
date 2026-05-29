---
name: era-map
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-05-28
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

## COHORT_N9 — 2026-04-22 → 2026-04-25 (superseded by COHORT_N10)

- **Carrier:** IMCOH_ABS + 9-patient locked cohort. Members: Pat_02,
  03, 05, 06, 07, 08, 10, 13, 15. Excluded: Pat_14 (corrupt
  `task_test`).
- **Trigger:** data normalization + quality pass on 2026-04-22.
- **Superseded:** 2026-04-25 when Pat_14's `task_test.mat` was
  vendor-replaced, restoring n=10 at the cross-phase level. Cite COHORT_N9
  only for the n=9 snapshot in `.agents/reports/2026-04-24_pipeline-status.md`;
  the live cohort tag is COHORT_N10.
- **Quirks (historical framing):**
  - Pat_03 at 1024 Hz was framed as "negative control, flagged
    distinctly". This was retired 2026-05-18 — Pat_03 is now a full
    cohort member at analysis layer; sampling rate is handled at the
    config layer (`FS_OVERRIDES`, `nperseg_for_fs(fs)`).
  - Pat_10 task rows [53, 54, 55] dropped at load.
- **Landmark docs:**
  - `.agents/guides/03_implementation/data-layout.md`
  - `.agents/reports/2026-04-22_migration-dryrun.md`

## COHORT_N10 — 2026-04-25 → present

- **Carrier:** IMCOH_ABS + 10-patient locked cohort. Members:
  Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15.
- **Trigger:** Pat_14 `task_test.mat` vendor replacement on
  2026-04-25; cross-phase cohort returned to n=10.
- **State:** the live cohort. Every current verdict and preprint claim
  references n=10. The only legitimate per-cell drop is Pat_15 in
  LRG-native β contexts (right-hemisphere-only implant, biology-driven
  anti-alignment); the older n=8 "pro-cohort restriction" was retired
  the same day.
- **Sign convention locked 2026-05-26:**
  `T_d := d(rest_pre, task) − d(task, rest_post)`; positive = TRACE
  at every layer (raw FC, D_coph, KC, Grassmann). Wilcoxon
  `alternative='greater'`. Per-patient `(T > 0).sum()`. Surrogate
  upper-tail p = `mean(s ≥ obs)`.
- **Landmark docs:**
  - `.agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md`
  - `.agents/reports/2026-05-07_epileptic-n10-revisit.md`
  - `data/reports/section_5_lrg_trace/README.md`
  - `.agents/preprint/WRITING_GUIDE.md` (preprint routing)

## CROSS_ERA

- **Meaning:** methodology, tooling, rules — not bound to an epoch.
- **Examples:** coding rules, naming conventions, frontmatter schema,
  renormalization style, never/always list.

## Sub-episode (not a separate era)

**April 2026 scalar-testing gauntlet** — ran under IMCOH_ABS +
COHORT_N9, produced nothing publishable. Documented in
`.agents/reports/2026-04-24_post-mortem-scalar-session.md`. Not a
separate era; a lesson inside COHORT_N9.
