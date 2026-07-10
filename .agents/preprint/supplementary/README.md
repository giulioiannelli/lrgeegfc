---
name: supp-index
era: IMCOH_ABS_COHORT_N10
status: current
kind: supplementary
scope: Supplementary negative-control analyses — the pairwise-blind ladder and the within-recording drift null. Establishes what simple pairwise descriptors and session drift CAN and CANNOT reproduce, and thereby what the multiscale (cophenetic) comparison uniquely contributes.
updated: 2026-07-10
---

# Supplementary — pairwise & drift controls

**Head.** Two negative-control analyses that fence the multiscale trace from
below (classical pairwise scalars) and from the temporal-confound side
(session drift). Together they establish that (a) the *usual* low-level graph
descriptors do not reproduce the trace, (b) what the cophenetic comparison
uniquely adds is **band-selectivity**, and (c) **session drift is an orthogonal
confound to matched-strength** — raw-FC *whole-task* traces fail it entirely; the
whole-task β trace survives it (fair locked C2 + audit_167). **Caution:** a windowed
drift null is *not* valid for the four-phase *conditional* functionals — a brief
"inference-specific is drift-confounded" reading was **withdrawn** (the sham returns
spurious positives on a no-signal pre-task arc); the inference-specific trace is
validated by matched-strength, not drift. The β→OFC whole-task flagship survives.

## Contents
- [`S1_pairwise_descriptor_ladder.md`](S1_pairwise_descriptor_ladder.md) — classical
  pairwise descriptors (strength, weighted clustering, centralities, closeness)
  run through the identical ρ_sym + matched-strength pipeline. Verdict: local
  scalars blind (strength degenerate under MS by construction); band-selectivity
  unique to cophenetic. audit_166.
- [`S2_drift_controls.md`](S2_drift_controls.md) — drift nulls. Raw-FC *whole-task*
  traces are drift; cophenetic β whole-task survives (fair locked C2 + audit_167); α
  whole-task construction-dependent. **Methodological cautionary record:** a windowed
  drift null is INVALID for the *conditional* inference-specific functional (returns
  spurious positives on a no-signal pre-task arc); a brief "drift-confounded" reading is
  **withdrawn** — N2's inference-specific claim rests on matched-strength, not a drift
  null. C2 + audit_167 + audit_168.

## Provenance
- Scope reports: `.agents/guides/task-persistence-investigation/2026-07-09_pairwise-descriptor-ladder.md`
- Scripts: `scripts/01_compute/audit/audit_166_pairwise_descriptor_ladder.py`,
  `audit_167_drift_null_ladder.py`, `audit_168_duration_matched_drift_infspec.py`
- Data: `data/audit/pairwise_descriptor_ladder/`, `data/audit/drift_null_ladder/`,
  `data/audit/duration_matched_drift_infspec/`
- Figure: `data/outputs/figures/drift_infspec_duration_matched.pdf`
- Control battery: also logged in `locked/CONTROLS.md` (drift = new mandatory control).
