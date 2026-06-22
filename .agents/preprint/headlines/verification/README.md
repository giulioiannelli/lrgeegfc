---
name: headlines-verification-index
era: IMCOH_ABS_COHORT_N10
status: current
kind: verification-index
scope: one self-contained verification brief per investigation flagged as most interesting for neurophysiology / epilepsy. Hand a brief to its owning agent; it returns a verdict + updated brief. This (crystallization) chat does not run them.
updated: 2026-06-22
---

# Verification briefs — hand these to the per-headline agents

**Head.** Each file here is a **self-contained verification brief**: an agent can
take it, run it, and return a pass/fail verdict without further context. Every
brief follows the project's **mandatory 5-point critical preamble** (claim · null ·
strongest alternative · does-the-null-control-it-by-mechanism · falsification +
residual limits), then exact steps, data/scripts, pass/fail thresholds, caveats,
and **what to return**. No brief is run in the crystallization chat.

## How to use

1. Pick the brief for the investigation.
2. Hand the file to the owning agent (routing in `../README.md` §6).
3. The agent runs it on **repo data + cached CSVs** (no behavioral data exists),
   honors the global rules (matched-strength is the referee; no-hardcoded-tables;
   report LOO-max; X-epi figures yes / C4 no), and **returns a verdict + the
   updated brief** (fill the "Result" stub at the foot).

## The set

- **`verify_beta_spares_alpha_recruits.md`** — N1.6 bridge: β spares / α recruits
  the epileptic core. *(Mostly evidenced; verify it's a genuine cross-band
  difference, not node-count.)* Owner: white-matter / localization + null-model.
- **`verify_virtual_resection.md`** — epilepsy: does removing the predicted
  epileptogenic community fragment the network beyond a matched node-removal null?
  *(New; mechanistic validation of N3.)* Owner: epi-marker-trace-analysis.
- **`verify_stage2_arc_duration.md`** — gates N2's inference component: does the
  whole-brain inference-specific arc survive a length-matched null? *(audit_103b is
  live.)* Owner: inference chat.
- **`verify_seedfree_epi_and_rigidity.md`** — N3.3 + a candidate sub-result: can
  the SOZ community be found seed-free, and is epileptic routing rigid across
  phases? *(audit_118/119/120/121 live.)* Owner: epi-marker-trace-analysis.

**Replay (N4)** has its execution+verification brief at the canonical task-trace
scope location:
`../../../guides/task-persistence-investigation/2026-06-22_replay-states-multiscale-reinstatement.md`
(task-trace tooling scope lives there by project rule).

## Two standing rules for every brief

- **Matched-strength (or a node-count decimation control where a subset is
  removed/rebuilt) is mandatory** — geometry/count baselines are not controls
  (the lesson that killed the taxonomy dissociation and anatomy v1).
- **No brain–behavior test** — TI performance data is unavailable and will not be
  obtained (PI 2026-06-22).
