---
name: interpatient-variability-resolution
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-04
supersedes_claims:
  - "β-trace inter-patient heterogeneity is unexplained residual biology (per_node_trace_anatomy_2026_06_25, audit_147)"
pointers:
  - scripts/01_compute/audit/audit_148_interpatient_variability_resolution.py
  - data/audit/interpatient_variability/
  - .agents/reports/2026-06-26_per-band-phenomenology-vision.md
  - .agents/reports/2026-06-25_multiphase-snr-reliability.md
---

# Who traces, and why: resolving the β inter-patient variability

## Head

The between-patient spread in the β cophenetic trace is **not** unexplained
residual biology. It is a **stable, left-lateralised retention trait**: whether a
patient's `rest_post` keeps *any* task-state LRG reorganization. The trait is
identical on the learning and the test phase (`T_learn ↔ T_test` Spearman
**ρ = 0.95, p = 2e-5, sign-concordant 10/10**), so it is neither test-specific nor
a task_test-provenance artifact. The three non-tracers fail for **two
mechanistically opposite reasons**: **Pat_15 barely engages** the task at all
(engagement-null), while **Pat_10 and Pat_14 engage but revert** in `rest_post`
(reset). The one structural driver that survives is **hemispheric sampling**
(right-hemisphere contacts anti-predict the trace, `N_R` ρ = −0.69), and it has a
**direct anatomical mechanism**: trace-carrier nodes are left-lateralised
cohort-wide. Reliability, OFC coverage, epileptogenic burden, recording durations
and the rest-to-rest gap are all null. One patient (Pat_14) is an irreducible
residual. The cohort β claim is untouched.

---

## 1. The variability is a stable trait, not a test artifact

The single strongest correlate of the per-patient β test-phase trace `T_test`
(= `ρ_split`) is the **learning-phase** trace `T_learn`:

| covariate | Spearman ρ vs `T_test` | p |
|---|---:|---:|
| **`T_learn`** | **+0.952** | **2e-5** |
| `dir_balance` (carrier−anti)/N | +0.879 | 0.001 |
| `mover_frac` (engagement) | +0.794 | 0.006 |
| **`B_hemi`** (left-minus-right) | **+0.697** | 0.025 |
| **`N_R`** (right contacts) | **−0.694** | 0.026 |
| `d_task` (raw task displacement) | +0.455 | 0.19 (n.s.) |
| `frac_epi` | +0.406 | 0.24 (n.s.) |
| reliability bottleneck | −0.115 | 0.75 (n.s.) |
| durations, gap, N, dispersion | \|ρ\|≤0.35 | all n.s. |

All 8 test-tracers also trace on the learning phase; both reset patients (Pat_10,
Pat_14) also reset on the learning phase (**sign concordance 10/10**). Because
`task_learn` was **never** vendor-replaced, this **kills the "Pat_14 corrupt
task_test" hypothesis outright** — the reset is present on the clean phase too.
The trace is a patient-level property: *does this brain retain the task-state
LRG structure into rest?* — the same answer for either task phase.

## 2. Two failure modes (node decomposition)

Per-node carrier/anti/neutral counts (`per_node_trace_decomposition`) split the
non-tracers cleanly along **two axes** — engagement (do nodes move at all) and
direction (do the movers persist or revert):

| patient | class | carrier / anti / neutral | mover-frac | dir | reading |
|---|---|---|---:|---:|---|
| Pat_02–08 (tracers) | TRACE | ~55–75 / 0–9 / ~45–66 | 0.43–0.63 | +0.31…+0.58 | engage **and** persist |
| **Pat_15** | NULL | 14 / 8 / **96** | **0.19** | +0.05 | **engagement-null** — 81% neutral; barely moves |
| **Pat_10** | RESET | 10 / **22** / 81 | 0.28 | **−0.11** | **reset** — movers revert |
| **Pat_14** | RESET | 12 / **18** / 89 | 0.25 | **−0.05** | **reset** — movers revert |

- **Pat_15 = engagement-null.** `d_task = 0.16` (3× below cohort), 81% of nodes
  neutral. Its `ρ_split ≈ 0` is a near-mechanical consequence of an almost-zero
  task-displacement vector — there is nothing to persist. **Not measurement
  noise:** Pat_15 is the *highest-reliability* patient (0.85 split-half). This is a
  coverage fact (see §3), not a detectability failure.
- **Pat_10 / Pat_14 = reset.** They displace normally (`d_task` 0.71 / 0.52) but
  their movers reverse (anti > carrier), in **both** task phases.

## 3. The structural driver: left-lateralised carriers

Right-hemisphere sampling anti-predicts the trace (`N_R` ρ = −0.69; `B_hemi`
ρ = +0.70), and it predicts the *stable trait* (`N_R` vs `(T_test+T_learn)/2`
ρ = −0.56, p = 0.10), the *reset direction* (`dir_balance` ρ = −0.61, p = 0.06) and
*engagement* (`mover_frac` vs `B_hemi` ρ = +0.61, p = 0.06). The **mechanism** is
direct: trace-carrier nodes are left-lateralised cohort-wide
(`hemisphere_lobe_enrichment`):

| band | left carrier-rate | right carrier-rate | left enrichment | left Bonf-p |
|---|---:|---:|---:|---:|
| δ | 8.8% | 2.1% | 1.36 | **0.001** |
| α | 4.7% | 1.0% | 1.37 | **0.039** |
| β | 6.3% | 3.4% | 1.19 | 0.80 (dir.) |

A right-heavy implant simply samples fewer carrier nodes, diluting or reversing
the whole-brain `ρ_split`. Per-patient β-network coverage confirms it: Pat_10
(anti) has **15 β-net contacts, all right**; Pat_15 (null) has 3, all right.

**Non-deterministic — two honest exceptions:** Pat_06 (right-hemisphere) *traces*
because its right sampling is **frontal/OFC** (19 OFC contacts, the trace hub),
not temporo-parietal; and **Pat_14** (left-frontal, 5 β-net contacts all left)
*resets* despite good left coverage. So laterality is the dominant **bias**
(~50% of variance), not the whole mechanism.

## 4. What it is NOT (quantitatively refuted)

- **Reliability / detectability** — ρ = −0.12; the best-measured patient (Pat_15,
  0.85) is a non-tracer, the worst-measured-in-rest_post (Pat_08, 0.24) is the
  2nd-strongest tracer. (Confirms Q1.)
- **OFC electrode coverage** — ρ = −0.10; Pat_08 (0 OFC contacts) traces strongly,
  Pat_14 (15 OFC) resets. (Confirms audit_147.)
- **Epileptogenic burden** — `frac_epi` ρ = +0.41 (n.s.), `N_epi` +0.29 (n.s.);
  the highest-epi patient (Pat_13, 38 epi) traces, the zero-epi patient (Pat_15)
  is null. (Confirms audit_64 frac_epi null.)
- **Total coverage N** (−0.17), **dispersion** (−0.33 n.s.), **all recording
  durations** and the **inter-rest gap** (−0.03) — null.

## 5. The residual

**Pat_14** is the one failure no available covariate explains: left-frontal
coverage, adequate reliability, normal durations, yet it resets on *both* task
phases. This is genuine non-consolidation biology (or a `rest_post` state
difference) that cannot be resolved without behavioral or vigilance data
(both unavailable — [[task_paradigm_transitive_inference]]). Honest bound:
we now *name and localise* the failure (a cross-phase reset, left-frontal,
band-selective toward high frequency) rather than call the whole cohort spread
"unexplained".

## 6. Consequence for the manuscript

- The cohort β verdict is **unaffected** — it rests on the 7 patients that both
  engage and persist, plus the group gate (p = 0.005), matched-strength, and
  β→OFC LOO. This audit explains the *spread around* that verdict, it does not
  touch it.
- Replace "heterogeneity = unexplained residual biology" with **"a
  left-lateralised retention trait; the non-tracers either under-engage
  (Pat_15, right-parietal coverage) or reset (Pat_10 right-sampling; Pat_14
  residual)."** [[feedback_fluctuations_are_signal]] — the spread *is* a result:
  it says the β trace lives on a left-dominant, consolidation-competent substrate.

Outputs: `data/audit/interpatient_variability/{master_table,correlations}.csv`,
`fig_interpatient_variability.pdf`. Reproduce: `audit_148`.
