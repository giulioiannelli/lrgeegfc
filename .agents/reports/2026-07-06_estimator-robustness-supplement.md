---
name: estimator-robustness-supplement
type: report
era: "IMCOH_ABS × COHORT_N10"
status: current
created: 2026-07-06
audience: referee-facing (supplementary methods draft)
pointers:
  - scripts/01_compute/audit/audit_149_estimator_robustness_regate.py
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - scripts/01_compute/audit/audit_151_localization_rhosym.py
  - scripts/01_compute/audit/audit_152_consolidation_arc_rhosym.py
  - scripts/01_compute/audit/audit_153_cross_phase_taxonomy_rhosym.py
  - scripts/01_compute/audit/audit_154_per_node_trace_decomposition_rhosym.py
---

# Supplement — the cophenetic trace is not an arbitrary-half artifact (ρ_split → ρ_sym)

## Head

The cophenetic trace is estimated on a **split baseline**: the two rest_pre halves A
and B are used as independent baselines for the task-side and rest-side contrasts, so
shared baseline measurement error cannot inflate the trace. Bare ρ_split hard-codes
*which* half feeds *which* arm — an arbitrary choice. We show (i) this choice is
ill-conditioned for near-zero patients but (ii) leaves every cohort verdict intact,
and we adopt the **symmetric estimator ρ_sym** that averages both arm assignments and
removes the artifact at zero data cost. **All flagship results are estimator-invariant.**

## 1. The split-baseline design and its one free choice

For per-pair LRG cophenetic distances `D_x` in phases x ∈ {rest_pre_A, rest_pre_B,
task_test, rest_post}, the trace is

    ρ_split = Spearman(D_task − D_preA,  D_post − D_preB).

Using two *independent* halves A ≠ B for the two contrasts is deliberate: a shared
baseline `D_preA` on both axes would correlate their measurement errors and inflate ρ.
The cost is that the assignment (A→task-side, B→rest-side) is arbitrary — the mirror
assignment (B→task-side, A→rest-side) is equally valid.

## 2. ρ_split is ill-conditioned for near-zero patients (audit_149)

The two halves disagree on ~45 % of pairs (half the data each). Swapping A↔B changes ρ
by up to 0.48 and **flips the sign of 20 of 60 (patient × band) cells** — but only for
low-reorganizability / near-zero patients (e.g. β: Pat_15 +0.083↔−0.250, Pat_13
+0.208↔−0.077). Strong tracers are rock-stable (split-uncertainty ½|ρ_AB−ρ_BA| ≤ 0.035;
Pat_02/05/08). The instability tracks **proximity to zero**, not the trace/no-trace
label: it is estimator conditioning, not biology.

## 3. The fix: ρ_sym (symmetric estimator, adopted)

    ρ_sym = ½ [ Spearman(D_task − D_preA, D_post − D_preB)
              + Spearman(D_task − D_preB, D_post − D_preA) ].

Averaging the two equally-valid assignments makes the estimator invariant to the
half-labelling. It uses the same data and the same matched-strength surrogate; strong
tracers are unchanged. Per-patient we now report **ρ_sym ± ½|ρ_AB−ρ_BA|** and label
`|ρ_sym| < 1 SE` **"undetermined"** — the borderline third of the cohort whose sign is
estimator noise rather than a biological verdict.

**A rejected alternative — full (shared) baseline.** Collapsing to a single baseline
removes the half-choice but reintroduces exactly the shared-error inflation the split
design prevents: for Pat_02 the full-baseline observed statistic 0.775 versus surrogate
median 0.767 leaves a de-inflated signal of 0.008 — catastrophic cancellation. The
split design is vindicated; ρ_sym is the *free* refinement of it, not a return to the
shared baseline.

## 4. Estimator-invariance is proven, not assumed

`audit_149` re-ran the cohort matched-strength gate under **both** estimators with
**independent per-estimator surrogates** (so agreement is not a shared-draw artifact):

| estimator | α gate p | β gate p | δ/θ/low-γ/high-γ | verdict flips |
|---|---|---|---|---|
| ρ_split | 0.002 | 0.007 | all fail | — |
| **ρ_sym** | **0.024** | **0.032** | all fail | **0 / 6 bands** |

Both rows are computed at the **same** R=200 with independent surrogates, so the larger
ρ_sym p-values are **not** a grid effect — they are the estimator doing its job. ρ_sym is
a less *extreme* statistic than ρ_split by construction: it pulls the ill-conditioned
near-zero patients toward zero (rather than letting a lucky half-assignment push them
positive), so fewer patients sit spuriously above their surrogate and the cohort Wilcoxon
is less extreme. Strong tracers are unchanged. This is the correct price for an unbiased
test: α and β still clear 0.05 with margin, the other four still fail, and **zero verdicts
flip**. Trading a little significance for the removal of an arbitrary-half artifact is the
right trade — a slightly larger but *arbitrary-free* p beats a smaller but *arbitrary* one.

## 5. Every flagship result reproduces under ρ_sym

New audits only; the locked ρ_split audits (63/83/103/105/144) were **not edited**.

| result | audit | ρ_sym outcome |
|---|---|---|
| cohort gate | 150 | **α p=.024, β p=.032 CLEAR**; rest fail |
| β → OFC localization | 151 | **OFC carrier q=.025** (both epi); sensorimotor + PFC depleted q≤.033 |
| consolidation arc | 152 | **β T_infspec_pe p=.0098, β-ONLY** (others min p=.246); α/β T_learn p=.014/.032 |
| anchor/trace/reset taxonomy | 153 | β trace-dominant (comp_trace .263 vs anchor-dominated elsewhere); anchor ≠ hubness |
| per-node carrier/anti | 154 | node-level ρ_sym↔ρ_split ρ=.96 (β); carrier/anti split + anti-node property preserved |

Cross-checks: audit_152's arm1 reproduces the locked audit_83 full-graph statistic to
1e-16; audit_153's arm1 ρ_split matches the locked audit_63 gate bit-exactly (60/60 cells).

## 6. Reproducibility / performance

The strength-preserving 4-cycle ±δ surrogate shuffle is numba-JIT-compiled with the RNG
draws generated outside the compiled loop, making it **bit-identical** to the reference
pure-Python implementation (max|Δ| = 0.0) while running **98× faster** (882 → 9 ms per
shuffle). This made the R=200 six-band gate a ~5-minute run and is what allowed the full
estimator re-analysis to be done exhaustively rather than sampled.

## Bottom line

The trace is a property of the data, not of the arbitrary half-labelling. ρ_sym is the
estimator of record; ρ_split is retained only as a supplementary column. Grassmann
`d_G(k)` and raw |ImCoh| are separate measures, unaffected by this change.
