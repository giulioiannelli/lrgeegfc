---
name: lane_e_00_preregistration
kind: preregistration
era: PAPER_FINALIZATION (Wave 1, lane E — encoding vs inference)
status: frozen_2026-08-31_before_any_number
created: 2026-08-31
scope: What lane E will test, how it will be gated, what would make encoding-vs-inference a real result and what would falsify it. Frozen before any lane-E number was computed. Deviations from this document are recorded in the report as deviations, never edited into this file.
pointers:
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - .agents/reports/2026-08-25_w0b-null-ladder.md
  - .agents/reports/2026-08-25_w0c-cohort-gate-and-tau.md
---

# Lane E — pre-registration

## Head

Lane E asks whether the resting brain after a transitive-inference task holds **the premises it was shown** (encoding) or **the ordering it had to compute** (inference), and whether those two are separable at all. The load-bearing statistic for the second half of that question, `T_infspec_pe`, is currently **uncalibrated**: it returns a large positive value on a no-signal input, so nobody knows what zero means for it. This lane's first job is therefore not to find a p-value but to build a formulation whose no-signal value is known. Only after that is demonstrated does anything get gated.

---

## 0. Notation, fixed here

Five phases per (patient, band): `A`, `B` (contiguous split halves of `rest_pre`), `task_learn` (ordered premises presented), `task_test` (novel non-adjacent pairs judged), `rest_post`. `D_x` is the condensed cophenetic distance vector of the UPGMA tree of the diffusion communication distance at scale `s = τ·λ_max`, on the graph fixed by `PIPELINE_CONTRACT.md`.

```
e  = D_learn − D_A     e2 = D_learn − D_B      encoding change
f  = D_test  − D_learn                          inference-specific change (arm-invariant)
g  = D_test  − D_A     g2 = D_test  − D_B       whole-task change
p  = D_post  − D_B     p2 = D_post  − D_A       persistence
```

```
T_test       = ½[ρ(g,p)     + ρ(g2,p2)   ]
T_learn      = ½[ρ(e,p)     + ρ(e2,p2)   ]
T_infspec    = ½[ρ(f,p)     + ρ(f, p2)   ]
T_infspec_pe = ½[pr(f,p|e)  + pr(f,p2|e2)]
```

`ρ` is Spearman; `pr` is the first-order partial Spearman. Every functional is symmetrised over the arbitrary A/B arm assignment and every pairing is cross-baseline.

---

## 1. Substrate, gate and family — declared before any number

**Substrate.** `imcoh_abs`, mst-union backbone, **knob-integrated over `f ∈ {0.07, 0.10, 0.14, 0.20}`** via `canonical_graph_ensemble`. No single-fraction number is reportable. The knob-integrated readout is the **median over fractions** of the per-patient margin, with the across-fraction spread reported alongside; the per-fraction surfaces are reported in full.

**Scale grid.** The locked 16-point `s = τ·λ_max ∈ logspace(0, log10 180, 16)`. Reported per scale, never best-scale.

**Null.** Matched-strength on the dense FC, then sparsified — the same draw reused across the four fractions inside a cell, so fraction-to-fraction differences are attributable to the backbone and not to surrogate noise. R = 200.

**Cohort statistic.** The margin `m_k = obs_k − median_r surr_{k,r}` across all n = 10 patients. No patient is dropped. Nothing is tested against zero.

**Multiplicity family, chosen deliberately and declared now.** W0-A measured the 16-scale axis at `n_eff = 1.51` independent tests and W0-C measured a 28-point axis at `n_eff = 1.2–1.9`; per-scale BH over the grid therefore over-charges by roughly an order of magnitude and is the wrong family. The **primary** gate is one sign-flip cluster-mass test (`axis_cluster_gate`) per (band, functional), which collapses the scale axis to a single test. The **primary family is the six bands within one functional**, because the scientific question is per-functional ("in which band does encoding persist?" is one question; "in which band does inference persist beyond encoding?" is another) and the paper reports them separately. The **24-cell family** (4 functionals × 6 bands) is reported as the conservative secondary in the same table, and the per-scale grid is reported in full because the shape of the curve is itself a result. Neither family is chosen after seeing which one clears.

**Effect size and robustness.** Every gated cell carries a bootstrap CI on the mean margin, a rank-biserial effect size, and the leave-one-patient-out swing of the verdict (the whole gate re-run per drop, not the statistic re-computed).

---

## 2. The 5-point critical preamble

**(1) The claim.** Two claims, separately gated.

- **C1 (persistence).** The encoding change and the inference-specific change each leave a trace in `rest_post` that a matched-strength surrogate does not reproduce.
- **C2 (separability).** Encoding and inference are *distinguishable* — the resting trace depends on which task block played which cognitive role, and/or the two persist with different band, scale or spatial signatures.

C2 is the claim that makes this about cognition. C1 without C2 is a generic two-block task effect.

**(2) The null.** Three, doing different jobs.

- **Matched-strength** (incumbent, for C1): each phase's dense FC shuffled strength-preservingly, all five pushed through the identical pipeline. Controls for node strength and for the construction, because the surrogate arc is built the same way as the real one.
- **Learn/test role swap** (new, for C2): the five *observed* graphs, with `task_learn` and `task_test` exchanged, per patient. The construction, the split-half baseline structure, the durations and the graphs are all preserved exactly; only the semantic assignment changes. The realization set is the full `2^10 = 1024` cohort sign pattern, so the cohort test is exact rather than Monte-Carlo.
- **No-signal sham arc** (new here as a *gate calibration*, not as a probe): five pseudo-phases carved from a single resting recording — the ordered variant and the block-shuffled (drift-free) variant — pushed through the *entire gate*, matched-strength surrogates included. W0-B ran the sham against **zero** and found every functional positive at 16/16 scales; nobody has run the sham against **its own surrogate**, which is what the gate actually does.

**(3) The strongest plausible alternative each null must control for.**

For C1: that the arc's construction manufactures a positive value out of shared structure between five phases of one subject, independently of any task. This is documented, quantified, and larger than the effect in at least one cell (β `T_infspec_pe`: sham +0.119 vs real +0.090).

For C2: that `task_learn` and `task_test` differ only in ways that have nothing to do with cognition — duration, recording order, fatigue, electrode drift, and the number of Welch segments available to the estimator. A swap null cannot separate "premises vs inference" from "first task block vs second task block". This is the single most important limitation of the lane and is stated wherever C2 is stated.

**(4) Whether the nulls actually control for it, by mechanism.**

Matched-strength **does** control for the construction offset, by mechanism: the surrogate arc has the same five phases, the same durations, the same estimator and the same shared-subject structure — it differs only in the strength-preserving shuffle. So the margin, unlike the raw value, has no reason to inherit the construction baseline. **But that is an argument, not a measurement**, and the measurement has never been made. It is E2's first deliverable: run the gate on the no-signal sham and check the margin sits at zero and the false-positive rate is nominal. What matched-strength *cannot* reach: the coherency estimator, the band split, session nonstationarity, drift, artifact epochs (W0-B's ladder), and it cannot reject that the trace is carried by coherence magnitude rather than lag (N1, which the headline β trace also fails, at a cohort median margin of −0.003 — enc/inf and the headline trace are in the same evidential position on that rung and neither may be described as evidence about lagged coupling).

The role swap **does** control for the construction offset by an exactly different mechanism — the construction is not merely matched, it is *identical*, because the same five graphs are used. It **cannot** control for anything that distinguishes the two task blocks other than their cognitive content, as stated in (3).

**(5) What would falsify each claim, and what remains.**

- If the gate's margin on the **no-signal sham** is significantly positive, or its false-positive rate exceeds nominal, then the gate does not calibrate for these functionals and **no p-value from it is reportable** — including W0-C's δ/α/β q = 0.043. That is a lane-killing outcome and will be reported as the headline if it occurs.
- If `T_learn` and `T_infspec_pe` fail the cohort gate under the declared family on the knob-integrated substrate, C1 is negative.
- If the real arc is indistinguishable from its learn/test-swapped counterpart, C2 is negative: the readout does not know which block was which, and encoding-versus-inference is not separable by this method however significant C1 is.
- If C2's positive rests on band, scale or magnitude differences whose directions disagree across bands and which do not survive the declared family, it is suggestive and must be written as suggestive.

---

## 3. What would make this a real result

All four, together:

1. **The gate is demonstrated calibrated for these functionals** — measured on no-signal input, not argued.
2. **C1 clears** on the knob-integrated substrate under the declared family, with LOO robustness reported.
3. **C2 clears** on the construction-preserving swap null, i.e. the resting trace demonstrably depends on which task block was the premises block.
4. **The confound in (3) is addressed** or explicitly conceded — either a block-order/duration control exists, or the claim is stated as "the two task blocks are distinguishable and the cognitive reading is one of several", not as "inference is separable from encoding".

Anything less is reported as what it is.

## 4. What would make this a negative, and that is an acceptable outcome

A negative here is: the gate calibrates, C1 clears for the whole-task and possibly the encoding functional, and **the swap test finds no dependence on which block was which**. That would say the resting trace is a two-block task effect that the LRG readout cannot resolve into encoding and inference. It would be reported as the headline, in the first line, and no further constructions would be searched over to escape it.

Every variant tried is listed in the report, including the ones that failed. No construction is adopted after the fact because it produced significance.

---

## 5. Pre-declared analysis list

| step | what | pre-declared decision rule |
|---|---|---|
| E1 | knob-integrated 5-phase grid, 4 fractions × 6 bands × 16 scales × 10 patients, R = 200 matched-strength | primary gate = per-band axis-cluster, family = 6 bands within functional; secondary = 24 cells |
| E2a | full gate on the no-signal sham arc (ordered + block-shuffled), matched-strength included | calibrated iff the cohort sham margin is not significantly positive and the measured FPR at α = 0.05 is ≤ 0.10 per cell |
| E2b | learn/test role-swap exact permutation test, 2^10 cohort sign patterns | reportable iff its FPR on the no-signal sham is ≤ 0.10 |
| E2c | the antisymmetry of `T_infspec` under the swap, derived and verified numerically | reported as a property, whether or not it helps |
| E3 | separability: swap test per band; paired `T_learn` vs `T_infspec_pe` contrasts; per-pair / anatomical distribution of the two persistences | a dissociation counts only if it is tested paired within patient and survives the declared family |
| E4 | scale: `n_eff` per functional, scale-position centroids paired within patient, re-run on Lane S's readout if it delivers one | directions must agree across bands or the result is suggestive |

Deviations from this table are recorded as deviations in `.agents/reports/2026-08-31_lane-e-encoding-vs-inference.md`.
