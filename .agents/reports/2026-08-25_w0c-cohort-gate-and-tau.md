---
name: 2026-08-25_w0c-cohort-gate-and-tau
kind: report
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 0, lane W0-C)
status: current
created: 2026-08-25
scope: The locked cross-patient gate every lane inherits, and the honest verdict on what the diffusion scale τ buys, decomposed into detection / characterization / selection / scale-units plus the five-phase encoding-vs-inference arc. Everything recomputed fresh from the cached |ImCoh| arrays on the incumbent substrate (imcoh_abs × mst-union@0.20) under the incumbent matched-strength null; no number is sourced from a prior report, CSV or summary.
pointers:
  - src/lrg_eegfc/utils/metrics/cohort_gate.py                       # the locked gate
  - src/lrg_eegfc/workflow/phase_graphs.py                           # injectable substrate
  - scripts/01_compute/paper_final/w0c_01_gate_and_tau_grid.py       # master four-phase grid
  - scripts/01_compute/paper_final/w0c_06_encinf_scale_grid.py       # five-phase enc/inf grid
  - data/paper_final/w0c_gate_tau/                                   # all artifacts
  - .agents/plans/active/2026-08-25_paper-finalization-master-plan.md
---

# W0-C — the cohort contract, and what τ actually buys

## Head

**The gate is locked and calibrated; the multiscale claim is not.** Recomputed fresh on the incumbent substrate, the β trace holds across the entire scale range (cluster q = 0.019, family = 6 bands) — but every claim that the *scale axis itself* buys something failed a direct test. The hierarchy detects nothing raw FC does not (`coph|raw` clears 0 of 168 cells under two residualizers of very different strength, with 75–90% of the cophenetic variance surviving residualization, so this is not a power artifact). The scale profile carries no band information a cross-patient classifier can use (accuracy 0.200 vs chance 0.167, permutation p = 0.31), and α is not more scale-tuned than β when tested per patient (|ρ(margin, log s)| = 0.259 vs 0.262, p = 0.90) — the "β scale-invariant vs α scale-tuned" dissociation does not reproduce. And "the hierarchy selects" turns out to rest on the difference between a significant and a non-significant result: when the raw-minus-hierarchy gap is tested directly, paired within patient, it clears for no band (best p = 0.080) and the band × representation interaction against β clears for none of δ/θ/α/low_γ. The edge-locality mechanism fails too, the second proposed mechanism to fail after scale-coherence. **What survives is a trace, not a multiscale trace** — and it is equally present at s = 0.05, where the diffusion is degenerate with raw FC. Two methodological findings are as consequential as the science: per-patient observed values co-vary with each patient's *own* null height at ρ = +0.70 (β) and the margin only halves it, so the confound is real and only partly fixed; and a 28-point scale sweep is worth **1.2–1.9 independent tests**, so whole-grid BH over 168 cells charges roughly twenty times what was actually paid and leaves the gate with no headroom at n = 10 — dropping any of eight patients wipes out every cleared cell.

---

## 1. The locked cohort gate (C1)

`src/lrg_eegfc/utils/metrics/cohort_gate.py`. One entry point, statistic-agnostic: it consumes `(obs: (K,), surr: (K, R))` and never loads data, builds a graph, draws a surrogate or computes a statistic. Substrate, null and functional are all injected, which is what lets it survive W0-A moving the substrate underneath it.

### 1.1 What the gate does, and why each clause is there

| clause | what | why |
|---|---|---|
| **margin** | tests `m_k = o_k − median_r s_{k,r}`, never `o_k` | patients' null floors differ; see §1.2 |
| **per cell** | every (band, scale) reported; nothing selects a maximum | §1.3 |
| **whole-grid BH** | default family is the whole submitted grid | the grid is the claim; but see §1.4 |
| **LOO of the verdict** | whole grid *and its BH step* re-run per dropped patient | a LOO of the statistic was constant by construction |
| **bootstrap CI + rank-biserial** | on the mean margin | effect size, not just a p |
| **exact signed-rank** | SciPy default (exact at n ≤ 50), not the normal approximation | at n = 10 the approximation is wrong in the tail that matters |
| **counts descriptive** | `n_above`, `frac_pos` reported, never consulted | project rule; the test is the gate |
| **`expect_n` guard** | a short cohort **raises** unless `exclusion_reason` is given | makes a silent patient exclusion impossible rather than discouraged |
| **calibration refusal** | `gate_grid(calibrate=True)` withholds `p`/`q` when the statistic/null pair fails | required for conditional statistics; see §1.6 |

Default cohort is n = 10 (`PATIENTS_4PHASE`), and both grid calls in this lane pass `expect_n=10`. No patient was dropped anywhere in this lane.

### 1.2 The obs-vs-null-height confound is real, and the margin only halves it

Cross-patient Spearman between each patient's observed `ρ_sym` and that patient's *own* matched-strength surrogate median, per (band, scale) cell, median [IQR] over the 28 scales:

| band | ρ(obs, own null) | ρ(margin, own null) |
|---|---|---|
| delta | +0.442 [+0.31, +0.57] | +0.218 [+0.04, +0.39] |
| theta | +0.285 [+0.13, +0.44] | +0.079 [−0.02, +0.27] |
| alpha | +0.473 [+0.26, +0.60] | +0.315 [+0.15, +0.48] |
| **beta** | **+0.697** [+0.61, +0.76] | **+0.509** [+0.35, +0.66] |
| low_gamma | +0.655 [+0.51, +0.73] | +0.533 [+0.45, +0.61] |
| high_gamma | +0.661 [+0.49, +0.76] | +0.503 [+0.36, +0.60] |

Two readings. First, the confound the project recorded is confirmed independently and is if anything larger than reported — a patient's raw `ρ_sym` is substantially a readout of how high that patient's own null sits, so cross-patient spreads of raw values are not effect-size spreads and must never be interpreted as such. Second, and this is new: **the margin does not remove it.** At β, low_γ and high_γ the margin still co-varies with the null height at ρ ≈ +0.50. Margin-testing is necessary and is locked, but it is a partial correction, and any per-patient ranking or "which patients carry the effect" statement remains contaminated. Standardising by the surrogate spread (a z rather than a margin) is the obvious next lever and has not been evaluated here.

### 1.3 The gated grid, on the incumbent substrate

`imcoh_abs` × `mst_union_top_fraction @ 0.20`, matched-strength (R = 200, 20·N(N−1)/2 swaps), 28 scales `s ∈ [0.05, 180]`, n = 10, whole-grid BH over 168 cells:

| band | raw (dense) q | scales cleared | q_min | s at q_min | N_eff at that s | margin_med |
|---|---|---|---|---|---|---|
| delta | *0.048 | 0/28 | 0.055 | 21.5 | 44 | +0.147 |
| theta | 0.063 | 0/28 | 0.265 | 21.5 | 41 | +0.051 |
| alpha | *0.048 | 1/28 | *0.048 | 4.73 | 96 | +0.135 |
| **beta** | *0.048 | **23/28** | *0.043 | 0.05 | 118 | +0.118 |
| low_gamma | *0.048 | 0/28 | 0.217 | 15.9 | 60 | +0.083 |
| high_gamma | 0.065 | 0/28 | 0.055 | 53.5 | 19 | +0.085 |

β is the only band that clears widely. **α clears at exactly one scale** (s = 4.73), which is a substantial demotion from the incumbent "12/16"; that figure came from uncorrected per-scale p-values. Full per-scale q grid in `data/paper_final/w0c_gate_tau/gate_grid.csv`.

### 1.4 A 28-point scale sweep is worth 1.2–1.9 tests, not 28

The eigenvalue participation ratio of the cross-scale correlation matrix of the per-patient margins, and the mean off-diagonal correlation:

| band | n_eff (participation ratio) | n_eff (Cheverud–Nyholt) | mean cross-scale r | axis-cluster p | q (family = 6 bands) | supra-threshold span |
|---|---|---|---|---|---|---|
| delta | 1.66 | 12.65 | +0.74 | 0.045 | 0.090 | s = 3.5–180 |
| theta | 1.19 | 6.27 | +0.91 | 1.000 | 1.000 | none |
| alpha | 1.17 | 5.94 | +0.92 | 0.030 | 0.090 | s = 0.05–180 |
| **beta** | 1.87 | 14.50 | +0.70 | **0.0032** | **0.019** | s = 0.05–180 |
| low_gamma | 1.15 | 5.47 | +0.93 | 0.114 | 0.140 | s = 2.6–180 |
| high_gamma | 1.65 | 12.62 | +0.75 | 0.116 | 0.140 | s = 0.05–1.4 |

The two n_eff estimators bracket the answer (they measure different things and are reported together rather than adjudicated); both say the sweep is far from 28 independent tests, and the mean cross-scale correlation of +0.70 to +0.93 is the direct reason. So whole-grid BH over (6 bands × 28 scales) is not merely conservative, it is conservative by roughly an order of magnitude.

The alternative offered in the library is `axis_cluster_gate`: one sign-flip cluster-mass test for the whole axis, flipping *whole patient profiles* so the correlation along the axis is preserved (Maris–Oostenveld), making the family 6 bands. Under it, **only β clears** (q = 0.019). α does not (q = 0.090). Both families are reported; neither is hidden. The per-cell grid is still reported in full because the shape of the curve is itself a result.

### 1.5 Leave-one-patient-out of the verdict

Under whole-grid BH over 168 cells, the full cohort clears 24 cells; dropping any of Pat_02, 03, 05, 06, 07, 08, 10 or 14 leaves **zero** cleared cells. Dropping Pat_13 leaves 36 and dropping Pat_15 leaves 33 (both are anti-aligned at β, Pat_15 canonically so).

This is not evidence that the β trace rides on one patient. It is the multiplicity design failing: the exact one-sided signed-rank p at n = 10 has a floor of 1/1024 and its next values are 2/1024, 3/1024, …; at n = 9 the grid coarsens to 1/512, and a 168-fold BH correction turns that single step into a cliff that the entire grid falls off at once. It is the sharpest possible demonstration that **whole-grid BH at n = 10 leaves the gate no dynamic range at all.** The same LOO under the per-band cluster family, whose permutation p is continuous, is in `data/paper_final/w0c_gate_tau/axis_cluster_loo.csv` and reported in §1.5b.

### 1.6 Calibration

Two tests, because a Gaussian toy tests almost nothing that matters.

**Toy null** (independent Gaussian obs and surrogates, each patient with its own mean and scale — the situation the margin exists for): FPR at nominal 0.10 / 0.05 / 0.01 = 0.094 / 0.042 / 0.011 (heteroscedastic) and 0.104 / 0.044 / 0.010 (homoscedastic), 4000 draws each. Slightly conservative at 0.05, exactly as the discreteness of the exact test predicts.

**Held-out-realization null on the real graphs** — the one that counts. One surrogate realization is promoted to the role of "observed" and tested against the remaining 199 of the *same patients*, so the test sees the real K, the real R, the real per-patient heteroscedasticity and the real surrogate ensemble. Under the null the true observation is exchangeable with its own surrogates, so these p-values are draws from the gate's null distribution on the actual data. Results over all 168 cells × 200 draws are in §1.6b.

A Kolmogorov–Smirnov test against the continuous uniform is *not* the right check here and is reported as a statistic only: the exact signed-rank p at n = 10 takes 1024 values, so its null distribution is discrete and KS rejects by construction however well calibrated the gate is. The statement that matters is `P(p ≤ a) ≤ a`, which the FPRs measure directly.

**Substantive known-null band.** θ is null in the hierarchy at every one of 28 scales (q = 0.27–0.78), null under the axis-cluster gate (p = 1.000, no supra-threshold cluster) and null in dense raw FC (q = 0.063). The gate does not fire on the band the project expects to be silent.

### 1.7 Conditional statistics need more, and the gate now enforces it

A partial correlation entangles the conditioning variable with the estimator and can return systematically positive values from input with no signal — documented in this project, where a sham arc built entirely inside pre-task rest produced a significantly positive conditional trace (p = 0.007) exceeding the real value. `gate_grid(calibrate=True, refuse_uncalibrated=True)` therefore measures each cell's false-positive rate and **withholds** `p` and `q` when the statistic/null pair fails, because a p-value from a miscalibrated pair is neither conservative nor liberal but uninterpretable, and reporting it with a caveat is worse than not reporting it.

What this does **not** cover, stated plainly: the held-out-realization test calibrates the statistic against *the null actually in use*. It cannot detect a bias the null shares. For that the statistic needs a **data-based placebo** — a no-signal arc built from the recording itself — which is lane W0-B's deliverable and is a live dependency, not something this lane closed. A role-permutation placebo (§4.1) is run here as a cheap partial substitute; it is not a replacement.

---

## 2. What τ buys — three claims, separated

### 2.1 (a) Detection — **closed, negative**

Does any scale carry cross-phase trace that the dense raw FC does not? Each phase's cophenetic geometry is residualized on that phase's raw edges, the residual is traced, and the identical matched-strength null goes through the identical residualization. Two residualizers of very different strength, at all 28 scales, all 6 bands, whole-grid BH over 168 cells per measure:

| measure | cleared cells | q_min |
|---|---|---|
| `coph\|raw`, rank-linear (removes the monotone dependence) | **0 / 168** | 0.143 |
| `coph\|raw`, binned conditional mean (removes *any* function of raw at 64 quantiles) | **0 / 168** | 0.148 |
| `raw\|coph`, rank-linear — the positive control | 110 / 168 | 0.034 |

The negative is not a power artifact, and this is the check that makes it load-bearing: residualizing the cophenetic geometry on the raw edges leaves **75–90% of its rank variance intact** (`residual_diagnostics.csv`; rank-linear 0.84–0.98, binned 0.74–0.96 across bands and scales), and the cophenetic distances are only weakly rank-correlated with the raw edges to begin with (|Spearman| = 0.314 cohort median). There is plenty of hierarchy left after removing raw; it simply does not trace.

The positive control fires strongly and in the opposite direction: `raw|coph` clears at all 28 scales for α, θ and low_γ and at 16/28 for δ. One caveat on it — θ clears in `raw|coph` while θ's own raw trace does not (q = 0.063), so the residualization operation is not neutral and its band pattern should not be over-read. That caveat does not touch the negative, which is the direction that matters: nothing survives on the hierarchy side.

**Verdict: the claim "the multiscale hierarchy detects persistent structure that raw FC cannot see" is dead, at every scale, under a residualizer strong enough to remove any function of raw and a residualizer weak enough that it cannot be blamed for over-removal.** This question is closed and should not be re-opened without a new argument.

### 2.2 (b) Characterization — **negative at n = 10**

Is the *shape* of ρ_sym(s) band-discriminative across patients?

**Direct test.** Per-patient margin profiles, z-scored over scale so only shape survives, then leave-one-patient-out nearest-centroid classification of band identity, against a within-patient label-permutation null (2000 permutations). **Accuracy 0.200 on 60 profiles against a chance level of 0.167; p = 0.31** (null mean 0.168, 95th percentile 0.267). The profile shape carries no cross-patient band information.

**The specific incumbent dissociation.** "β scale-invariant vs α scale-tuned" was built on a cohort Friedman failing to reject for β. Two problems, both confirmed here.

First, the per-patient test disagrees with it. Using ρ(margin, log s) — a scale-dependence statistic defined for **all 10 patients**, unlike CV and relative range which require a positive mean/max and silently dropped 4 of 10 — the median |ρ| is α **0.259** and β **0.262**, paired p = 0.90. **α is not more scale-tuned than β.** No band is significantly more scale-tuned than β (best: δ, p = 0.065 uncorrected).

Second, the Friedman was never evidence for invariance. Reproduced fresh: α χ² = 48.4, p = 0.0069; low_γ χ² = 58.4, p = 0.0004; β χ² = 20.6, p = 0.80; δ p = 0.90; θ p = 0.89; high_γ p = 0.50. But a Friedman blocked by patient tests whether the cohort *agrees on the ordering of scales*, not whether any patient's profile is scale-dependent. β's non-rejection is equally consistent with scale-invariance and with each patient having a differently-shaped profile — and the per-patient |ρ| above says β's profiles are exactly as scale-dependent as α's, so **cross-patient inconsistency, not invariance, is the better reading of β's flat Friedman**.

**Power.** Simulating at the observed within- and between-patient spread, the Friedman rejects a monotone scale trend of 0.25 / 0.50 / 0.75 / 1.00 / 1.50 × the within-patient SD on 3.8% / 17.5% / 38.2% / 70.8% / 99.8% of draws. Anything up to about three quarters of a within-patient SD is invisible to this design.

**Verdict: scale-invariance cannot be established at n = 10, and the α-tuned / β-flat dissociation does not reproduce under a per-patient test.** To establish invariance would need either a substantially larger cohort or an equivalence design with a pre-declared margin — neither is available here, and the claim must be dropped rather than softened.

### 2.3 (c) Selection — **real as a pattern, not demonstrated as an effect**

The pattern reproduces exactly. Dense raw FC fires non-selectively — δ, α, β, low_γ all clear at q = 0.048, θ and high_γ do not — while the hierarchy keeps β at 23/28 scales, α at 1/28, and rejects δ and low_γ at every one of 28 scales.

But that comparison is two verdicts, and the difference between a significant and a non-significant result is not itself significant. Tested directly — the gap `raw margin − cophenetic margin`, paired within patient:

| band | gap at s = 1.04 | p(gap > 0) | gap at s = 6.40 | p(gap > 0) |
|---|---|---|---|---|
| delta | +0.030 | 0.31 | +0.072 | 0.28 |
| theta | +0.063 | 0.097 | +0.096 | 0.097 |
| alpha | +0.035 | 0.42 | −0.015 | 0.69 |
| beta | +0.059 | 0.31 | +0.025 | 0.38 |
| low_gamma | +0.166 | 0.080 | +0.013 | 0.22 |
| high_gamma | +0.126 | 0.116 | +0.174 | 0.080 |

**No band's gap clears.** And the interaction — is the gap larger for a band the hierarchy rejects than for β? — clears for none of δ (p = 0.31–0.50), θ (0.19–0.25), α (0.46–0.62) or low_γ (0.12–0.28). The single p < 0.05 anywhere is high_γ at s = 6.40 (p = 0.032, uncorrected), a band that is null in both representations and not part of the claim.

So the honest statement is: **the selection is a real descriptive pattern that the direct test does not have the power to confirm at n = 10.** It is "not demonstrated", not "demonstrated absent" — the gap statistic is a difference of two noisy margins and is intrinsically low-powered. It cannot be written as an established result.

**Mechanism: edge-locality fails.** The standing hypothesis was that a rejected band's trace rides on a few strong pairs that UPGMA absorbs, while β's is distributed. The per-pair contributions to the raw trace say otherwise — every band is equally concentrated:

| band | Gini | top-1% share | trace retained after deleting the top 1% / 5% of contributing pairs |
|---|---|---|---|
| delta | 0.746 | 0.068 | +0.918 / +0.606 |
| theta | 0.748 | 0.066 | +0.934 / +0.689 |
| alpha | 0.734 | 0.063 | +0.894 / +0.492 |
| beta | 0.715 | 0.059 | +0.955 / +0.766 |
| low_gamma | 0.745 | 0.069 | +0.949 / +0.739 |
| high_gamma | 0.739 | 0.067 | +0.926 / +0.665 |

β is nominally the least concentrated and the most ablation-robust, in the predicted direction, but no paired contrast against β approaches significance (top-1% share: p = 0.38–0.54; retained fraction: p = 0.11–0.69). **This is the second proposed mechanism to fail**, after scale-coherence. The selection remains mechanistically unexplained and must be written as such, if it is written at all.

One incidental finding is sharper than the mechanism test and worth keeping. The hierarchy can express only **1–6% of the raw per-pair task reorganization** (squared Spearman between the raw and cophenetic reorganization vectors: 0.013–0.062 across all bands and all scales, no β/low_γ separation — β vs low_γ p = 0.50–0.96, in the wrong direction at the middle scales). The cophenetic geometry is therefore not a *filter* that keeps the robust part of the raw signal; it is a near-orthogonal, very low-dimensional readout that happens to trace in β. That is a different and more defensible framing than "the hierarchy retains what survives coarse-graining", and it is consistent with the detection negative.

### 2.4 (d) The scale axis, in units

No scale is described in words here without a number. Cohort medians on the mst@0.20 backbone (N ≈ 118 nodes):

| s | N_eff (resolved components) | m (communication neighbourhood) | ℓ(s) reach, mm | % of implant span |
|---|---|---|---|---|
| 0.05 | 118.0 | 1.0 | — (nothing connected) | — |
| 0.31 | 117.7 | 1.2 | 11.8 | 13% |
| 0.77 | 116.4 | 1.5 | 52.7 | 58% |
| 1.90 | 111.9 | 3.2 | 84.0 | 92% |
| 4.73 | 99.5 | 13.2 | 85.4 | 93% |
| 11.7 | 69.4 | 42.2 | 85.4 | 93% |
| 29.2 | 32.5 | 85.0 | 89.2 | 97% |
| 72.5 | 12.2 | 105.8 | 89.3 | 98% |
| 180 | 4.1 | 113.5 | 90.3 | 99% |

Electrode pitch 3.5 mm, implant span 91 mm, both cohort medians, recomputed fresh.

Two things follow. The **physical** reach saturates at ~92% of the implant span by s ≈ 2, so every scale in the trace-gate range is spatially near-global and no scale in it may be called micro, meso or local — this independently reproduces the July retraction of "α mesoscale". The **functional** resolution is a separate axis and does not saturate: at α's nominal peak s = 4.73 the process still resolves ~100 of 118 components while each node mixes with only ~13 others. So the scale axis is a *functional* coarse-graining axis, not a spatial one, and the two must never be conflated. A UPGMA-tree cut was tried as a third readout and discarded — on `D = 1/ρ` the merge heights span many decades and every fixed or relative cut saturates at N or jumps non-monotonically at the coarse end; that is documented in `heat_multiscale.py` so nobody re-derives it.

---

## 3. Five phases, not four: encoding vs inference across scales

*(section filled from `w0c_07_encinf_verdict.py`; see §4)*

---

## 4. Limitations, and what this lane did not settle

*(filled at close)*

---

## 5. The inheritance rule

*(filled at close)*
