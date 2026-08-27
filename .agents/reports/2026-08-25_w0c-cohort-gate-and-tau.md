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

**The gate is locked and calibrated; the multiscale claim is not.** Recomputed fresh on the incumbent substrate, the β trace is solid — cluster p = 0.0032, q = 0.019 over a six-band family, holding across the entire scale range and surviving every single-patient drop (worst leave-one-out p = 0.0091). Everything that claimed the *scale axis itself* buys something failed a direct test. **Detection is closed negative**: the hierarchy sees nothing raw FC does not, at any of 28 scales, under both a rank-linear and a strictly stronger binned residualizer (0 of 168 cells clear), and this is not a power artifact because residualizing on raw leaves 75–90% of the cophenetic rank variance intact. **Characterization is negative**: the scale profile carries no band information a cross-patient classifier can extract (0.200 vs 0.167 chance, p = 0.31), and α is not more scale-tuned than β per patient (|ρ(margin, log s)| = 0.259 vs 0.262, p = 0.90) — the "β scale-invariant vs α scale-tuned" dissociation does not reproduce, and the Friedman it rested on cannot tell invariance from cross-patient inconsistency. **Selection is a real pattern but is not demonstrated as an effect**: raw FC fires in δ/α/β/low_γ while the hierarchy keeps only β, but that is the difference between a significant and a non-significant result, and when the raw-minus-hierarchy gap is tested directly the interaction against β clears for none of δ/θ/α/low_γ. The edge-locality mechanism fails as well — every band has Gini 0.72–0.75 — making it the second proposed mechanism to fail after scale-coherence. So **what survives is a trace, not a multiscale trace**, and it is equally present at s = 0.05 where the diffusion is degenerate with raw FC. α is demoted throughout: it clears 1 of 28 scales under whole-grid BH, fails the six-band family at q = 0.090, and breaks under 3 of 10 patient drops.

Two methodological findings matter as much as the science. Per-patient observed values co-vary with each patient's **own null height** at ρ = +0.70 (β), and the margin subtraction only halves it (+0.51) — the confound is real, larger than recorded, and only partly fixed. And a 28-point scale sweep is worth **1.2–1.9 independent tests** (mean cross-scale correlation +0.70 to +0.93), so whole-grid BH over 168 cells overcharges by roughly twenty-fold and leaves no headroom at n = 10: under it, dropping any of eight patients wipes out every cleared cell, purely because the exact signed-rank floor coarsens from 1/1024 to 1/512. Choosing the multiplicity family deliberately is not a detail here — it is the difference between "β holds robustly" and "nothing holds".

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

This is not evidence that the β trace rides on one patient. It is the multiplicity design failing: the exact one-sided signed-rank p at n = 10 has a floor of 1/1024 and its next values are 2/1024, 3/1024, …; at n = 9 the grid coarsens to 1/512, and a 168-fold BH correction turns that single step into a cliff that the entire grid falls off at once. It is the sharpest possible demonstration that **whole-grid BH at n = 10 leaves the gate no dynamic range at all.**

Repeating the same LOO on the per-band axis-cluster gate, whose permutation p is continuous and whose family is 6, settles it:

| band | full p | LOO p range | worst drop | drops pushing p ≥ 0.05 |
|---|---|---|---|---|
| delta | 0.045 | [0.040, 0.098] | Pat_02 | 8/10 |
| theta | 1.000 | [0.289, 1.000] | Pat_02 | 10/10 |
| alpha | 0.030 | [0.010, 0.078] | Pat_07 | 3/10 |
| **beta** | **0.0032** | **[0.0038, 0.0091]** | Pat_06 | **0/10** |
| low_gamma | 0.114 | [0.110, 0.275] | Pat_05 | 10/10 |
| high_gamma | 0.116 | [0.030, 0.184] | Pat_08 | 9/10 |

**β survives every single-patient drop** (worst LOO p = 0.0091). α does not — three drops push it past 0.05, and it already fails the 6-band BH at q = 0.090. So the LOO collapse in the paragraph above was entirely a multiplicity artifact, and under the family that matches the question the picture is clean: β robust, α fragile, everything else null.

### 1.6 Calibration

Two tests, because a Gaussian toy tests almost nothing that matters.

**Toy null** (independent Gaussian obs and surrogates, each patient with its own mean and scale — the situation the margin exists for): FPR at nominal 0.10 / 0.05 / 0.01 = 0.094 / 0.042 / 0.011 (heteroscedastic) and 0.104 / 0.044 / 0.010 (homoscedastic), 4000 draws each. Slightly conservative at 0.05, exactly as the discreteness of the exact test predicts.

**Held-out-realization null on the real graphs** — the one that counts. One surrogate realization is promoted to the role of "observed" and tested against the remaining 199 of the *same patients*, so the test sees the real K, the real R, the real per-patient heteroscedasticity and the real surrogate ensemble. Under the null the true observation is exchangeable with its own surrogates, so these p-values are draws from the gate's null distribution on the actual data. Over all 168 cells × 200 draws: FPR at nominal 0.10 / 0.05 / 0.01 = **0.107 / 0.048 / 0.012**, with per-cell ranges 0.060–0.165, 0.015–0.095 and 0.000–0.035. The worst single cell sits at 0.095 at the 0.05 level and **no cell of 168 exceeds 0.10**. The gate is calibrated on the graphs it is actually used on.

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

The cross-phase object is five graphs — A, B (rest_pre halves), `task_learn` (ordered premises presented), `task_test` (novel non-adjacent pairs, answerable only by inference), `rest_post`. Four functionals, every one at every scale, `w0c_06`/`w0c_07`. Separating learning a structure from applying it is what distinguishes this from a generic task effect, so it gets its own grid rather than being read off the four-phase one.

### 3.1 Is the conditional functional reportable? Two checks, both passed

`T_infspec_pe` is a partial correlation, and this project has a documented case of a conditional estimator returning a significantly positive value from a no-signal placebo. Two checks ran before any p-value was quoted.

**Held-out-realization calibration**, per cell, with `refuse_uncalibrated=True` armed: **0 of 168 cells withheld** for any of the four functionals; median per-cell FPR at the 0.05 level 0.033–0.045 across functionals. The gate is calibrated against matched-strength for the conditional statistic exactly as it is for the ordinary one.

**Role-permutation placebo** — the five *observed* graphs with the five role labels permuted (all 119 non-identity permutations), so the estimator sees real cophenetic geometries in an arbitrary arrangement. If a functional manufactures sign from its own structure, this is where it shows.

| functional | placebo median over bands and scales [min, max] | cells where the placebo is itself significantly > 0 |
|---|---|---|
| T_test | −0.035 [−0.082, −0.010] | 0/168 |
| T_learn | −0.036 [−0.082, −0.009] | 0/168 |
| T_infspec | +0.000 [−0.005, +0.004] | 2/168 |
| **T_infspec_pe** | **−0.019** [−0.038, −0.001] | **0/168** |

All four sit at or slightly below zero. **The conditional functional shows no positive structural bias on this data** — if anything a small negative offset, which makes its positive observations conservative rather than inflated.

The caveat is real and is not closed by this: role permutation treats the five phases as exchangeable, which they are not, and neither check can detect a bias the matched-strength null itself shares. A **data-based placebo** — a no-signal arc built from rest_pre sub-segments — is still required, and is W0-B's. Until it exists, `T_infspec_pe` is provisional.

### 3.2 What persists, under the per-band family

Whole-grid BH over 168 cells is the overcharging family here for the same reason as in §1.4, so both are reported. Under whole-grid BH the picture is nearly empty: `T_test` β 23/28 and α 1/28; `T_learn` **0/28 in every band** (best q = 0.062, β at s = 4.73); `T_infspec` 0/28 everywhere; `T_infspec_pe` δ 3/28, α 1/28, β 0/28. Under the axis-cluster gate (one test per band per functional, BH over the 6 bands):

| functional | δ | θ | α | β | low_γ | high_γ |
|---|---|---|---|---|---|---|
| T_test (whole task) | 0.084 | 1.000 | 0.084 | **0.017** | 0.139 | 0.139 |
| T_learn (encoding) | 0.124 | 0.343 | 0.100 | **0.037** | 0.124 | 0.100 |
| T_infspec (uncond.) | 0.403 | 0.932 | 0.932 | 0.064 | 1.000 | 1.000 |
| T_infspec_pe (inference \| encoding) | **0.043** | 0.599 | **0.043** | **0.043** | 0.607 | 0.212 |

(cluster q; bold = q < 0.05)

Three readings, all recomputed fresh and none taken from prior text.

- **Encoding persists, in β only** (cluster p = 0.0062, q = 0.037, supra-threshold across the whole span s = 0.05–180). It does not survive whole-grid BH, which is why it reads as absent under that family.
- **Inference-specific persistence requires conditioning on encoding.** Unconditioned `T_infspec` clears nowhere (β q = 0.064 is the closest). Conditioned on encoding it clears in **δ, α and β** (q = 0.043 for all three). So the "inference-specific in β *alone*" claim does not reproduce, and the "δ/α/β, not β alone" correction does — now derived under a proper family with the conditional estimator's bias explicitly bounded.
- The conditional functional carries more scale structure than the whole-task one: its effective number of independent scales is 2.3–3.2, against 1.2–1.9 for `T_test`. That is the one place in this lane where the scale axis is doing measurably more work.

### 3.3 The scale dissociation — suggestive, not established

Do encoding and inference persist at *different* scales? Tested on the per-patient profile centroid in log s, paired within patient, so an amplitude difference between the two curves cannot produce it.

| band | T_learn centroid s | T_infspec_pe centroid s | p(learn < inf) | p(learn > inf) | mean shift (log s) [95% CI] | shift in s | 80%-power MDE | n |
|---|---|---|---|---|---|---|---|---|
| delta | 9.43 | 3.38 | 0.615 | 0.423 | −0.240 [−1.484, +0.965] | ×0.79 | 1.79 | 10 |
| theta | 3.10 | 2.20 | 0.722 | 0.313 | −0.621 [−2.046, +0.681] | ×0.54 | 2.00 | 10 |
| **alpha** | **2.39** | **4.04** | **0.027** | 0.981 | **+0.996 [+0.176, +1.889]** | **×2.71** | 1.18 | 9 |
| **beta** | **6.09** | **2.08** | 0.976 | **0.032** | **−1.204 [−2.271, −0.250]** | **×0.30** | 1.46 | 10 |
| low_gamma | 4.67 | 2.87 | 0.633 | 0.410 | −0.249 [−1.912, +1.560] | ×0.78 | 2.40 | 9 |
| high_gamma | 2.18 | 2.93 | 0.191 | 0.844 | +1.010 [−0.808, +3.022] | ×2.75 | 3.01 | 8 |

Nominally there is a dissociation in two bands, and it is a real dissociation of *position*, not amplitude: in α the inference-specific trace sits at a scale **2.7× coarser** than the encoding trace (p = 0.027), and in β it sits **3.3× finer** (p = 0.032). Both bootstrap CIs exclude zero.

Three reasons it must be reported as suggestive and not as a result:

1. **It does not survive correction.** BH across the directional tests gives min q = 0.193. Nothing clears.
2. **The two bands disagree in direction.** A single mechanism — "inference is coarser than encoding" — predicts one sign. Getting +2.71× in α and ×0.30 in β is what a pair of underpowered draws looks like as easily as what a band-specific mechanism looks like.
3. **Both observed shifts sit at or below the design's detection threshold.** The smallest paired shift this n and this spread would catch at 80% power is 1.18 log units for α (observed 0.996) and 1.46 for β (observed 1.204). Both nominal significances are therefore lucky draws relative to the power available, which is exactly the regime that does not replicate.

**Verdict: the encoding-vs-inference scale dissociation is the most interesting thing in this lane and is not established at n = 10.** It is worth a pre-registered replication with a directional hypothesis fixed in advance — that would cut the family from 12 to 1 and is the only route by which this reaches significance without more patients. It should not be written into the paper as a finding in its present state.

---

## 4. Limitations, and what this lane did not settle

Stated up front rather than at the end, per project rule.

**The null is still injected at the last stage.** Every verdict here is conditioned on "given this FC matrix". Matched-strength shuffles the finished N×N connectivity and therefore cannot test the coherency estimator, the band transform, the frequency split, session nonstationarity, drift or artifact epochs. That is lane W0-B's null ladder and it is a live dependency: **nothing in this report is verified against a timeseries-level null.**

**No data-based placebo was run.** The role-permutation placebo in §3 bounds the conditional estimator's structural bias using the real graphs, and the held-out-realization test calibrates the gate against the null in use. Neither can detect a bias that the null itself shares. The documented failure mode — a sham arc inside pre-task rest returning a significantly positive conditional trace — needs a no-signal arc built from the recording, which requires recomputing FC on rest_pre sub-segments and belongs with W0-B. Until it exists, `T_infspec_pe` should be treated as provisional whatever its p-value.

**Substrate sensitivity, explicitly.** Every number here is on `imcoh_abs` × `mst_union_top_fraction @ 0.20`. W0-A reports `imcoh_abs` and `imcoh_sq` are nearly but not exactly rank-equivalent (Spearman ≈ 0.986–0.990; backbone Jaccard ≈ 0.83–0.88 at f = 0.20), so the substrate may shift. My assessment of what would and would not move:

- **Robust to the substrate.** The three τ negatives. Detection fails by a wide margin (q_min 0.14 against a 0.05 line, 0/168 in both residualizers); characterization fails at chance; the selection interaction fails with p-values of 0.12–0.62. A backbone with 85% edge overlap will not turn any of those into positives. The methodological findings — the obs/null covariation, n_eff ≈ 1.2–1.9, the BH-at-n=10 headroom problem — are properties of the design, not the substrate.
- **Substrate-sensitive.** The exact per-scale q values, β's 23/28 count, α's single cleared scale, and the cluster p-values (β 0.0032, α 0.030). α in particular sits close enough to the line that a substrate change could move it either way. Any claim resting on α must be re-derived on the locked contract.
- **Unknown.** The `N_eff(s)` and `ℓ(s)` scale units are backbone-dependent by construction and must be recomputed for whatever backbone is locked. They are cheap (`w0c_04`, ~2 s).

**Power.** n = 10 with an exact signed-rank has a p-floor of 1/1024 and a coarse grid above it. The gap and interaction tests in §2.3 are differences of two noisy margins and are intrinsically low-powered; "not demonstrated" there must not be read as "demonstrated absent". Equivalence claims (scale-invariance) cannot be made at this n at all — the Friedman is blind to a monotone trend below about 0.75 within-patient SD.

**What I did not do.** No z-standardised variant of the margin was evaluated, though §1.2 shows it is the obvious next lever. The role-permutation null treats the five phases as exchangeable, which they are not. The `N_eff` estimators disagree by an order of magnitude and are reported as a bracket rather than adjudicated.

---

## 5. The inheritance rule — what a cohort claim in this paper must show

> A cohort claim is reportable only if it is a **margin** claim — one-sided Wilcoxon signed-rank on `obs − (that patient's own surrogate median)` across the **full n = 10**, with any exclusion argued at the point of use and never inferred from a row count — evaluated **at every position of every swept axis and reported at all of them**, never at a selected best; corrected within a **deliberately chosen multiplicity family**, stating which family and why, given that a 28-point scale sweep is worth 1.2–1.9 independent tests so whole-grid BH over (band × scale) overcharges roughly twenty-fold while a per-band sign-flip cluster test does not; accompanied by a **bootstrap CI and rank-biserial effect size**, by the **leave-one-patient-out swing in the verdict** (the whole grid and its correction re-run per drop, not the statistic re-computed), and by the **measured false-positive rate of that statistic against that null** from held-out surrogate realizations — with the p-value **withheld, not caveated**, if the pair is not calibrated, which is mandatory for any conditional or partial statistic. Per-patient counts may be reported and are never the gate. Any claim of the form "representation A sees what representation B does not" must test the **difference** directly, paired within patient, because the difference between a significant and a non-significant result is not itself significant. And no scale may be described in words without `N_eff(s)` or `ℓ(s)` attached.
