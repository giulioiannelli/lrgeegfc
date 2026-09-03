---
name: 2026-08-31_lane-s-scale-variability
kind: report
era: PAPER_FINALIZATION (Wave 0, lane W0-S)
status: current
created: 2026-08-31
updated: 2026-09-03
scope: Does the cross-phase trace vary across diffusion scales, or is the incumbent readout blind to scale? Pre-registered before any number, run on the locked substrate knob-integrated over the plateau f in [0.07, 0.20], every readout referenced to a matched-strength null and gated by the locked cohort gate. Verdict is NEGATIVE and that is the headline: the effect is the same at every scale, the readout is not the reason, and the mechanism is measured rather than asserted.
pointers:
  - .agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md
  - scripts/01_compute/paper_final/w0s_01_scale_locality_grid.py
  - scripts/01_compute/paper_final/w0s_06_quantization_ceiling.py
  - data/paper_final/lane_s_scale/
---

# W0-S — is the scale axis flat, or is the readout blind to it?

## Head

**The scale axis is genuinely flat, the readout is not the reason, and the project should drop the multiscale framing.** Knob-integrated over the locked plateau, the trace clears the per-band axis-cluster gate in delta, alpha and beta (cluster q = 0.034 each) with the supra-threshold cluster spanning the **entire** axis, s = 0.02 to 180, and clearing separately in both halves of it; theta is null (p = 1.000) exactly as the known-null band should be. The 32-point sweep is worth **1.06-1.26 independent tests**, while the same statistic's own held-out-surrogate noise floor would give **2.9-3.4** — so the axis is *more* redundant than chance, in every band, at p = 1.000. Nothing survives correction for scale variation (slope min q = 0.078; cross-patient shape agreement min q = 0.088). And the suspect readout is exonerated: **none of five scale-local constructions passes the pre-registered criteria**, every one of them sitting *below* its own noise floor on the very statistic that was supposed to show it carrying more information.

**Two limitations first, per project rule.** Everything here is conditioned on "given this FC matrix" — matched-strength shuffles the finished connectivity, so nothing in this lane is verified against a timeseries-level null, and W0-B's ladder remains a live dependency. And the surrogate ensemble is R = 50, reduced from a pre-registered 100 during a period of machine contention and stated rather than tuned; it is adequate for margin medians, for calibration at tolerance 0.10, and for a null-referenced p at 1/51 resolution, but a borderline calibration number would deserve a deeper ensemble.

**The mechanism, measured.** Moving four decades along the scale axis perturbs the readout's input **about as much as splitting one rest recording in half does** — the cophenetic vector at s = 0.02 versus s = 180 is Spearman 0.46-0.59 similar, while the same vector at the same scale on the two halves of pre-task rest is 0.44-0.59 similar. Underneath that, a UPGMA tree's merge groups partition ~6 900 contact pairs into 117 cells, and that coarse-graining's capacity is **scale-invariant** (its measured ceiling on the raw per-pair reorganisation moves only 0.110 to 0.130 across the whole axis) even though the tree's granularity falls 30-fold over the same range. This also closes an unrelated open question: W0-C's unexplained "the hierarchy expresses only 1-6% of raw reorganisation" **is** that coarse-graining, and needs no filtering or orthogonality story.

---

## 1. What was pre-registered, and what the criteria were for

The scope report — written and committed before any number, and not edited since — is `.agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md`. It fixed the notation, two stratifications of the pair set, three candidate readouts, and three criteria a readout had to satisfy to replace the incumbent: **(a)** reproduce the incumbent when aggregated over scales, so it is the same phenomenon and not a new one; **(b)** carry more independent cross-scale information *than its own noise produces*; **(c)** be calibrated.

Criterion (b)'s second clause is the reason this lane did not reach a false positive, and it is worth stating plainly because it is a trap any repeat of this work will hit. The effective-number-of-independent-tests statistic **rises automatically with noise** — independent noise decorrelates the columns of the per-patient x per-scale margin matrix. A criterion built on n_eff alone is therefore satisfied by making the statistic *worse*. It has to be referenced to held-out surrogate realizations of the same readout, which carry the identical noise level and no signal. §4 shows exactly how much that mattered.

Nothing in this lane was ever tested against zero. W0-B established that all four cross-phase functionals are significantly positive at 16/16 scales on a no-task sham arc with temporal order destroyed, so zero is not the reference; every number below is a margin against a matched-strength surrogate that went through the identical sparsify -> Laplacian -> heat kernel -> UPGMA -> stratify -> rank-correlate pipeline.

**Equivalence was verified on real FC before anything was gated**, rather than asserted: the graphs are bit-identical to `canonical_graph_ensemble` at every plateau fraction; the lane's `T` equals `rho_sym_over_scales` and `cross_phase_functionals`' `T_probe` to 2.2e-16; and both stratified decompositions sum back to `T` to 2.2e-16 (`w0s_verify_readouts.py`).

---

## 2. Two sub-questions closed before the main result

**The grid-sampling explanation is falsified.** The premise that the informative regime is undersampled is true only of the *locked contract grid*, which has 16 points starting at s = 1 — four below s = 3 and **none** below s = 1. The 28-point sweep whose n_eff was actually measured at 1.2-1.9 already had **14 of its 28 points below s = 3** and 10 below s = 1. The flatness was never a consequence of skipping the interesting region. This lane's own grid puts 18 of 32 points below s = 3 and 14 below s = 1, and finds the same flatness there (§3).

**The scale axis has two degenerate ends, and they are not symmetric.** Below s ~ 0.3 the kernel is `I + tau*W` to leading order, so the communication distance is a monotone transform of the raw FC and the tree is effectively UPGMA on raw connectivity — the fine end is the raw-FC limit *by construction*. Above s ~ 2 the communication reach saturates near the implant span. The genuinely "hierarchical" window is therefore narrow, and this is a property of the operator, not a sampling choice.

---

## 3. Is the trace scale-dependent? No — it is present at every scale, equally

The primary family is the per-band sign-flip axis-cluster test, with BH over the 6 bands, because a 32-point sweep is not 32 tests (§4) and whole-grid BH over 192 cells overcharges by roughly an order of magnitude. Knob-integrated over `f in {0.07, 0.10, 0.14, 0.20}`, n = 10, no patient dropped:

| band | margin (cohort median) | cluster p | cluster q | supra-threshold span | fine half s < 3, q | coarse half s >= 3, q |
|---|---|---|---|---|---|---|
| delta | +0.089 | 0.014 | **0.034** | s = 0.02-180 | 0.053 | 0.055 |
| theta | -0.035 | 1.000 | 1.000 | none | 1.000 | 1.000 |
| alpha | +0.122 | 0.013 | **0.034** | s = 0.02-180 | 0.053 | **0.043** |
| beta | +0.092 | 0.017 | **0.034** | s = 0.02-180 | 0.053 | **0.043** |
| low_gamma | +0.095 | 0.042 | 0.063 | s = 0.02-180 | 0.084 | 0.055 |
| high_gamma | +0.029 | 0.184 | 0.221 | s = 0.02-0.28 | 0.150 | 0.187 |

Three readings. The trace is real and it clears in **delta, alpha and beta** under knob integration — a different band pattern from W0-C's single-fraction f = 0.20 result, and the direction the W0-A contract predicts, since the plateau window is where alpha and beta are both positive. Theta is null at every scale and under every family, so the gate is not firing on the band the project expects to be silent. And the supra-threshold cluster spans the **whole axis** in all three clearing bands: there is no sub-range of scale where the effect lives.

Under whole-grid BH over the 192 (band x scale) cells, **no band clears a single cell** (q_min = 0.088). That is the multiplicity family failing, not the effect vanishing — the same phenomenon W0-C documented — and both families are reported here rather than the favourable one being chosen.

**The margin is not literally constant, and this is the one place a reader could be misled.** The beta cohort margin rises about 6-fold along the axis, from +0.035 at s = 0.02 to +0.221 at s = 134, and the per-patient monotone trend is consistent in sign (mean within-patient Spearman of margin against log s = **+0.489**). But the *standardized* effect does not move — the per-scale p stays in 0.019-0.12 throughout — because the surrogate spread grows with the margin. So the effect size in margin units grows with scale while its detectability does not, and no cohort-level test separates the profile from flat (§4).

---

## 4. The scale axis carries less independent information than noise would

`effective_tests` on the per-patient margin matrix, exactly as W0-C measured it, with the null being held-out surrogate realizations of the **same** readout:

| band | n_eff observed | n_eff, own noise floor | p (observed above null) |
|---|---|---|---|
| delta | 1.18 | 2.88 | 1.000 |
| theta | 1.08 | 3.11 | 1.000 |
| alpha | 1.06 | 3.07 | 1.000 |
| beta | 1.22 | 3.18 | 1.000 |
| low_gamma | 1.26 | 3.44 | 1.000 |
| high_gamma | 1.20 | 3.41 | 1.000 |

This reproduces W0-C's 1.2-1.9 independently and sharpens it: the observed value is not merely low, it is **below what pure noise gives at this cohort size and this axis length**, in every band, without exception. The scale axis is delivering one thing repeated, and the redundancy is far too strong to be a sampling accident.

Two direct tests of variation, both null-referenced and both corrected over the 6-band family:

| statistic | best band | observed | null median | p | min q over 6 bands |
|---|---|---|---|---|---|
| within-patient monotone scale dependence, mean abs Spearman(margin, log s) | beta | 0.531 | 0.398 | 0.020 | 0.078 |
| cross-patient agreement on profile **shape** (mean pairwise r of z-scored profiles) | beta | 0.200 | -0.017 | 0.020 | 0.088 |

The shape statistic is the one that matters, because it separates the two readings a flat cohort curve admits — genuine within-patient invariance, versus patients each having a differently-shaped profile that averages flat — which a Friedman cannot do (W0-C §2.2). Beta's +0.200 against a -0.017 null says patients do *partially* agree on where along the axis the effect is larger. **It does not survive correction (q = 0.088), and it is reported as a directional hint, not a result.** Together with the +0.489 slope it is the single thing in this lane worth a pre-registered directional replication; it is not something to write into the paper.

---

## 5. Every construction tried, including the failures

Five scale-local readouts, all pre-registered, all reported:

| readout | what it does | rho vs incumbent (beta) | (a) same phenomenon | n_eff obs / own null | (b) more info | (c) calibrated | passes |
|---|---|---|---|---|---|---|---|
| `T` incumbent | Spearman over all ~6 900 pairs | 1.000 | yes | 1.19 / 3.15 | no | yes (0.995) | — |
| `Qcon_diag` | additive contribution of the equal-count stratum the diffusion is resolving | 0.891 | **yes** | 1.31 / 3.70 | no | yes (1.000) | **no** |
| `Ccon_diag` | same, octave strata | 0.758 | no | 1.66 / 4.07 | no | yes (0.984) | **no** |
| `Tloc_diag` | within-octave trace, re-ranked inside the stratum | 0.624 | no | 1.76 / 4.28 | no | yes (0.995) | **no** |
| `Qloc_diag` | within-quintile trace, re-ranked | 0.661 | no | 1.51 / 4.10 | no | yes (0.990) | **no** |
| `Thei` | trace on the log merge-height profile — not pairwise at all | 0.236 | no | 1.85 / 2.53 | no | yes (1.000) | **no** |

**This table is the reason criterion (b) was written with a null reference.** Read the absolute n_eff column alone and four of the five candidates look like an improvement on the incumbent's 1.19 — `Thei` at 1.85 is 55% higher, `Tloc_diag` at 1.76 is 48% higher. Every one of them is *below its own noise floor*. The apparent gain is entirely the extra variance of a statistic computed on fewer pairs, which is exactly what A3.1 of the pre-registration predicted and exactly what would have been reported as "scale structure the incumbent hides" had the criterion been written the obvious way.

Calibration is not the problem: held-out-realization FPR at nominal 0.05 is within tolerance in 98.4-100% of cells for every readout, so no p-value here is being withheld as uninterpretable.

Only `Qcon_diag` also satisfies criterion (a) — it is the same phenomenon, which is expected since it is an exact additive part of the incumbent. `Thei` fails (a) badly at rho = 0.236: it measures something genuinely different, and that difference is not a trace.

---

## 6. Why the readout cannot see scale, measured rather than asserted

**The input barely moves relative to its own measurement noise.** The fair yardstick for "how much does four decades of scale change the cophenetic vector" is how much the *same* vector changes between two contiguous halves of one pre-task rest recording — same scale, same object, only measurement differing:

| band | cophenetic vector, s = 0.02 vs s = 180 | same vector, split-half of one rest recording |
|---|---|---|
| delta | 0.565 | 0.492 |
| theta | 0.532 | 0.515 |
| alpha | 0.462 | 0.493 |
| beta | 0.587 | 0.593 |
| low_gamma | 0.519 | 0.569 |
| high_gamma | 0.535 | 0.443 |

They are the same magnitude. Traversing the entire scale axis perturbs the readout's input no more than splitting one rest recording in half does. That is the cleanest single statement of why the axis is flat, and it needs no appeal to the estimator at all.

**The representation's capacity is scale-invariant even though the tree's granularity is not.** A UPGMA tree on N leaves has exactly K = N-1 merges, and the merge groups **partition** the M = N(N-1)/2 pairs, so the cophenetic vector lies in an orthogonal coordinate subspace of dimension K inside R^M and can carry only the between-group component of any per-pair signal. Measured on the cohort (N = 118, M = 6 903, K = 117):

| quantity | value |
|---|---|
| chance ceiling for a tree-unrelated signal, (K-1)/(M-1) | 0.0168 — empirically 0.0167-0.0168 |
| **measured ceiling on the raw per-pair reorganisation** | **0.092-0.139** |
| observed expressible share, recomputed fresh | 0.010-0.032 |
| the same quantity as reported by W0-C | 0.013-0.062 |

So **W0-C's unexplained "the hierarchy expresses only 1-6% of raw reorganisation" is the coarse-graining, not a filter and not orthogonality**: a 117-cell summary of a 6 903-cell object retains 1.7% of an arbitrary signal, and at most 9-14% even of this particular, tree-aligned one. Two separately-reported project numbers become one.

And the capacity does not move with scale. Across s = 0.02 to 134 the measured ceiling goes 0.110 -> 0.130 and the pair-mass-weighted effective number of levels goes 15.8 -> 21.7, while `N_eff` — the tree's actual granularity — falls from 118 to about 5. **The hierarchy's resolution changes 30-fold; the readout's information capacity does not change at all.**

**One correction to the hypothesis I was handed.** The mechanism is *not* rank-tie saturation. All 117 merge heights are distinct, and the exact tie ceiling `sqrt(1 - sum(n^3-n)/(M^3-M))` on the cophenetic vector is **0.986-0.993** — for the reorganisation vector it is 0.9996. Ties leave essentially all the rank correlation available. The binding constraint is dimensional, not tie-based, and the distinction matters because it predicts the opposite remedy: reweighting the quantised values cannot help, only escaping the K-cell partition can.

**That prediction was made before the candidates were scored, and it held.** Because a stratum is *defined* by merge level, restricting to one leaves only ~6 distinct cophenetic values in place of 117, and its split-half reliability collapses from 0.44-0.59 to **0.10-0.28**. The band-pass makes the quantisation about twenty times worse, which is why `Tloc`/`Qloc` were always going to fail.

**A flaw in my own diagnostics, reported rather than buried.** I intended to report a split-half reliability for the merge-height readout too. The proxy I computed is degenerate: Spearman between two *sorted* height vectors is +1 by construction, so the 1.000 in `carrier_overlap_summary.csv` carries no information and must not be cited. `Thei` itself is unaffected — it correlates *differences* of log heights, not sorted heights — but its reliability is simply unmeasured here.

**The carriers do turn over, so the flatness is not a trivial identity.** The per-pair contribution vectors at the two ends of the axis overlap at cosine 0.19-0.30, against a matched-strength reference of 0.19-0.36 — i.e. at the extremes, essentially different pairs are carrying the trace. The mean off-diagonal overlap is 0.51-0.59 against a null of 0.19-0.36, so there is real shared structure, but it is far from the same pairs everywhere. **The scale axis does change which contacts build the trace; what it does not change is anything the cohort statistic can see.** That is a more precise and more interesting negative than "every scale is the same measurement".

---

## 7. Answer to the lane's question

**The effect is the same at every scale.** The trace is present across the whole axis in delta, alpha and beta and in both halves of it separately; the sweep carries fewer independent tests than its own noise floor in all six bands; no test of scale variation survives correction; and no scale-local readout is a better instrument than the incumbent by the criteria fixed in advance. The multiscale framing is not supported, and the honest description of the finding is a **trace that is scale-pervasive rather than multiscale** — it exists at every diffusion scale from the raw-FC limit to near-global, with no scale doing distinguishable work.

## 8. What I could not settle

**Whether beta's profile really has a reproducible shape.** Slope +0.489 and shape agreement +0.200 against a -0.017 null are both nominally significant and both die under a 6-band correction (q = 0.078, 0.088). At n = 10 this design cannot resolve it; a pre-registered directional test in one band is the only route that does not need more patients.

**Whether the flatness survives a timeseries-level null.** Every verdict is conditioned on the finished FC matrix. W0-B's ladder is the live dependency and nothing here is verified against it.

**Equivalence is not established, only failure to reject.** "No scale structure detectable at n = 10" is not "scale-invariant". Establishing invariance needs an equivalence design with a pre-declared margin, which this cohort size does not support — the same limitation W0-C recorded.

**Surrogate depth.** R = 50, not the pre-registered 100, for machine-contention reasons stated at the time. Adequate for everything used here; a borderline calibration verdict would deserve more.

**The merge-height readout is under-characterised.** `Thei` fails criterion (a) at rho = 0.236 and so was not pursued, but it is the only construction that escapes the quantisation, and its reliability was not validly measured (§6). If anyone revisits the scale question, that is where to look — not at another band-pass in tree level.

---

## 9. Reproduce

```
./scripts/01_compute/paper_final/run_lane_s.sh          # everything, in order
```

Detached runs are mandatory for the grid (`setsid ... < /dev/null &`), which is resumable per cell; the verdict stage reuses its stage-1 gate grid and checkpoints stage 2 per readout, so an interruption costs one unit. Injection points: `W0S_R` (surrogate depth), `W0S_NPERM` (sign-flip permutations, default 2 000 — p-resolution 5e-4, far finer than a 6-band family needs), `W0S_WORKERS`, `W0S_OUT`.
