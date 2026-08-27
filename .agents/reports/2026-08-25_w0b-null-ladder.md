---
name: 2026-08-25_w0b-null-ladder
kind: report
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 0, lane W0-B)
status: current
created: 2026-08-25
scope: The timeseries-level null ladder for the cross-phase trace. Builds and runs the nulls that act UPSTREAM of the finished FC matrix — where matched-strength cannot reach — and reports, per band and per scale, which rungs the α/β trace survives and which it does not, plus an explicit statement per rung of what that rung cannot reject. Includes the null-calibration of the five-phase encoding/inference functionals.
pointers:
  - src/lrg_eegfc/utils/surrogate/timeseries_nulls.py          # the library deliverable
  - src/lrg_eegfc/utils/fc/heat_multiscale.py                  # cross_phase_functionals, partial_spearman
  - scripts/01_compute/null_ladder/                            # runners 00–08
  - data/paper_final/w0b_nulls/                                # all computed artifacts
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py   # the incumbent gate reproduced exactly
  - .agents/preprint/supplementary/S2_drift_controls.md        # the broken-null precedent this lane re-tests
---

# W0-B — the null ladder

## Head

**The β and α traces do not survive the one null that preserves coherence magnitude and destroys only lag, and they do not survive it under the control that was designed to rescue them either.** Under N1 (per-channel circular shift, realized as its exact frequency-domain phase ramp) β clears **0 of 16 scales**, best cohort p = 0.080; α clears **0 of 16**, best p = 0.138. Holding the backbone edge set fixed at the observed one — which removes the objection that the surrogate is being scored on a different graph — the negative becomes total: β and α fail at **0 of 16 scales** with cohort median margins of **−0.003 and +0.004**, against +0.214 and +0.130 under matched-strength. That is not a weakened effect; it is no separation whatsoever between the observed trace and a surrogate built from the same graph with the lag structure destroyed. They survive everything else on the ladder: the incumbent matched-strength null (β 16/16, α 12/16) and the segment-lattice independence null N1b (β 16/16, α 13/16) — so the trace is not noise, and genuine cross-channel coupling is required to produce it. What fails is the *attribution*: a surrogate that keeps each phase's own coherence-magnitude structure `|C_ij(f)|` intact and randomizes only the phase reproduces the trace as well as the observed `|ImCoh|` does. The effect size is not what collapses — β's median cohort margin under N1 is +0.174, against +0.214 under matched-strength — the **null floor rises and becomes patient-specific**, by up to an order of magnitude (Pat_06 α: observed +0.886, N1 null +0.865). Stated precisely, because the distinction matters: this is **not** a finding that the measurement is contaminated by volume conduction. `|ImCoh|` remains volume-conduction-immune by construction — zero-phase mixing contributes exactly zero, and that is a property of the estimator that no null can take away. What fails is narrower and still serious: the *discriminative content* of the cross-phase trace is **not lag-specific**. A surrogate that discards the observed lag structure entirely, keeping only each phase's coherence-magnitude geometry, reproduces the trace at cohort level. So the result cannot be presented as evidence about time-lagged interaction, and the choice of `|ImCoh|` over ordinary coherence cannot be justified by pointing at this result — the result does not distinguish them.

---

## 1. What was built

`src/lrg_eegfc/utils/surrogate/timeseries_nulls.py` — surrogates that act on the per-segment Fourier coefficients or on the complex band coherency `C(f)`, i.e. **upstream** of the `N × N` adjacency that matched-strength shuffles. No matrix-level rotation null is rebuilt (the archived CRC/SB-CRC family stays archived).

| rung | construction | preserved exactly | destroyed |
|---|---|---|---|
| **N1** | per-channel circular shift as a frequency-domain phase ramp | `\|C_ij(f)\|`, every per-channel spectrum | cross-spectral phase = **lag** |
| **N2** | iid uniform phase per (channel, frequency) | `\|C_ij(f)\|`, per-channel spectra | lag + cross-frequency phase consistency |
| **N1b** | circular shift of each channel's Welch **segment index** | each channel's Welch PSD | cross-channel temporal correspondence (magnitude included) |
| **N3** | block phase-label permutation, free + order-preserving rotation | session content, per-phase durations | the placement of the phase boundaries |
| **N4** | within-rest placebo, duration-matched | — | the task itself |

Literature anchors: Theiler et al. (1992) for phase randomization; Perkel, Gerstein & Moore (1967) for the shift predictor that N1b generalizes; Nolte et al. (2004) and Ewald et al. (2012) for the estimator; Welch (1967) for the segmentation.

## 2. Three preflight facts that shaped the design

**The observed statistic is reproduced.** ρ_sym recomputed from raw timeseries through the surrogate code path matches the value script 13 gates **exactly (bit-identical) on 791 of 800 cells**, with a worst-case deviation of `6.0e-05` on the remaining 9 — the residue of the float32 FC cache, four orders of magnitude below any effect discussed here. The library CSD reproduces the cached `imcoh_abs` adjacency to `5.9e-08`. Every rung is therefore attached to the right observed value, which was the stated sanity requirement.

**The circular-shift identity is exact only in the single-segment DFT sense, and the brief's premise needed correcting.** For one un-windowed segment the shift theorem holds to `2e-15` and `|C_ij(f)|` is preserved to `7e-16`. Under the **Welch estimator the pipeline actually uses**, a genuine whole-recording roll does *not* preserve magnitude: mean `|C|` collapses from **0.241 to 0.036**, because rolling channel `c` makes its segment `k` pair with an entirely different epoch of channel `j`, so the cross-spectral terms average incoherently. A whole-recording circular shift under Welch is therefore an **independence null**, not a lag null. The two are separated in this ladder rather than conflated: N1 is the frequency-domain, magnitude-preserving lag null; N1b is the honest time-domain version, labelled as the independence null it is.

**N1/N2 inflate the edge weights, and that is licensed.** Randomizing the phase drives `⟨|Im C|⟩_f → (2/π)⟨|C|⟩_f`: measured ×4.6 on Pat_05 β (0.0332 → 0.1538, against the analytic prediction 0.1535). The surrogate therefore sits **above** the data in raw edge weight. This is harmless *only* because the readout is invariant to a global positive rescale of `W` — the `mst_union_top_fraction` backbone is rank-based and `s = τλ_max` makes `τL` scale-free. That was verified numerically, not assumed: `max |ρ_sym(W) − ρ_sym(6.5W)|` over the 16 scales is **0.0e+00**. The one-sided direction was verified independently by re-running the identical gate code on the incumbent matched-strength rung, which still clears (β 16/16, α 12/16); a sign error would have broken that too.

## 3. The verdict — ρ_sym, per band, per scale

Scales cleared out of 16, cohort one-sided Wilcoxon on the margin `obs − surr_p50` across the 10 patients. `high_gamma` was **not run** on the new rungs (see §8).

| band | ms (incumbent) | **N1 (lag)** | **N1 fixed-backbone** | **N2 (lag)** | N1b (independence) |
|---|---|---|---|---|---|
| delta | 8 | **0** | **0** | **0** | 8 |
| theta | 0 | **0** | **0** | **0** | 0 |
| **alpha** | 12 | **0** | **0** | **0** | 13 |
| **beta** | 16 | **0** | **0** | **0** | 16 |
| low_gamma | 0 | **0** | **0** | **0** | 0 |

Best (smallest) cohort p over the 16 scales — best-scale selection is already optimistic, and every lag rung still fails:

| band | ms | **N1** | **N1 fixed-bb** | **N2** | N1b |
|---|---|---|---|---|---|
| delta | 0.0137 | 0.7842 | 0.2158 | 0.4229 | 0.0029 |
| theta | 0.2158 | 0.5771 | 0.2783 | 0.6875 | 0.1875 |
| **alpha** | 0.0068 | **0.1377** | **0.1377** | **0.1611** | 0.0068 |
| **beta** | 0.0010 | **0.0801** | **0.2158** | **0.0654** | 0.0010 |
| low_gamma | 0.1162 | 0.2461 | 0.0801 | 0.2158 | 0.1162 |

**The fixed-backbone control removes the strongest objection to this negative.** N1 drives `W` toward `(2/π)|C|`, and the `|C|` edge ranking is not the `|ImCoh|` ranking — on sEEG, `|C|` is dominated by same-probe pairs (a 2–8× bias, CLAUDE.md invariant 5). A plain N1 surrogate therefore lives on a *different* graph from the observed, so one could argue its high floor is the cross-phase stability of that anatomical scaffold rather than anything about lag. `n1_fixedbb` pins the backbone edge set to the observed `|ImCoh|` one for surrogate and observed alike, so only the weights' lag content varies. The trace still fails, and far more decisively than under plain N1. The cohort **median margin** — the gate statistic, taken across all 16 scales — is **−0.003 for β and +0.004 for α**, against +0.214 and +0.130 under matched-strength. That is not a reduced effect; it is **no separation at all**. Once the graph is held fixed, a surrogate that has discarded every trace of the observed lag structure reproduces the α and β cross-phase traces essentially exactly.

Under BH-FDR across the common (5 bands × 16 scales) family, within rung: ms β 16/16 and N1b β 16/16 survive; N1 and N2 reach a minimum q of 0.905 and 0.819. Under leave-one-patient-out of the *verdict*, ms β 16/16 and N1b β 16/16 hold with any single patient dropped; N1 and N2 hold at 0 scales for every band.

**N1b reproduces the incumbent almost exactly** (β 16/16, α 13/16, δ 8/16 against ms 16/12/8). This is not an impression from the counts: across all 800 shared cells the two rungs' margins correlate at **Pearson r = 0.983**, with a median absolute difference of **0.0115** and null medians of `+0.0026` (ms) against `−0.0005` (N1b). Two nulls built on completely different mechanisms — one shuffling a finished matrix under a strength constraint, one destroying cross-channel temporal correspondence in the raw Welch segments — land on the same floor cell for cell. That is a genuine cross-validation of the incumbent gate: the existing β result is *not* an artifact of the matched-strength construction. It also means matched-strength has been behaving as an **independence null** all along, which is precisely why it never reached the magnitude-versus-lag question.

## 4. Why N1 fails — it is the null floor, not the effect

The effect does not disappear under N1. Cohort median margin over scales: β **+0.174** under N1 against **+0.214** under matched-strength; α +0.016 against +0.130.

What changes is the null, in two ways that compound. Its **spread** is 3–4× wider than matched-strength's (β 0.375 vs 0.096; α 0.285 vs 0.100) and 7–8× wider than N1b's (0.050, 0.045). And its **height varies across patients**, which matched-strength's and N1b's do not: the interquartile range of the per-patient null median, taken across the cohort and pooled over scales, is **0.232 (β) and 0.534 (α) under N1**, against **0.0024 and 0.0029 under N1b**. Matched-strength and N1b supply an essentially constant floor near zero for every patient; N1 supplies a floor that is a property of that patient's own coherence geometry. At β's best scale:

| patient | observed | N1 null p50 | N1 margin | MS null p50 | MS margin |
|---|---|---|---|---|---|
| Pat_02 | +0.548 | +0.225 | +0.322 | +0.156 | +0.392 |
| Pat_03 | +0.359 | +0.424 | **−0.065** | +0.037 | +0.322 |
| Pat_05 | +0.269 | −0.121 | +0.390 | +0.071 | +0.198 |
| Pat_06 | +0.227 | +0.929 | **−0.702** | +0.050 | +0.177 |
| Pat_07 | +0.177 | −0.007 | +0.184 | −0.007 | +0.185 |
| Pat_08 | +0.482 | −0.037 | +0.519 | +0.061 | +0.422 |
| Pat_10 | +0.172 | −0.011 | +0.183 | +0.002 | +0.170 |
| Pat_13 | +0.039 | +0.103 | **−0.064** | +0.013 | +0.026 |
| Pat_14 | +0.079 | −0.018 | +0.097 | +0.031 | +0.048 |
| Pat_15 | +0.015 | −0.112 | +0.128 | −0.046 | +0.061 |

7/10 positive, median +0.155, but Pat_06's null sits at +0.929 — for that patient a lag-randomized surrogate reproduces essentially the entire trace. The α table is worse: Pat_06 observed +0.886 against an N1 null of +0.865. **In a subset of patients the coherence-magnitude geometry alone produces a trace as large as the observed one.** The cohort splits, and the Wilcoxon fails. This is a real property of the data, not a harness defect: N1 is a genuine stochastic null (its spread is wider than matched-strength's, so it is not a deterministic substitution in disguise), and the same gate code clears matched-strength.

## 5. What each rung CANNOT reject

**N1 / N2** cannot reject an account in which the trace lives in cross-phase changes of coherence **magnitude** — they hold each phase's own `|C_ij(f)|` fixed by construction. They cannot reach session nonstationarity, drift, artifact epochs, or the phase segmentation, because all phases keep their true boundaries. N2 is additionally **not an independent rung**, and the measured reason is not the one theory predicts. Both send the cross-spectral phase to uniform, so they share a first moment; the expectation was that N2 would be the *low-variance* limit, because its per-bin phases average down while N1's single coherent ramp does not. In ρ_sym that difference does **not** materialize: the null spread `p95 − p50` is 0.335 (N1) against 0.322 (N2) pooled, and 0.375 against 0.418 for β. With only 7–35 in-band bins the frequency averaging is too weak to separate them, and the variance that matters for ρ_sym is structural (which edges get up-weighted) rather than in the band mean. N1 and N2 are therefore near-duplicates on this pipeline: report one, and never treat the second as corroboration of the first.

**N1b** destroys coupling magnitude along with lag (mean `|C|` 0.241 → 0.036), so clearing it proves only that *some* genuine simultaneous cross-channel coupling is required. It cannot separate lag from magnitude, and it must never be described as "the circular-shift null" without that qualification.

**N3** changes the phase content, so it says nothing about lag versus magnitude. It cannot reject a nonstationarity whose change point happens to coincide with the true task onset — rotation moves the boundary away from such an event by design, so no within-session null can separate them. Block quantization (30 s) means a ±1-block rotation is nearly the true partition, which inflates the null's upper tail and makes the test conservative.

**N4** is a placebo, not a gate: its pseudo-phases are carved from one ~600 s resting recording and hold roughly a quarter of the segments of their real counterparts, so it is reported against a **duration-matched** real comparator; without that, "no task" is confounded with "less data".

**Matched-strength**, the incumbent, cannot test the coherency estimator, the band transform, the frequency-band split, session nonstationarity, drift, artifact epochs, or the phase segmentation. Its empirical behaviour here (cell-for-cell agreement with N1b) shows it is functioning as an independence null.

## 6. Null calibration of the five-phase functionals — none of them is zero-centred

**Every one of the four cross-phase functionals returns a significantly positive value on an input that contains no task at all, and this is not confined to the conditional statistic.** Five pseudo-phases {A, B, task_learn, task_test, rest_post} were carved out of a *single resting recording*, sized in proportion to the true phase durations, in two variants: **ordered** (contiguous windows in true temporal order — the audit_168 construction, drift retained) and **shuffled** (same sizes, randomly reassigned 30 s blocks — drift destroyed). Both were run on `rest_pre` and, independently, on `rest_post`. Cohort medians over the 16 scales, with the count of scales at which the value is significantly above zero (one-sided Wilcoxon, n = 10):

| functional | band | real (full) | real (dur-matched) | sham ordered (rest_pre) | sham shuffled (rest_pre) | sham shuffled (rest_post) |
|---|---|---|---|---|---|---|
| **T_test** | alpha | +0.127 · 12/16 | +0.004 · 1/16 | +0.110 · 9/16 | +0.062 · **16/16** | +0.039 · **16/16** |
| **T_test** | beta | +0.237 · 16/16 | +0.046 · 0/16 | +0.167 · 14/16 | +0.081 · **16/16** | +0.090 · **16/16** |
| **T_learn** | beta | +0.208 · 16/16 | +0.007 · 0/16 | +0.132 · 15/16 | +0.077 · **16/16** | +0.054 · 15/16 |
| **T_infspec** | beta | +0.043 · 6/16 | +0.069 · 5/16 | +0.089 · 9/16 | +0.037 · **16/16** | +0.058 · **16/16** |
| **T_infspec_pe** | beta | +0.090 · 9/16 | +0.066 · 2/16 | **+0.119 · 16/16** | +0.053 · **16/16** | +0.077 · **16/16** |

Three things follow, and the third is the actionable one.

**The offset is not only drift.** The ordered sham runs roughly twice the shuffled sham (β `T_infspec_pe` +0.119 vs +0.053; β `T_test` +0.167 vs +0.081), so session order does inflate these statistics. But the **shuffled** sham — in which temporal order, and therefore drift, has been destroyed — is still significantly positive at 16 of 16 scales for α and β on both source recordings. The construction itself manufactures a positive value: five pseudo-phases carved from one recording share structure, and these difference-correlations convert shared structure into a positive score. This is the mechanism separation the earlier work could not perform, and it rules out "just fix the drift null" as a remedy.

**`T_infspec_pe` reproduces the S2 precedent and generalizes it.** The ordered sham (+0.119, 16/16) **exceeds** the real value (+0.090, 9/16) for β, exactly as `S2_drift_controls.md` recorded (+0.243 sham vs +0.091 real). The withdrawal of the "inference is drift-confounded" verdict was correct — the null was invalid — but the deeper problem is now visible: the statistic is not zero-centred under *any* of these constructions, drift-free ones included.

**What this does and does not invalidate.** It is not the algebra: fed five independent white-noise vectors the four functionals return −0.031 to +0.030, so the estimator is unbiased on independent input and the offset comes from shared structure in real cophenetic distances. And critically, **the project's gate has never been "T > 0"** — it is "observed > surrogate", and a per-phase surrogate (matched-strength, N1) inherits the same construction and therefore carries the same baseline. So this calibration does **not** retroactively void the surrogate-referenced results. What it does void is any statement of the form "`T_infspec_pe` is positive, therefore there is an inference-specific trace", and it imposes a standing requirement: **before quoting a p-value for a conditional functional against any null, show that the null reproduces the no-signal baseline** rather than sitting at zero.

**Limits of this calibration, stated plainly.** The sham is a *probe, not a gate*, for two reasons. Its five pseudo-phases come from one recording and so share more structure than genuinely distinct recordings do, which biases it **upward** relative to a fair no-task arc. And it necessarily uses ~20 % of each real phase's duration, which matters enormously: truncating the real arc to the same duration collapses β `T_test` from +0.237 (16/16) to +0.046 (0/16). These functionals are strongly data-hungry, so "real full-duration exceeds sham" and "real duration-matched falls below sham" are both true and neither is the decisive comparison. The decisive comparison remains observed-versus-surrogate at full duration. One further inconsistency to note: the sham estimates all five pseudo-phases at the full `nperseg`, while the real arc uses `nperseg // 2` for A and B (see §7).

## 7. A pipeline asymmetry found in passing — the A/B `nperseg`

Not a null, but it surfaced while building one and it affects the observed statistic every lane uses. The production pipeline estimates the two `rest_pre` halves A and B at `nperseg_for_fs(fs) // 2` (the audit_63 pre-flight choice, made to keep the segment count comparable on a half-length recording) while `task_test` and `rest_post` use the full `nperseg_for_fs(fs)`. A different Welch segment length means a different frequency resolution and hence a different estimation bias. Now note where those matrices land: in `ρ_sym = ½[ρ(D_task − D_A, D_post − D_B) + ρ(D_task − D_B, D_post − D_A)]`, **both** arguments of each Spearman have the form (full-`nperseg` estimate) − (half-`nperseg` estimate). Any systematic part of that estimator difference is common to the two arguments, and a shared additive component inflates a correlation. Script 10 measures the size of this by recomputing the observed statistic with A and B matched to the full `nperseg`.

**It is not a material artifact — this one is clean.** Across α and β at all 16 scales (n = 10), the median shift from matching `nperseg` is between −0.056 and +0.065, with **1 of 32 cells** nominally significant (β at s = 45.1, p = 0.042) — exactly what 32 tests produce by chance. For β the shift is not even consistently signed: the matched construction is *higher* than the canonical one at 9 of 16 scales (e.g. s = 4.0: +0.310 matched vs +0.267 canonical). The cohort medians track each other closely throughout (β s = 1: +0.253 canonical vs +0.207 matched; α s = 4: +0.173 vs +0.173).

So the shared full-minus-half estimator-bias route, which is real in principle, does not carry measurable weight here, and the A/B construction the incumbent gate rests on does not need respecifying. Worth noting because the mechanism was plausible enough that it had to be checked rather than argued away. The caveat still applies to how strongly one can read the null result: matching `nperseg` also roughly halves the number of Welch segments in A and B, so the comparison mixes the removed bias with added estimator noise — but that noise would push the matched value *down*, and it did not go down, which makes the negative the clean direction of failure.

## 8. What was not run, and why

Everything cut is listed here explicitly; nothing was silently truncated.

**`high_gamma`** (80–300 Hz, 441 in-band bins) was excluded from every new rung: its per-segment coefficient array is 284 MB per phase per worker, and the box was already OOM-killing processes across the three Wave-0 lanes. The incumbent matched-strength gate gives it 4/16 scales at raw p and 0/16 under BH, so it is not load-bearing — but its status under the lag null is **unknown**, not null.

**`n3_free`** (the free block permutation) was dropped after the mid-run resource throttle capped this lane at 4 workers. It is the *looser* of the two N3 brackets — it destroys temporal order entirely, so it cannot distinguish "the boundaries matter" from "the session has any temporal structure at all". `n3_order`, the exact rotation test, is the scientifically sharper variant and was kept. A partial `n3_free` run before the throttle showed per-patient behaviour consistent with `n3_order` but was not checkpointed and is not reported.

**N2 on the five-phase functionals** was not run, on the grounds established in §5 that N2 is a near-duplicate of N1 on this pipeline rather than an independent rung.

**Surrogate count** was held at R = 200 throughout; no rung had its R reduced.

## 9. The inheritance rule — what must accompany every cohort claim

**Matched-strength alone is no longer sufficient, and the reason is now empirical rather than theoretical: it and the independent segment-lattice null N1b agree cell for cell, which shows matched-strength has been functioning as an *independence* null — it tests whether any genuine cross-channel coupling is required, and nothing more.** From here, every cohort-level claim in this paper must be reported against **two** nulls, not one: matched-strength (or equivalently N1b) for "is there coupling structure at all", **and N1, the lag-randomizing null that preserves each phase's coherence magnitude `|C_ij(f)|` exactly**, for "is the effect carried by the time-lagged interaction `|ImCoh|` was chosen to isolate". Where N1 is cleared, the claim may be stated as being about lagged interaction. Where it is not — which is currently the case for α and β at every scale — the claim must be stated as being about **coherence structure**, a property `|ImCoh|` shares with ordinary coherence, and the manuscript must not use volume-conduction immunity as its justification for the result. Three riders. First, N2 is not corroboration of N1: it has the same first moment and merely lower variance, so report it as N1's low-variance limit or omit it. Second, any claim that the trace is tied to the task rather than to the session requires **N3 order-preserving**, because no rung that keeps the true phase boundaries can address session nonstationarity. Third, and non-negotiable: **no cross-phase functional may ever be tested against zero.** §6 shows that all four — `T_test` included, not merely the conditional `T_infspec_pe` — return significantly positive values at 16 of 16 scales for α and β on a five-phase arc carved from a single resting recording with temporal order destroyed. Zero is not the null for any of them; the surrogate is. Every reported effect must therefore be a **margin against a per-phase surrogate that inherits the same construction**, never a raw value, never a sign count, and never a "T > 0 hence a trace" statement. For conditional statistics the requirement is stronger still: show that the chosen null reproduces the no-signal baseline before quoting its p-value at all, because for β `T_infspec_pe` the no-signal sham (+0.119) *exceeds* the real value (+0.090).

## 10. Where this leaves the paper

Three of the four Wave-0 fracture-6 questions are now answered, and the answers are not the ones the plan hoped for.

**Is the trace real?** Yes. It clears matched-strength and, independently, the segment-lattice independence null N1b, which agree at r = 0.983 across 800 cells. β holds at 16/16 scales under both, LOO-robust, and survives BH across the band × scale family. Nothing in this lane weakens that.

**Is it about lagged interaction?** No — and this is the load-bearing negative. Every lag-destroying, magnitude-preserving rung fails for every band at every scale, and the fixed-backbone control (the one designed to give the trace its best chance) shows *zero* separation: β median margin −0.003. `|ImCoh|` is still volume-conduction-immune as an estimator, but the trace it reports is reproduced by the coherence-magnitude geometry alone. The manuscript cannot present this result as evidence about time-lagged coupling, and cannot cite it to justify `|ImCoh|` over ordinary coherence.

**Are the encoding/inference functionals safe to gate on?** Only against a surrogate, never against zero. All four are significantly positive at 16/16 scales for α and β on a no-signal arc with drift destroyed, and for β `T_infspec_pe` the no-signal sham exceeds the real value.

**What is genuinely unsettled.** Whether the trace survives session repartitioning (N3, running at the time of writing — §11); `high_gamma` under any new rung; and, most importantly, the constructive question this lane can pose but not answer — **if the discriminative content is coherence magnitude rather than lag, is ordinary coherence the better substrate for this paper?** That belongs to W0-A, which owns the transform decision, and it should be told that the lag-specificity argument for `|ImCoh|` is not supported by the trace result.

## 11. Still running at the time of writing

Two rungs were still computing when this report was finalized. Both checkpoint every cell to disk as it completes, so neither can be lost, and both write their own cohort-gate table on exit. **Whoever picks this up should read these two before treating the ladder as complete** — neither can overturn §3 (they test different alternatives), but §6's and §10's open items depend on them.

| what | script | where the result lands | log |
|---|---|---|---|
| **N1 on the five-phase functionals** (`T_test`, `T_learn`, `T_infspec`, `T_infspec_pe`) | `scripts/01_compute/null_ladder/08_run_functionals_n1.py` | `data/paper_final/w0b_nulls/func_n1/{per_patient_scale,cohort_gate}.csv` | `logs/func_n1.log` |
| **N3 order-preserving (exact rotation) + N4 placebo** | `scripts/01_compute/null_ladder/04_run_n3_n4.py` | `data/paper_final/w0b_nulls/{n3_order,n4}/` | `logs/n3n4_final.log` |

Partial evidence already in hand for each. For the **functional N1**, the per-cell observed values printed during the run reproduce the canonical arc (e.g. Pat_05 β `T_test` +0.449, `T_learn` +0.393, `T_infspec` +0.099, `T_infspec_pe` +0.200), and the pattern visible so far in α is that `T_infspec` and `T_infspec_pe` are frequently *negative* per patient (Pat_02 −0.059/−0.053, Pat_07 −0.087/−0.056, Pat_10 −0.050/−0.031) — consistent with §6's finding that these functionals carry a construction baseline rather than a robust task signal, but not yet a cohort verdict. For **N3 order-preserving**, the pre-throttle partial run showed observed values sitting clearly above the rotation null for β in several patients (Pat_07 +0.299 vs +0.068; Pat_08 +0.575 vs +0.072) and not in others (Pat_10 +0.072 vs +0.052), i.e. plausibly a real but not unanimous effect. **Neither of these partial readings should be quoted as a result.**

Re-run `scripts/01_compute/null_ladder/05_summarize_ladder.py` once they land; it picks up `n3_order` and `n4` automatically and re-derives the whole ladder table with BH over the common band × scale family.

## 12. Reproducing any of this

```
ln -s <repo>/data data                        # data/ is gitignored and not in the worktree
scripts/01_compute/null_ladder/run.sh <script>   # pins PYTHONPATH to THIS worktree's src/
```

`00_validate_identity_and_observed.py` (the preflight in §2) and `02_validate_library.py` (primitive-level invariants, the ×4.6 magnitude calibration, and the global-rescale-invariance proof) should be run first on any machine — they are the checks that make the rest interpretable.

N3 (free + order-preserving), N4, the functional calibration, the functional N1, and the `n1_fixedbb` control were running when this report was first written; their sections are filled in as they complete.
