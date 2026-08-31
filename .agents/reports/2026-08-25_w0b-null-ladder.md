---
name: 2026-08-25_w0b-null-ladder
kind: report
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 0, lane W0-B)
status: current
created: 2026-08-25
updated: 2026-08-31
scope: The timeseries-level null ladder for the cross-phase trace. Builds and runs the nulls that act UPSTREAM of the finished FC matrix — where matched-strength cannot reach — and reports, per band and per scale, which rungs the α/β trace survives and which it does not, plus an explicit statement per rung of what that rung cannot reject. Includes the null-calibration of the five-phase encoding/inference functionals, which closes the placebo gap W0-C left open on T_infspec_pe.
pointers:
  - src/lrg_eegfc/utils/surrogate/timeseries_nulls.py          # the library deliverable
  - src/lrg_eegfc/utils/fc/heat_multiscale.py                  # cross_phase_functionals, partial_spearman
  - scripts/01_compute/null_ladder/                            # runners 00–10
  - data/paper_final/w0b_nulls/                                # all computed artifacts
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py   # the incumbent gate reproduced exactly
  - .agents/preprint/supplementary/S2_drift_controls.md        # the broken-null precedent this lane re-tests
---

# W0-B — the null ladder

## Head

**The α and β traces are real and are tied to the task, but their content is not lag-specific — and that last part is the finding that changes what this paper can claim.** Three results, in the order they constrain each other. (i) The trace survives every *independence* null: the incumbent matched-strength gate (β 16/16 scales, α 12/16) and, independently, the segment-lattice shift N1b (β 16/16, α 13/16), which agree with matched-strength at **r = 0.983** across 800 cells — so the β result is not an artifact of the matched-strength construction. (ii) It survives **session repartitioning**: under the exact order-preserving rotation test, which holds the session, its drift and the phase durations fixed and moves only where the phase boundaries fall, β clears **10 of 16 scales** (best p = 0.0049, 8/16 leave-one-patient-out robust) — so the geometry is tied to the real task boundaries, not to the session's own structure. (iii) But **every lag-destroying, coherence-magnitude-preserving null fails, for every band, at every scale.** N1 clears 0/16 (β best p = 0.080); the fixed-backbone control `n1_fixedbb`, built specifically to give the trace its best chance by holding the graph identical, shows **no separation at all** — cohort median margins of **−0.003 (β)** and **+0.004 (α)**, against +0.214 and +0.130 under matched-strength. A surrogate that discards every trace of the observed lag structure reproduces the α and β traces exactly. The same holds for the whole five-phase encoding/inference decomposition: across 320 (functional × band × scale) tests under N1, **not one cell survives BH** (minimum q = 0.801), `T_infspec_pe` included.

Stated precisely, because the distinction matters and is easy to get wrong: this is **not** a finding that the measurement is contaminated by volume conduction. `|ImCoh|` remains volume-conduction-immune by construction — zero-phase mixing contributes exactly zero, and no null takes that away. Nor is it a finding that the trace is an artifact; (i) and (ii) rule that out. It is narrower and still serious: **the discriminative content of the trace is coherence-magnitude geometry, not time-lagged interaction.** The paper therefore cannot present these results as evidence about lagged coupling, and cannot cite them to justify `|ImCoh|` over ordinary coherence — on this statistic the two are not distinguishable.

---

## 1. What was built

`src/lrg_eegfc/utils/surrogate/timeseries_nulls.py` — surrogates that act on the per-segment Fourier coefficients or on the complex band coherency `C(f)`, i.e. **upstream** of the `N × N` adjacency that matched-strength shuffles. No matrix-level rotation null is rebuilt (the archived CRC/SB-CRC family stays archived).

| rung | construction | preserved exactly | destroyed |
|---|---|---|---|
| **N1** | per-channel circular shift as a frequency-domain phase ramp | `\|C_ij(f)\|`, every per-channel spectrum | cross-spectral phase = **lag** |
| **N1 fixed-bb** | N1, with the backbone edge set pinned to the observed one | as N1, **plus the graph topology** | lag only |
| **N2** | iid uniform phase per (channel, frequency) | `\|C_ij(f)\|`, per-channel spectra | lag + cross-frequency phase consistency |
| **N1b** | circular shift of each channel's Welch **segment index** | each channel's Welch PSD | cross-channel temporal correspondence (magnitude included) |
| **N3 order** | exact rotation of the phase labels along the session ring | session content, drift, durations, contiguity | where the boundaries fall |
| **N4** | within-rest placebo, duration-matched | durations, degrees of freedom | the task itself |

Literature anchors: Theiler et al. (1992) for phase randomization; Perkel, Gerstein & Moore (1967) for the shift predictor that N1b generalizes; Nolte et al. (2004) and Ewald et al. (2012) for the estimator; Welch (1967) for the segmentation.

## 2. Three preflight facts that shaped the design

**The observed statistic is reproduced.** ρ_sym recomputed from raw timeseries through the surrogate code path matches the value script 13 gates **exactly (bit-identical) on 791 of 800 cells**, worst-case deviation `6.0e-05` on the remaining 9 — the residue of the float32 FC cache, four orders of magnitude below any effect here. The library CSD reproduces the cached `imcoh_abs` adjacency to `5.9e-08`.

**The circular-shift identity is exact only in the single-segment DFT sense, and the brief's premise needed correcting.** For one un-windowed segment the shift theorem holds to `2e-15` and `|C_ij(f)|` is preserved to `7e-16`. Under the **Welch estimator the pipeline actually uses**, a genuine whole-recording roll does *not* preserve magnitude: mean `|C|` collapses from **0.241 to 0.036**, because rolling channel `c` makes its segment `k` pair with an entirely different epoch of channel `j`, so cross-spectral terms average incoherently. A whole-recording circular shift under Welch is therefore an **independence null**, not a lag null. The two are separated here rather than conflated: N1 is the frequency-domain, magnitude-preserving lag null; N1b is the honest time-domain version, labelled as the independence null it is.

**N1/N2 inflate the edge weights, and that is licensed.** Randomizing the phase drives `⟨|Im C|⟩_f → (2/π)⟨|C|⟩_f`: measured ×4.6 on Pat_05 β (0.0332 → 0.1538, analytic prediction 0.1535). The surrogate sits **above** the data in raw edge weight. Harmless *only* because the readout is invariant to a global positive rescale of `W` — the backbone is rank-based and `s = τλ_max` makes `τL` scale-free. Verified numerically, not assumed: `max |ρ_sym(W) − ρ_sym(6.5W)|` over the 16 scales is **0.0e+00**. The one-sided direction was verified by re-running the identical gate code on matched-strength, which still clears; a sign error would have broken that too.

## 3. The verdict — ρ_sym, per band, per scale

Scales cleared of 16, cohort one-sided Wilcoxon on the margin `obs − surr_p50` across the 10 patients. `high_gamma` was not run on the new rungs (§8).

| band | ms | N1b | **N1** | **N1 fixed-bb** | **N2** | **N3 order** | N4 |
|---|---|---|---|---|---|---|---|
| delta | 8 | 8 | **0** | **0** | **0** | 0 | 0 |
| theta | 0 | 0 | **0** | **0** | **0** | 0 | 0 |
| **alpha** | 12 | 13 | **0** | **0** | **0** | 1 | 0 |
| **beta** | 16 | 16 | **0** | **0** | **0** | **10** | 0 |
| low_gamma | 0 | 0 | **0** | **0** | **0** | 0 | 0 |

Best (smallest) cohort p over the 16 scales — best-scale selection is already optimistic, and every lag rung still fails:

| band | ms | N1b | **N1** | **N1 fixed-bb** | **N2** | **N3 order** | N4 |
|---|---|---|---|---|---|---|---|
| delta | 0.0137 | 0.0029 | 0.7842 | 0.2158 | 0.4229 | 0.0527 | 0.5000 |
| theta | 0.2158 | 0.1875 | 0.5771 | 0.2783 | 0.6875 | 0.7539 | 0.3125 |
| **alpha** | 0.0068 | 0.0068 | **0.1377** | **0.1377** | **0.1611** | 0.0186 | 0.7842 |
| **beta** | 0.0010 | 0.0010 | **0.0801** | **0.2158** | **0.0654** | **0.0049** | 0.3125 |
| low_gamma | 0.1162 | 0.1162 | 0.2461 | 0.0801 | 0.2158 | 0.0967 | 0.1611 |

Cohort median margin `obs − surr_p50`, taken across all 16 scales — the effect size the gate is built on:

| band | ms | N1b | **N1** | **N1 fixed-bb** | **N3 order** | N4 |
|---|---|---|---|---|---|---|
| **alpha** | +0.130 | +0.127 | +0.016 | **+0.004** | +0.048 | −0.079 |
| **beta** | +0.214 | +0.234 | +0.174 | **−0.003** | +0.105 | −0.026 |

Under BH-FDR across the common (5 bands × 16 scales) family, within rung: **ms β 16/16 (min q 0.042) and N1b β 16/16 (min q 0.030) survive; every lag rung reaches a minimum q of 0.81–0.99; N3 order reaches min q 0.186.** Under leave-one-patient-out of the verdict, ms β 16/16 and N1b β 16/16 hold with any single patient dropped, N3 order holds at 8/16 for β, and the lag rungs hold at 0 for every band.

**The fixed-backbone control removes the strongest objection to this negative.** N1 drives `W` toward `(2/π)|C|`, and the `|C|` edge ranking is not the `|ImCoh|` ranking — on sEEG, `|C|` is dominated by same-probe pairs (a 2–8× bias, CLAUDE.md invariant 5). A plain N1 surrogate therefore lives on a *different* graph, so one could argue its high floor is the cross-phase stability of that anatomical scaffold rather than anything about lag. `n1_fixedbb` pins the backbone edge set to the observed `|ImCoh|` one for surrogate and observed alike, so only the weights' lag content varies. The observed statistic is unchanged by this (`obs_dev = 0.00e+00` on every cell — masking `Wobs` with its own backbone returns the backbone). The trace still fails, and more decisively than under plain N1: the cohort median margin is **−0.003 for β, +0.004 for α**. Not a reduced effect — no separation.

**This negative is attributable to the science, not to a construction artifact.** One candidate artifact was found and cleared (§7, the A/B `nperseg` asymmetry — not material). The gate direction was verified against a rung known to clear. The observed statistic reproduces the incumbent bit-identically. The surrogate's weight inflation is provably absorbed by a scale-free readout. And the topology objection is answered by `n1_fixedbb`. The remaining honest caveat is stated in §5: N1 holds each phase's own `|C|` fixed, so it cannot distinguish "the trace is magnitude-carried" from "the trace is magnitude-carried *and* lag-carried, with magnitude sufficient on its own" — but either way the lag-specific claim is unsupported.

## 4. Why the lag nulls fail — the null floor, not the effect

Under plain N1 the effect does not disappear: β's cohort median margin is +0.174 against +0.214 under matched-strength. What changes is the null, in two ways that compound. Its **spread** is 3–4× wider than matched-strength's (β 0.375 vs 0.096) and 7–8× wider than N1b's (0.050). And its **height varies across patients**, which the independence nulls' does not: the interquartile range of the per-patient null median, pooled over scales, is **0.232 (β) and 0.534 (α) under N1** against **0.0024 and 0.0029 under N1b**. Matched-strength and N1b give every patient the same near-zero floor; N1 gives each patient a floor set by their own coherence geometry. At β's best scale:

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

7/10 positive, median +0.155, but Pat_06's null sits at +0.929 — for that patient a lag-randomized surrogate reproduces essentially the entire trace. The α table is worse (Pat_06 observed +0.886 against an N1 null of +0.865). In a subset of patients the coherence-magnitude geometry alone produces a trace as large as the observed one. Once the graph is also held fixed (`n1_fixedbb`) that becomes true of the cohort as a whole.

## 5. Session repartitioning — the one rung the trace passes

**β's trace is tied to the real task boundaries.** The order-preserving N3 treats the session as a ring, rotates the phase labels along it, and cuts at the true boundaries: every pseudo-phase stays a contiguous interval of the true length, local temporal order survives everywhere but one seam, and any smooth session trend is *retained* rather than removed. Only the placement of the cut points moves. The realization set is the rotation group (78–91 rotations per patient), so all rotations are enumerated and the per-patient test is **exact** rather than Monte-Carlo.

β clears **10 of 16 scales** (best p = 0.0049, cohort median margin +0.105, observed median +0.211 against a rotation null of +0.059), and **8 of 16 survive leave-one-patient-out**. α clears 1/16 (best p = 0.0186); δ, θ and low_γ clear none.

Two qualifications, both real. First, β's N3 result **does not survive BH** across the band × scale family (min q = 0.186), where matched-strength and N1b both do — so this is suggestive rather than decisive. Second, the test is conservative by construction: a rotation of ±1 block is nearly the true partition, which inflates the null's upper tail. That conservatism cuts in the trace's favour for interpreting the raw-p result and against it for interpreting the BH failure, so the fair summary is that **the session-structure explanation is disfavoured for β but not excluded at corrected significance.**

**N4, the within-rest placebo, is underpowered by construction and should not be read as a negative.** All bands clear 0/16, and the margins are *negative* for α (−0.079) and β (−0.026) — the duration-matched real arc sits slightly *below* the within-rest placebo. That is the direct consequence of the data-hungriness measured in §6: truncating the real arc to the ~20 % of each phase that a single resting recording can supply collapses β `T_test` from +0.237 (16/16 significant) to +0.046 (0/16). At that duration neither arm has the resolution to separate. N4 therefore establishes the *baseline level* of these functionals on a no-task session — which is its load-bearing use in §6 — but it cannot adjudicate the trace.

## 6. Null calibration of the five-phase functionals — none of them is zero-centred

**Every one of the four cross-phase functionals returns a significantly positive value on an input that contains no task at all, and this is not confined to the conditional statistic.** Five pseudo-phases {A, B, task_learn, task_test, rest_post} were carved out of a *single resting recording*, sized in proportion to the true phase durations, in two variants: **ordered** (contiguous windows in true temporal order — the audit_168 construction, drift retained) and **shuffled** (same sizes, randomly reassigned 30 s blocks — drift destroyed). Both were run on `rest_pre` and, independently, on `rest_post`. Cohort medians over the 16 scales, with the count of scales significantly above zero (one-sided Wilcoxon, n = 10):

| functional | band | real (full) | real (dur-matched) | sham ordered (rest_pre) | sham shuffled (rest_pre) | sham shuffled (rest_post) |
|---|---|---|---|---|---|---|
| **T_test** | alpha | +0.127 · 12/16 | +0.004 · 1/16 | +0.110 · 9/16 | +0.062 · **16/16** | +0.039 · **16/16** |
| **T_test** | beta | +0.237 · 16/16 | +0.046 · 0/16 | +0.167 · 14/16 | +0.081 · **16/16** | +0.090 · **16/16** |
| **T_learn** | beta | +0.208 · 16/16 | +0.007 · 0/16 | +0.132 · 15/16 | +0.077 · **16/16** | +0.054 · 15/16 |
| **T_infspec** | beta | +0.043 · 6/16 | +0.069 · 5/16 | +0.089 · 9/16 | +0.037 · **16/16** | +0.058 · **16/16** |
| **T_infspec_pe** | beta | +0.090 · 9/16 | +0.066 · 2/16 | **+0.119 · 16/16** | +0.053 · **16/16** | +0.077 · **16/16** |

**The offset is not only drift.** The ordered sham runs roughly twice the shuffled sham (β `T_infspec_pe` +0.119 vs +0.053; β `T_test` +0.167 vs +0.081), so session order does inflate these statistics. But the **shuffled** sham — temporal order destroyed, hence drift destroyed — is still significantly positive at 16 of 16 scales for α and β on both source recordings. The construction itself manufactures a positive value: five pseudo-phases carved from one recording share structure, and these difference-correlations convert shared structure into a positive score. This is the mechanism separation the earlier work could not perform, and it rules out "just fix the drift null" as a remedy.

**`T_infspec_pe` reproduces the S2 precedent and generalizes it.** The ordered sham (+0.119, 16/16) **exceeds** the real value (+0.090, 9/16), exactly as `S2_drift_controls.md` recorded (+0.243 sham vs +0.091 real). The withdrawal of the "inference is drift-confounded" verdict was correct — the null was invalid — but the deeper problem is now visible: the statistic is not zero-centred under *any* of these constructions, drift-free ones included. **This closes the placebo gap W0-C left open**: that lane could not build a data-based placebo for the conditional statistic and marked `T_infspec_pe` provisional pending this ladder. The placebo now exists, and it says the conditional functional carries a construction baseline that a comparison against zero would misread as signal.

**It is not the algebra.** Fed five independent white-noise vectors the four functionals return −0.031 to +0.030, so the estimator is unbiased on independent input; the offset comes from shared structure in real cophenetic distances, not from the formula.

**What this does and does not invalidate.** The project's gate has never been "T > 0" — it is "observed > surrogate", and a per-phase surrogate inherits the same construction and therefore the same baseline. So this calibration does **not** retroactively void surrogate-referenced results, W0-C's included. What it voids is any statement of the form "`T_infspec_pe` is positive, therefore there is an inference-specific trace", and it imposes a standing requirement: before quoting a p-value for a conditional functional, show that the null reproduces the no-signal baseline rather than sitting at zero.

**Limits, stated plainly.** The sham is a *probe, not a gate*. Its pseudo-phases come from one recording and so share more structure than genuinely distinct recordings do, biasing it **upward**. And it uses ~20 % of each real phase's duration, which matters enormously (β `T_test` +0.237 → +0.046). So "real full-duration exceeds sham" and "real duration-matched falls below sham" are both true and neither is decisive; the decisive comparison remains observed-versus-surrogate at full duration. One inconsistency to note: the sham estimates all five pseudo-phases at the full `nperseg`, while the real arc uses `nperseg // 2` for A and B (§7).

## 6a. The lag null on the five-phase functionals — nothing survives

N1 was run over the full five-phase arc, scoring all four functionals per scale, R = 200, 50 cells. **Across the whole 320-test family (4 functionals × 5 bands × 16 scales), not one cell survives BH; the minimum q is 0.801.**

| functional | band | raw p < 0.05 | min p | min q |
|---|---|---|---|---|
| T_test | beta | 0/16 | 0.0654 | 0.801 |
| T_learn | alpha | 0/16 | 0.0527 | 0.801 |
| T_learn | beta | 0/16 | 0.0654 | 0.801 |
| **T_infspec** | delta | 5/16 | 0.0137 | 0.801 |
| **T_infspec** | beta | 2/16 | 0.0244 | 0.801 |
| T_infspec_pe | beta | 0/16 | 0.0801 | 0.801 |
| T_infspec_pe | delta | 0/16 | 0.0527 | 0.801 |

The only cells reaching raw significance are `T_infspec` in δ (5/16) and β (2/16) — the one contrast whose first argument, `D_test − D_learn`, never references `rest_pre` and is therefore structurally immune to split-half baseline artifacts. That is a coherent place for residual signal to sit, and it is worth a targeted follow-up, but **at 2–5 of 16 scales with q = 0.80 it is not a result.** Taken with §6, the position on the encoding/inference decomposition is: it is not zero-centred, and what remains after conditioning does not survive a lag null.

## 6b. What N1's null actually is — the pipeline on ordinary coherence

The N1 surrogate is not an abstract construct. Because randomizing the phase sends `⟨|Im C|⟩_f → (2/π)⟨|C|⟩_f`, N1's null should be, near enough, *the incumbent pipeline run on ordinary coherence magnitude*. That was tested directly: `W = ⟨|C|⟩_f` was pushed through the identical backbone → Laplacian → heat-kernel → UPGMA → ρ_sym path for 5 patients × {α, β} × 16 scales, and compared against the N1 surrogate medians from the production run.

They track at **Pearson r = 0.879** (mean absolute difference 0.156) — so yes, N1's floor is essentially "what ordinary coherence would have given you".

That makes the negative concrete, and it carries a consequence the transform lane needs:

| band | median ρ_sym on `\|ImCoh\|` | median ρ_sym on `⟨\|C\|⟩` |
|---|---|---|
| alpha | +0.212 | **+0.277** |
| beta | +0.281 | **+0.355** |

**Ordinary coherence produces a *larger* cross-phase trace than `|ImCoh|` does** (5 patients, all 16 scales). The estimator this project selected for its volume-conduction immunity is the weaker of the two on the statistic the paper is built on. Per-patient the gap is not uniform — Pat_05 β runs the other way (`|ImCoh|` +0.449 vs `|C|` −0.329), and Pat_06 β is the extreme case in the opposite direction (+0.219 vs +0.949, which is exactly why Pat_06 drives the N1 floor in §4) — so this is a 5-patient indication, not a cohort verdict. It is enough to define the experiment W0-A should run, and it is stated in §11.

## 7. A pipeline asymmetry found in passing — the A/B `nperseg`

Not a null, but it surfaced while building one and it affects the observed statistic every lane uses. The production pipeline estimates the two `rest_pre` halves at `nperseg_for_fs(fs) // 2` while `task_test` and `rest_post` use the full `nperseg_for_fs(fs)`. In `ρ_sym = ½[ρ(D_task − D_A, D_post − D_B) + ρ(D_task − D_B, D_post − D_A)]`, **both** arguments of each Spearman then have the form (full-`nperseg`) − (half-`nperseg`). Any systematic part of that estimator difference is common to the two arguments, and a shared additive component inflates a correlation.

**It is not a material artifact — this one is clean.** Across α and β at all 16 scales (n = 10), the median shift from matching `nperseg` is between −0.056 and +0.065, with **1 of 32 cells** nominally significant (β at s = 45.1, p = 0.042) — what 32 tests produce by chance. For β the shift is not even consistently signed: the matched construction is *higher* at 9 of 16 scales. The A/B construction does not need respecifying. Worth recording because the mechanism was plausible enough that it had to be checked rather than argued away, and because it is the artifact that was found and cleared before attributing the §3 negative to the science.

## 8. What was not run, and why

Everything cut is listed here; nothing was silently truncated.

**`high_gamma`** (80–300 Hz, 441 in-band bins) was excluded from every new rung: its per-segment coefficient array is 284 MB per phase per worker, and the box was OOM-killing processes across the three Wave-0 lanes at the time. Matched-strength gives it 4/16 scales at raw p and 0/16 under BH, so it is not load-bearing — but its status under the lag null is **unknown**, not null.

**`n3_free`** (the free block permutation) was dropped at the resource throttle. It is the *looser* of the two N3 brackets — destroying temporal order entirely, it cannot distinguish "the boundaries matter" from "the session has any temporal structure at all". `n3_order`, the exact rotation test, is the sharper variant and was kept.

**N2 on the five-phase functionals** was not run: N2 is a near-duplicate of N1 on this pipeline (§9), not an independent rung.

**Surrogate count** was held at R = 200 throughout; no rung had its R reduced.

## 9. What each rung CANNOT reject

**N1 / N1 fixed-bb / N2** cannot reject an account in which the trace lives in cross-phase changes of coherence **magnitude** — they hold each phase's own `|C_ij(f)|` fixed by construction. They cannot reach session nonstationarity, drift, artifact epochs, or the phase segmentation, because all phases keep their true boundaries (that is N3's job). And they cannot distinguish "magnitude-carried" from "magnitude *and* lag carried, with magnitude sufficient alone" — only that the lag-specific claim is unsupported.

**N2 is not an independent rung**, and the measured reason is not the one theory predicts. Both send the phase to uniform, so they share a first moment; N2 was expected to be the low-variance limit because its per-bin phases average down. In ρ_sym that does not materialize: null spread 0.335 (N1) vs 0.322 (N2) pooled, 0.375 vs 0.418 for β. With 7–35 in-band bins the frequency averaging is too weak, and the variance that matters is structural, not in the band mean. Report one; never treat the second as corroboration.

**N1b** destroys coupling magnitude along with lag (mean `|C|` 0.241 → 0.036), so clearing it proves only that *some* genuine simultaneous coupling is required. It cannot separate lag from magnitude, and must never be called "the circular-shift null" without that qualification.

**N3 order** says nothing about lag versus magnitude — it changes the phase content. It cannot reject a nonstationarity whose change point coincides with the true task onset: rotation moves the boundary away from such an event by design, so no within-session null can separate them. Block quantization (30 s) makes ±1-block rotations nearly the true partition, inflating the null's upper tail.

**N4** is a placebo, not a gate. Its pseudo-phases hold ~¼ the segments of their real counterparts and come from a single recording, so it is underpowered and biased upward; it establishes the no-task baseline level (§6) and cannot adjudicate the trace.

**Matched-strength**, the incumbent, cannot test the coherency estimator, the band transform, the frequency-band split, session nonstationarity, drift, artifact epochs, or the phase segmentation. Its measured behaviour here — cell-for-cell agreement with N1b at r = 0.983 — shows it functions as an independence null.

## 10. The inheritance rule — what must accompany every cohort claim

**Matched-strength alone is no longer sufficient, and the reason is now empirical rather than theoretical: it and the mechanically unrelated N1b agree at r = 0.983 across 800 cells, which shows it has been functioning as an *independence* null — it tests whether any genuine cross-channel coupling is required, and nothing more.** From here, every cohort-level claim in this paper must be reported against **three** rungs, not one: matched-strength (or equivalently N1b) for "is there coupling structure at all"; **N1 with the backbone held fixed** for "is the effect carried by the time-lagged interaction `|ImCoh|` was chosen to isolate"; and **N3 order-preserving** for "is it tied to the task rather than to the session", since no rung that keeps the true phase boundaries can address session nonstationarity. Where N1 is cleared, a claim may be stated as being about lagged interaction; where it is not — currently α and β at every scale, and all four five-phase functionals — the claim must be stated as being about **coherence structure**, a property `|ImCoh|` shares with ordinary coherence, and volume-conduction immunity must not be offered as its justification. Two riders. N2 is not corroboration of N1; it is a near-duplicate, so report one. And **no cross-phase functional may ever be tested against zero** — all four return significantly positive values at 16 of 16 scales for α and β on a no-signal arc with drift destroyed, so every reported effect must be a margin against a per-phase surrogate that inherits the same construction, never a raw value, a sign count, or a "T > 0 hence a trace" statement.

## 11. Where this leaves the paper

**Is the trace real?** Yes. β holds at 16/16 scales under matched-strength and under the independent N1b, LOO-robust, surviving BH. Nothing here weakens that.

**Is it tied to the task rather than the session?** For β, probably yes — 10/16 scales and 8/16 LOO-robust under the exact rotation test, best p = 0.0049, though it does not survive BH (min q 0.186). For α, marginal (1/16). This is the rung the trace passes, and it is worth having.

**Is it about lagged interaction?** No. Every lag-destroying, magnitude-preserving rung fails for every band at every scale, and the fixed-backbone control shows zero separation (β margin −0.003). This is the load-bearing negative.

**Are the encoding/inference functionals safe to gate on?** Only against a surrogate, never zero — and under the lag null none of them survives BH anywhere (min q 0.801), `T_infspec_pe` included.

**What is genuinely unsettled.** `high_gamma` under any new rung. Whether β's N3 result would survive BH with a less conservative session null.

**And one question this lane can now pose sharply but not settle, which belongs to W0-A.** §6b shows N1's null is, at r = 0.879, just the incumbent pipeline run on ordinary coherence — and that on 5 patients across all 16 scales, `⟨|C|⟩` yields a *larger* median ρ_sym than `|ImCoh|` in both trace bands (α +0.277 vs +0.212; β +0.355 vs +0.281). So the transform chosen for volume-conduction immunity appears to be the weaker one on the statistic the paper rests on. W0-A owns the transform decision and should be told two things: the lag-specificity argument for `|ImCoh|` is **not supported** by the trace result, and the head-to-head experiment is cheap and well-defined — run the incumbent gate on `⟨|C|⟩_f` across the full cohort and compare effect size, null floor and band selectivity against `|ImCoh|`. Note the trade this would force into the open: `⟨|C|⟩` is *not* volume-conduction-immune, so a larger trace on it is not automatically a better result — it may be a larger trace on a partly artefactual substrate. That is precisely the comparison the paper currently makes by assertion and has never made by evidence.

## 12. Reproducing any of this

```
ln -s <repo>/data data                           # data/ is gitignored, not in the worktree
scripts/01_compute/null_ladder/run.sh <script>   # pins PYTHONPATH to THIS worktree's src/
```

Run `00_validate_identity_and_observed.py` (§2 preflight) and `02_validate_library.py` (primitive invariants, the ×4.6 magnitude calibration, the global-rescale-invariance proof) first on any machine — they are the checks that make the rest interpretable. `05_summarize_ladder.py` regenerates the §3 tables from all rungs on disk.
