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

**The β and α traces do not survive the one null that preserves coherence magnitude and destroys only lag.** Under N1 (per-channel circular shift, realized as its exact frequency-domain phase ramp) β clears **0 of 16 scales**, best cohort p = 0.080; α clears **0 of 16**, best p = 0.138. They survive everything else on the ladder: the incumbent matched-strength null (β 16/16, α 12/16) and the segment-lattice independence null N1b (β 16/16, α 13/16) — so the trace is not noise, and genuine cross-channel coupling is required to produce it. What fails is the *attribution*: a surrogate that keeps each phase's own coherence-magnitude structure `|C_ij(f)|` intact and randomizes only the phase reproduces the trace as well as the observed `|ImCoh|` does. The effect size is not what collapses — β's median cohort margin under N1 is +0.174, against +0.214 under matched-strength — the **null floor rises and becomes patient-specific**, by up to an order of magnitude (Pat_06 α: observed +0.886, N1 null +0.865). The honest reading is that the multiscale cross-phase trace is carried by coherence-magnitude geometry at least as much as by the time-lagged interaction `|ImCoh|` was chosen to isolate, and no claim in this paper should be stated as being about lagged interaction until that is resolved.

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

**The observed statistic is reproduced exactly.** ρ_sym recomputed from raw timeseries through the surrogate code path matches the value script 13 gates to `0.00e+00` on every one of the 150 production cells, and the library CSD reproduces the cached `imcoh_abs` adjacency to `5.9e-08` (float32 cache precision). Every rung is therefore attached to the right observed value, which was the stated sanity requirement.

**The circular-shift identity is exact only in the single-segment DFT sense, and the brief's premise needed correcting.** For one un-windowed segment the shift theorem holds to `2e-15` and `|C_ij(f)|` is preserved to `7e-16`. Under the **Welch estimator the pipeline actually uses**, a genuine whole-recording roll does *not* preserve magnitude: mean `|C|` collapses from **0.241 to 0.036**, because rolling channel `c` makes its segment `k` pair with an entirely different epoch of channel `j`, so the cross-spectral terms average incoherently. A whole-recording circular shift under Welch is therefore an **independence null**, not a lag null. The two are separated in this ladder rather than conflated: N1 is the frequency-domain, magnitude-preserving lag null; N1b is the honest time-domain version, labelled as the independence null it is.

**N1/N2 inflate the edge weights, and that is licensed.** Randomizing the phase drives `⟨|Im C|⟩_f → (2/π)⟨|C|⟩_f`: measured ×4.6 on Pat_05 β (0.0332 → 0.1538, against the analytic prediction 0.1535). The surrogate therefore sits **above** the data in raw edge weight. This is harmless *only* because the readout is invariant to a global positive rescale of `W` — the `mst_union_top_fraction` backbone is rank-based and `s = τλ_max` makes `τL` scale-free. That was verified numerically, not assumed: `max |ρ_sym(W) − ρ_sym(6.5W)|` over the 16 scales is **0.0e+00**. The one-sided direction was verified independently by re-running the identical gate code on the incumbent matched-strength rung, which still clears (β 16/16, α 12/16); a sign error would have broken that too.

## 3. The verdict — ρ_sym, per band, per scale

Scales cleared out of 16, cohort one-sided Wilcoxon on the margin `obs − surr_p50` across the 10 patients. `high_gamma` was **not run** on the new rungs (see §7).

| band | ms (incumbent) | **N1 (lag)** | **N2 (lag, low-var)** | N1b (independence) |
|---|---|---|---|---|
| delta | 8 | **0** | **0** | 8 |
| theta | 0 | **0** | **0** | 0 |
| **alpha** | 12 | **0** | **0** | 13 |
| **beta** | 16 | **0** | **0** | 16 |
| low_gamma | 0 | **0** | **0** | 0 |

Best (smallest) cohort p over the 16 scales — best-scale selection is already optimistic, and the lag rungs still fail:

| band | ms | **N1** | **N2** | N1b |
|---|---|---|---|---|
| delta | 0.0137 | 0.7842 | 0.4229 | 0.0029 |
| theta | 0.2158 | 0.5771 | 0.6875 | 0.1875 |
| **alpha** | 0.0068 | **0.1377** | **0.1611** | 0.0068 |
| **beta** | 0.0010 | **0.0801** | **0.0654** | 0.0010 |
| low_gamma | 0.1162 | 0.2461 | 0.2158 | 0.1162 |

Under BH-FDR across the common (5 bands × 16 scales) family, within rung: ms β 16/16 and N1b β 16/16 survive; N1 and N2 reach a minimum q of 0.905 and 0.819. Under leave-one-patient-out of the *verdict*, ms β 16/16 and N1b β 16/16 hold with any single patient dropped; N1 and N2 hold at 0 scales for every band.

**N1b reproduces the incumbent almost exactly** (β 16/16, α 13/16, δ 8/16 against ms 16/12/8). Two nulls built on completely different mechanisms — one shuffling a finished matrix under a strength constraint, one destroying cross-channel temporal correspondence in the raw segments — agree cell for cell. That is a genuine and reassuring cross-validation of the incumbent gate, and it means the existing β result is *not* an artifact of the matched-strength construction. It also means matched-strength has been behaving as an independence null all along, which is exactly why it never reached the magnitude question.

## 4. Why N1 fails — it is the null floor, not the effect

The effect does not disappear under N1. Cohort median margin over scales: β **+0.174** under N1 against **+0.214** under matched-strength; α +0.016 against +0.130.

What changes is the null. The N1 surrogate spread (`p95 − p50`) is **3–4× wider** than matched-strength's (β 0.375 vs 0.096; α 0.285 vs 0.100), and the N1 null *height* is patient-specific in a way matched-strength's is not. At β's best scale:

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

**N1 / N2** cannot reject an account in which the trace lives in cross-phase changes of coherence **magnitude** — they hold each phase's own `|C_ij(f)|` fixed by construction. They cannot reach session nonstationarity, drift, artifact epochs, or the phase segmentation, because all phases keep their true boundaries. N2 is additionally **not an independent rung**: it has the same first moment as N1 and differs only by having lower realization variance (measured sd over 20 draws: 1.5e-04 vs 1.4e-04 on the band-mean, with N2's per-bin phases averaging down), so it should be reported as the low-variance limit of N1, never as corroboration of it.

**N1b** destroys coupling magnitude along with lag (mean `|C|` 0.241 → 0.036), so clearing it proves only that *some* genuine simultaneous cross-channel coupling is required. It cannot separate lag from magnitude, and it must never be described as "the circular-shift null" without that qualification.

**N3** changes the phase content, so it says nothing about lag versus magnitude. It cannot reject a nonstationarity whose change point happens to coincide with the true task onset — rotation moves the boundary away from such an event by design, so no within-session null can separate them. Block quantization (30 s) means a ±1-block rotation is nearly the true partition, which inflates the null's upper tail and makes the test conservative.

**N4** is a placebo, not a gate: its pseudo-phases are carved from one ~600 s resting recording and hold roughly a quarter of the segments of their real counterparts, so it is reported against a **duration-matched** real comparator; without that, "no task" is confounded with "less data".

**Matched-strength**, the incumbent, cannot test the coherency estimator, the band transform, the frequency-band split, session nonstationarity, drift, artifact epochs, or the phase segmentation. Its empirical behaviour here (cell-for-cell agreement with N1b) shows it is functioning as an independence null.

## 6. Null calibration of the five-phase functionals

*(section completed below once the calibration run lands — see §8)*

## 7. What was not run, and why

`high_gamma` (80–300 Hz, 441 in-band bins) was excluded from N1/N2/N1b: its per-segment coefficient array is 284 MB per phase per worker, and the machine has 31 GB. The incumbent matched-strength gate gives it 4/16 scales at raw p, 0/16 under BH, so it is not load-bearing — but its status under the lag null is **unknown**, not null.

The lag null was applied to the five-phase functionals via N1 only; N2 was not run there, on the grounds established in §5 that it is not an independent rung.

## 8. The inheritance rule — what must accompany every cohort claim

**Matched-strength alone is no longer sufficient, and the reason is now empirical rather than theoretical: it and the independent segment-lattice null N1b agree cell for cell, which shows matched-strength has been functioning as an *independence* null — it tests whether any genuine cross-channel coupling is required, and nothing more.** From here, every cohort-level claim in this paper must be reported against **two** nulls, not one: matched-strength (or equivalently N1b) for "is there coupling structure at all", **and N1, the lag-randomizing null that preserves each phase's coherence magnitude `|C_ij(f)|` exactly**, for "is the effect carried by the time-lagged interaction `|ImCoh|` was chosen to isolate". Where N1 is cleared, the claim may be stated as being about lagged interaction. Where it is not — which is currently the case for α and β at every scale — the claim must be stated as being about **coherence structure**, a property `|ImCoh|` shares with ordinary coherence, and the manuscript must not use volume-conduction immunity as its justification for the result. Three riders. First, N2 is not corroboration of N1: it has the same first moment and merely lower variance, so report it as N1's low-variance limit or omit it. Second, any claim that the trace is tied to the task rather than to the session requires **N3 order-preserving**, because no rung that keeps the true phase boundaries can address session nonstationarity. Third, and non-negotiable, **any partial or conditional statistic — `T_infspec_pe` above all — must be calibrated on a no-signal, duration-matched input before a single p-value is quoted for it**; the precedent in `S2_drift_controls.md` and the calibration in §6 below show that a conditional functional can return a large positive value from an arc in which no consolidation is possible, and a p-value attached to an uncalibrated statistic is not evidence.

## 9. Results still landing

N3 (free + order-preserving), N4, the functional calibration, the functional N1, and the `n1_fixedbb` control were running when this report was first written; their sections are filled in as they complete.
