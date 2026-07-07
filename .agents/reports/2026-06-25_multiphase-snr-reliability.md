---
name: multiphase-snr-reliability
type: report
era: "IMCOH_ABS × COHORT_N10"
status: current
created: 2026-06-25
updated: 2026-06-25
pointers:
  - scripts/01_compute/audit/audit_145_multiphase_snr.py
  - data/audit/multiphase_snr/per_patient_per_band_reliability.csv
  - data/audit/multiphase_snr/per_band_summary.csv
  - data/audit/multiphase_snr/patient_ranking.csv
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
---

# Full per-phase reliability refinement of the detectability ("SNR") axis

**HEAD (the juice).** Correcting the cophenetic trace `ρ_split` for the
reliability of **all four** phase-distance vectors (not just rest_pre)
does **not** make the between-patient heterogeneity go away — so the
heterogeneity is **NOT entirely measurement reliability; a residual
cross-patient axis survives**. The PI's objection was half right and
half wrong: rest_post really is the dominant noise source the rest_pre-only
floor was blind to (it is the least-reliable phase in **38/60**
patient×band cells, median split-half ρ ≈ 0.42 vs rest_pre ≈ 0.56), so the
refined SNR is a better instrument; but after disattenuation the patient-level
trace ordering is essentially **unchanged**, the reliability composite
explains only **~6% of the pooled variance** in `ρ_split`
(R² = 0.064, residual-variance fraction 0.94), and the single most
striking datum is that **Pat_15 — the highest-reliability patient
(mean composite 0.82) — is the lowest tracer (mean ρ_split ≈ −0.11)**, the
exact opposite of the pure-reliability prediction. Verdict:
**universal-trace-detectability-limited is REJECTED as the complete story;
a real biological axis of who-traces survives reliability correction.**
The within-band SNR↔ρ correlation stays positive (reliability matters
*within* a band), but it cannot explain *which patients trace across many
bands*. Honest caveat up front: at n=10 every per-band CI is wide; the
verdict rests on point estimates + the cross-band patient axis, not on a
powered test.

---

## 1. What was refined and why

The session-current detectability axis is
`snr = d_task / d_noise` with `d_noise = 1 − Spearman(D_preA, D_preB)` —
**rest_pre split-half only**. `ρ_split = Spearman(Δ_task, Δ_rest)` where
`Δ_task = D_tt − D_preA` and `Δ_rest = D_post − D_preB`, so its attenuation
is governed by the reliability of **four** cophenetic vectors
(D_preA, D_preB, D_tt, D_post). A patient can look like a non-tracer purely
because **task_test or rest_post is noisy** — which `d_noise` never sees.

Classical correction-for-attenuation:
`ρ_obs ≈ ρ_true · √(rel(Δ_task)·rel(Δ_rest))`. We measured each phase's
split-half cophenetic reliability and propagated it to the two difference
vectors, then disattenuated.

## 2. Data added (task_test split-half FCs)

rest_pre and rest_post halves were already cached (audit_63 et al.).
The missing **task_test** half-FCs were computed the same way
(`compute_imcoh_abs_halves`, Welch `nperseg = max(256, nperseg_for_fs(fs)//2)`)
and cached at
`data/cache/imcoh_halves_fc/{pat}/{band}_task_test_{A|B}_imcoh_abs.npy`
for all 10 patients × 6 bands (120 files).

**Pat_10 quirk — resolved + flagged.** The brief said "apply the Pat_10
[53,54,55] row drop." In fact `load_timeseries('Pat_10','task_test',…)`
**already returns 113 channels with those rows removed** (it equals the
canonical full FC N=113). Re-dropping gave N=110 — caught on the first run
by a hard N-vs-full-FC guard (`ValueError`), which is exactly what the
guard is for. The manual drop was removed; the guard stays unconditional.
Every cached task_test half-FC N now equals the corresponding full FC N for
all 10 patients (verified in-script before writing).

## 3. The reliability estimator (chosen + justified)

- **Per-phase reliability**: `rel_phase = Spearman(canonical_cophenet(W_A),
  canonical_cophenet(W_B))` for phase ∈ {rest_pre, task_test, rest_post},
  using the **canonical** LRG cophenetic (`canonical_cophenet`, the exact
  function the locked analyses and matched-strength surrogate use).
- **Difference-vector reliability**: under the classical true-score model
  `X = T + E`, the reliability of a difference `D_a − D_b` is
  `(rel_a + rel_b − 2 r_TT) / (2 − 2 r_EE)`. We have no third repeat, so we
  cannot estimate the cross-phase true-score correlation `r_TT` or error
  correlation `r_EE`. We adopt the **equal-variance, independent-error,
  orthogonal-true-score** approximation (`r_TT = r_EE = 0`):
  **`rel(Δ) ≈ (rel_a + rel_b)/2`** — the standard "difference of two equally
  reliable, uncorrelated measurements" result. This is deliberately
  **conservative**: difference scores usually have *positively* correlated
  true scores (`r_TT > 0`), which would make the true `rel(Δ)` *higher* and
  the disattenuation *larger*; our choice therefore **under-corrects** and
  biases **against** finding residual biology — the safe direction for our
  verdict.
- **Composite + disattenuation**:
  `R_comp = √(rel(Δ_task)·rel(Δ_rest))`,
  `ρ_corrected = ρ_obs / R_comp` (guarded: composite < 0.05 → reported as
  capped/NaN, never divided through ~0).
- **Spearman-Brown** full-length versions of every reliability are also
  reported (`*_sb`), since split-half halves the data; they tell the same
  story.

`ρ_obs` is read from the cached headline
(`matched_strength_surrogate_split_baseline/per_patient_per_band.csv`) — the
trace itself is **not recomputed**.

## 4. Per-phase reliability — the PI's objection, confirmed in part

Median split-half cophenetic reliability across all 60 cells:

| phase | median rel (half) | role |
|---|---|---|
| rest_pre | 0.555 | the only phase the OLD `d_noise` saw |
| task_test | 0.585 | invisible to old `d_noise` |
| **rest_post** | **0.417** | **invisible to old `d_noise`, and the WORST** |

**rest_post is the reliability bottleneck in 38/60 cells** (rest_pre 18,
task_test 4). The rest_pre-only floor was systematically blind to the
**dominant** source of attenuation. So the refinement is a genuinely better
instrument — the PI was right that the old axis under-measured noise.

## 5. THE deliverable — is the heterogeneity entirely reliability?

**No. A residual survives.** Three independent cuts agree:

1. **Reliability explains almost none of the pooled trace variance.**
   Pooled over 60 cells, `Spearman(ρ_split, R_comp) = +0.199`,
   OLS `R² = 0.064`, **residual-variance fraction = 0.936**. Per band the
   reliability-only `R²` point estimates are tiny — β **0.0006**, γ_l 0.037,
   γ_h 0.033, α 0.048, δ 0.123, θ 0.186 — with bootstrap CIs running from
   ~0 to ~0.95 (uninformative at n=10, but the centre of mass is near zero,
   not near one).

2. **The cross-band patient "who-traces" axis is orthogonal to reliability.**
   Per-patient means across the 6 bands:

   | patient | mean ρ_split | mean ρ_corrected | mean R_comp |
   |---|---|---|---|
   | Pat_06 | +0.620 | +0.885 | 0.688 |
   | Pat_05 | +0.412 | +0.767 | 0.433 |
   | Pat_08 | +0.219 | +0.525 | 0.456 |
   | Pat_03 | +0.256 | +0.477 | 0.476 |
   | Pat_02 | +0.246 | +0.186 | 0.532 |
   | Pat_07 | +0.060 | +0.133 | 0.585 |
   | Pat_13 | −0.008 | −0.017 | 0.423 |
   | Pat_10 | −0.007 | −0.037 | 0.424 |
   | Pat_14 | −0.039 | −0.117 | 0.466 |
   | Pat_15 | −0.105 | −0.140 | **0.816** |

   `Spearman(patient-mean ρ_corrected, patient-mean R_comp) = +0.006`;
   for observed, `+0.079`. **Flat.** The patient axis of trace strength is
   **not** a reliability axis. **Pat_15 is the highest-reliability patient
   and the lowest tracer** — the cleanest possible falsifier of "non-tracers
   are just noisy." Pat_05 is the mirror image: a strong tracer on
   *below-median* reliability. Disattenuation *widens* the top of the
   distribution (Pat_05/06 rise) without rescuing any low tracer.

3. **The β trace ordering is intact after correction.** Pat_08/05/02/03 stay
   the top-4 β tracers (corrected ρ 1.18/0.80/0.77/0.73 — the 1.18 is a
   reliability-floor artifact, see §7), Pat_14/10 stay the only negatives
   (−0.11/−0.17). Correction rescales, it does not reorder.

**Verdict: RESIDUAL BIOLOGY.** The bare trace is not "universal,
detectability-limited." Reliability is a real *within-band* modulator (high
SNR helps β clear, see below) but a **biological who-traces axis survives**
full per-phase reliability correction. This is consistent with — and
strengthens — the detectability-as-one-factor reading in the SNR/band-taxonomy
handoff, but it forecloses the strong "entirely reliability" claim.

## 6. Corrected SNR ↔ ρ_split, per band

Within a band the detectability story holds (more reliable patient ⇒ larger
observed trace), but the correlation **weakens** under the corrected SNR
(because correcting partly *removes* the reliability that drove it):

| band | corrected ρ(SNR, ρ_split) [95% CI] | current ρ(SNR, ρ_split) [95% CI] |
|---|---|---|
| delta | **+0.770** [+0.170, +0.988] | +0.867 [+0.371, +1.000] |
| theta | +0.224 [−0.558, +0.847] | +0.370 [−0.444, +0.923] |
| alpha | +0.236 [−0.456, +0.752] | +0.333 [−0.488, +0.896] |
| **beta** | **+0.491** [−0.275, +0.961] | +0.782 [+0.200, +0.987] |
| low_gamma | +0.442 [−0.359, +0.883] | +0.721 [+0.019, +1.000] |
| high_gamma | +0.733 [+0.149, +0.987] | +0.964 [+0.779, +1.000] |

"Corrected SNR" here is `d_task / d_noise_full` with
`d_noise_full = 1 − min(rel_pre, rel_tt, rel_post)` (the **bottleneck**
phase, not rest_pre). The drop from current→corrected (β +0.78→+0.49) is the
expected signature of the SNR↔ρ link being **partly** a reliability
tautology — but only partly: δ and γ_h stay strongly positive even after
correction, and none of the within-band correlations is the cross-band
who-traces axis (which is flat, §5).

## 7. Patient re-ranking (current vs corrected SNR)

The current high-SNR set {06, 05, 02, 03, 08} is **mostly preserved but
reshuffled** once rest_post/task_test reliability enters. β example (sorted
by corrected rank):

| patient | snr_current | snr_corrected | rank_current → corrected |
|---|---|---|---|
| Pat_02 | 2.08 | 1.47 | 2 → **1** |
| Pat_03 | 1.71 | 1.38 | 4 → **2** |
| Pat_10 | 1.34 | 1.34 | 6 → **3** |
| Pat_05 | 2.85 | 1.29 | **1 → 4** |
| Pat_08 | 1.94 | 1.25 | 3 → 5 |
| Pat_06 | 1.40 | 0.89 | 5 → 6 |
| Pat_15 | 1.09 | 0.60 | 8 → 10 |

Key moves: **Pat_05 falls 1→4** (its rest_pre was unusually clean — `d_noise`
0.247 — but its bottleneck phase is much noisier, 0.547, so the old SNR
*over-rated* it); **Pat_10 rises 6→3** (clean across phases, was penalised by
nothing in particular). The corrected top-5 per band vs the current top-5
differ by 2–4 patients in every band (e.g. β swaps P06↔P10; α swaps in
P05/P08/P10 for P06/P14). So **the ranking does change** — the rest_pre-only
SNR mis-ranks several patients — but the *set* of strong tracers is stable
because the residual biology (§5) dominates the ranking, not the reliability
correction.

## 8. Honesty — what this estimator can and cannot separate

**Can:** phase-specific unreliability *outside* rest_pre now enters the
attenuation composite (task_test and rest_post reliabilities are measured
from their own halves), which is the failure mode the PI named and the old
`d_noise` was blind to.

**Cannot / caveats (load-bearing):**
- **Split-half halves the data**, so every `rel_phase` is a downward-biased,
  noisy estimate of full-length reliability (fewer Welch segments). The
  *absolute* disattenuation magnitude is model-dependent; Spearman-Brown
  full-length values are reported but rest on the parallel-halves assumption.
  The **rank** verdict (who-traces axis ⟂ reliability) is robust to this; the
  *level* of `ρ_corrected` is not.
- **Spearman on a cophenetic vector that is ~98% tied** at the N−1 merge
  heights is depressed by ties, so `rel_phase` is **conservatively low**,
  which **inflates** `ρ_corrected`. Corrected ρ > 1 (e.g. β Pat_08 = 1.18,
  γ_h Pat_02 = −1.47 from a near-zero composite) are **reliability-floor
  artifacts**, not super-reliability — read them as "ceiling/edge," not as
  literal disattenuated correlations. They do not affect the residual-biology
  verdict, which is driven by the flat patient-mean correlation.
- The `rel(Δ) ≈ (rel_a + rel_b)/2` propagation assumes **equal signal
  variance and independent error** across the two phases. Task and rest share
  electrodes/montage; if their measurement errors are **correlated**, the
  difference reliability is mis-estimated and the estimator **cannot separate
  correlated-error from shared true signal**. Our orthogonal-true-score choice
  under-corrects, biasing against residual biology — so the residual we *do*
  find is, if anything, an underestimate.
- **3/60 cells have a negative phase reliability** (all high_gamma:
  Pat_02/10/14 — the noisiest band). Two yield NaN composites (correctly
  excluded); one (Pat_02 γ_h) yields a small positive composite and an
  out-of-range ρ_corrected, flagged above. High-γ disattenuation should be
  treated as unreliable.
- **n = 10.** Every per-band CI is wide and no per-band test is powered. The
  verdict rests on (a) tiny pooled `R²`, (b) the *flat* cross-band
  patient-mean correlation with a clean falsifier (Pat_15), and (c) the
  preserved trace ordering — a convergent, not a single-test, argument.

## 9. Files

- `scripts/01_compute/audit/audit_145_multiphase_snr.py` — compute (5-point
  preamble docstring; caches task_test halves with N-guard; reliabilities;
  disattenuation; residual regression + bootstraps; ranking).
- `data/audit/multiphase_snr/per_patient_per_band_reliability.csv` — 60 rows:
  per-phase rels (half + Spearman-Brown), difference-vector rels, composite,
  ρ_corrected (+ capped flags), bottleneck noise.
- `data/audit/multiphase_snr/per_band_summary.csv` — per-band reliability-only
  R² (+ CI), residual-variance fraction, SNR↔ρ correlation, capped counts.
- `data/audit/multiphase_snr/patient_ranking.csv` — per-(patient,band)
  current vs corrected SNR and ranks.
- `data/cache/imcoh_halves_fc/{pat}/{band}_task_test_{A|B}_imcoh_abs.npy` —
  NEW task_test half-FCs (120 files).
