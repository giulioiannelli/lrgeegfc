---
name: 2026-09-03_lane-s-contrast-scale-structure
kind: report
era: PAPER_FINALIZATION (Wave 0, lane W0-S, part B)
status: current
created: 2026-09-03
updated: 2026-09-03
scope: Part B of lane W0-S. Part A tested scale-locality of one functional at a time; this tests a CONTRAST between functionals (T_learn - T_test, and T_infspec alone), which the Part-A verdict did not cover, since two flat quantities can have a non-flat difference. Also settles the low_gamma band-selectivity near-miss that Part A tabulated and never discussed. Both verdicts are NEGATIVE and that is the headline.
pointers:
  - .agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md
  - .agents/reports/2026-08-31_lane-s-scale-variability.md
  - .agents/reports/2026-08-25_w0a-substrate-contract.md
  - scripts/01_compute/paper_final/w0s_07_contrast_scale_structure.py
  - scripts/01_compute/paper_final/w0s_08_band_separation_audit.py
  - scripts/01_compute/paper_final/w0s_09_cohort_median_slope_null.py
  - data/paper_final/lane_s_scale/contrast/
---

# W0-S part B — contrasts are flat too, and the band-selectivity claim is not currently supported

## Head

**The contrast is flat as well, and this is the cleaner, more final negative the coordinator asked for.** The sign change in `T_learn − T_test` along the scale axis is an artefact of reading a slope off a cohort-median curve: the same slope statistic, computed on pure surrogate data with no task information at all, reaches `|rho| = 0.69–0.82` at its 95th percentile, so the quoted `rho = +0.87 / −0.83 / +0.64` sit **inside the noise band of the statistic that produced them**. Recomputed per patient and referenced to each object's own noise floor, **nothing passes**: the contrast is worth 1.50–2.84 independent tests against its own floor of 2.97–3.68 (`q = 1.000` in all six bands), no band shows the two clusters a reversal requires, and per-patient zero-crossing scales are no more concentrated than chance. **Separately: `low_gamma`'s near-miss is a gate artefact, but the band-selectivity claim is still not supported** — `low_gamma`'s own clearing collapses when one patient is dropped and exists only at the sparse end of the plateau, yet the direct paired test of `beta − low_gamma` is null in *both* directions, so at n = 10 we can neither establish `low_gamma`'s trace nor demonstrate that `beta` differs from it.

**Three limitations first, per project rule.** Nothing here is verified against a timeseries-level null — matched-strength shuffles the finished connectivity matrix, so W0-B's ladder remains a live dependency for every number below. Lane E's grid samples `s ∈ [1, 180]` only, so this part speaks for the locked contract axis and not for the `s < 1` region Part A covered. And the ordered sham exists for four bands (`delta`, `theta`, `alpha`, `beta`) and supplies only **two** cohort-level realizations, one per source recording, so at cohort level it can support a magnitude comparison but never a p-value; the per-patient pairing described in §4 is what makes it a test at all.

**One correction to my own working diagnostic, made before publishing it.** My first pass reported a "2.5× inflation" of the cohort-median slope over the per-patient slope. That number was wrong — it compared an absolute value against a *signed* mean, and signed per-patient slopes partly cancel. Recomputed against the mean *absolute* per-patient slope the ratio is 0.85, and it is 0.87 under the null, i.e. **the median-smoothing story is not the mechanism**. The real mechanism is §3 and it is a different and simpler one.

---

## 1. The four answers, in one table

| Question | Answer | The number that decides it |
|---|---|---|
| Does `T_learn − T_test` carry more independent scale information than its own noise? | **No** | `n_eff` 1.50–2.84 vs own noise floor 2.97–3.68; `p = 0.66–1.00`, `q = 1.000` in every band |
| Does it genuinely reverse sign along `tau`? | **No** | 0 of 36 (object x band) cells show both a negative and a positive cluster; crossing scales no tighter than null (`p >= 0.44`) |
| Does `T_infspec` have a real scale slope? | **No** | cohort-median `rho = −0.83` in delta, but empirical `p = 0.035`, `q = 0.209`; 0 of 30 cells clear |
| Is `low_gamma` a genuine near-miss or a gate artefact? | **Gate artefact — and the bands still cannot be separated** | LOO `p` rises 0.045 → 0.121; clears only at `f <= 0.10`; paired `beta − low_gamma` null both ways (`p = 1.000`) |

---

## 2. What was pre-registered, and when

Amendment B was appended to the scope report `.agents/guides/task-persistence-investigation/2026-08-31_scale-local-trace-readouts.md` and **committed before the first contrast statistic was computed** (commit `3a0e8ee`; the first analysis commit is `61afd38`). It fixed the paired-contrast definition, the two nulls, three criteria B1–B3, the reporting rule, and — in §12.5 — the `low_gamma` call, written down before it was known which answer would be convenient.

The criteria, verbatim in effect: **B1**, the contrast must exceed the 95th percentile of *its own* held-out null `n_eff`; **B2**, a reversal must be tested as a reversal (a negative cluster and a positive cluster, in the right order along `s`, with per-patient crossings more concentrated than both nulls); **B3**, every statistic reported as a pair `(p_vs_surrogate, p_vs_sham)`, with a statistic that clears the first and not the second declared to be within-recording drift. All three were required. **None of the six objects passes any pairing of them, and none passes B1 or B2 at all.**

---

## 3. What actually killed the lead: a cohort-median slope has almost no degrees of freedom

This is the transferable result of part B and it applies well beyond this lane.

The quoted evidence for a scale-dependent contrast was `rho(log s)` computed on a **cohort-median profile** — median the margin across the ten patients at each scale, then correlate the resulting 16-point curve with `log s`. Rather than argue about the correct degrees of freedom, I built the null distribution of exactly that statistic: promote a held-out matched-strength realization to the observed slot, form the same cohort-median profile, take the same Spearman, 200 times per cell.

**Under that null the cohort-median `rho` has a median absolute value of 0.27–0.40 and a 95th percentile of 0.69–0.82.** A `|rho|` of 0.8 on this axis is roughly a one-in-twenty event on data containing no task information whatsoever. The consequence for the specific leads, with the naive Spearman p that a 16-point correlation would advertise beside the empirical p from this null:

| Object | Band | `rho` (cohort median) | naive p | empirical p | BH q | p inflated by |
|---|---|---|---|---|---|---|
| `T_learn` | beta | +0.841 | 4.4e-05 | 0.055 | 0.313 | **1 200×** |
| `T_infspec` | delta | −0.832 | 6.3e-05 | 0.035 | 0.209 | **550×** |
| `T_infspec` | alpha | −0.635 | 8.2e-03 | 0.199 | 0.597 | 24× |
| `C = T_learn − T_test` | delta | +0.597 | 1.5e-02 | 0.234 | 0.507 | 16× |
| `C` | beta | +0.453 | 7.8e-02 | 0.423 | 0.507 | 5× |
| `C` | theta | +0.429 | 9.7e-02 | 0.378 | 0.507 | 4× |

**0 of 30 (object × band) cells clear at `q < 0.05`; the smallest q anywhere is 0.060.** The coordinator's estimate that the implied p was "inflated by an order of magnitude" was correct in direction and conservative by two orders in the cells that matter most — precisely because the inflation grows with the apparent strength of the correlation, so the most persuasive-looking numbers are the most inflated ones.

The mechanism is not that the median smooths the data (I checked, and the smoothing ratio is the same under the null — see the correction in the Head). It is that Part A already measured this axis as worth about one independent test, so a 16-point Spearman on it has effectively no degrees of freedom, and its null distribution is nearly as wide as `[−1, +1]` allows.

---

## 4. The contrast itself, per patient, per scale

Formed the honest way: the difference is taken at raw level, before any surrogate subtraction, so that the surrogate of the difference is the difference of the surrogates (Lane E's realization index is shared across the functional axis, so the pairing survives) and the locked `patient_margin` contract stays exact. The unpaired variant `margin[T_learn] − margin[T_test]` is carried alongside as a sensitivity check and agrees throughout.

**B1 — independent information: fails everywhere, and not marginally.**

| Object | `n_eff` observed | own held-out noise floor | `p` | `q` |
|---|---|---|---|---|
| `T_test` | 1.06–1.26 | 2.87–3.75 | 1.000 | 1.000 |
| `T_learn` | 1.07–1.30 | 2.96–3.59 | 1.000 | 1.000 |
| `T_infspec` | 1.46–2.75 | 2.87–4.31 | 0.58–1.00 | 1.000 |
| `C` (paired) | 1.50–2.84 | 2.97–3.68 | 0.66–1.00 | 1.000 |
| `C` (unpaired) | 1.54–2.66 | 2.97–3.68 | 0.71–1.00 | 1.000 |

The contrast **does** decorrelate its scales more than the raw functionals do — 1.50–2.84 against `T_test`'s 1.06–1.26, a genuine increase of roughly two-fold. Quoted alone that reads like the hypothesis confirmed. It is not: the contrast's own noise floor is 2.97–3.68, so the increase is smaller than what noise alone delivers. This is the same trap that killed four of Part A's five candidates, and it is why criterion (b) is null-referenced rather than absolute.

**B2 — reversal: fails everywhere.** Zero bands show both a supra-threshold negative cluster and a supra-threshold positive cluster. In `beta`, the band with the largest quoted contrast slope, there is a significant *negative* cluster (`p = 0.035`) and **no supra-threshold positive cluster anywhere** (`p = 1.000`, cluster mass exactly 0): the cohort median does creep from `-0.016` at `s = 1` to `+0.007` at `s = 180`, but that final positive excursion never reaches the cluster-forming threshold, so there is no positive cluster for the negative one to pair with. Per-patient zero crossings are found in 4–7 of 10 patients against a null median of 5–7, and their scatter in `log10 s` (MAD 0.30–0.94) is *wider* than the null's (0.25–0.36) in every band. Patients do not agree on where a crossing happens because, on this evidence, there is no crossing to agree about.

**The noise-inflation caveat does not apply here, and I measured that rather than assuming it either way.** Lane E was right that the contrast is not simply a sum of two independent errors: `T_learn` and `T_test` share the vectors they are built from, and their surrogate ensembles correlate at `rho = 0.33–0.56` per band. The consequence, measured on the same cells: the contrast's null SD is **0.52–1.08×** its inputs' — at most 8 % noisier than either part in any band, and roughly half the noise in the gamma bands. **That strengthens the negative rather than weakening it.** The contrast is a cleaner instrument than either of its parts and still finds nothing.

---

## 5. The ordered sham changes what "clears matched-strength" is worth

This is the one place where a statistic *did* clear the matched-strength null, and the sham is what stops it being reported as a result.

The direction-agnostic scale-dependence statistic `slope_abs` — the cohort mean of each patient's own `|Spearman(margin, log s)|` — exceeds its matched-strength floor in several cells (`T_test` delta `p = 0.020`, beta `p = 0.030`; `T_infspec` delta `p = 0.030`; the contrast in delta, low_gamma and high_gamma at `p = 0.0099`). Read against matched-strength alone, that is a scale-structure finding.

Read against the ordered sham it is not. The sham's `slope_abs` is **0.38–0.70**, systematically *above* the matched-strength null's 0.37–0.45 and squarely inside the real data's 0.32–0.80. In plain terms: **matched-strength surrogates understate how much a real recording's scale profile wobbles when nothing is happening.** A real graph carries community structure that makes the diffusion cascade lumpy along `tau`; a strength-matched shuffle does not. Every cell above stops clearing when the comparison is made against the same patient's own no-task, drift-carrying profile.

This is a methodological finding with reach beyond part B, and it is the concrete justification for criterion B3: **matched-strength is the minimum null, not a sufficient one, for any claim about the *shape* of a curve along a swept axis.** For a claim about the *height* of the margin it remains the locked standard; for a claim about how that height varies with `tau`, it is too generous.

---

## 6. The one live thread, which belongs to Lane E

`T_learn` is the single object whose scale-dependence exceeds the ordered sham. Paired per patient — is patient k's real profile more scale-dependent than patient k's own sham profile? — it gives `p = 0.0020` (delta, rest_pre), `0.0137` (delta, rest_post), `0.0244` and `0.0244` (beta, both sources). Those are the **top four of 40** comparisons in that family, and the two bands each replicate across two independent source recordings.

**It is not a result and I am not reporting it as one.** BH over the 40-test family gives a minimum `q = 0.078`; nothing clears. And `T_learn` fails B1 hard (`n_eff` 1.07–1.30 against a floor of 2.96–3.59), which has a coherent reading rather than being a contradiction: a profile that is strongly *monotone* is maximally *redundant*, because a monotone curve over 16 points is essentially one number — its trend — not sixteen. So the most that can be said is that `T_learn`'s margin may carry a single smooth trend along `tau` that the sham does not reproduce.

That is an encoding-side claim about `T_learn`, not a scale-structure claim about the contrast, so it belongs to Lane E. Handing it over with the pre-registration it would need: **one** directional test of `slope_signed` for `T_learn` in delta and beta against the ordered sham, direction fixed in advance from these numbers, on the re-derived tensor Lane E is producing.

---

## 7. low_gamma: a gate artefact, and the band-selectivity claim is still unsupported

The call was pre-registered in §12.5 before any of these numbers existed. It resolves to **gate artefact**, but the fuller answer is worse for the paper than that label alone suggests, so both halves are stated.

**Why `low_gamma`'s own clearing is an artefact — two independent grounds.**

1. **It does not survive leave-one-patient-out.** `low_gamma`'s axis-cluster `p = 0.045` rises to `0.121` when the most influential single patient (Pat_02) is dropped, flipping the verdict in 6 of 10 drops. It is the only clearing band that does this: `alpha` `p_max = 0.030`, `delta` `0.042`, `beta` `0.044`, all with 0 flips. Calibration is not the problem — `low_gamma` is calibrated (`FPR` at nominal 0.05: median 0.04, worst scale 0.10) and 7 of 10 patients are positive, so the pre-registered "minority-carried" and "uncalibrated" routes to the artefact verdict are both **not** what triggered it.
2. **It exists only at the sparse end of the plateau, which independently reproduces W0-A.** Resolving the knob instead of integrating it: `low_gamma`'s cluster `p` is `0.0073` at `f = 0.07`, `0.0234` at `0.10`, then `0.141` at `0.14` and `0.125` at `0.20`, with its median margin collapsing `+0.147 → +0.174 → +0.042 → −0.033`. W0-A reached the same conclusion from a different statistic on the same substrate ("low_γ leaks at f = 0.07 and f = 0.10 and only goes null at f ≥ 0.14"). Two lanes, two statistics, one answer.

**Why the band-selectivity claim is nevertheless not supported.** Reading two marginal q-values side by side (`beta` 0.042, `low_gamma` 0.067) is not a test of whether the bands differ. Tested directly — per patient, per scale, `margin_beta − margin_low_gamma` through the locked sign-flip cluster gate — the answer is null **in both directions**: `p = 1.000` for `beta > low_gamma` and `p = 1.000` for `low_gamma > beta`, with no supra-threshold cluster of either sign anywhere on the axis, per-scale minimum `q = 0.59` and `0.88`, and a median cohort difference of 0.026 with 4 of 10 patients on the `beta` side. Their knob-integrated median margins are `+0.092` (`beta`) and `+0.095` (`low_gamma`) — `low_gamma`'s is the *larger* of the two.

I also tested the one thing that visibly separates them, and it too fails at n = 10. `beta`'s margin *rises* with the sparsification fraction (per-patient slope `+0.022` per log-unit of `f`) while `low_gamma`'s *falls* (`−0.064`), which is why integrating over the knob averages the difference away. Paired across patients that difference gives `p = 0.278` with 6 of 10 patients in the expected direction, and the paired cluster gate run separately at each fraction gives a minimum `q = 0.459`. (Exploratory, and labelled as such: it composes the two pre-registered tests rather than adding a third object.)

**The honest summary.** At n = 10, on this substrate, no test I ran can separate `beta` from `low_gamma`. That does *not* say `low_gamma` has a trace — its own clearing fails LOO and is confined to `f <= 0.10`. It does *not* impugn `beta` — `beta` clears, survives LOO with 0 flips, and is calibrated. It says the claim "`beta` yes, `low_gamma` no" currently rests on two marginal q-values that were never compared, and when they are compared the comparison is null. **Any paper text asserting band-selectivity against `low_gamma` must either carry this non-separation and W0-A's 0.51-octave window, or be softened to a statement about `beta` alone.**

---

## 8. Every construction tried, including the failures

Six objects, all reported: `T_test`, `T_learn`, `T_infspec`, `T_infspec_pe`, the paired contrast `C = T_learn − T_test`, and the unpaired `margin[T_learn] − margin[T_test]`. Seven scale-structure statistics on each, all reported: `n_eff` (participation ratio of the cross-scale correlation matrix), `slope_abs`, `slope_signed`, `shape_agreement` (cross-patient agreement on profile shape), the axis-cluster gate in both signs, the count of per-patient zero crossings, and their concentration. Two nulls on each, both reported: held-out matched-strength realizations of the same object, and the ordered sham paired per patient.

Nothing was dropped for being inconvenient, and two things were kept specifically because they were: the `slope_abs` cells that clear matched-strength (§5, killed by the sham) and `T_learn`'s sham excess (§6, killed by multiplicity and handed on anyway).

One library helper was added rather than inlined, with the trap it avoids written into its docstring: `cohort_gate.axis_reversal_gate` tests a reversal as a reversal, and documents why the sign-flip null that licenses `axis_cluster_gate` is **invalid** for a crossing location — flipping a patient's whole profile negates both the intercept and the slope and leaves the crossing `exp(−a/b)` *exactly* unchanged, so a sign-flip null there would return `p = 1` by construction while looking like a fair test. The null has to come from draws that break cross-patient alignment instead, which is what the held-out realizations and the sham provide.

---

## 9. What I could not settle

- **No timeseries-level null.** Every number is conditioned on "given this FC matrix". W0-B's ladder is the live dependency, unchanged from part A.
- **`s < 1` is unsampled** on Lane E's grid, so part B speaks only for the locked contract axis `s ∈ [1, 180]`. Part A covered `s ∈ [0.02, 180]` for the incumbent and found the flatness extends across all of it, which makes it unlikely but not impossible that a contrast behaves differently below `s = 1`.
- **Two bands have no sham.** `low_gamma` and `high_gamma` have no ordered-sham cells, so in those bands only the matched-strength null applies — and §5 shows that null is too generous for shape claims. Lane E is computing the missing two; the §5 conclusion should be re-checked there when they land.
- **Power.** With `K = 10` and an axis worth roughly one independent test, "flat" means *not detectable at n = 10*, never *proven constant*. The `beta`/`low_gamma` non-separation in §7 is the sharpest instance: it is a statement about the instrument as much as about the bands.
- **These cells are about to be superseded.** Lane E is re-deriving all four functionals from an 8-vector rank-correlation tensor on the same substrate, seeds and scale grid. Nothing here depends on the collapsed 4-functional storage, so the numbers should reproduce exactly; §5 and §7 are the two places worth re-running first.

---

## 10. Reproduce

```
bash scripts/01_compute/paper_final/run_py.sh scripts/01_compute/paper_final/w0s_07_contrast_scale_structure.py
bash scripts/01_compute/paper_final/run_py.sh scripts/01_compute/paper_final/w0s_08_band_separation_audit.py
bash scripts/01_compute/paper_final/run_py.sh scripts/01_compute/paper_final/w0s_09_cohort_median_slope_null.py
```

Artifacts, all under `data/paper_final/lane_s_scale/`: `contrast/contrast_scale_structure.csv` (every object × band × statistic × both nulls), `contrast/contrast_criteria.csv` (B1–B3 verdicts), `contrast/contrast_noise_levels.csv` (§4's SD ratios and null correlations), `contrast/contrast_per_scale_profiles.csv`, `contrast/contrast_gate_grid.csv`, `contrast/cohort_median_slope_null.csv` (§3), `band_separation_bands.csv`, `band_separation_verdict.csv`, `band_separation_by_fraction.csv`, `band_separation_paired_per_scale.csv`, `band_separation_paired_by_fraction.csv`, `band_separation_knob_slopes.csv` (§7).
