---
date: 2026-05-11
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-batch
scope: §5.2 KC null verification + §5.3/§6 cohort heterogeneity vs implant geometry
inputs:
  - data/audit/implant_geometry/per_patient_features.csv
  - data/audit/implant_geometry/correlation_matrix.csv
  - data/audit/implant_geometry/anti_alignment_count.csv
  - data/audit/kc_matched_strength_surrogate/cohort_summary.csv
  - data/audit/kc_matched_strength_surrogate/per_patient_per_band.csv
  - data/audit/kc_matched_strength_surrogate/joint_signature.csv
scripts:
  - scripts/01_compute/audit/audit_64_implant_geometry_test.py
  - scripts/01_compute/audit/audit_65_kc_matched_strength_surrogate.py
---

# 2026-05-11 — implant geometry test + KC matched-strength null verification

**Head.** Two parallel bias / null-verification investigations.
`audit_64` (implant geometry, complete) tests whether per-patient
trace-direction heterogeneity at α / β / γ_l correlates with implant
features. β returns moderate uncorrected association (B_hemi vs obs_ρ
Spearman ρ = +0.697, p = 0.025, BH-FDR q = 0.994 fails); α / γ_l null.
`audit_65` (KC tree-distance matched-strength surrogate, **complete**)
puts §5.2 KC λ ∈ {0, 1} T_KC under the same R=200 4-cycle ±δ surrogate
as audit_63. **Major null result**: NO band-λ cell achieves "separated"
verdict; α λ=0 lands as "also_negative" (surrogate median -0.122 more
negative than observed -0.015); β λ=0 surrogate gives ~50% of observed
trace negativity; only 1-2/10 patients per cell pass individual
lower-tail at p<0.05; cohort-paired Wilcoxon p ranges 0.22–0.72 across
all six cells. The §5.2 KC 10/10 cohort claim against the within-baseline
null does **not** robustify against the wider matched-strength null.

The reading: the §5 / §6 manuscript needs to reframe the §5.2 KC β
result as **narrow-null detection of a tree-distance reorganization
that is recoverable from matched-strength edge-identity randomization
at the cohort level**. The within-baseline KC result is not invalidated
(the per-patient tree-distance does shift more under task than under
half-baseline drift), but matched-strength noise propagates similar
shifts through the average-linkage construction. The §5.3 ρ_split
result (audit_63: β cohort Wilcoxon p=0.005, α p=0.002) becomes the
load-bearing LRG-layer claim; §5.2 KC becomes a scale-resolved
companion that confirms reorganization geometry without isolating a
shape-driven signature beyond per-node strength heterogeneity.

## Companion to: 2026-05-10 matched-strength surrogate batch

This file pairs with `2026-05-10_matched-strength-surrogate-batch.md`,
which covers audit_61 (epi-fraction) + audit_62 (shared-baseline ρ) +
audit_63 (split-baseline ρ_split). Read the May-10 file first for the
surrogate algorithm context and the §5.3 verdicts. The two files
together brief the writing agent on every probe of the §5.3 / §5.2
trace claim landed since 2026-05-08.

## audit_64 — implant geometry test (complete)

**Question.** Does per-patient implant geometry account for the cohort
heterogeneity in trace direction observed at α / β / γ_l?

**Method.** Pure correlational analysis on pre-existing per-patient
outputs (no LRG, no FC, no surrogate). For each patient, extract
9 implant features from `implant_pat_NN.csv`:

- N (channel count), N_L / N_R (per-hemisphere counts),
  B_hemi = (N_L − N_R) / (N_L + N_R) (hemispheric balance),
  N_temporal / N_frontal (lobar counts),
  frac_lobe_max (max lobar fraction; coverage concentration proxy),
  sigma_disp (spatial dispersion),
  frac_epi (epileptogenic-zone contact fraction).

For each band b ∈ {α, β, γ_l}, extract 5 per-patient trace measures
(load-bearing scalars from §5.2 / §5.3 / §4):

- obs_ρ (split-baseline ρ_split from audit_63)
- obs_z (split-baseline z from audit_63)
- T_KC λ=0 (from kc_triangle, sign-flipped so positive = trace)
- T_KC λ=1 (sign-flipped)
- T_d (from raw_fc_phase_distance, sign-flipped)

Spearman correlation across 9 features × 5 measures × 3 bands = 135
cells. BH-FDR over m = 135 at q = 0.05.

**Per-band verdicts:**

| band | top feature | Spearman ρ | p uncorrected | BH-FDR q (m=135) | verdict |
|---|---|---|---|---|---|
| α | frac_epi | -0.406 | 0.244 | 0.994 | **null** |
| β | B_hemi | +0.697 | **0.025** | 0.994 | **moderate** |
| γ_l | frac_lobe_max | +0.358 | 0.31 | 0.994 | **null** |

**β top three by |ρ|:**

| feature | Spearman ρ | p | BH-FDR q |
|---|---|---|---|
| B_hemi | +0.697 | 0.025 | 0.994 |
| N_R | -0.694 | 0.026 | 0.994 |
| N_L | +0.590 | 0.073 | 0.994 |

**Anti-alignment count per patient (across 12 = 4 probes × 3 bands):**

| patient | n_anti / 12 | class |
|---|---|---|
| Pat_02 | 2 | consistently_pro |
| Pat_03 | 2 | consistently_pro |
| Pat_05 | 2 | consistently_pro |
| Pat_06 | 0 | consistently_pro |
| Pat_07 | 6 | consistently_anti |
| Pat_08 | 1 | consistently_pro |
| Pat_10 | 3 | mixed |
| Pat_13 | 3 | mixed |
| Pat_14 | 7 | consistently_anti |
| Pat_15 | 10 | consistently_anti |

The three consistently-anti patients (Pat_07, Pat_14, Pat_15) are the
ones to watch in §6 cohort-heterogeneity prose. Pat_15 is anti at
10/12 trace measures.

**Reading.**

- 1/3 trace bands admit moderate implant-geometry correlation (β,
  uncorrected only). The β B_hemi association is directionally consistent
  with the hypothesis that left/right contact balance contributes to the
  per-patient trace direction, but does not survive m=135 BH correction.
- α / γ_l show no implant-geometry correlation at the load-bearing
  scalars.
- **The §6.3 outlook framing should be band-resolved**: for β, disclose
  the moderate uncorrected B_hemi correlation as a candidate
  attribution mechanism (with the BH-correction caveat); for α / γ_l,
  state the cohort heterogeneity is documented but the mechanism is
  not identified in this 10-patient sample.

## audit_65 — KC matched-strength surrogate (COMPLETE)

**Question.** The §5.2 β KC 10/10 cohort claim (q=0.006 at m=12, ten
patients individually below their own within-baseline split-half null
at both λ=0 and λ=1) currently rests on the *narrow* within-baseline
null. audit_63 returned an *intermediate* verdict for §5.3 ρ_split
against a wider matched-strength null. Does the §5.2 KC tree-distance
result survive when tested against the same wider null?

**Method.** For each (patient, band ∈ {α, β, γ_l}, phase ∈ {rsPre_A,
taskTest, rsPost}), generate R=200 strength-preserving 4-cycle ±δ
surrogate Laplacians (algorithm matches audit_63 exactly; SWAP_FACTOR =
20). For each surrogate r and λ ∈ {0, 1}:

    T_KC_surr_r(λ) = d_KC(λ; D^tt_surr, D^post_surr)
                    - d_KC(λ; D^pre_A_surr, D^tt_surr)

with `kc_distance(...,normalize=True)` matching audit_36 / §5.2. The
observed split-baseline T_KC (using D^pre_A as the rsPre leg) is
benchmarked at the start to confirm baseline-shift is small —
verified at Pat_02 β λ=0 (split-A: +0.503 vs full: +0.508).

**Verdict labels (spec):**

- **separated**: cohort Wilcoxon p < 0.05 (one-sided, T_obs < T_surr)
  AND |median surr T_KC across patients| < 0.05 AND ≥ 8/10 patients
  individually below their own surrogate at p<0.05 lower-tail.
- **also_negative**: median surr T_KC < 0.5 × median obs T_KC AND
  median obs T_KC < 0 (surrogate substantially recovers trace
  direction).
- **intermediate**: anything else.

**Run.** R=200 × 30 cells (10 patients × 3 bands × 2 λ), 3 phases per
surrogate. Wall-clock 5773.3 s (~96 min). Algorithm matches audit_63
exactly (independent rng realization).

### Cohort verdicts

| band | λ | obs median T_KC | surr median (per-pat med) | n below own surr | Wilcoxon p | verdict |
|---|---|---|---|---|---|---|
| α | 0 | -0.015 | -0.122 | 1/10 | 0.539 | **also_negative** |
| α | 1 | -0.396 | +0.007 | 1/10 | 0.722 | intermediate |
| β | 0 | -0.890 | -0.403 | 1/10 | 0.216 | intermediate |
| β | 1 | -0.651 | +0.739 | 1/10 | 0.539 | intermediate |
| γ_l | 0 | -0.695 | -0.101 | 2/10 | 0.313 | intermediate |
| γ_l | 1 | +1.422 | +1.489 | 1/10 | 0.500 | intermediate |

**Reading.**

1. **No "separated" cell.** Across all six (band, λ) cells, no cohort
   passes the strict separation gate. Cohort-paired Wilcoxon p ranges
   0.22–0.72; per-patient n_below ranges 1–2/10.
2. **α λ=0 also_negative.** Observed median T_KC (-0.015) is barely
   negative; matched-strength surrogate median (-0.122) is *more*
   negative. The α λ=0 trace direction is reproduced (and exceeded)
   by matched-strength noise — the observation is consistent with
   strength-driven noise, not a structural shape signature.
3. **β λ=0 leaning-contaminated.** Observed median -0.890; surrogate
   median -0.403 (≈ 50% of observed in the trace direction). The
   verdict is "intermediate" by the strict 0.5× cut, but the
   matched-strength null produces substantial trace-direction signal
   on its own — the β λ=0 observed margin over surrogate is small at
   the cohort scale.
4. **β λ=1 sign-divergent surrogate, no cohort separation.** Observed
   median -0.651; surrogate median +0.739 (opposite sign). At first
   glance this is "separated"-looking, but per-patient distributions
   are wide and only 1/10 individually pass at p<0.05. Cohort
   Wilcoxon p = 0.54.
5. **γ_l λ=0 leaning-null.** Surrogate median (-0.10) close to zero,
   observed (-0.69) more negative, but only 2/10 individually pass.
   Cohort p = 0.31.
6. **γ_l λ=1 absolute null.** Both observed and surrogate medians
   positive (+1.42, +1.49). The §5.2 γ_l λ=1 8/10 negative-T_KC count
   reflects sign-only majority; the matched-strength null produces the
   same sign distribution.

### Joint signature across three probes

Per (patient, band), the `n_probes_significant_at_p05` count combines:

- ρ_split z (audit_63, upper-tail at p<0.05)
- T_KC λ=0 z (this run, lower-tail at p<0.05)
- T_KC λ=1 z (this run, lower-tail at p<0.05)

| n_probes_sig | n cells (patient × band) | examples |
|---|---|---|
| 3/3 | 2 | Pat_06/α, Pat_05/γ_l |
| 2/3 | 3 | Pat_03/β, Pat_06/γ_l, Pat_08/β |
| 1/3 | 11 | Pat_02/β (ρ_split only), Pat_03/α, Pat_03/γ_l, … |
| 0/3 | 14 | Pat_02/α, Pat_05/α, Pat_07/α, Pat_07/γ_l, Pat_10/β, Pat_10/γ_l, Pat_13/α, Pat_13/γ_l, Pat_14/β, Pat_14/γ_l, Pat_15/α, Pat_15/β, Pat_15/γ_l, Pat_07/α |

At β specifically: only Pat_03 and Pat_08 fire on 2/3 probes; 5
patients fire on 1/3 (rho_split only); 3 patients fire on 0/3
(Pat_10, Pat_14, Pat_15). The β cohort claim is **not** built on
co-confirming evidence at the per-patient level — most patients show
the §5.3 ρ_split signature but not the §5.2 KC tree-distance signature
under matched-strength surrogacy.

### Manuscript-ready disclosure (matches ticket's "also_negative" template, adapted)

> The β KC tree-distance signature is not separated from matched-strength
> edge-identity randomization at the cohort level (cohort-paired
> Wilcoxon p = 0.22 at λ=0; p = 0.54 at λ=1; only 1/10 patients
> individually below their own surrogate at p<0.05 lower-tail at either
> λ). Surrogate median T_KC at λ=0 (-0.40) reaches roughly 50% of the
> observed median (-0.89); at λ=1 the surrogate median is sign-divergent
> (+0.74 vs observed -0.65) but per-patient distributions are wide. The
> §5.2 β KC 10/10 result against the within-baseline split-half null
> reflects narrow-null detection of a tree-distance reorganization
> whose tree-distance signature is, however, recoverable from
> matched-strength randomization at the cohort scale. The within-baseline
> KC result is not invalidated, but it is reframed: the dendrogram
> reorganization at β is shape-driven only to the extent that the
> per-node edge-strength sequence carries the cross-phase shift through
> the average-linkage construction. The §5.3 split-baseline ρ_split
> result (audit_63: β cohort-paired Wilcoxon p = 0.005) becomes the
> load-bearing LRG-layer claim; §5.2 KC becomes a scale-resolved
> geometric companion that confirms tree reorganization without
> isolating a strength-independent signature.

## §5 / §6 manuscript implications (joint, post-audit_65)

**β cohort claim — final epistemic level after three probes:**

1. **Substrate non-stationarity (§4 / §5.3)** — confirmed: split-baseline
   ρ_split cohort-paired Wilcoxon p = 0.005, n=7/10 individually above
   own surrogate (audit_63).
2. **LRG-layer tree-distance trace (§5.2)** — **does not survive
   matched-strength surrogacy**. Cohort Wilcoxon p > 0.2 at both λ;
   only 1/10 patients individually below own surrogate at either λ. The
   §5.2 within-baseline 10/10 result reflects narrow-null detection.
3. **Cohort heterogeneity attribution to implant geometry (§6 outlook)** —
   moderate uncorrected B_hemi correlation +0.697 (p=0.025, fails BH
   q=0.05 over m=135). Disclose with caveat.

The recommended reframing: §5.3 ρ_split is the load-bearing LRG-layer
claim; §5.2 KC is a scale-resolved geometric companion whose
tree-distance signature does NOT separate from matched-strength noise
at the cohort scale.

**α / γ_l cohort claims:**

- α: substrate non-stationarity at cohort scale (audit_63 p=0.002,
  per-patient 5/10). KC λ=0 lands as "also_negative" — α KC trace
  direction is reproduced by matched-strength noise. KC layer should
  not be cited as supporting α in §5.2.
- γ_l: substrate cohort claim does not survive (audit_63 p=0.116);
  KC layer also intermediate-leaning-null. The γ_l finding stays at
  the per-pair / per-leaf descriptive level (§5.4 / §5.5).

**Cohort-heterogeneity §6 prose:** the three consistently-anti patients
(Pat_07, Pat_14, Pat_15 at n_anti ≥ 6/12 from audit_64) account for
most per-patient mixedness. Implant geometry is a candidate mechanism
for β (uncorrected only) and not for α / γ_l in this sample.

## Acceptance gate (per audit_65 ticket)

> "After this run lands, the §5.2 prose finalizes its phrasing of the
> cohort universality claim at the verdict label returned: If
> **also_negative**: the manuscript reframes the β KC result as
> within-baseline-detection of a tree-distance reorganization that is
> not robust to matched-strength surrogacy at the per-pair eigenvector
> level."

The returned verdict is **also_negative** at α λ=0 explicitly, and
**intermediate-leaning-contaminated** at β λ=0 / λ=1 (50% surrogate
recovery at λ=0; sign-divergent surrogate at λ=1 but no cohort
separation). The substantive epistemic level is the same as
"also_negative" for the β claim: the §5.2 KC tree-distance shift at β
is recoverable from matched-strength randomization at the cohort scale.

The writing agent picks the disclosure paragraph drafted above (β
section) for §5.2 / §6.

## Files

### audit_64 (complete)

- `data/audit/implant_geometry/per_patient_features.csv` — per-patient
  9-feature implant table.
- `data/audit/implant_geometry/per_patient_trace_measures.csv` — per-band
  5-measure trace scalars.
- `data/audit/implant_geometry/correlation_matrix.csv` — full
  9×5×3 = 135 Spearman cells with raw p and BH-FDR q.
- `data/audit/implant_geometry/group_comparison.csv` — pro vs anti
  Mann-Whitney comparison per feature.
- `data/audit/implant_geometry/anti_alignment_count.csv` — per-patient
  anti-count across the 12 measure × band cells.
- `data/audit/implant_geometry/figures/{correlation_heatmap,pro_vs_anti_panel}.pdf`
- `data/audit/implant_geometry/README.md`

### audit_65 (in flight)

- `data/audit/kc_matched_strength_surrogate/cohort_summary.csv` —
  per-(band, λ) cohort verdict.
- `data/audit/kc_matched_strength_surrogate/per_patient_per_band.csv` —
  per-(patient, band, λ) observed + surrogate stats.
- `data/audit/kc_matched_strength_surrogate/joint_signature.csv` —
  three-probe z-score table per (patient, band).
- `data/audit/kc_matched_strength_surrogate/figures/{kc_cohort_distribution,kc_per_patient_panel}.pdf`
- `data/audit/kc_matched_strength_surrogate/README.md`
- `data/audit/kc_matched_strength_surrogate/run.log` — full per-cell
  trace from the background run.

## Provenance

- Cohort: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15 (n=10, COHORT_N10)
- Bands: alpha, beta, low_gamma
- FC method: imcoh_abs
- LRG: τ = 1/λ_max, ultrametric via average linkage on Trho.
- audit_64 wall-clock: ~30 s (correlational only)
- audit_65 wall-clock target: ~95 min (R=200 × 30 cells × 3 phases)
- KC `kc_distance(..., normalize=True)` matches audit_36 / §5.2 default.
- audit_65 surrogate algorithm matches audit_63 exactly; rng seeds differ
  (independent realization of the surrogate ensemble).
