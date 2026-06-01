---
name: pair-trace-measures-methodology-audit
type: report
era: COHORT_N10
status: current
created: 2026-05-30
updated: 2026-05-30
kind: methodology-comparison
fc_method: imcoh_abs
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
pointers:
  - .agents/guides/task-persistence-investigation/2026-05-30_asymmetric-pair-trace-regression.md
  - .agents/preprint/bands/01_beta.md
  - scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py
  - data/audit/pair_trace_measures_comparison/cohort_summary.csv
sources:
  - data/audit/pair_trace_measures_comparison/per_patient_per_band.csv
  - data/audit/pair_trace_measures_comparison/cohort_summary.csv
  - data/audit/pair_trace_measures_comparison/figures/comparison_heatmap.pdf
---

# Asymmetric pair-trace measures vs the symmetric `ρ_split^coph` — methodology audit

> **NOT FOR PREPRINT INCLUSION (user directive 2026-06-01).** These findings are
> kept as an internal methodology record only. The asymmetric slope `s_TR`, the
> `ρ_Pearson` comparison, and the per-patient heterogeneity discussion below are
> **not** to be added to the manuscript for now. The current preprint cycle
> continues with `ρ_split^coph` (Spearman) unchanged. Bear the results in mind;
> do not surface them in `.agents/preprint/` band files or LaTeX.

## Head

**Switching the §5.3 headline from the rank-symmetric `ρ_split^coph`
(Spearman) to the asymmetric regression slope `s_TR = ⟨Δ_task, Δ_rest⟩ /
‖Δ_task‖²` would LOSE the β-band trace, not strengthen it.** On the cophenet
primitive `D_coph` — the one primitive that carries band-selectivity — the
matched-strength cohort verdict goes from **β: 7/10 patients, p = 0.007**
(Spearman) to **β: 2/10, p = 0.161** (slope). The mechanism is clean and
decisive: `s_TR` is magnitude-weighted, and on `D_coph` the per-pair
magnitudes are dendrogram merge heights that are themselves strength-driven —
exactly the structure the mandatory matched-strength null nullifies. The
surrogate slope reproduces the observed slope almost exactly per patient
(e.g. Pat_03 obs +0.612 vs its own surrogate median +0.626), so the margin
collapses. Spearman survives because ranking strips that strength-driven
magnitude. **Recommendation: keep `ρ_split^coph` as the primary measure.** If
a magnitude-aware companion is wanted, the symmetric `ρ_Pearson` on `D_coph`
is the defensible choice (it preserves the α+β selectivity: β 8/10, p = 0.005);
the asymmetric slope and `R²` are both unsuitable (slope loses the signal; `R²`
fires on all six bands including the *anti*-trace θ band).

This is a comparative methodology audit, not a manuscript change. Nothing in
the locked preprint is touched. The current cycle continues with `ρ_split^coph`.

## What was run

`audit_73_pair_trace_measures_comparison.py` computes **four cross-phase
measures** — `ρ_Spearman` (current), `ρ_Pearson`, the asymmetric slope `s_TR`,
and `R² = ρ_Pearson²` — on **three LRG primitives** — raw density `ρ̂(τ_max)`,
raw distance `D(τ_max) = 1/ρ̂`, and cophenet `D_coph = cophenet(UPGMA(D))` —
over the full n=10 cohort × 6 bands, each against the **same R=200
matched-strength surrogate ensemble** already on disk
(`data/cache/matched_strength_surrogate_lrg/`, seed 20260511). Split-baseline
convention identical to audit_63: `Δ_task = X(task) − X(rsPre_A)`,
`Δ_rest = X(rsPost) − X(rsPre_B)` on independent rsPre halves.

Per `feedback_matched_strength_mandatory`, the matched-strength null was re-run
for `s_TR` — ρ_split's clearance does not transfer to a magnitude-aware
measure. Per `feedback_no_pre_registered_acceptance`, no acceptance gate was
pre-registered; the cohort-paired one-sided Wilcoxon (obs > surrogate-median)
is the gate. Per `feedback_no_single_patient_p_driven`, LOO sensitivity is
reported per cell.

**Reproduction check (harness validity).** The `(D_coph, ρ_Spearman)` cell
reproduces the published `ρ_split^coph`: Pat_05 β obs = 0.490999 vs published
0.490999 (Δ = 4.6e-07); cohort-median β = +0.2206 (published +0.2206). The raw
`D` rows reproduce `preprint_05` **bit-exactly** (Pat_05 β obs 0.461800,
surr_p50 0.033965, obs_z 37.5476 — identical). The harness is sound.

## The verdict — `D_coph` (the band-selective primitive)

Cohort-paired one-sided Wilcoxon p (obs > surrogate-median), matched-strength:

| measure | δ | θ | α | β | γ_l | γ_h | β n_above | β ratio (med) |
|---|---|---|---|---|---|---|---|---|
| **ρ_Spearman** (current) | 0.216 | 0.688 (anti) | **0.005** | **0.007** | 0.116 | 0.278 | 7/10 | 18.2× |
| **ρ_Pearson** | 0.188 | 0.385 (anti) | **0.003** | **0.005** | **0.042** | 0.097 | 8/10 | 66× |
| **s_TR (slope)** | 0.042 | 0.539 | 0.216 | 0.161 | 0.138 | 0.246 | **2/10** | 8.0× |
| **R²** | **0.005** | **0.003** (anti) | **0.010** | **0.001** | **0.001** | **0.002** | 8/10 | 242× |

Bold = p < 0.05. "anti" = cohort obs-median sits *below* the surrogate median
(the one-sided "greater" test correctly returns high p there). Figure:
`figures/comparison_heatmap.pdf`.

Three readings, top to bottom of the `D_coph` block:

1. **`ρ_Spearman` (current §5.3) reproduces its known signature:** α and β
   light up (p = 0.005, 0.007), everything else dark; θ is anti-trace. This is
   the band-selective cophenet result the manuscript rests on.
2. **`ρ_Pearson` agrees and is slightly sharper:** same α+β core, β strengthens
   to 8/10 and ratio 66× (vs 18× for Spearman), γ_l becomes marginal (0.042).
   A magnitude-aware *symmetric* measure keeps the selectivity.
3. **`s_TR` (the proposed asymmetric headline) does NOT recover α or β.** β
   p = 0.161 (2/10), α p = 0.216. The only "significant" `s_TR` cell is δ
   (0.042) — spurious, and not where the trace lives. **The slope fails the
   gate at exactly the band the manuscript depends on.**
4. **`R²` fires on all six bands** — including θ, which is *anti*-trace on
   Spearman/Pearson (p = 0.69 / 0.39). A measure that flags the anti-trace band
   at p = 0.003 cannot distinguish trace from drift. `R² = Pearson²` discards
   sign, so any linear association (trace, anti-trace, or shared drift) reads
   as "significant". Unusable as a trace statistic on its own.

## Why `s_TR` fails — the mechanism

The slope is `s_TR = ⟨Δ_task, Δ_rest⟩ / ‖Δ_task‖²`: a magnitude-weighted
quantity dominated by the highest-magnitude pairs. On `D_coph` those are the
deepest dendrogram merge-height shifts — and merge-height magnitude is governed
by the node-strength distribution, which is precisely what the matched-strength
surrogate preserves. So the surrogate reproduces the observed slope nearly
pair-for-pair. Per-patient evidence at β `D_coph`:

| patient | `s_TR` obs | own-surrogate `s_TR` median | obs above own surrogate? |
|---|---|---|---|
| Pat_02 | +0.528 | +0.443 | p = 0.050 (just) |
| Pat_03 | +0.612 | **+0.626** | **no** (obs *below* surrogate) |
| Pat_05 | +0.912 | +0.802 | p = 0.205 |
| Pat_06 | +0.712 | +0.704 | p = 0.380 |
| Pat_08 | +0.375 | −0.179 | p = 0.005 (yes) |

Pat_03 is the smoking gun: an impressive-looking 61% "recovery slope" that is
*entirely* reproduced by a strength-matched random graph (surrogate 63%). Even
the strongest patient (Pat_05, 91%) sits inside its surrogate distribution
(p = 0.21). Only 2/10 patients clear their own surrogate.

By contrast, `ρ_Spearman` strips magnitude to ranks: Pat_05 β obs Spearman
0.491 vs surrogate median 0.051 — a 9.6× margin that the surrogate cannot
reproduce, because rank co-movement of *which* pairs shift together is not a
function of the strength distribution alone.

Two secondary strikes against `s_TR` as a cohort statistic:
- **Between-patient instability.** β `D_coph` `s_TR` ranges −0.18 … +0.91
  across patients (coefficient of variation 1.42) vs −0.09 … +0.51 for
  Spearman (CV 0.90). The unbounded magnitude-weighted slope is far noisier
  patient-to-patient.
- **Units are not cross-primitive comparable** (scope §5.2): `s_TR` on ρ̂, D,
  and `D_coph` live on different dynamic ranges, so the slope cannot be read
  across the primitive axis the way `R²` or a ratio-over-null can.

## Raw `ρ̂` and raw `D` — the permissive layer (context)

On the raw primitives, `ρ_Spearman`/`ρ_Pearson` fire in several bands
(δ, β, γ_l on raw D; β, γ_l on raw ρ̂) — the "trace-everywhere" permissiveness
that `cophenet_methodology_rationale` documents. The cophenet wrap is what
buys band-selectivity. This audit confirms that story holds for the *symmetric*
measures. It does **not** rescue `s_TR`: the slope is non-selective and
sub-threshold on raw D too (β p = 0.080), and the `R²` over-firing is universal
across primitives. So the cophenet projection cannot fix the asymmetric slope —
the problem is the magnitude weighting, not the primitive.

## Recommendation

1. **Keep `ρ_split^coph` (Spearman on `D_coph`) as the primary §5.3 measure.**
   It is principled *for this substrate* precisely because rank-invariance
   discards the strength-driven magnitude that the mandatory null nullifies.
   The "more principled, units-bearing" appeal of the asymmetric slope is
   real in the abstract but self-defeating here: the units it bears (fraction
   of merge-height shift recovered) are dominated by strength structure.
2. **If a magnitude-aware companion is desired, use the symmetric
   `ρ_Pearson` on `D_coph`**, not the slope. It preserves the α+β selectivity
   (β 8/10, p = 0.005) and is reportable as "the trace survives when pair
   magnitudes are retained, in a sign-and-direction sense". This is a
   legitimate robustness sentence; the slope is not.
3. **Do not adopt `s_TR` as a headline and do not report `R²` as a trace
   statistic.** `s_TR` can appear as a descriptive per-leaf decomposition
   (`figures/per_leaf_slope_beta/`) with the explicit caveat that the
   cohort-level slope fails matched-strength at β; `R²` should not be cited as
   evidence of trace at all (it flags the anti-trace θ band).
4. **No manuscript change.** Per `feedback_preprint_folder_scope`, if any of
   this enters the methods, it does so as a sensitivity paragraph
   ("we considered the asymmetric regression slope and found it loses the β
   trace under matched-strength because…"), authored under
   `.agents/preprint/methods/`, only on explicit user approval.

## Addendum (2026-06-01) — cross-band `ρ_Pearson`, patient heterogeneity, and Grassmann coherence

### `ρ_Pearson` across all six bands (`D_coph`)

The magnitude-aware *symmetric* Pearson does not collapse like `s_TR`; it is
sharper than Spearman on the bands the manuscript already claims, but slightly
more permissive at the margin.

```
              SPEARMAN (current)            PEARSON
band       p_MS  n_above  LOO_max     p_MS  n_above  LOO_max
delta     0.216   4/10    0.367      0.188   5/10    0.326   both null
theta     0.688*  2/10    0.875      0.385*  4/10    0.590   both anti
alpha     0.005   5/10    0.010      0.003   7/10    0.006   both SIG; Pear sharper
beta      0.007   7/10    0.014      0.005   8/10    0.010   both SIG; Pear sharper
low_gamma 0.116   5/10    0.213      0.042   6/10    0.082   <-- the only flip
high_gamma 0.278  4/10    0.455      0.097   6/10    0.180   both non-sig
```
`*` = anti (obs below surrogate). The α/β core agrees under both measures
(robustness across the magnitude axis). The single divergence, γ_low, is
LOO-fragile (`LOO_max = 0.082` → removing Pat_02 alone pushes it back above
0.05) and a 6/10 split — a genuine 5–6-patient signal in both measures, not a
clean new band. **Recommendation: report Spearman+Pearson as a robustness pair
on α/β; do not elevate γ_low. Do not switch primary.**

### The per-patient heterogeneity problem (the central concern)

No band is *immune* (clean null in all 10 patients) and no band is *universal*
(effect in all 10). Every band is a cohort split, and the **same patients
recur as dissenters** (Pat_10, Pat_14, Pat_15 at β cophenet per-pair). So the
band-selectivity claim rests on the cohort Wilcoxon crossing threshold, which
is itself a function of cohort composition. Per-patient β `D_coph` Spearman:
7/10 trace (Pat_10 −0.091, Pat_14 −0.049 anti; Pat_15 +0.083 below own
surrogate). This is the honest weak point the matched-strength null does not
resolve: it validates the *cohort median* is above a strength-matched null, but
it does not make the effect *present in every patient we claim it for*.

### Grassmann vs cophenet per-pair — cross-patient coherence (β)

The Grassmann subspace probe (audit_66) is **more patient-coherent** than the
cophenet per-pair probe at β. Convention: `T_G > 0 = trace`, obs **above**
surrogate (locked 2026-05-26; `audit_66` code line 248; see the convention
correction logged below).

```
probe                       β trace patients   dissenters
cophenet per-pair (Spearman)     7/10          Pat_10, Pat_14, Pat_15
Grassmann subspace (k=27..55)    9/10          Pat_15 only
```

Pat_10 and Pat_14 carry the β reorganization in the **leading-eigenmode
subspace rotation** (Grassmann z = +5.07, +5.11, ~100% of the k-window above
their own surrogate) even though they are flat/anti in the per-pair
merge-height distance. Only Pat_15 (right-hemisphere-only implant, the
biology-driven anti-anchor) dissents on Grassmann. So under the subspace probe,
"where is the β trace in Pat_10?" has an answer (in the slow-mode subspace);
under the per-pair cophenet probe it is a genuine dissent. This is the
empirical basis for the user's intuition that Grassmann is the more coherent
probe — confirmed, at β. (Whether Grassmann is *also* more band-selective —
coherent AND selective — is the open question for the "where to go" decision;
not yet checked across all bands here.)

### Convention correction logged (2026-06-01)

While pulling the Grassmann per-patient data, the `01_beta.md` §3.3 lab-note was
found to present `T_G` in the **pre-lock "negative = trace" convention** (obs
−0.489 etc.), contradicting the locked 2026-05-26 `T_G > 0 = trace` convention
that the `audit_66` code and `cohort_summary.csv` already use. The manuscript
table (`beta_per_patient.tex`) was already correct (positive throughout, user
directive 2026-05-19). Fixed `01_beta.md` §3.3 (formula order, both tables,
Wilcoxon direction, `n_below`→`n_above`, prose) to positive=trace. **Still
open** (flagged, not auto-fixed): (a) the substrate §3.1 `T_d^(d_S)` is *also*
in the old convention AND shows −0.038 where the current CSV says +0.088 at
8/10 — a sign **and** staleness mismatch needing provenance review; (b) the
Grassmann verification reports (`2026-05-11`, `2026-05-14`) and the
`2026-05-19_beta_red_paragraphs` writing directive still use "negative =
trace / below own surrogate" and should be swept or marked historical.

## Honest limitations of this audit

- The matched-strength null is the only null run here; the drift-floor null for
  `s_TR` (scope §10 Q4) was deferred. It would not change the verdict — `s_TR`
  already fails the *stronger* matched-strength gate — but it is owed before
  any positive `s_TR` claim, which we are not making.
- The third-common-cause confound (scope §5.5) is unfalsifiable in this
  four-phase design for every measure here, `ρ_split^coph` included. This audit
  compares measures under a shared null; it does not certify the underlying
  trace beyond what §5.3 already establishes.
- `R²` "fires everywhere" is a property of comparing a non-negative statistic's
  observed value against a near-zero surrogate median with a one-sided
  "greater" test. A two-sided or sign-aware framing would demote θ, but `R²`
  still cannot carry directional trace information by construction — the
  recommendation stands.

## Provenance

| Artifact | Path |
|---|---|
| Scope report | `.agents/guides/task-persistence-investigation/2026-05-30_asymmetric-pair-trace-regression.md` |
| Compute script | `scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py` |
| Library helper | `src/lrg_eegfc/utils/metrics/hypothesis.py::regression_slope_through_origin` (+ test in `tests/test_hypothesis.py`) |
| Per-patient CSV (720 rows) | `data/audit/pair_trace_measures_comparison/per_patient_per_band.csv` |
| Cohort CSV (72 rows) | `data/audit/pair_trace_measures_comparison/cohort_summary.csv` |
| Verdict heatmap | `data/audit/pair_trace_measures_comparison/figures/comparison_heatmap.pdf` |
| Per-pair scatter (β backbone) | `data/audit/pair_trace_measures_comparison/figures/scatter_slope_beta/{Pat_02,05,06,08}.pdf` |
| Per-leaf `s_TR(c)` dendrogram | `data/audit/pair_trace_measures_comparison/figures/per_leaf_slope_beta/{Pat_02,05,06,08}.pdf` |
| Surrogate cache (reused) | `data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` |
| Reproduction reference | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` (audit_63), `data/preprint/rho_split_raw_D/beta_matched_strength.csv` (preprint_05) |
