---
name: arc-mesoscale-inference-null
type: report
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-06-22
headline: N2.2/N2.5 + N2 §G mesoscale (.agents/preprint/headlines/02_encoding_vs_inference.md)
pointers:
  - data/audit/consolidation_arc/arc_mesoscale_null_R200.csv
  - scripts/01_compute/audit/audit_103d_arc_mesoscale_null.py
  - .agents/reports/2026-06-22_tau-sweep-arc-result.md
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md
---

# The β inference bump at the mesoscale is REAL — it clears matched-strength at τ≈2.6 and 6.8 while every control band stays null

**Head.** The τ-sweep (`audit_103c`, observed-only) saw a mild bump: β
inference-specific consolidation looked slightly *stronger* at coarse diffusion
times (τ≈2.6, 6.8 × 1/λ_max) than at the fine scale τ=1 — but with no null, it was
only a hint. This is the deferred referee: the **mandatory matched-strength null at
the peak τ**, reusing the canonical ensemble (seed 20260511) with its cophenetic
reformed at the mesoscale τ (zero regeneration). **Verdict: the bump survives.** β
`T_infspec·e` clears the strict null at **τ≈2.6 (p=0.014)** and **τ≈6.8 (p=0.010)**,
the effect size is modestly *larger* than fine-scale, it is LO-Pat_15 robust, and —
the decisive check — **every negative-control band stays null** at the mesoscale.
Inference-specific consolidation is therefore a genuinely **multiscale** phenomenon,
mildly favouring the mesoscale — not a fine-scale artifact, and (crucially) not the
coarse-graining artifact that sank the truncation null.

## What was run

`audit_103d` runs the matched-strength surrogate null on the arc functionals at
**τ_mult ∈ {1.0, 2.61, 6.81}** (the fine-scale anchor + the two β bumps, taken
exactly from the `audit_103c` grid). The canonical surrogate eigendecompositions are
cache hits (the same ensemble that validated the fine-scale arc, `audit_103 --null`,
R=200); only the cophenetic is reformed at the new τ. Per-patient p = mean(surr ≥
obs); cohort = one-sided paired Wilcoxon(obs > surr-median) — audit_103's exact
logic. All 6 bands run; **δ/θ/low-γ/high-γ are the negative controls** (the locked
`audit_103b` lesson: a new τ is a new test surface, so run a known-null band through
it). **Correctness anchor PASS:** at τ=1 the per-patient `T_infspec·e` reproduces
`audit_103 --null` to max|Δ| = 9.7e-17.

## Result — β inference-specific (`T_infspec·e`)

| τ / λ_max | obs median | obs−null gap | obs>0 | p<.05 | Wilcoxon p | LO-P15 p |
|----------:|-----------:|-------------:|:-----:|:-----:|-----------:|---------:|
| **1.00** (fine) | +0.122 | +0.109 | 9/10 | 7/10 | **0.0068** | 0.0137 |
| **2.61** (meso) | +0.145 | +0.135 | 8/10 | 7/10 | **0.0137** | 0.0273 |
| **6.81** (meso) | +0.153 | +0.148 | 8/10 | 6/10 | **0.0098** | 0.0195 |

- **β clears at every scale**, and the **effect size (obs−null gap) widens
  monotonically** fine→meso (+0.109 → +0.135 → +0.148). The mesoscale signal is
  *real beyond strength geometry* — the strength-matched null coarse-grains
  identically, so if the bump were generic coarse-graining the surrogate would
  reproduce it and β would not clear. It does.
- **Significance is comparable across scales** (τ=1 is marginally the most
  significant by p, with the tightest per-patient agreement, 7/10). So the honest
  claim is **mesoscale-robust with a mild effect-size tilt toward coarser scale** —
  not a dramatic mesoscale-only peak.

## The negative controls held (the decisive part)

`T_infspec·e` for every known-null band, min Wilcoxon p across the mesoscale τ:

| band | p @ τ=1 | p @ 2.61 | p @ 6.81 | min |
|------|--------:|---------:|---------:|----:|
| α (encoding-only) | 0.080 | 0.278 | 0.216 | 0.080 |
| δ | 0.539 | 0.246 | 0.188 | 0.188 |
| θ | 0.812 | 0.754 | 0.784 | 0.754 |
| low-γ | 0.995 | 0.903 | 0.947 | 0.903 |
| high-γ | 0.348 | 0.246 | 0.188 | 0.188 |

**None clears at any mesoscale τ.** This is the opposite of the retired truncation
null (`audit_103b`), where δ false-positived at p≈0.014. There, the test surface
manufactured significance; here the mesoscale surface does not — only β passes, and
the frequency-specificity is preserved at coarse scale. The mesoscale bump is a
property of β's geometry, not of coarse-graining.

> Note vs the sibling N1 τ-sensitivity result: there, the *single-phase* trace's
> coarse-τ "gains" were a **collapse artifact** (the placebo rose to match). Here
> the statistic is the four-phase *difference-of-differences* partial concordance,
> and the matched-strength null **is** that placebo — β beats it at the mesoscale
> while the controls don't. Different statistic, genuinely different (clean) verdict.

## What it means

The relations the brain *figured out* (inference) leave a persistent β trace that is
**at least as strong when the network is read in larger, multi-step chunks** as at
the finest scale — modestly stronger, in fact. That is the scale signature you would
expect if inference is **integration over multi-step relational paths** (slower
diffusion mixes longer paths), as opposed to a purely local edge change. It is a
correlational scale signature, not a proof of mechanism — but it is now a
*verified* one, upgrading the τ-sweep's observed-only hint.

## Honest ceiling

- **Mesoscale-robust, not mesoscale-exclusive.** Effect size tilts coarser;
  significance is comparable across scales. Do not sell it as "inference lives at the
  mesoscale" — sell it as "inference-specific consolidation is multiscale and
  survives the strict null out to ~7× coarser diffusion, mildly favouring coarse."
- **Independent-per-phase null.** Like all arc results, this uses the
  independent-per-phase matched-strength ensemble; the coordinated-cross-phase null
  (audit_63 gap) is a standing caveat across the whole arc, not specific here.
- **No behavioral anchor** (none obtainable, PI 2026-06-22).
- **Secondary observation (not headlined):** the *encoding* echo `T_learn` broadens
  at the mesoscale — δ goes from n.s. (p=0.097) at τ=1 to clearing at τ=6.8
  (p=0.024), and low-γ trends the same way. So memorising echoes become more
  broadband at coarse scale, while inference stays β-exclusive. Logged, not claimed.

## Provenance

- Null → `data/audit/consolidation_arc/arc_mesoscale_null_R200.csv` (180 rows:
  10 pat × 6 band × 3 τ) · `audit_103d_arc_mesoscale_null.py` · 2026-06-22.
- Anchor: τ=1 reproduces `arc_null_per_patient.csv` (audit_103 --null) to 1e-16.
- Upstream hint: `.agents/reports/2026-06-22_tau-sweep-arc-result.md` (`audit_103c`).
- Scope (5-point preamble): `audit_103d` docstring + the τ-sweep scope
  `.agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md`
  (P2 scale-of-inference). Memory `arc_inference_consolidation_2026_06_18`.
