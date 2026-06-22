---
name: tau-sweep-arc-result
type: report
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-06-22
headline: N2.2/N2.5 (.agents/preprint/headlines/02_encoding_vs_inference.md)
pointers:
  - scripts/01_compute/audit/audit_103c_tau_sweep_arc.py
  - data/audit/consolidation_arc/arc_tau_sweep.csv
  - data/audit/consolidation_arc/figures/arc_tau_sweep.pdf
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md
---

# The encoding/inference decomposition is SCALE-ROBUST — β inference holds across diffusion times, not just the finest one

**Head.** Every N2 number was computed at a single diffusion time τ = 1/λ_max (the
finest scale). Sweeping τ from there out to 10× shows the decomposition does **not**
depend on that choice: the **β inference-specific consolidation stays positive and
significant across the whole window**, the **α/β dissociation holds at every scale**,
and the other four bands stay null throughout. So N2.2 graduates from "robust at one
scale" to **scale-robust** — a strengthening, no new claim. (Observed-statistic
diagnostic; the matched-strength verdict still lives at τ = 1/λ_max.)

## What was run

`audit_103c` lifts the diffusion time in the arc cophenetic to a parameter
(τ = τ_mult/λ_max per phase) and recomputes the encoding echo `T_learn(τ)` and the
inference-specific `T_infspec·e(τ)` across a log grid τ_mult ∈ [1, 10] (13 points),
for all 10 patients × 6 bands. Each phase Laplacian is eigendecomposed once; the
heat kernel is reformed per τ. **Observed-only** (no surrogate — that is mandatory
only for a new positive claim at some τ ≠ 1, per the scope). Significance shown is a
cohort one-sided Wilcoxon of the observed concordance against 0 — a **shape
diagnostic**, explicitly not the matched-strength referee.

**Correctness anchor: PASS** — at τ_mult = 1 the functionals reproduce `audit_103`
to max|Δ| ≈ 1e-16. **Numerical floor: clean** — the cophenetic is fully finite/
non-degenerate (finite fraction = 1.000) at every τ, so nothing below is a
spectral-extreme artifact.

## Result

- **β inference is scale-robust.** `T_infspec·e(τ)` for β is positive at every τ in
  [1, 10], cohort median ≈ +0.09 to +0.15, observed-significant (p < 0.05) across
  almost the entire window. It is **not** a fine-scale artifact of τ = 1/λ_max.
- **The α/β dissociation holds at every scale.** β stays clearly above α throughout;
  α's inference concordance is weak (≈ +0.06 to +0.11) and mostly non-significant at
  every τ — never robustly positive. low-γ, θ, δ, high-γ inference stay null/negative
  across the whole sweep. The frequency-specificity is a property of the geometry,
  not of the chosen scale.
- **Encoding echo `T_learn` is scale-robust too** — positive and mostly significant
  for α, β (and δ at larger τ) across the window; duration-immune at every τ by
  construction (no `task_test` term).
- **A mesoscale bump — now VERIFIED (`audit_103d`, 2026-06-22).** β inference shows
  shallow local maxima at τ ≈ 2.6 (+0.145) and τ ≈ 6.8 (+0.153) — above the
  fine-scale +0.122. This observed-only bump was deferred to the mandatory referee;
  that referee has now run: β `T_infspec·e` **clears the matched-strength null at the
  mesoscale** (τ≈2.6 p=0.014, τ≈6.8 p=0.010, LO-P15 robust), effect size modestly
  larger than fine-scale, **while every control band stays null** (δ did not
  false-positive). So inference-specific consolidation is multiscale, mildly
  mesoscale-favouring — *consistent with* integration over multi-step relational
  paths. Result `.agents/reports/2026-06-22_arc-mesoscale-inference-null.md`.

## Why it matters

`τ = 1/λ_max` is borrowed from the LRG papers' sparse-spectrum regime, which does
not strictly apply to our dense, weight-heterogeneous graphs (a standing caveat on
every cohort claim). This sweep discharges that caveat for N2: the inference
decomposition is the same story from the finest scale out to 10× coarser. The
headline does not hinge on a scale choice.

## Critical issues

- **Observed-only.** The p-values here are vs 0, not vs the strength-matched null.
  The *validated* verdict remains at τ = 1/λ_max (β `T_infspec·e` p ≈ 0.007,
  matched-strength). The sweep establishes **shape/robustness**, not a new
  significance at any other τ.
- **Mesoscale bump is a hint only** — needs a surrogate at the peak τ before it can
  be reported as "inference is mesoscale."
- **Single dataset / n = 10.** The τ-grid stops at 10× (the SOZ-marker informative
  regime); far-coarser τ (toward τ → ∞, uniform kernel) was not probed and would be
  degenerate anyway.

## Provenance

- Sweep → `data/audit/consolidation_arc/arc_tau_sweep.csv` (780 rows) ·
  `audit_103c_tau_sweep_arc.py` · 2026-06-22.
- Figure → `data/audit/consolidation_arc/figures/arc_tau_sweep.pdf`
  (`--figure`; τ-response curves, β/α emphasized).
- Scope (5-point preamble) →
  `.agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md`.
- Headline home: N2.2/N2.5 in `.agents/preprint/headlines/02_encoding_vs_inference.md`.
