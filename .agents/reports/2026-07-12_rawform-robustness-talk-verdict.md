---
name: 2026-07-12_rawform-robustness-talk-verdict
type: report
era: IMCOH_ABS × COHORT_N10 × rho_sym
status: current
created: 2026-07-12
updated: 2026-07-12
pointers:
  - .agents/reports/2026-07-12_brutal-review-failure-directions.md
  - scripts/01_compute/audit/audit_176_rawform_robustness.py
  - scripts/01_compute/figures_embedded/fig_rawform_robustness.py
  - data/audit/rawform_robustness/cohort_summary.csv
---

# Talk verdict — the cross-phase trace is a CONNECTIVITY property, not a diffusion artifact

**Head (the juice).** The scary discovery — `1/K(τ_min)` ≈ `1/A` — is not a
crisis and does not force a propagator-vs-adjacency choice. The apparent
α/low_γ *disagreement* between them is an artifact of comparing the propagator
against **`1/A`, the one pathological distance** (the reciprocal blows up on
weak/noisy edges; it flips Pat_02's α negative). Replace it with a *well-behaved*
raw distance — `−log A`, the textbook similarity→distance — and it **agrees with
the propagator**: **α and β clear the matched-strength gate under both; low_γ and
all nulls fail under both.** So the trace is a property of the connectivity, and
the diffusion machinery is not load-bearing *for the trace*. This is the
single-scale, honest result to present Thursday; it directly operationalizes the
brutal review's "not multiscale" finding rather than fighting it.

## Decision (locked 2026-07-12)

**Option 1 — propagator stays the spine, add ONE robustness slide.** Zero pipeline
redo. The LRG propagator is kept because it is needed anyway for (a) the epi §3
marker (coarse-τ heat *spreading* — genuine diffusion, not reproducible by any raw
distance) and (b) all β→OFC localization figures. For the cross-phase **trace**
specifically, we add one slide showing it is reproduced by plain hierarchical
clustering on a raw connectivity distance under the same matched-strength null.
The "1/K ≈ 1/A" worry becomes the talk's **strongest methods point: robustness.**

Rejected: leading with the raw `−log A` HAC as the primary result (cleaner
"not overselling LRG" story, but 3 days out it needs a β→OFC localization
recompute on `−log A`; unnecessary since the science is identical either way).
Also rejected: migrating the pipeline to any single raw distance (the earlier
audit_176 "migrate to −logA" framing is **withdrawn** — do not introduce a second
functional competing with the propagator; `A` is the τ→0 limit of the one
propagator operator).

## Verified gate table (ρ_sym matched-strength, R=200, reproduced bit-for-bit)

Cohort Wilcoxon gate p on `ρ_sym_obs − ρ_sym_surr_p50`, one-sided. Clean raw forms
= `−log A`, `1−A`; the reciprocal `1/A` is the pathological form.

| band | `1/K` propagator | `−log A` | `1−A` | `1/A` (pathological) | verdict |
|---|---|---|---|---|---|
| **β** | .032 ✓ | **.019 ✓** | .042 ✓ | .032 ✓ | **ROBUST — all four distances** |
| **α** | .024 ✓ | **.0098 ✓** | .032 ✓ | .053 ✗ | **robust — fails only `1/A`** |
| low_γ | .080 ✗ | .053 ✗ | .019 ✓ | .019 ✓ | form-dependent → **do not claim** |
| δ | .080 ✗ | .053 ✗ | .116 ✗ | .065 ✗ | null |
| θ | .784 ✗ | .348 ✗ | .423 ✗ | .065 ✗ | null |
| high_γ | .080 ✗ | .278 ✗ | .278 ✗ | .116 ✗ | null |

Reading: `−log A` reproduces the propagator's **entire** portrait, *stronger*
(α .0098 vs .024; β .019 vs .032). low_γ clears only the two *bounded* forms and
fails both *log/diffusion* forms — the fluctuation-heavy band the taxonomy already
flags as patient-specific. δ/θ/high_γ null under everything.

## Why discarding `1/A` is principled (not cherry-picking)

`A_ij = ⟨|ImCoh|⟩ ∈ (0,1]` is a similarity. `1/A → ∞` as `A → 0`, so the weakest,
noisiest edges dominate the distance — the opposite of what clustering needs.
`−log A` and `1−A` compress that tail. This is a standard property, statable
*before* looking at which bands each yields, so dropping `1/A` is a design call,
not a result-driven one. The `−log(weight)` transform is the canonical way to turn
a coherence network into an additive distance.

## What to say on the slide (talk language)

> "One might worry the diffusion machinery is doing the heavy lifting. It is not.
> At the finest scale the heat-kernel distance reduces to a re-weighted connectivity
> distance, so the cross-phase trace is a property of the **connectivity**, not the
> diffusion. We verify this directly: plain hierarchical clustering on a raw
> connectivity distance recovers the **same α and β cohort trace** under the same
> matched-strength null. The diffusion propagator is what buys the epileptic-focus
> marker and the anatomical localization — the trace itself is connectivity-level."

Frame ρ^coph as comparing the **hierarchical, nested-community** organization
across phases (real in any agglomerative tree) at a **single scale**. Do **not**
claim diffusion-*multiscale* resolution — that is degenerate on the dense graph and
lives in the sparsification work (separate track). Nothing extra is computed for
this; audit_176 already is the evidence.

## Scope / caveats (state these, don't bury)

- **Dense-graph regime.** All at τ_min on the fully-connected FC graph (the
  Villegas-2025 degenerate regime). The MS null is passed; the "is LRG *necessary*"
  answer for the trace is honestly **no** — and that is fine, because LRG earns its
  place via epi + localization + the (ongoing) sparsified multiscale.
- **Epi marker is NOT in scope of "propagator dispensable."** It uses coarse-τ
  spreading = real heat dynamics; keep it propagator-only.
- **α has no clean anatomical home** (band-taxonomy). Present α as a *robust
  cohort-level* trace (mesoscale-tuned, audit_172); β remains the flagship
  (robust **and** localized → OFC).
- **low_γ is not a result.** Consistent with prior "γ_l patient-specific / high
  fluctuation." Do not call it artifact or anatomy either — just don't headline it.

## Provenance

- Gate: `audit_176_rawform_robustness.py` — 4 cophenetic arms fed identical
  matched-strength surrogate draws (R=200, audit_150 seed/order; `diff` arm
  reproduces audit_150). Reproduced bit-for-bit on re-run (285 s).
- Figures: `fig_rawform_robustness.py` →
  `data/outputs/figures/robustness/fig_robustness_propagator_vs_raw.pdf` (main
  slide) and `…/fig_robustness_distance_matrix.pdf` (4-distance backup).
- Converges with `.agents/reports/2026-07-12_brutal-review-failure-directions.md`
  (propagator not distinctive) and the sparse-backbone track (audit_175).
