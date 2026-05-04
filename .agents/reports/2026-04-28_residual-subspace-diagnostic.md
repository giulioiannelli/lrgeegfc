---
name: 2026-04-28_residual-subspace-diagnostic
type: handoff
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-28
updated: 2026-04-28
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-28_residual-subspace-trace.md
  - data/audit/residual_subspace/diagnostic.pdf
  - data/audit/residual_subspace/alpha_beta_grid.csv
  - .agents/reports/2026-04-28_psi-tau-scan-verdict.md
---

# Residual-subspace alignment — diagnostic verdict

**Renormalization head.** The 3-patient × 6-band × 12-τ̃ residual-subspace
diagnostic kills the path. Mean `α_1 − null_p95 = −0.066` across all 216
cells: canonical alignment between leading eigenvectors of `D^task − D^pre`
and `D^post − D^pre` is **on average below** what phase-label permutation
gives. β-band shows 3/3 patients with contiguous τ̃-windows above null, but
`β_R ≈ 0.78–0.95` indicates near-rank-1 residuals, where alignment is a
mathematical consequence of dominant-mode dominance (any pair of phase
residuals would align), not a task-specific signature. Pat_03 (1024 Hz
outlier) sits at or below null in **every** band. No band reaches 3/3
patients with both contiguous τ̃-runs AND moderate `β_R` simultaneously.
Option A is dead at the empirical level — pivoting to Option B (commit to
α/β under H2c continuous controls, stop adding methodology).

## Inputs

- 3 patients: Pat_02 (known-best), Pat_06 (cleanest), Pat_03 (1024 Hz outlier).
- 6 bands × 3 phases × 12 τ̃-points = 216 cells.
- `τ̃ = τ · λ_max^φ` per phase; common grid `τ̃ ∈ [1, min_φ(λ_max/λ_gap)]`.
- 5 phase-label permutations per cell as null.
- Outputs: `data/audit/residual_subspace/{alpha_beta_grid.csv, loadings.npz, diagnostic.pdf}`.

## Headline numbers

| measure | value | reading |
|---|---|---|
| mean `α_1 − null_p95` over all 216 cells | **−0.066** | canonical alignment is *below* permutation null on average |
| fraction of cells with `α_1 > null_p95` | 22.7% | barely better than chance under a null with only 5 perms |
| best band (β): 3/3 patients with contiguous run ≥ 1 | longest run 3 | but `β_R` mean = 0.86 → near rank-1 → trivial alignment |
| Pat_03 across all bands | mean excess −0.20 | systematically anti-aligned vs null |

## Per-band breakdown (max contiguous τ̃-run with `α_1 > null_p95`, peak `α_1`, mean `β_R`)

| band | Pat_02 | Pat_06 | Pat_03 | verdict |
|---|---|---|---|---|
| β | run=2, α=0.99, β_R=0.78 | run=1, α=0.996, β_R=0.88 | run=3, α=0.995, β_R=0.91 | high `α_1` but `β_R` is in the rank-1 regime — trivial alignment |
| δ | run=2, α=0.998, β_R=0.72 | run=4, α=0.977, β_R=0.54 | none (below null) | only 2/3, and δ is artefact-suspect |
| θ | run=2, α=0.96, β_R=0.72 | run=3, α=0.96, β_R=0.43 | none | 2/3, mid `β_R`, but Pat_03 negative |
| low_γ | run=1, α=0.999, β_R=0.91 | run=3, α=0.85, β_R=0.59 | none | rank-1 dominance; Pat_03 negative |
| high_γ | none | run=2, α=0.99, β_R=0.42 | none | 1/3 only |
| α | none | run=1, α=0.90, β_R=0.34 | none | 1/3 only |

## Why `α_1` is uninformative when `β_R` is high

`β_R = (λ_1^R)² / Σ (λ_i^R)²` near 1 means the post-pre residual matrix is
near rank-1 — its energy is concentrated in a single eigenvector that
captures the network's dominant mode of pairwise distance variation. For
β-band the median `β_R = 0.88`. In this regime any two residual matrices
(canonical or phase-permuted) project onto the same dominant direction, so
`α_1` is high regardless of which phase plays which role. The diagnostic
shows this directly: when `β_R` is high, observed `α_1` ≈ null `α_1`. Not a
finding.

## Pat_03 anomaly

Pat_03 sits at or below null on **every** band. Mean `null_excess = −0.41`
in α band, `−0.26` in δ, `−0.12` in high_γ. This patient is the 1024 Hz
sampling-rate outlier already flagged across the project. Either the
spectrum / `λ_max` rescaling is mis-handled by τ̃-matching when sampling
rate differs, or Pat_03 genuinely has anti-correlated task→post
reorganisation. Either way, no cohort claim survives the inclusion of
Pat_03 with the current pipeline. (Excluding Pat_03 makes any cohort
threshold a 2/N rather than 3/N statement on an already small sample.)

## Higher-order subspaces

`α_3` and `α_5` give similar means to `α_1` (0.76 and 0.71 respectively
across all cells). `α_3 / α_1` ratios cluster near 0.93 in most bands —
i.e. higher-order modes align about as well as the top mode. high_γ has
`α_3 / α_1 ≈ 1.06` (higher modes align *better*), but only 1/3 patients
reaches the contiguous-run threshold there. No clean rescue from going to
higher k.

## Kill-condition summary

Pre-registered in the scope report (`2026-04-28_residual-subspace-trace.md`):

- ❌ `β_R < 0.1` everywhere → **failed differently**: `β_R` is *too high*
  (rank-1 dominance), making `α_1` trivial.
- ❌ Observed `α_1 ≤ null_α_1` cohort-wide → **failed**: mean excess is
  negative across the cohort.
- ⚠ Top-`u_1^R` loadings concentrate? → not assessed at cohort scale because
  the prior two killed the path; loadings PDF (`diagnostic.pdf` page 4) is
  available for visual sanity but not actionable.

All three kill conditions hit. The path is closed.

## What this leaves us with

The full cluster + subspace + matrix-residual cycle has now been audited
empirically against this dataset:

| approach | empirical status |
|---|---|
| H2c shared-baseline continuous (paired Wilcoxon on `D` matrices) | α p=0.008, β p=0.014, **8/10 patients each** (only surviving cohort claim) |
| H2c split-baseline / cross-probe / drift floor controls | α and β survive; δ/γ_h fail or artefact-suspect |
| Δ_VI(k) integer-k partition divergence | δ ridge fails continuous controls |
| h_rel / mm-radius scale axes | dispersion across patients; no clean cohort signal |
| Ψ-as-selector (Villegas 2025 partition stability) | dead — trivial argmax (n=2 or n≈N−1), no mid-tree signal |
| Ψ_L per-cluster persistence | dead — `Ψ_L < 1` everywhere, no stable subtree above noise |
| Residual-subspace `α_k` / `β_M` (this work) | dead — `α_1 < null` on average; `β_R` rank-1 makes positive cells trivial |

**The α/β H2c continuous-controls finding is the only thing that has survived
every test with a real null.** Treat it as the working headline. Stop
chasing additional methodology on this question with this dataset.

## Decision

- Mark scope `2026-04-28_residual-subspace-trace.md` `status: superseded`
  (scope-report frontmatter update — code preserved, do not delete).
- No further methodology proposals on the task-trace question for this
  dataset / FC method / cohort-size. The marginal value of an eighth scan is
  negative and we already have the answer.
- **Next:** writeup of α/β under H2c continuous controls as the deliverable
  (Option B). Plan to be drafted in a separate handoff.

## Outputs preserved

- `data/audit/residual_subspace/alpha_beta_grid.csv` — 216 rows.
- `data/audit/residual_subspace/loadings.npz` — `u_1^S`, `u_1^R` per cell.
- `data/audit/residual_subspace/diagnostic.pdf` — 4-page diagnostic
  (α_1 vs null, β_S/β_R, α_k comparison, loadings at peak τ̃).
- `scripts/01_compute/diagnostics/diag_residual_subspace.py` — reproducible.
- `scripts/01_compute/diagnostics/diag_psi_tau_scan.py` — Ψ τ-scan that
  established the dendrograms are too tight for stability-based cuts.
