---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: audit-report
supersedes: psi_argmax_histogram (3-cell tolerance band counts)
---

# Ψ argmax pinning — strict 3-class cohort sweep (top-down convention)

**Head.** Under the new top-down Ψ convention (n=0 = root cut, n=N−3 = leaf-pair cut) and a strict no-tolerance 3-class taxonomy `{top_cut, bottom_cut, interior}`, **89/180 cells (49%) pin at exactly one of the two Ψ extrema** (21 top-cut, 68 bottom-cut) and 91/180 are strict-interior. **However**, 88% of those strict-interior cells sit within 5 indices of a boundary, so the histogram is still overwhelmingly bimodal pile-up — the strict count alone *understates* the visual story.

The earlier 137/180 (76%) figure-caption number was under an asymmetric 3-cell tolerance band `(n ≤ 2 OR n ≥ N−3)` in the old ascending convention. The earlier "180/180 generic" coding-agent claim was rhetorical (mechanism in the limit), not a count.

## Convention

- **Top-down (descending).** `heights = sort(Z[:, 2])[::-1]` so heights[0] = root (largest), heights[N−2] = leaf-pair (smallest).
- **Ψ(n) = N · [log₁₀ tₙ − log₁₀ tₙ₊₁]**, indexed n ∈ [0, N−3].
- **n=0 ⇒ root cut singularity** (gap between root and first sub-root merge).
- **n=N−3 ⇒ leaf-pair cut singularity** (gap between second-smallest and smallest merge).
- Matches `redo1_psi_irrelevance` (psi_irrelevance.pdf).

## Headline counts (strict, no tolerance)

| class | n cells | % |
|---|---|---|
| top_cut (n = 0) | 21 | 12% |
| bottom_cut (n = N−3) | 68 | 38% |
| interior (0 < n < N−3) | 91 | 51% |
| **pinned (top + bottom)** | **89** | **49%** |
| total | 180 | 100% |

Bottom-pinning dominates top-pinning ~3:1, consistent with the "Ψ amplifies smallest-magnitude tail" mechanism: under linearly-spaced merge heights, Δlog₁₀ is largest where t is smallest.

## Distance-from-boundary distribution

Min-distance of argmax to the nearer boundary, across all 180 cells:

| min_dist | n cells | cumulative % |
|---|---|---|
| 0 (strict pinned) | 89 | 49% |
| 1 | 33 | 68% |
| 2 | 22 | 80% |
| 3 | 14 | 88% |
| 4 | 5 | 90% |
| 5 | 6 | 94% |
| ≥ 6 | 11 | 100% |
| **> 10 (deep interior)** | **2** | **1%** |

**Reading.** Only 2/180 cells have argmax more than 10 indices into the interior. The "interior" 91-cell bucket is essentially a near-boundary collar, not a distribution of meaningful Ψ peaks. Whatever scale Ψ picks, it is overwhelmingly within ±5 indices of one extremum.

## Per-band heterogeneity

| band | top_cut | bottom_cut | interior | n | % pinned (strict) |
|---|---|---|---|---|---|
| δ (delta) | 3 | 5 | 22 | 30 | 27% |
| α (alpha) | 6 | 9 | 15 | 30 | 50% |
| β (beta) | 1 | 12 | 17 | 30 | 43% |
| θ (theta) | 9 | 10 | 11 | 30 | 63% |
| γ_low | 1 | 12 | 17 | 30 | 43% |
| γ_high | 1 | 20 | 9 | 30 | 70% |

**Real band gradient.** δ shows the most interior structure (only 27% strict-pinned); γ_high is most boundary-pinned (70%). High-frequency bands favor bottom-pinning more strongly. θ has an unusually high top-cut count (9), driven by patient outliers.

## Per-patient

| patient | top | bot | int | n | % pinned |
|---|---|---|---|---|---|
| Pat_02 | 1 | 11 | 6 | 18 | 67% |
| Pat_03 (1024 Hz outlier) | 3 | 4 | 11 | 18 | 39% |
| Pat_05 | 1 | 6 | 11 | 18 | 39% |
| Pat_06 | 0 | 6 | 12 | 18 | 33% |
| Pat_07 | 2 | 8 | 8 | 18 | 56% |
| Pat_08 | 1 | 7 | 10 | 18 | 44% |
| Pat_10 | 1 | 6 | 11 | 18 | 39% |
| Pat_13 | 3 | 7 | 8 | 18 | 56% |
| Pat_14 | 1 | 8 | 9 | 18 | 50% |
| Pat_15 | 8 | 5 | 5 | 18 | 72% |

Pat_15 is anomalous: the only patient with more top-pinning than bottom-pinning (8 vs 5).

## What changed at the repo level

- `redo2_psi_boundary_audit()` now uses **descending** heights and the **3-class** taxonomy. Output renamed:
  - `tables/psi_argmax_pinning.csv` (was `psi_boundary_audit.csv` — old file path kept-as-stale, will be cleaned up next pass)
  - `tables/psi_argmax_pinning_by_band.csv`
  - `tables/psi_argmax_pinning_by_patient.csv`
  - `figures/psi_argmax_histogram.pdf` regenerated; x-axis is now top-down (0 = root cut, 1 = leaf-pair cut). Tolerance-band shading removed.
- `redo1_psi_irrelevance` already in top-down (done earlier today).
- Shared rcParams block (`plt.rcParams.update`) at top of script controls all font sizes.

## Reconciliation with prior writing

| number | source | meaning | status |
|---|---|---|---|
| 180/180 | coding-agent integration note | mechanism-in-the-limit rhetoric | **drop, not a count** |
| 137/180 (76%) | psi_argmax_histogram.pdf caption | asymmetric 3-cell tolerance band, old ascending convention | superseded |
| 89/180 (49%) | this audit | strict no-tolerance, top-down convention | **canonical** |
| 169/180 (94%) | derived from this audit | within 5 indices of boundary | useful "weak" framing |

Recommendation for the §5.1 figure caption and prose: cite **89/180 (49%) strict** as the headline, paired with **94% within 5 indices of a boundary** as the visual-pile-up qualifier. The "boundary-pinning is universal" framing was overclaim; the honest reading is "boundary-favoring with 1% deep-interior outliers and a real band gradient (δ 27% → γ_high 70%)".

## Files

- `data/audit/section5_v2_round3_redo/tables/psi_argmax_pinning.csv`
- `data/audit/section5_v2_round3_redo/tables/psi_argmax_pinning_by_band.csv`
- `data/audit/section5_v2_round3_redo/tables/psi_argmax_pinning_by_patient.csv`
- `data/audit/section5_v2_round3_redo/figures/psi_argmax_histogram.pdf` (regenerated)
- `scripts/01_compute/audit/audit_round3_section5_redo.py:234-310` (new redo2 body)
