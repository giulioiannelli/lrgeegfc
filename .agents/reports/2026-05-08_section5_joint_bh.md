---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-report
scope: §5.3 closing joint-BH cell count and smallest-q citation
inputs:
  - data/audit/section5_v2_round2/tables/joint_fdr_table.csv
  - scripts/01_compute/audit/audit_round2_section5_v2.py (joint_fdr_table builder, lines 122–142)
---

# §5.3 joint-BH cell count — verification

**Head.** The §5.3 prose has an internal mismatch between **how it lists** the joint family (sums to 48 cells) and **the smallest-q figure it cites** (0.287, which is the 84-cell number from the cached granular table). Both numbers are real but they belong to different family definitions. Anatomy cells from §5.5 are correctly **excluded** from the joint LRG-probe family; they live under a band × region Bonferroni layer that is structurally different.

The cached `joint_fdr_table.csv` was built with a **granular family of 84 cells** (5 KC λ-values per band, 4 Grassmann k-values per band, 3 d_rank distances per band, plus 1 cell each for CTM and VI(k)-best). At m = 84, smallest joint q = 0.287, achieved by four cells: CTM α, CTM low-γ, CTM β, and Grassmann low-γ k=13. **Zero cells survive q ≤ 0.05.**

If the §5.3 prose family of 48 cells (KC restricted to the two structural axes, Grassmann collapsed to one summary per band) is the intended family, the joint BH re-runs to **smallest q = 0.164** at the same four cells; still zero cells survive q ≤ 0.05.

So the headline "no joint survival at q = 0.05" is robust across both family definitions, but the cited smallest-q value (0.287 vs 0.164) is a function of cell count and must be made consistent with the chosen family.

## Family definitions

### Family A: cached granular (m = 84)

What the file `joint_fdr_table.csv` actually contains:

| probe | sub-param scan | cells |
|---|---|---:|
| d_rank | {d_S, d_P, d_F} × 6 bands | 18 |
| kc | {λ ∈ 0, 0.25, 0.5, 0.75, 1} × 6 bands | 30 |
| grassmann | {k ∈ 3, 5, 8, 13} × 6 bands | 24 |
| ctm | {split_vs_drift} × 6 bands | 6 |
| vi_k | {argmin-k summary} × 6 bands | 6 |
| **TOTAL** | | **84** |

Built by `joint_fdr_table()` in `scripts/01_compute/audit/audit_round2_section5_v2.py:122–142`, applying `false_discovery_control(..., method="bh")` over all 84 raw p's.

### Family B: §5.3 prose breakdown (m = 48)

What the §5.3 prose says it adds up:

| probe | sub-param scan | cells |
|---|---|---:|
| KC at λ=0 and λ=1 | 6 bands × 2 axes | 12 |
| CTM | 6 bands | 6 |
| Matrix distances on D(τ) at d_S, d_P, d_F | 6 bands × 3 distances | 18 |
| VI(k) summary scalar | 6 bands × 1 | 6 |
| Grassmann summary scalar | 6 bands × 1 | 6 |
| **TOTAL** | | **48** |

Grassmann summary: I used the per-band best-k (smallest p among k ∈ {3, 5, 8, 13}), matching the spirit of the existing VI(k) `best_k argmin` cell. This is the sharpest legitimate scalar per band. Other choices (k = 13 only, or median over k) only loosen the headline; the dominant low-γ k=13 raw p of 0.0137 controls the family either way.

### Family C: alternative compact (m = 36)

Drop VI(k) (it is a partition-similarity probe rather than an LRG distance probe) and drop Grassmann (it is an eigenmode probe rather than a distance/topology probe), leaving the 3 distance-family probes only:

| probe | cells |
|---|---:|
| d_rank | 18 |
| kc (λ=0, λ=1) | 12 |
| ctm | 6 |
| **TOTAL** | **36** |

## Joint BH outcomes

| family | m | smallest joint q | cells achieving smallest q | n surviving q ≤ 0.05 |
|---|---:|---:|---|---:|
| A (cached, granular) | 84 | **0.287** | CTM α, CTM low-γ, CTM β, Grassmann low-γ k=13 | 0 |
| B (§5.3 prose, 48-cell) | 48 | **0.164** | CTM α, CTM low-γ, CTM β, Grassmann low-γ k=13 | 0 |
| C (compact, distance-only) | 36 | **0.164** | CTM α, CTM low-γ, CTM β | 0 |

The four-cell tie at the top is preserved across all definitions because the four smallest raw p-values are far enough below the rest of the distribution that the BH step-up bottoms out on them under any of the three masks.

## Top-10 cells by raw p (joint family B; q under m=48)

| rank | probe | band | sub-param | raw p | q_joint (m=48) |
|---:|---|---|---|---:|---:|
| 1 | ctm | α | split_vs_drift | 0.006836 | 0.1641 |
| 2 | ctm | low-γ | split_vs_drift | 0.009766 | 0.1641 |
| 3 | grassmann | low-γ | k=13 | 0.013672 | 0.1641 |
| 4 | ctm | β | split_vs_drift | 0.013672 | 0.1641 |
| 5 | kc | β | λ=0 | 0.018555 | 0.1674 |
| 6 | vi_k | low-γ | best_k=68 | 0.024414 | 0.1674 |
| 7 | vi_k | α | best_k=3 | 0.024414 | 0.1674 |
| 8 | d_rank | β | d_P | 0.041992 | 0.1680 |
| 9 | kc | low-γ | λ=1 | 0.041992 | 0.1680 |
| 10 | vi_k | β | best_k=15 | 0.041992 | 0.1680 |

The same ranking holds under family A (84 cells); only the q-values rescale (e.g. rank 1 q goes 0.164 → 0.287).

## Anatomy cells: excluded, correctly

§5.5 anatomy survival uses a Bonferroni correction at m = 48 over (KC trace cluster × DK region) cells from a topology-only KC trace audit. That family is **structurally orthogonal** to the §5.3 LRG-probe joint family:

- §5.3 family: rows are (probe, band, sub-param) cells where each cell is a *cohort one-sided Wilcoxon* on a per-patient probe statistic.
- §5.5 anatomy: rows are (band, DK region) cells where each cell is a *region-level enrichment count* against a sampling-corrected null for KC trace residency.

Mixing the two would (i) double-count β/low-γ cells that already appear in the joint LRG family and (ii) merge a continuous-statistic Wilcoxon family with a count-statistic enrichment family. The right reporting structure is what §5 already does: report joint BH within the LRG probe family in §5.3 and Bonferroni at m = 48 within the anatomy family in §5.5, separately.

## Verdict

The "84 (probe, band) cells" figure in the §5.3 closing paragraph corresponds to the cached granular family and is correct **for that family**. The probe-family breakdown immediately preceding the count (12 + 6 + 18 + 6 + 6 = 48) does **not** sum to 84 — it leaves out the within-probe sub-parameter scans that the cached table exercises (5 KC λ-levels including {0.25, 0.5, 0.75}; 4 Grassmann k-levels at k ∈ {3, 5, 8}).

The author has **two consistent options**:

1. **Keep "84 cells, smallest q = 0.287"** and rewrite the breakdown to expose the granular sub-parameter scans:
   > "we assemble eighty-four cells: KC × {λ = 0, 0.25, 0.5, 0.75, 1} × 6 bands (30); Grassmann × {k = 3, 5, 8, 13} × 6 bands (24); d_rank × {d_S, d_P, d_F} × 6 bands (18); CTM × 6 bands (6); VI(k) summary × 6 bands (6)"

2. **Restrict to the 48-cell family** (KC at the two structural axes, Grassmann collapsed to a single summary scalar per band) and update the cited smallest-q to 0.164 — which is the joint q one gets when the granular sub-parameter scans are pulled out and replaced with one summary per band.

In either case, the qualitative claim — *no cell survives q ≤ 0.05 under any joint-BH scope* — is unchanged; this is the load-bearing statement for the §5.3 closing paragraph and reproduces under all three family definitions tested above. The CTM trace bands lead the rank ordering at the head of the joint family in every definition.

## Sufficient numbers for §5.3 prose

- Family A (cached, m = 84): smallest joint q = **0.287** at 4 cells (CTM α, low-γ, β; Grassmann low-γ k=13). 0 survivors at q ≤ 0.05.
- Family B (m = 48, §5.3 prose breakdown): smallest joint q = **0.164** at the same 4 cells. 0 survivors at q ≤ 0.05.
- Family C (m = 36, distance-only): smallest joint q = **0.164** at 3 CTM cells. 0 survivors at q ≤ 0.05.
- Anatomy cells (§5.5 m = 48 Bonferroni) are correctly excluded — distinct family.
- The top-of-family cells are CTM α (p = 0.0068), CTM low-γ (p = 0.0098), CTM β (p = 0.0137), and Grassmann low-γ k = 13 (p = 0.0137). KC β λ=0 (p = 0.0186) is the next-strongest cell but does not enter the four-way tie at the smallest joint q.
