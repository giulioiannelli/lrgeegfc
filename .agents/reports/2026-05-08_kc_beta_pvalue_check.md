---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-report
scope: §5.2 KC decomposition β-band p-value ordering (λ=0 topology vs λ=1 heights)
inputs:
  - data/audit/kc_triangle/Td_per_patient_per_band_lambda.csv
  - data/audit/section5_v2_round3_redo/tables/kc_decomposition_points.csv
  - data/audit/section5_v2_round2/tables/kc_decomposition_summary.csv
---

# KC β p-value ordering — verification

**Head.** The smaller-count-smaller-p apparent paradox is real and reproduces from the cached per-patient T_KC values. It resolves cleanly into a **magnitude effect**: at λ=0 (topology) the 3 positive (counter-imprint) patients carry tiny |T_KC| (mean 1.02, ranks {1, 2, 4}), so the negative-direction signed-rank sum claims W⁻=48/55 (W⁺=7). At λ=1 (heights) the 2 positive patients carry larger |T_KC| (mean 1.83, ranks {4, 6}), so W⁻=45/55 (W⁺=10) despite one extra negative patient. Wilcoxon rewards the cleaner sign-by-magnitude alignment at λ=0, hence the smaller p (0.0186 vs 0.0420). The hypothesis the writing agent posed is **confirmed**.

## Per-patient T_KC at β, τ=1/λ_max

Source: `data/audit/kc_triangle/Td_per_patient_per_band_lambda.csv` (rows where `band=beta`, `lam ∈ {0.0, 1.0}`).

| patient | T_KC (λ=0, topology) | T_KC (λ=1, heights) | sgn₀ | sgn₁ |
|---|---:|---:|:-:|:-:|
| Pat_02 |   0.5081 |  −7.6450 | + | − |
| Pat_03 | −19.7298 |  −2.7111 | − | − |
| Pat_05 |   0.1749 |  −0.3326 | + | − |
| Pat_06 |  −5.0902 |  −2.2197 | − | − |
| Pat_07 |  −2.4326 |  −1.3483 | − | − |
| Pat_08 |  −5.0689 |  −6.2626 | − | − |
| Pat_10 |  −2.4519 |  −1.4171 | − | − |
| Pat_13 |  −3.1590 |  −1.0314 | − | − |
| Pat_14 |  −1.3667 |   1.2500 | − | + |
| Pat_15 |   2.3806 |   2.4069 | + | + |

## Cohort statistics

| quantity | λ=0 (topology) | λ=1 (heights) |
|---|---:|---:|
| n_neg / 10 | **7** | **8** |
| median \|T_KC\| (cohort) | 2.4422 | 1.8184 |
| mean \|T_KC\| (cohort) | 4.2363 | 2.6625 |
| median \|T_KC\| over negatives | **3.1590** | **1.8184** |
| mean \|T_KC\| over negatives | **5.6142** | **2.8710** |
| median \|T_KC\| over positives | 0.5081 | 1.8285 |
| mean \|T_KC\| over positives | **1.0212** | **1.8285** |
| sum of signed ranks W⁺ (anti-H1) | **7.00** | **10.00** |
| sum of signed ranks W⁻ | 48.00 | 45.00 |
| Wilcoxon W (= min) | 7.0000 | 10.0000 |
| one-sided p (H₁: T_KC < 0) | **0.018555** | **0.041992** |

Numbers match the canonical summary `kc_decomposition_summary.csv` (0.0185546875 / 0.0419921875) to all printed digits — these are exact rational p-values from the n=10 null permutation distribution.

## Diagnosis

1. **At λ=0 the positive patients are weakly positive.** Pat_02 (+0.51), Pat_05 (+0.17), Pat_15 (+2.38) carry |T| ranks 1, 2, 4. They contribute only W⁺ = 7 of the total 55 rank mass.
2. **At λ=1 the positive patients are middling.** Pat_14 (+1.25), Pat_15 (+2.41) carry |T| ranks 4 and 6. They contribute W⁺ = 10 of 55.
3. **At λ=0 a single dominant patient anchors the negative side.** Pat_03 has |T_KC| = 19.73 — the largest single value in either column — and it points to imprint, claiming rank 10.
4. **Median |T_KC| over the negative subgroup is 1.74× larger at λ=0 than λ=1** (3.16 vs 1.82). The negatives are not just more numerous at λ=1; they are individually smaller.

These four facts compound: λ=0 has fewer signs in the imprint direction, but those signs are bigger and the dissenting signs are tiny. Wilcoxon’s signed-rank statistic measures exactly that compound, not the sign count.

## Verdict

The published numbers (7/10 at p=0.0186 for λ=0; 8/10 at p=0.0420 for λ=1) are correct and reproduce from the cached LRG outputs. The hypothesis that magnitude rather than sign count drives the ordering is confirmed: the topology axis is "carried by louder patients with quieter dissenters", which is the precise asymmetry the signed-rank test rewards.

## Sufficient numbers for §5.2 prose

- λ=0: 7/10 imprint, W=7, p=0.0186; W⁺=7 vs W⁻=48; positives have mean |T|=1.02 vs negatives 5.61 (5.5× louder).
- λ=1: 8/10 imprint, W=10, p=0.0420; W⁺=10 vs W⁻=45; positives have mean |T|=1.83 vs negatives 2.87 (1.6× louder).
- Pat_03 at λ=0 (|T|=19.73) is the dominant single contributor; the topology effect is not Pat_03-only — Pat_06 (5.09), Pat_08 (5.07), Pat_13 (3.16), Pat_07 (2.43), Pat_10 (2.45) all line up — but Pat_03 is the rank-10 anchor.
