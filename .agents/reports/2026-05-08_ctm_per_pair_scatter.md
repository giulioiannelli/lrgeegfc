---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: figure-build-report
scope: section_5_3_ctm_per_pair_scatter
---

# §5.3 per-pair scatter — selection, verification, sanity

**Head.** Three-panel scatter of Δ_task(i,j) vs Δ_rest(i,j) over contact pairs at three (patient, band) cells: **β @ Pat_06** (ρ_split = +0.259, |d| to cohort median +0.222 = 0.0367); **α @ Pat_05** (ρ_split = +0.071, |d| = 0.0431); **θ @ Pat_08** (ρ_split = -0.051, |d| = 0.0018). Same-probe pairs in light gray, cross-probe pairs in navy on top. Axes clipped to ±0.345 (98th percentile of |Δ|).

## Selection table

| panel | patient | band | ρ_split | cohort med | |d| | ρ_xprobe | N | same | cross |
|---|---|---|---|---|---|---|---|---|---|
| β | Pat_06 | beta | +0.259 | +0.222 | 0.0367 | +0.249 | 6555 | 485 | 6070 |
| α | Pat_05 | alpha | +0.071 | +0.115 | 0.0431 | +0.057 | 6903 | 500 | 6403 |
| θ | Pat_08 | theta | -0.051 | -0.049 | 0.0018 | -0.048 | 7140 | 608 | 6532 |

## Null-band rationale

Cohort medians at candidate null bands: θ -0.0492, γ_h -0.0135; δ +0.0313 (for reference). **θ chosen** even though |median θ| = 0.0492 is larger than |median γ_h| = 0.0135: ρ_split measures the rank agreement (diagonal tilt) between Δ_task and Δ_rest, not the magnitude of either. γ_h cohort-median exemplars produce a pathologically tight cloud near the origin, which reads as "no dynamics" rather than "no rank agreement". θ has Δ spread comparable to α and β at the cohort-median exemplar, so the visual contrast against the β / α panels isolates the absence of diagonal tilt at comparable cloud size — which is the property ρ_split actually tests.

## ρ verification (Spearman recomputed from the npz)

| panel | ρ_split (CSV) | ρ_split (npz) | ρ_xprobe (CSV) | ρ_xprobe (npz) | match? |
|---|---|---|---|---|---|
| β | +0.2590 | +0.2590 | +0.2486 | +0.2486 | ✓ |
| α | +0.0715 | +0.0715 | +0.0566 | +0.0566 | ✓ |
| θ | -0.0510 | -0.0510 | -0.0479 | -0.0479 | ✓ |

## Sanity checks

1. **β panel diagonal tilt** (full): ρ = +0.259. Expected positive ✓
2. **α panel diagonal tilt** (full): ρ = +0.071. Expected positive ✓
3. **θ panel cloud shape** (full): ρ = -0.051. Expected near-zero ✓
4. **β cross-probe subset tilt**: ρ_xprobe = +0.249. Expected positive ✓
5. **α cross-probe subset tilt**: ρ_xprobe = +0.057. Expected positive ✓
6. **θ cross-probe subset shape**: ρ_xprobe = -0.048. Expected near-zero ✓

## Files

- Figure: `data/audit/ctm_per_pair_scatter/figures/manuscript_ctm_per_pair_scatter.pdf`
- Selection: `data/audit/ctm_per_pair_scatter/tables/selection.csv`
- Source npz: `data/reports/imcoh_continuous_trace/per_pair_split/{Pat}_{band}.npz`
- Cohort table: `data/audit/ctm_triangle/Td_per_patient_per_band.csv`
- Build script: `scripts/01_compute/audit/audit_55_ctm_per_pair_scatter.py`
