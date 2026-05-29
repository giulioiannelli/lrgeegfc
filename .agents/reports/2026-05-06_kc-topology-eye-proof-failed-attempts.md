---
name: 2026-05-06_kc-topology-eye-proof-failed-attempts
type: report
era: IMCOH_ABS × COHORT_N10
status: dead-end
created: 2026-05-06
pointers:
  - .agents/reports/archive/2026-05/2026-05-06_section-5-critical-review.md
  - data/audit/section5_v2_round3_redo/figures/kc_topology_eye_proof.pdf
  - data/audit/section5_v2_round3_redo/figures/kc_topology_per_pair_signed.pdf
  - data/audit/section5_v2_round3_redo/tables/kc_topology_eye_proof_panel_summary.csv
  - data/audit/section5_v2_round3_redo/tables/kc_topology_per_pair_signed_summary.csv
  - scripts/01_compute/audit/audit_round3_section5_redo.py  (function `redo6_kc_topology_eye_proof`)
---

# KC λ=0 topology eye-proof — failed visualization attempts

**Head.** A per-pair "eye-proof" of the KC λ=0 cohort trace at β was attempted
in four iterations and **none of them is publishable**. The numerical content
(mean s = +1.38 at β, +1.35 at low_γ, −0.70 at δ) tracks the cohort Wilcoxon
ordering correctly, but the per-pair distributions are sufficiently
diagonal-symmetric and integer-quantized that no rendering of the
`(|m_RPre − m_TT|, |m_TT − m_RPost|)` data carries a reader-visible signal.
The signal that drives `T_KC(λ=0; β) < 0` cohort-wide is real but lives in a
<10% asymmetry of a heavy-mass-at-zero distribution; it is statistically
detectable but visually flat. **Do not use any of these figures in the
manuscript.** The cohort-level KC λ=0 finding stands on the Wilcoxon p-value
and the within-baseline-null control (see
`2026-05-06_section-5-kc-controls.md`); a per-pair eye-proof is not
achievable with this metric.

## What was tried (in order)

### Attempt 1 — patient-resolved scatter of `(|Δ1|, |Δ2|)`

- Layout: 2×3 grid, one band per panel, 412 k points (10 patients × 6 bands ×
  ~6.9 k pairs each), patient-resolved colors with α = 0.18 + jitter 0.30.
- Result: integer-valued depths (m ∈ {0..15}) plus jitter produced a
  checkerboard mess. Overplotting saturated the small-magnitude region;
  larger-magnitude pairs were sparse and uninformative.
- Verdict: unreadable.

### Attempt 2 — signed `(Δ1, Δ2)` restricted to top-decile by max-magnitude

- Goal: spec's fallback formulation if absolute-value scatter fails.
- Result: integer m values + sign-folding produced four perpendicular bands
  forming a "plus sign" pattern in each panel. The corridor shading
  (|y| < |x| = trace) was overwhelmed by overplotting on the integer grid.
- Verdict: unreadable. Deleted.

### Attempt 3 — 2D density heatmap (Blues) with cohort + per-patient medians

- Layout: same 2×3 grid; integer-bin 2D histogram of `(|Δ1|, |Δ2|)`; white
  dashed diagonal; gold cohort star + colored per-patient circles.
- Result: each panel's density is essentially symmetric around the diagonal
  to the eye. The cohort median shifts by 1 unit at most (β at (4, 3) vs δ
  at (4, 4)) on a 0–25 axis — invisible without zooming. The per-patient
  circles cluster on/near the diagonal.
- Verdict: technically clean, geometrically faithful, but the trace direction
  is not visible. The user explicitly flagged: "how on earth does this
  figure prove our point — I can't see this from the picture."

### Attempt 4 — asymmetry heatmap `H[x,y] − H[y,x]` (PiYG diverging)

- Goal: directly visualize which side of the diagonal carries excess
  pair-mass — green = trace excess (more `|Δ1| > |Δ2|`), red = anti excess.
- Result: the asymmetry IS visible if the reader knows what to look for —
  β panel has slightly greener cells below the diagonal and slightly redder
  cells above; δ shows the opposite. But the cell-by-cell magnitudes are
  small relative to the colormap range, so the panels look mostly white
  with scattered low-saturation patches. The trace direction is not legible
  at a glance.
- Verdict: unreadable. The user's exact words: "completely uninformative,
  or at least the differences are minimal — can't really claim any result
  from here."

### Attempt 5 (companion) — 1D histogram of `s = |Δ1| − |Δ2|` per band

- Layout: 2×3 grid; histogram of `s` colored green/red by sign; gold mean line.
- Result: this one DOES show the asymmetry — visibly longer right tail at
  β and low_γ (green mass), longer left tail at δ (red mass). The
  cohort-mean annotations match the cohort Wilcoxon ordering.
- Verdict: the cleanest of the five attempts, but the asymmetry magnitude
  (mean s ≈ 1 unit shift on a ±15 axis with mass concentrated near 0) is
  still modest. Useful as a supplementary plot, not as a headline figure.

## Why all of these fail

The cohort-level finding `T_KC(λ=0; β) < 0` (uncorrected p = 0.019, n = 7/10)
corresponds to a per-pair distribution where:

- Cohort mean of `s = |Δ1| − |Δ2|` is **+1.38** at β (out of an IQR of ~6).
- Cohort median is **0** because of integer quantization with heavy mass at
  zero (~10–14% of pairs have `s = 0`).
- The fraction of pairs with `s > 0` (trace direction) is **48.7%** — under
  50%, but still the highest across bands.

So the trace direction lives in a small (sub-unit on integer m's) shift of
the per-pair distribution. The L2 norm at the cohort level amplifies this
small per-pair shift into a detectable cohort-level Wilcoxon, but at the
per-pair level there is no visible cloud asymmetry to plot. The signal-to-
noise of "asymmetry per pair" is ~0.1; the signal-to-noise of "asymmetry
cumulated into cohort L2" via 6.9 k pairs × 10 patients is what the
Wilcoxon picks up.

## Recommendation

- **Do not** include any per-pair (Δ1, Δ2) scatter, density, or asymmetry
  figure in the manuscript. The KC λ=0 cohort claim stands on the Wilcoxon
  p-value and the within-baseline-null control, not on a visualizable per-
  pair pattern.
- **If a §5.2 visual companion is needed** for the KC λ=0 trace, the
  acceptable options are:
  1. The existing `kc_decomposition.pdf` (cohort-median (T_KC(λ=0),
     T_KC(λ=1)) scatter per band — already in the round-3 redo bundle and
     is the load-bearing figure for §5.2).
  2. A per-patient bar chart of T_KC(λ=0; p, β) showing 7/10 below zero
     directly. This is a cohort-level visualization, not per-pair.
  3. The new `kc_topology_per_pair_signed.pdf` 1D histogram, AT MOST as
     supplementary material with a caveat that the asymmetry is real but
     small.
- **Do not retry per-pair eye-proofs for KC λ=1** without first checking
  whether the asymmetry is bigger there. The `m`-vector at λ=0 is integer-
  valued; the `M`-vector (heights) is continuous and might support a
  cleaner per-pair view, but I have not tried it.

## Files left on disk

These were generated during the failed iterations. Keep for record but do
not cite:

- `data/audit/section5_v2_round3_redo/figures/kc_topology_eye_proof.pdf`
  (current rendering = asymmetry heatmap, attempt 4)
- `data/audit/section5_v2_round3_redo/figures/kc_topology_per_pair_signed.pdf`
  (1D histogram, attempt 5; the cleanest of the lot but still modest)
- `data/audit/section5_v2_round3_redo/tables/kc_pair_depth_shifts.csv.gz`
  (412 k rows of per-pair m-shifts; useful raw data, not a figure)
- `data/audit/section5_v2_round3_redo/tables/kc_topology_eye_proof_panel_summary.csv`
- `data/audit/section5_v2_round3_redo/tables/kc_topology_eye_proof_per_patient_medians.csv`
- `data/audit/section5_v2_round3_redo/tables/kc_topology_per_pair_signed_summary.csv`

## Numerical verdict (still valid, independent of the figures)

| band | mean s | verdict | matches Wilcoxon? |
|:--|:--:|:--|:--|
| δ | −0.70 | anti | ✓ (cohort 2/10, p = 0.96) |
| θ | −0.10 | null | ✓ (4/10, p = 0.75) |
| α | −0.40 | null | ✓ (6/10, p = 0.46) |
| **β** | **+1.38** | **trace** | ✓ (7/10, p = 0.019) |
| low_γ | +1.35 | trace | divergence with cohort λ=0 (4/10) but consistent with the cohort λ=1 trace finding |
| γ_h | +0.52 | trace (mild) | borderline (5/10, p = 0.35) |

The numbers are sound. The figures are not.
