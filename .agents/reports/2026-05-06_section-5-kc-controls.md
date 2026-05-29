---
name: 2026-05-06_section-5-kc-controls
type: report
era: IMCOH_ABS × COHORT_N10
status: locked
created: 2026-05-06
updated: 2026-05-06
pointers:
  - .agents/reports/archive/2026-05/2026-05-06_section-5-critical-review.md
  - .agents/reports/2026-05-06_section-5-writing-brief.md
  - data/audit/section5_v2_kc_controls/cohort_summary.csv
  - data/audit/section5_v2_kc_controls/per_patient_table.csv
  - data/audit/section5_v2_kc_controls/pat03_dropout.csv
  - data/audit/lrg_global_probe_controls/cohort_controls_summary.csv
  - scripts/01_compute/audit/audit_46_kc_section5_controls.py
  - data/outputs/figures/section_5_lrg_trace/lrg_controls/kc_controls_summary.pdf
---

# Section 5.2 — KC controls (split-baseline null, cross-probe, Pat_03 dropout)

**Head.** Under the same control regime CTM carries in §5.3, β at KC λ=0
(topology) and β at KC λ=1 (heights) **both acquire a fully controlled
cohort backing**: within-baseline-null Wilcoxon p = 0.001 each, within-probe
BH-FDR (m = 12) q = 0.006 each, 10/10 patients below their own
half-baseline null, **Pat_03-robust** under the controlled test (p_drop =
0.002 for both). The cross-probe restriction agrees in direction (β λ=0
p = 0.024, β λ=1 p = 0.042) but does not pass BH (q = 0.146 / 0.168);
it is reported as a sanity-check that same-probe geometry does not drive
the finding. low-γ at λ=1 sits below uncorrected p = 0.05 under all
three control layers (absolute 0.042, null 0.014, cross-probe 0.024)
and the cross-probe restriction *strengthens* it (8/10 → 9/10), but
within-probe BH at m=12 is missed (q_null = 0.055, q_xprobe = 0.146).
Section 5.2 can therefore promote β to a controlled cohort headline
parallel to §5.3 CTM α/β/low-γ; low-γ heights stays as a directional
companion, not a controlled headline; α / δ / θ / γ_h remain silent on
KC under every control. This pass changes no number in the §5.2 prose
but lets the manuscript drop the methodological-asymmetry caveat
between §5.2 and §5.3.

## What was done

### The three control layers

Mirrors the CTM control regime of §5.3:

1. **Absolute Wilcoxon** (existing baseline; the §5.2 round-2 table).
   `T_KC(p, b, λ) = d_KC(λ; T^TT, T^RPost) − d_KC(λ; T^RPre, T^TT)`,
   one-sided paired Wilcoxon `T_KC < 0` on n = 10. Source:
   `data/reports/section_5_lrg_trace/03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv`.

2. **Within-baseline null** (CTM ρ_drift analog, **already computed**
   in `audit_41_lrg_global_probe_controls.py`). Build dendrograms
   `T^RPreA`, `T^RPreB`, `T^RPostA` on the half-cache, define
   `T_KC_null(λ; p, b) = d_KC(λ; T^RPreB, T^RPostA) − d_KC(λ; T^RPreA, T^RPreB)`,
   and run paired one-sided Wilcoxon `T_KC < T_KC_null`. Source:
   `data/audit/lrg_global_probe_controls/{kc_null,cohort_controls_summary}.csv`.

3. **Cross-probe KC** (CTM ρ_xprobe analog, **new in this pass**, computed
   by `scripts/01_compute/audit/audit_46_kc_section5_controls.py`). Per
   patient, leaves are tagged with their sEEG probe via
   `parse_seeg_label(channel_labels.csv[i])`. The KC `(m, M)` per-leaf-pair
   vectors of each phase are computed on the full dendrogram, then the L2
   norm of the lambda-blended difference is restricted to cross-probe
   leaf pairs:

   ~~~
   d_KC_xprobe(λ; T_a, T_b) = || ((1−λ)·m_norm + λ·M_norm)|_xprobe ||₂
   ~~~

   with the same pooled-max normalization the library `kc_distance`
   uses (so the cross-probe number lives in the same scale as the
   unrestricted KC). Triangle and Wilcoxon as in (1) and (2).

Cross-probe pair counts per patient:

| patient | n_leaves | n_pairs | n_xprobe_pairs | xprobe % |
|:--|:--:|:--:|:--:|:--:|
| Pat_02 | 117 | 6786 | 6169 | 90.9% |
| Pat_03 | 122 | 7381 | 6780 | 91.9% |
| Pat_05 | 118 | 6903 | 6403 | 92.8% |
| Pat_06 | 115 | 6555 | 6070 | 92.6% |
| Pat_07 | 116 | 6670 | 6108 | 91.6% |
| Pat_08 | 120 | 7140 | 6532 | 91.5% |
| Pat_10 | 113 | 6328 | 5780 | 91.3% |
| Pat_13 | 119 | 7021 | 6489 | 92.4% |
| Pat_14 | 119 | 7021 | 6473 | 92.2% |
| Pat_15 | 118 | 6903 | 6277 | 90.9% |

Same-probe pairs are ~8–9% of all pairs (the diagonal-like blocks of
contacts on the same shaft).

### Multiple-comparisons framework

Within-probe Benjamini-Hochberg FDR at **m = 12** (six bands × two
λ values, λ = 0 and λ = 1) per control layer. λ = 0.5 is reported
uncorrected for completeness; the §5.2 prose only cites λ = 0 and λ = 1.
This mirrors CTM's within-probe BH at m = 6 in §5.3.

The n = 10 paired Wilcoxon one-sided floor is **1/2¹⁰ = 0.000977**;
when this is hit, n_real_below_null = 10/10 (or 0/10 in the wrong
direction). It is not differentiable from any other 10/10 outcome
under the same test. The tables flag this where it occurs.

## Cohort verdict — 12 cells × 3 layers

Per-cell numbers (from `cohort_summary.csv`); rank-biserial `r_rb` ∈
[−1, +1] alongside each p-value, BH-FDR q within m = 12 per layer.
Symbols: `↓` = direction is trace (p one-sided "less"); bold = p ≤ 0.05;
**bold red-equivalent** marks BH-FDR q ≤ 0.05.

### λ = 0 (pure topology)

| band | n_abs↓ | p_abs | q_abs | n_null↓ | p_null | **q_null** | n_xpr↓ | p_xpr | q_xpr |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| δ | 2 | 0.958 | 0.958 | 4 | 0.539 | 0.539 | 2 | 0.958 | 0.958 |
| θ | 4 | 0.754 | 0.886 | 6 | 0.161 | 0.259 | 3 | 0.754 | 0.886 |
| α | 6 | 0.461 | 0.857 | 5 | 0.138 | 0.259 | 6 | 0.500 | 0.886 |
| **β** | 7 | **0.019** | 0.168 | **10** | **0.001** | **0.006** ✓ | 7 | **0.024** | 0.146 |
| low-γ | 4 | 0.577 | 0.866 | 6 | 0.216 | 0.259 | 4 | 0.652 | 0.886 |
| γ_h | 5 | 0.348 | 0.834 | 6 | 0.188 | 0.259 | 5 | 0.313 | 0.750 |

### λ = 1 (pure heights)

| band | n_abs↓ | p_abs | q_abs | n_null↓ | p_null | **q_null** | n_xpr↓ | p_xpr | q_xpr |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| δ | 4 | 0.812 | 0.886 | 6 | 0.423 | 0.461 | 4 | 0.812 | 0.886 |
| θ | 5 | 0.754 | 0.886 | 5 | 0.161 | 0.259 | 5 | 0.754 | 0.886 |
| α | 6 | 0.500 | 0.857 | 6 | 0.216 | 0.259 | 6 | 0.539 | 0.886 |
| **β** | 8 | **0.042** | 0.168 | **10** | **0.001** | **0.006** ✓ | 8 | **0.042** | 0.168 |
| **low-γ** | 8 | **0.042** | 0.168 | 8 | **0.014** | 0.055 | **9** | **0.024** | 0.146 |
| γ_h | 5 | 0.216 | 0.647 | 6 | 0.138 | 0.259 | 5 | 0.216 | 0.647 |

### Lock-in per cell

| (band, λ) | within-baseline null at q ≤ 0.05 | uncorrected p ≤ 0.05 (any layer) | controlled-headline status |
|:--|:--:|:--:|:--|
| β, λ = 0 | ✓ q = 0.006 | ✓ all three layers | **Controlled headline; cross-probe and Pat_03-dropout-robust under null.** |
| β, λ = 1 | ✓ q = 0.006 | ✓ all three layers | **Controlled headline; Pat_03-robust under null.** |
| low-γ, λ = 1 | q = 0.055 (just misses) | ✓ all three layers | Directional, uncorrected only. Cross-probe *strengthens* (9/10). |
| α, λ ∈ {0, 1} | — | — | Silent. |
| δ / θ / γ_h | — | — | Silent or anti-trace. |

The two cells that survive within-probe BH at m = 12 are exactly the two
cells that survive Bonferroni at m = 42 in `audit_41_lrg_global_probe_controls`
(the original within-baseline-null table; β λ = 0 and β λ = 1 are 2 of
the 5 cells passing m = 42 there). This is internal consistency, not
a new finding.

## Pat_03 dropout × 3 control layers — 18 cells

For each (band, λ): full cohort vs Pat_03-dropped (n = 9). Source:
`pat03_dropout.csv`.

| band | λ | abs p_full | abs p_drop | null p_full | null p_drop | xpr p_full | xpr p_drop |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| δ | 0 | 0.958 | 0.936 | 0.539 | 0.674 | 0.958 | 0.936 |
| δ | 0.5 | 0.278 | 0.410 | 0.278 | 0.285 | 0.278 | 0.410 |
| δ | 1 | 0.812 | 0.787 | 0.423 | 0.545 | 0.812 | 0.787 |
| θ | 0 | 0.754 | 0.674 | 0.161 | 0.180 | 0.754 | 0.674 |
| θ | 0.5 | 0.065 | 0.125 | 0.097 | 0.180 | 0.065 | 0.125 |
| θ | 1 | 0.754 | 0.850 | 0.161 | 0.285 | 0.754 | 0.850 |
| α | 0 | 0.461 | 0.367 | 0.138 | 0.102 | 0.500 | 0.410 |
| α | 0.5 | 0.688 | 0.500 | 0.577 | 0.410 | 0.722 | 0.545 |
| α | 1 | 0.500 | 0.285 | 0.216 | 0.125 | 0.539 | 0.326 |
| **β** | **0** | **0.019** | **0.037** ✓ | **0.001** | **0.002** ✓ | **0.024** | 0.049 ✓ |
| **β** | **0.5** | 0.188 | 0.326 | **0.001** | **0.002** ✓ | 0.138 | 0.248 |
| **β** | **1** | **0.042** | 0.082 ✗ | **0.001** | **0.002** ✓ | **0.042** | 0.082 ✗ |
| low-γ | 0 | 0.577 | 0.633 | 0.216 | 0.180 | 0.652 | 0.674 |
| low-γ | 0.5 | 0.903 | 0.875 | 0.722 | 0.633 | 0.920 | 0.875 |
| **low-γ** | **1** | **0.042** | 0.064 ✗ | **0.014** | **0.027** ✓ | **0.024** | **0.037** ✓ |
| γ_h | 0 | 0.348 | 0.410 | 0.188 | 0.326 | 0.313 | 0.367 |
| γ_h | 0.5 | 0.313 | 0.410 | 0.278 | 0.455 | 0.313 | 0.410 |
| γ_h | 1 | 0.216 | 0.326 | 0.138 | 0.248 | 0.216 | 0.326 |

(✓ = stays at uncorrected p ≤ 0.05 after Pat_03-dropout; ✗ = loses
uncorrected p ≤ 0.05.)

### Reading

- **β λ = 0**: robust under all three layers. The absolute and
  cross-probe Wilcoxons keep p ≤ 0.05 at n = 9 (0.037 and 0.049
  respectively); the within-baseline null amplifies to p = 0.002 (9/9
  patients still below their own null).
- **β λ = 1**: under absolute and cross-probe, fragile to Pat_03 (loses
  p ≤ 0.05 at n = 9). Under within-baseline null, **fully robust**
  (0.001 → 0.002), 9/9 below own null. This was the §5.2 manuscript's
  weak spot; the controlled test promotes it from "Pat_03-fragile" to
  "Pat_03-robust under control" because every patient — including
  Pat_03 — has a high own-null floor at λ = 1, so the paired contrast
  is consistent across the cohort.
- **low-γ λ = 1**: absolute fails Pat_03 dropout (0.042 → 0.064);
  within-baseline null and cross-probe both survive at uncorrected
  level (0.027 / 0.037). Companion finding, not promoted.

## What §5.2 should now say

Two prose updates land in the §5.2 manuscript paragraph and the
critical-review lock-in. Both are additive — the round-2 / round-3
numbers and figures stand.

### A. Add a controls paragraph after the paragraph that opens the
geometric decomposition and before "the geometric decomposition of
\FigRef{fig:kc_decomposition} is descriptive". Suggested text:

> Under a within-baseline drift-floor null built from rest-half
> dendrograms (`T^RPreB → T^RPostA` with `T^RPreA → T^RPreB`,
> mirroring the CTM `ρ_drift` of §5.3), the β trace at λ = 0 and at
> λ = 1 each carry 10/10 patients below their own half-baseline null,
> with paired one-sided Wilcoxon `T_KC < T_KC_null` at p = 0.001
> (the n = 10 floor 1/2¹⁰ = 0.000977). Within-probe BH-FDR at m = 12
> (six bands × two λ values) leaves both cells surviving at
> q = 0.006. Restricting the L2 norm to cross-probe leaf pairs only
> agrees in direction (β λ = 0 p = 0.024, β λ = 1 p = 0.042) but
> does not pass BH (q = 0.146 / 0.168); the same-probe geometry does
> not drive the β finding. low-γ heights (λ = 1) sits below
> uncorrected p = 0.05 under all three layers (absolute 0.042,
> null 0.014, cross-probe 0.024) and is *strengthened* by the
> cross-probe restriction (8/10 → 9/10), but within-probe BH at
> m = 12 is missed (q_null = 0.055, q_xprobe = 0.146); reported as
> a directional companion, not a controlled headline.

### B. Replace the absolute-Wilcoxon Pat_03 caveat for β λ = 1 with the
controlled-test result. Current §5.2 text:

> "(\(p = 0.042\); n = 9 dropout p = 0.082)"

becomes:

> "(\(p = 0.042\); under the within-baseline drift-floor null
> p = 0.001, BH-corrected q = 0.006, Pat_03-dropout p = 0.002 — robust)"

Section 5.3 already carries this structure for CTM. After this edit,
§5.2 and §5.3 carry parallel control regimes and the manuscript no
longer mixes controlled and uncontrolled cohort claims across
adjacent subsections.

## Files in this report

- `data/audit/section5_v2_kc_controls/per_patient_table.csv` — 180
  rows (10 patients × 6 bands × 3 λ) with real, null, xprobe T_KC and
  per-row indicator booleans.
- `data/audit/section5_v2_kc_controls/cohort_summary.csv` — 18 rows
  (band × λ) with all three layers' Wilcoxon p, n_below, rank-biserial,
  BH-corrected q (m = 12 over λ = 0/λ = 1).
- `data/audit/section5_v2_kc_controls/pat03_dropout.csv` — 18 rows
  with p_full and p_drop for all three control layers.
- `data/outputs/figures/section_5_lrg_trace/lrg_controls/kc_controls_summary.pdf`
  — companion figure, two-panel `−log₁₀ p` bar chart per band per
  control layer, λ = 0 and λ = 1.
- `scripts/01_compute/audit/audit_46_kc_section5_controls.py` — the
  cross-probe + Pat_03-dropout compute. Re-runnable in `lapbrain`.

## What this DOES and DOES NOT establish

**Establishes:**
- Section 5.2's β trace at topology (λ = 0) and heights (λ = 1) carry
  controlled cohort backing under within-probe BH-FDR m = 12, parallel
  to §5.3 CTM.
- The β λ = 1 weakness under Pat_03-dropout (absolute Wilcoxon 0.042
  → 0.082) is upgraded to robustness under the within-baseline-null
  control (0.001 → 0.002).
- Same-probe geometry does not drive β at either λ (cross-probe Wilcoxon
  agrees in direction with the unrestricted Wilcoxon).
- low-γ heights remains directional / companion-level. Cross-probe
  restriction strengthens it (n_xpr 9/10 vs n_full 8/10), but
  within-probe BH at m = 12 is missed.

**Does NOT establish:**
- A new statistically significant cohort finding beyond β. low-γ
  heights does not promote. α / δ / θ / γ_h are silent under every
  layer.
- Joint MTC across §5.2 + §5.3 + §5.4 is unchanged. The conservative
  joint-Bonferroni statement in `2026-05-06_section-5-critical-review.md`
  (m = 48 over 8 probes × 6 bands) still applies; this pass adds within-
  *probe* control parity, not joint-MTC promotion.
- Anatomical specificity, multi-scale localization, or per-leaf
  contributions are out of scope for §5.2 and unchanged.

## Reproducibility snippet

```python
import pandas as pd
df = pd.read_csv(
    "data/audit/section5_v2_kc_controls/cohort_summary.csv"
)
print(df[df["lam"].isin([0.0, 1.0])][[
    "band", "lam", "n_abs_below_zero", "p_abs", "q_abs_within_m12",
    "n_real_below_null", "p_null", "q_null_within_m12",
    "n_xprobe_below_zero", "p_xprobe", "q_xprobe_within_m12",
]].round(4).to_string(index=False))
```
