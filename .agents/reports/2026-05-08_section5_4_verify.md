---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-report
scope: §5.4 multiscale companion views — Grassmann formula convention, VI(k) cohort claims, Grassmann cohort claims, singleton mask robustness, Pat_03/07/15 dropout
inputs:
  - scripts/01_compute/audit/audit_46_grassmann_principal_angles.py
  - data/audit/section5_v2_round3_redo/tables/grassmann_principal_angles_per_patient.csv
  - data/audit/section5_v2_round3_redo/tables/grassmann_principal_angles_cohort.csv
  - data/audit/section5_v2_round2/tables/grassmann_k_sweep.csv
  - data/audit/vi_triangle_heatmap/cohort_n_trace_band_k.csv
  - data/audit/vi_triangle_heatmap/cohort_band_summary.csv
  - data/audit/vi_triangle_heatmap/T_VI_per_patient_per_band_per_k.csv
  - data/audit/vi_triangle_heatmap/singleton_share_band_k.csv
---

# §5.4 — multiscale companion views — verification

**Head.** The Grassmann pipeline computes **option (c)**: `d_G = sqrt(sum sin^2(theta_i))`. The §5.4 prose's pair of equations is mutually inconsistent — only `d_G^2 = sum sin^2(theta_i)` matches the code. The Frobenius identity in the prose would give `sqrt(2 · sum sin^2)`, off by a factor of √2. The principal-angle heatmap uses signed `Δθ_i(k) = θ_pre→tt − θ_tt→post` in **radians** (positive = trace direction); cohort-median is on the angles, not on sines.

The VI(k) prose is mostly correct on peak locations and rough span shapes but **the β "contiguous span k=17..35" claim is not strictly contiguous** — the run is k=17..32, drops below 7 at k=33, then k=34..35; the "k=17..35" framing should be qualified. The "low-γ concentrated band k=48..54" matches a longer 12-cell contiguous run k=45..56 (the longest contiguous n_trace≥7 span anywhere in any band).

**The Grassmann claim "low-γ k=13 with p=0.014 is the only (band,k) cell at uncorrected p<0.05" is wrong**. Even within the restricted k ∈ {13, 20, 30} k_sweep cache, low-γ k=30 is also at p=0.042 (<0.05). On the full continuous k ∈ {2..80} scan, ~35 cells fall under p=0.05 — dominated by β at k ≥ 40 and low-γ at k=10..28.

Singleton mask robustness: cohort N ranges 113 (Pat_10) to 122 (Pat_03), median 118; the longest contiguous trace spans (β k=17..32, low-γ k=45..56) survive thresholds 0.4/0.5/0.6 unchanged. Pat_03/07/15 dropouts do not collapse VI(k) trace bands at peak cells; Pat_15 dropout at β k=22 yields **9/9** unanimous trace.

---

## A. Grassmann convention

### Code (`scripts/01_compute/audit/audit_46_grassmann_principal_angles.py:73–74`)

```python
def chordal_from_angles(theta: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.sin(theta) ** 2)))
```

The chordal distance is computed as **option (c)**:

```
d_G(U, U') = sqrt(sum_{i=1..k} sin^2(theta_i))
```

This matches `d_G = sqrt(k − sum sigma_i^2)` (line 4 of script docstring), since `sigma_i = cos(theta_i)` and `sum cos^2 + sum sin^2 = k`.

### §5.4 prose — inconsistency

The prose currently states both:

```
d_G(U, U') = ||U U^T − U' U'^T||_F        ← Frobenius identity
d_G^2(U, U') = sum_{i=1..k} sin^2(theta_i) ← option (c)
```

These differ by a factor of √2: `||U U^T − U' U'^T||_F^2 = 2 · sum sin^2(theta_i)`.

**Correction:** Replace the Frobenius identity by either
- `d_G(U, U') = sqrt(sum sin^2(theta_i))` (matches code), or
- `d_G(U, U') = (1/√2) · ||U U^T − U' U'^T||_F` (equivalent to the code).

The "raw chordal" `d_G = ||P − P'||_F` (option a) is NOT what the code computes.

### Principal-angle heatmap convention

`grassmann_principal_angles.pdf` shows `delta_theta_median_radians`, i.e. cohort-median of `Δθ_i(k) = θ^{Pre→TT}_i − θ^{TT→Post}_i` in **radians**, with positive values indicating the trace direction (TT→Post subspace closer to TT than Pre→TT). The cohort aggregation is on **angles**, not on sines or sin². ✓

The diagonal `i ≈ k` reading in the prose is the correct way to read the heatmap.

---

## B. VI(k) cohort verification (mask threshold 0.5)

| band | masked cells | unmasked k range | peak n_trace | peak k | contiguous spans n_trace ≥ 7 |
|---|---:|---|---:|---:|---|
| δ        | 42 | 2..79 | **6** | 16 | **NONE** |
| θ        | 39 | 2..82 | 8 | 44 | k=34..36, 39..41, 44, 58 |
| α        | 40 | 2..81 | **9** | **3** | k=3, 13..15, 18..21, 23..24, 28, 34..37 |
| β        | 42 | 2..79 | **9** | **22** | **k=17..32**, 34..35, 60, 62, 67..68, 72..78 |
| low-γ    | 42 | 2..79 | **8** | **48** | k=25..27, 35..36, 40, **k=45..56**, 58, 60..64, 67..69, 72..77 |
| high-γ   | 42 | 2..79 | 8 | 79 | k=52, 54, 74..76, 78..79 |

The longest contiguous n_trace ≥ 7 run anywhere in any band is **low-γ k=45..56 (12 cells)**, followed by **β k=17..32 (16 cells)** — wait, β actually wins (16 vs 12). β is the band with the longest contiguous trace span.

### Verification of specific §5.4 prose claims

| prose claim | verdict |
|---|---|
| β contiguous span k=17..35 with n_trace ≥ 7 | **partial**. Strictly contiguous run is k=17..32; n_trace=6 at k=33; resumes at k=34..35. "k=17..35" reads as one block but has a one-cell dip at k=33. Recommend: "k=17..32 (with k=34..35 returning to ≥ 7 after a single-cell dip at k=33)". |
| β peak n_trace = 9 at k=22 | ✓ |
| β sparser pro-trace continues to k ~ 78 | ✓ (k=60, 62, 67..68, 72..78) |
| α pro-trace concentrated k=3..37; peak n_trace = 9 at k=3 | ✓ peak; "concentrated k=3..37" is fragmented (gaps at k=4..12, 16..17, 22, 25..27, 29..33, 38..). Recommend: "non-contiguous pro-trace cells across k=3..37 with peak at k=3 and additional plateaus at k=13..15, 18..21, 34..37". |
| α contiguity broken by isolated null cells | ✓ |
| low-γ pro-trace at coarser k=25..77 | ✓ rough framing |
| low-γ concentrated band k=48..54 reaching n_trace = 8 | ✓ — sits inside the larger contiguous span k=45..56 (peak at k=48). The "k=48..54" framing is a sub-window of the actual span. |
| δ anti-trace at fine k, null elsewhere, no contiguous pro-trace | ✓ |
| θ isolated peaks without contiguity | **slight overclaim**: there are short contiguous clusters k=34..36 and k=39..41, plus singletons at k=44 and k=58. Recommend: "two short contiguous clusters at k=34..36 and k=39..41 with peak n_trace=8 at k=44, but no extended span". |
| high-γ sparse pro-trace at the coarse edge of the unmasked range | ✓ (spans at k=52, 54, 74..76, 78..79) |

---

## C. Grassmann cohort verification

### C.3 — per-k Wilcoxon on chordal scalar (continuous k ∈ {2..80})

**Cells with raw p < 0.05 (one-sided H1: T_E1 < 0):**

| band | k | n_trace | raw p | T_median |
|---|---:|---:|---:|---:|
| low-γ | **13** | 8/10 | **0.0137** | −0.377 |
| low-γ | **14** | 9/10 | **0.0137** | −0.339 |
| low-γ | 12 | 7/10 | 0.0186 | −0.372 |
| α | 12 | 8/10 | 0.0322 | −0.207 |
| β | 70 | 8/10 | 0.0322 | −0.297 |
| β | 41, 42 | 7/10, 7/10 | 0.0322 | −0.46, −0.50 |
| β | 65 | 7/10 | 0.0322 | −0.319 |
| low-γ | 16 | 8/10 | 0.0322 | −0.278 |
| low-γ | 10, 15 | 8/10 | 0.0322 | −0.114, −0.320 |
| β | 78, 79 | 8/10, 8/10 | 0.0322 | −0.44, −0.40 |
| β | 43, 55, 57, 58, 63, 67, 68, 69, 71, 72, 73, 77, 80 | 6..8/10 | 0.0420 | various |
| α | 73, 74 | 8/10 | 0.0420 | −0.46, −0.45 |
| low-γ | 11, 24, 25, 26, 27, 28, 30 | 7..9/10 | 0.0420 | various |

**About 35 cells** fall under p < 0.05 across all six bands × continuous k ∈ {2..80}. The §5.4 prose claim that "low-γ k=13 (p=0.014) is the **only** Grassmann (band, k) cell at uncorrected p < 0.05" is **incorrect**. The most significant single cell tie is **low-γ k=13 and k=14, both at p=0.0137**.

If the prose intended the **restricted scan k ∈ {13, 20, 30}** that the cached `grassmann_k_sweep.csv` reports:

| band | k=13 | k=20 | k=30 |
|---|---:|---:|---:|
| δ        | 0.461 | 0.278 | 0.461 |
| θ        | 0.862 | 0.839 | 0.903 |
| α        | 0.053 | 0.188 | 0.116 |
| β        | 0.423 | 0.097 | 0.138 |
| **low-γ** | **0.014** | 0.053 | **0.042** |
| high-γ   | 0.313 | 0.313 | 0.615 |

Even under this 18-cell restricted scan, **two cells** fall below 0.05: low-γ k=13 (p=0.014) AND low-γ k=30 (p=0.042). The "only one cell" claim is false under either scan.

**Recommended fix:** state that low-γ k=13 (p=0.014, n_trace=8/10) and low-γ k=14 (p=0.014, n_trace=9/10) are the strongest single Grassmann cells; many additional β cells at k ≥ 40 and low-γ cells at k ∈ 10..28 fall under uncorrected p=0.05 but do not survive joint BH (covered in §5.3 verification).

### C.1, C.2 — diagonal-band signature (cohort-median Δθ_i(k))

For each band, top-mode (i = k) cohort-median Δθ in radians, and width of consecutive same-sign-as-top modes back from i = k at three reference k:

| band | top-mode median Δθ (rad) | top-mode median n_trace at k ≥ 12 | k cells at n_trace ≥ 7 (top mode) | width back at k=20 | k=40 | k=60 |
|---|---:|---:|---:|---:|---:|---:|
| δ        | +0.0003 | 6/10 | 35% | 1 | 4 | 51 |
| θ        | −0.0002 | 4/10 | 0% | 1 | 1 | 16 (sign flips) |
| α        | +0.0002 | 7/10 | 81% | **19** | 2 | **59** |
| β        | +0.0007 | 7/10 | **81%** | **13** | **26** | **49** |
| low-γ    | +0.0005 | 7/10 | **67%** | **20** | **39** | 37 |
| high-γ   | +0.0004 | 5/10 | 9% | 13 | 2 (sign flip) | 1 (sign flip) |

**Sustained positive diagonal band** (top mode systematically positive):
- **α**: starting k ~ 10, sustained through k=79; modal width back is variable (19 modes at k=20, 59 modes at k=60).
- **β**: starting k ~ 14, sustained through k=80; **most coherent diagonal in the cohort** (81% of k ≥ 12 cells at n_trace ≥ 7, deep widths back).
- **low-γ**: starting k ~ 10, sustained through k=80; broadest diagonal in absolute value (Δθ medians up to +0.0055 rad in the near-top window).

**Sustained negative diagonal:** none. θ's top mode is slightly negative at k ≥ 12 (median 4/10 trace = anti-trace) but width is short.

**Mixed / sign flips:**
- δ: weak positive Δθ; only 35% of k cells at top-mode n_trace ≥ 7.
- high-γ: sign flips along the diagonal across k (top-mode flip occurs around k ~ 30..40).
- θ: weak anti-trace tendency, no sustained band.

### Verification of §5.4 specific qualitative claims

| prose claim | verdict |
|---|---|
| α: n_trace ~ 7–8 throughout starting at k ~ 10 | ✓ (top-mode n_trace median 7/10 at k ≥ 12, 81% of cells at n_trace ≥ 7) |
| β: fine-k cold zone k ≲ 12 with n_trace ≤ 5 | ✓ (top-mode n_trace median ≤ 5 at k ≤ 13) |
| β: positive from k ~ 20 onward at n_trace ~ 7–8 sustained through k=79 | ✓ (median 7/10, fraction 0.81 at n_trace ≥ 7) |
| low-γ: broadest positive diagonal, 3–5 modes back from i=k | ✓ structurally; widths back are 20/39/37 at k=20/40/60 — the "3–5 modes back" framing under-counts. Actually nearer 20–40 modes back of same-sign-as-top. Recommend: "broad positive diagonal extending many modes back from i=k". |
| low-γ: n_trace ~ 7–8 across k ~ 10–79 | ✓ (top-mode median 7/10 at k ≥ 12) |
| θ: n_trace ~ 4 sustained for k ≳ 12 (anti-trace) | ✓ |
| δ: n_trace ~ 5–6 at coarse k without contiguous span | ✓ (top-mode median 6/10, 35% of cells at n_trace ≥ 7) |
| high-γ: sign-flips along the diagonal across k | ✓ |

---

## D. Singleton mask robustness

### D.1 — cohort N (channels)

| patient | N | | patient | N |
|---|---:|---|---|---:|
| Pat_02 | 117 | | Pat_10 | **113** |
| Pat_03 | **122** | | Pat_13 | 119 |
| Pat_05 | 118 | | Pat_14 | 119 |
| Pat_06 | 115 | | Pat_15 | 118 |
| Pat_07 | 116 | | | |
| Pat_08 | 120 | | | |

**min = 113 (Pat_10), median = 118, max = 122 (Pat_03)**. Pat_10 is the cohort-N floor and sets the singleton ceiling: at large k the partition collapses to mostly singletons earliest in Pat_10. The §5.4 "k ≳ 80" framing for the mask edge is consistent (mask 0.5 sets k_max_unmasked ∈ 79..82 across bands).

### D.2 — sensitivity to mask threshold (0.4, 0.5, 0.6)

| band | thr | k_max_unmasked | longest contig span n_trace ≥ 7 | n cells at n_trace ≥ 7 |
|---|---:|---:|---|---:|
| α | 0.4 | 72 | 4 (k=18..21 / k=34..37) | 15 |
| α | 0.5 | 81 | 4 | 15 |
| α | 0.6 | 89 | 4 | 15 |
| β | 0.4 | 69 | **16** (k=17..32) | 22 |
| β | 0.5 | 79 | **16** (k=17..32) | 29 |
| β | 0.6 | 88 | **16** (k=17..32) | 32 |
| low-γ | 0.4 | 70 | **12** (k=45..56) | 27 |
| low-γ | 0.5 | 79 | **12** (k=45..56) | 33 |
| low-γ | 0.6 | 88 | **12** (k=45..56) | 39 |

**Verdict per band:**
- α: longest contiguous span of 4 cells unchanged across thresholds; total n_trace ≥ 7 cells unchanged at 15.
- **β: 16-cell contiguous span k=17..32 unchanged across all three thresholds.**
- **low-γ: 12-cell contiguous span k=45..56 unchanged across all three thresholds.**

The contiguous-pro-trace spans reported in section B for the three trace bands **survive threshold 0.4 unchanged**. The β fine-k pile-up and low-γ concentrated band are robust to mask choice.

---

## E. Pat_03 / Pat_07 / Pat_15 dropout sensitivity at peak cells

| (band, k) | full cohort | pro-trace patients | anti-trace | drop Pat_03 | drop Pat_07 | drop Pat_15 |
|---|---:|---|---|---:|---:|---:|
| **α k=3 (peak)**     | 9/10 | Pat_02, 03, 05, 06, 07, 08, 10, 14, 15 | Pat_13 | **8/9** | **8/9** | **8/9** |
| **β k=22 (peak)**    | 9/10 | Pat_02, 03, 05, 06, 07, 08, 10, 13, 14 | Pat_15 | **8/9** | **8/9** | **9/9** |
| **low-γ k=48 (peak)**| 8/10 | Pat_02, 03, 05, 06, 08, 10, 13, 14 | Pat_07, Pat_15 | 7/9 | 8/9 | 8/9 |
| low-γ k=50 | 8/10 | Pat_02, 03, 05, 06, 08, 10, 13, 14 | Pat_07, Pat_15 | 7/9 | 8/9 | 8/9 |

(The §5.4 verification asked about low-γ k=50; the actual peak per `cohort_band_summary.csv` is k=48. Same pro/anti patient sets at both k.)

**No collapse under any of the three single-patient drops:**
- α k=3: every drop preserves 8/9 ≥ 7/9 majority.
- **β k=22: dropping Pat_15 promotes the cohort to unanimous 9/9.** Pat_15 is the lone β anti-trace patient at this k.
- low-γ k=48: dropping Pat_03 yields the lowest sub-cohort count (7/9) but still a clear majority. Dropping Pat_07 or Pat_15 (the anti-trace patients) preserves 8/9 trace.

**Pat_03 dropout does not collapse the VI(k) signal at any of the three peak cells.** The α and β peaks are in fact insensitive to Pat_03 (Pat_03 is a strong pro-trace patient there: T_VI = −0.168 at α k=3 and T_VI = −0.401 at β k=22).

---

## Sufficient numbers for §5.4 prose corrections

- **Grassmann formula:** code uses `d_G = sqrt(sum sin^2(theta_i))` (option c). Replace the Frobenius `||P − P'||_F` line; or scale it by `1/√2` to make the two equations consistent.
- **β VI(k) span:** "k=17..32 contiguous (with k=34..35 returning above 7 after a single-cell dip at k=33)" — not "k=17..35".
- **low-γ VI(k):** the prose's "k=48..54" sits inside a longer **k=45..56** contiguous run with peak at k=48. Optionally widen to k=45..56.
- **θ VI(k):** "two short clusters k=34..36 and k=39..41 with peak 8/10 at k=44" — not "isolated peaks without contiguity".
- **Grassmann uncorrected-p claim:** false. low-γ k=13 ties low-γ k=14 at p=0.014 (the cohort-strongest cells). Many additional β (k ≥ 40) and low-γ (k ∈ 10..28) cells fall under p < 0.05 — about 35 cells over the full {2..80} k grid, 2 cells over the restricted {13, 20, 30} grid. Replace "the only cell" with "the strongest cell" or "tied with low-γ k=14 (p=0.014, 9/10)".
- **low-γ Grassmann diagonal width:** "many modes back" (20–40+ at the cohort-median level) is more accurate than "3–5 modes back".
- **Mask robustness:** β k=17..32 and low-γ k=45..56 contiguous spans survive thresholds 0.4 and 0.6 unchanged.
- **Pat_03 dropout:** α k=3, β k=22, low-γ k=48 all preserve majority trace under Pat_03 drop. Pat_03 is a strong pro-trace patient at the α and β peaks; dropping it lowers the ratio by 1/9 only because the cohort shrinks.
- **Pat_15 dropout at β k=22:** promotes cohort to 9/9 (Pat_15 is the lone β anti-trace patient).
