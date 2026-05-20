---
date: 2026-05-08
era: COHORT_N10 / IMCOH_ABS
status: current
type: verification-report
scope: §5.2 KC cross-probe restriction at β (λ∈{0,1}) and low-γ (λ=1); α negative control
inputs:
  - data/audit/section5_v2_kc_controls/per_patient_table.csv
  - data/audit/section5_v2_kc_controls/cohort_summary.csv
  - scripts/01_compute/audit/audit_46_kc_section5_controls.py
---

# KC cross-probe restriction — verification

**Head.** All §5.2 cohort numbers reproduce exactly from the cached audit. The "low-γ behaves opposite to β" framing is a **prose-level apples-to-oranges artefact**, not a mechanism mismatch. §5.2 silently compares raw p (full) to BH-corrected q (cross-probe) for the β cells but raw-to-raw for low-γ. On a clean apples-to-apples comparison — either raw vs raw or BH vs BH within m=12 — the three cells are mutually consistent:

- under **BH(m=12)** on the cross-probe layer, all three cells *strengthen-or-unchanged* relative to BH(m=12) on the full layer (β λ=0: 0.168 → 0.146; β λ=1: 0.168 → 0.168; low-γ λ=1: 0.168 → 0.146).
- under **raw p**, β λ=0 weakens slightly (0.0186 → 0.0244, count unchanged at 7/10), β λ=1 is exactly unchanged (0.0420 → 0.0420, 8/10), and low-γ λ=1 strengthens via a single sign flip (0.0420 → 0.0244, 8/10 → 9/10) driven entirely by Pat_07 crossing zero (T_KC: +0.087 → −0.084, both within noise of zero).

So the "same-probe = real signal at β, removing it = remove signal" reading is intact at β. At low-γ heights the cohort is dominated by Pat_05 (T = −13.94) and Pat_02 (T = −9.71), with everyone except Pat_15 already piled near zero; the cross-probe boost is a **rank-test sensitivity artefact at the zero-crossing**, not evidence that same-probe pairs at low-γ contribute opposite-sign signal. Hypothesis **A** (same-probe pairs at low-γ are noisy/inconsistent and removing them denoises) is mildly supported by `std(T_xpr − T_full) = 0.268` at low-γ λ=1 vs `0.135–0.195` at the β cells, but the dominant mechanism is the proximity of Pat_07's full-T to zero rather than a systematic pattern.

α λ=0 and λ=1 are silent under both restrictions, as expected.

## Reproduction of §5.2 cohort numbers

Source: `cohort_summary.csv` matches `wilcoxon(... alternative='less')` on per-patient rows from `per_patient_table.csv` to all printed digits.

| cell | n_neg full | p_full | n_neg xpr | p_xpr | q_BH full (m=12) | q_BH xpr (m=12) |
|---|:-:|---:|:-:|---:|---:|---:|
| β λ=0   | 7/10 | 0.018555 | 7/10 | **0.024414** | 0.168 | **0.146** |
| β λ=1   | 8/10 | 0.041992 | 8/10 | **0.041992** | 0.168 | **0.168** |
| low-γ λ=1 | 8/10 | 0.041992 | 9/10 | **0.024414** | 0.168 | **0.146** |
| α λ=0   | 6/10 | 0.460938 | 6/10 | 0.500000 | 0.857 | 0.886 |
| α λ=1   | 6/10 | 0.500000 | 6/10 | 0.539062 | 0.857 | 0.886 |

The BH(m=12) layer was reproduced independently by ranking all 6 bands × {λ=0, λ=1} cross-probe raw p-values and applying the standard Benjamini–Hochberg step-up.

## Per-patient T_KC: full vs cross-probe

### β λ=0 (topology)

| patient | T_full | T_xpr | diff (xpr − full) | sgn_full | sgn_xpr | flip? |
|---|---:|---:|---:|:-:|:-:|:-:|
| Pat_02 |  +0.5081 |  +0.7917 | +0.2836 | + | + | |
| Pat_03 | −19.7298 | −19.1387 | +0.5911 | − | − | |
| Pat_05 |  +0.1749 |  +0.3089 | +0.1339 | + | + | |
| Pat_06 |  −5.0902 |  −4.8644 | +0.2258 | − | − | |
| Pat_07 |  −2.4326 |  −2.2320 | +0.2005 | − | − | |
| Pat_08 |  −5.0689 |  −4.6380 | +0.4309 | − | − | |
| Pat_10 |  −2.4519 |  −2.3165 | +0.1353 | − | − | |
| Pat_13 |  −3.1590 |  −3.0936 | +0.0654 | − | − | |
| Pat_14 |  −1.3667 |  −1.2841 | +0.0826 | − | − | |
| Pat_15 |  +2.3806 |  +2.2832 | −0.0974 | + | + | |

mean(xpr−full) = +0.205, std = 0.195. Same-probe pairs systematically push T_KC *down* (more imprinted) at β-topology, by ~0.2 per patient. Removing them moves all imprinted patients slightly closer to zero — the imprint-direction count is preserved (7/10) but the rank-test loses a little power, so raw p ticks up 0.0186 → 0.0244. **Consistent with same-probe = real signal at β-topology.**

### β λ=1 (heights)

| patient | T_full | T_xpr | diff | sgn_full | sgn_xpr | flip? |
|---|---:|---:|---:|:-:|:-:|:-:|
| Pat_02 | −7.6450 | −7.3547 | +0.2902 | − | − | |
| Pat_03 | −2.7111 | −2.7676 | −0.0566 | − | − | |
| Pat_05 | −0.3326 | −0.2069 | +0.1258 | − | − | |
| Pat_06 | −2.2197 | −2.0572 | +0.1625 | − | − | |
| Pat_07 | −1.3483 | −1.3737 | −0.0254 | − | − | |
| Pat_08 | −6.2626 | −6.1180 | +0.1446 | − | − | |
| Pat_10 | −1.4171 | −1.5152 | −0.0981 | − | − | |
| Pat_13 | −1.0314 | −1.1345 | −0.1031 | − | − | |
| Pat_14 | +1.2500 | +1.2900 | +0.0400 | + | + | |
| Pat_15 | +2.4069 | +2.3225 | −0.0845 | + | + | |

mean(xpr−full) = +0.040, std = 0.135. Mixed signs — same-probe contribution is small and inconsistent across patients. No sign flips, both raw p and BH-q are unchanged.

### low-γ λ=1 (heights)

| patient | T_full | T_xpr | diff | sgn_full | sgn_xpr | flip? |
|---|---:|---:|---:|:-:|:-:|:-:|
| Pat_02 |  −9.7130 |  −9.1384 | +0.5746 | − | − | |
| Pat_03 |  −2.1616 |  −2.0763 | +0.0852 | − | − | |
| Pat_05 | −13.9394 | −13.6859 | +0.2535 | − | − | |
| Pat_06 |  −6.4684 |  −6.5300 | −0.0616 | − | − | |
| Pat_07 |  **+0.0873** |  **−0.0840** | −0.1714 | + | − | **FLIP** |
| Pat_08 |  −0.3805 |  −0.3208 | +0.0597 | − | − | |
| Pat_10 |  −0.2588 |  −0.3153 | −0.0565 | − | − | |
| Pat_13 |  −0.3798 |  −0.4691 | −0.0893 | − | − | |
| Pat_14 |  −0.0410 |  −0.0554 | −0.0144 | − | − | |
| Pat_15 |  +9.1571 |  +8.7108 | −0.4463 | + | + | |

mean(xpr−full) = +0.013, std = **0.268** (largest of the three cells). Only **Pat_07 flips** sign — and it flips from +0.087 to −0.084, both within noise of zero. The "9/10 imprint" headline is built on a sub-0.1 magnitude crossing. Pat_15 (+9.16 → +8.71) stays robustly anti-imprint; Pat_05 and Pat_02 stay deeply imprint with magnitudes >9.

### α λ=0 / λ=1 (silent control)

α λ=0: 6/10 imprint full, 6/10 imprint xpr; raw p 0.461 → 0.500. Pat_03 (+1.99) and Pat_07 (+6.28) sit positive; Pat_06 (−7.16), Pat_08 (−6.19) sit negative. Both layers null.

α λ=1: 6/10 imprint full, 6/10 imprint xpr; raw p 0.500 → 0.539. Pat_03 (+7.15) and Pat_15 (+5.67) sit positive; Pat_06 (−3.60), Pat_08 (−3.84) sit negative. Both layers null.

**α confirms:** when the underlying cohort is null, cross-probe restriction does not invent a signal.

## Variance diagnosis (Hypothesis A vs B)

`std(T_xpr − T_full)` across the cohort, by cell:

| cell | std(xpr − full) | mean \|T_full\| | mean \|T_xpr\| |
|---|---:|---:|---:|
| β λ=0   | 0.195 | 4.236 | 4.095 |
| β λ=1   | 0.135 | 2.663 | 2.614 |
| low-γ λ=1 | **0.268** | 4.259 | 4.139 |

low-γ heights has the **most variable** same-probe contribution across patients — consistent with hypothesis A (same-probe pairs there are noisy/inconsistent rather than systematically signal-bearing). But:

1. The **size** of that variability (0.27) is small in absolute terms compared to the cohort's |T| range, which spans 0.04 to 13.9 at low-γ λ=1.
2. The "denoising" only matters for **Pat_07**, whose full T is so close to zero (+0.087) that any same-probe contribution can flip its sign.
3. Pat_15, the dominant counter-imprint patient at low-γ heights, *also* moves down by 0.45 under restriction — still strongly positive. The cohort headline remains "Pat_05/Pat_02 vs Pat_15".

So hypothesis A is **mildly supported**. The full mechanistic statement is: **at low-γ heights the cohort is bimodal (two strongly imprinted patients vs one strongly counter-imprinted), with everyone else piled within ±0.4 of zero. The cross-probe restriction does not change the bimodality but happens to nudge Pat_07's near-zero T from +0.087 to −0.084. That single sign flip is the entire 8/10 → 9/10 / p 0.042 → 0.024 effect.** The strengthening is real but operates on a single patient sitting at the noise floor, not on a structural property of same-probe vs cross-probe distinguishability at low-γ.

Hypothesis **B** (different probe convention) is **rejected**: the same `cross_probe_pair_mask(probe_ids)` from `audit_46_kc_section5_controls.py` is used at every band/axis cell.

## Verdict

The §5.2 numbers are correct. The "INCONSISTENT" framing for low-γ vs β is a comparison artefact: §5.2 mixes raw p and BH-q across cells. On a unified comparison (raw vs raw, or q vs q within m=12), the three cells are mutually consistent — all three cross-probe layers either match or slightly improve on their full counterparts after BH correction.

The mechanistic asymmetry §5.2 is trying to convey is real but more delicate than "low-γ contradicts same-probe = real signal":

- **At β-topology**, same-probe pairs add ~0.2 of imprint signal per patient on average; removing them reduces signal coherently across the cohort, costing the rank test some power.
- **At low-γ-heights**, same-probe pairs add zero on average but with the largest patient-to-patient variability of the three cells; removing them flips one patient (Pat_07) at the zero-crossing, which is enough to gain the rank test one positive vote.

Both readings are compatible with "same-probe pairs carry real but small lagged-coupling signal under |ImCoh|"; what differs is the cohort's *consistency* of that contribution — coherent at β, scattered at low-γ.

## Sufficient numbers for §5.2 prose

- β λ=0: full 7/10, p=0.019; xpr 7/10, p=0.024; q_BH 0.168 → 0.146.
- β λ=1: full 8/10, p=0.042; xpr 8/10, p=0.042; q_BH 0.168 → 0.168.
- low-γ λ=1: full 8/10, p=0.042; xpr 9/10, p=0.024; q_BH 0.168 → 0.146; the count change is **Pat_07 only** (+0.087 → −0.084, both in noise floor).
- α λ=0 / λ=1: 6/10 / 6/10 in both layers, all raw p ≥ 0.46. Negative control intact.
- std(T_xpr − T_full): β λ=0 = 0.195; β λ=1 = 0.135; low-γ λ=1 = 0.268 (largest of the three).
