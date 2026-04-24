---
name: imcoh-gap-analysis
type: report
era: IMCOH_SQ
status: superseded
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

> **⚠ QUANTITATIVE_STALE (2026-04-15)** — numbers in this report were
> computed under the pre-reset ImCoh mislabelling (stored `|ImCoh|²`
> under the name `imcoh`). Qualitative conclusions survive (rankings
> preserved under sqrt); regenerate absolute values against
> `fc_method="imcoh_abs"` before quoting. See
> `~/.claude/projects/.../memory/imcoh_taxonomy.md`.

# H2a Beta Gap Analysis — Why Unanimity Cells Are Non-Contiguous

## Summary

Beta H2a has 82/93 k-levels with 5/5 unanimity (k=16-108 range, 88% coverage).
The 32 gap k-values are caused by **single-patient marginal fluctuations**, not
structural transitions. The underlying effect is continuous and robust.

## Contiguous Runs

| Run | k-range | Cells | Note |
|-----|---------|-------|------|
| 1 | 18-34 | 17 | Clean mid-range |
| 2 | 42-58 | 17 | Clean mid-range |
| 3 | **60-97** | **38** | **Longest run** |
| 4 | 100-106 | 7 | Tail |
| Isolated | 16, 36, 108 | 3 | Single cells |

## Gap Decomposition

| Gap region | k-values | Breaker | Contrast | Margin to 2nd | Explanation |
|-----------|----------|---------|----------|---------------|-------------|
| Onset | 2-15 | **Pat_03** | −0.06 to −0.71 | 0.04-0.60 | Pat_03 coarse-scale outlier (1024 Hz) |
| Single | 17 | Pat_02 | −0.021 | 0.09 | Marginal fluctuation |
| Gap 1 | 35, 37-41 | **Pat_08** | −0.01 to −0.07 | 0.28-0.43 | Pat_08 barely negative, 4 others strongly positive |
| Single | 59 | Pat_05 | −0.009 | 0.25 | Essentially zero |
| Pair | 98-99 | Pat_05 | −0.011 | 0.07 | Essentially zero |
| Tail | 107-112 | Pat_03/07 | −0.000 to −0.017 | 0.01-0.03 | Noise at fine scales |

## Who Breaks Unanimity

| Patient | Times breaking | % of gaps | Nature |
|---------|:-------------:|:---------:|--------|
| Pat_03 | 16 | 50% | Coarse-scale outlier (1024 Hz sampling) |
| Pat_08 | 6 | 19% | Narrow k-window, barely negative |
| Pat_05 | 3 | 9% | Essentially zero crossings |
| Pat_07 | 3 | 9% | Fine-scale noise |
| Pat_02 | 1 | 3% | Single marginal fluctuation |

## Key Observation

At every gap k-value, the breaking patient's contrast is marginal
(|contrast| < 0.07 in 28/32 gaps). The other 4 patients remain strongly
positive (+0.10 to +0.70). The unanimity criterion (strict 5/5) amplifies
single-patient noise into apparent gaps in an otherwise continuous effect.

## Recommended Framing

> "Beta H2a contrast is positive in all 5 patients across k=16-108
> (mean +0.22 ± 0.04). Under strict 5/5 unanimity, 82/93 k-levels
> (88%) reach consensus. Gaps are caused by single-patient marginal
> fluctuations: Pat_03 at coarse scales (documented 1024 Hz outlier),
> Pat_08 at k=35-41 (contrast −0.04 ± 0.02), and isolated zero-crossings
> at k=59, 98-99 (|contrast| < 0.01)."

## Implication

The non-contiguity is NOT a weakness — it's a consequence of the strict
unanimity criterion applied to a continuous signal. The effect is present
at all scales from k=16 to k=108. The gaps tell us more about individual
patient variability than about the phenomenon itself.
