---
name: talk-slide-M4-the-null
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: M-4
part: II — Methods
duration: ~60 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# M-4 · How we test — and the null that matters

## Slide placeholder (copy into Canva)

**Main point.** Phase-to-phase similarity = **cophenetic correlation** (Spearman
on the tree distances). A trace must clear a **matched-strength surrogate null** —
the mandatory control that killed our earlier headline. The gate is the
**Wilcoxon vs that null**, not patient counts and not the ×null ratio.

**Concepts to land.**
- **Similarity of two phases = cophenetic correlation ρ_coph** — Spearman on the
  cophenetic (tree) distances. We use **ρ_sym**, a symmetric split-half estimator
  that fixes sign-flips on near-zero patients.
- **A trace must clear a matched-strength surrogate null** — surrogates that
  **preserve each node's strength** but destroy the specific structure. This is
  the **mandatory** control: an earlier within-baseline headline *collapsed* when
  matched-strength reproduced ~half the effect.
- **A second, complementary control — the drift null.** Matched-strength asks
  *"more than a strength change?"*; the **drift null** asks *"more than slow
  session non-stationarity?"* — the brain state drifting over time, unrelated to
  the task. Running both closes two different escape routes. (Matched-strength is
  the **verification gate**; the drift and within-baseline splits are
  **diagnostics**.)
- **The gate is the Wilcoxon vs matched-strength** — **not** patient counts,
  **not** the ×null ratio. Cross-phase drift alone ≠ trace.
- **How to read ρ_sym — two things make 0.20 big.** ρ_sym is a *correlation
  between the task-induced reorganization and the change that persists into rest*.
  **(1) The scale is set by the null, not by 1:** against the matched-strength
  null (~0.01), **β = +0.20 is 16.6× the null** (α = +0.10, 13.7×). **(2) It's
  diluted by the non-movers:** most of the structure doesn't move across phases
  (stable / anchored), so the global number averages the real effect over a stable
  majority — **condition on the pairs that *did* move, and the trace rises well
  above 0.20** ("of those that moved, this fraction traced"). Per-patient tracers
  reach **0.3–0.5**. Read ρ *with* its ×null ratio — don't read 0.20 as "only 20%."

**Figures / visuals.**
- Signal-vs-null **triangle**, per band:
  `data/preprint/figures/all_bands/fig_bands_null_triangle_coph.pdf`
  (optionally beside the raw-FC version for the raw-vs-multiscale preview:
  `…/fig_bands_null_triangle_rawfc.pdf`).

**References.** None external — this is our null methodology (matched-strength
surrogates; audit_63 is the internal record).

---

## Keep honest (content constraints, not styling)

- **Matched-strength is MANDATORY** before any FC cohort claim; it is the
  **verification gate**. The **drift** null and within-baseline / split-half nulls
  are **complementary diagnostics** — they rule out non-stationarity (passage of
  time), but they are **not verification** on their own (audit_63 — a
  within-baseline headline collapsed under matched-strength; that's the rule's
  origin).
- **Never rank bands by the ×null ratio.** γ_low's ratio is *larger* than β's yet
  it fails the gate — the matched-strength control certifies a band, not the ratio.
- Cross-phase similarity = **ρ_coph (Spearman on cophenetic distances)** — never
  ARI / NMI / fcluster(K).
- If numbers appear, they are **ρ_sym**: **β p = 0.032 (16.6× null)** ·
  **α p = 0.024 (13.7×)**. Verdicts are estimator-invariant (0/6 bands flip).
- **Per-patient counts are never the gate** (the Wilcoxon is); pair the cohort p
  with LOO. Trace = T_d > 0; surrogate upper-tail p = mean(s ≥ obs).
- **Distinguish "fraction moved" from "fraction of movers that traced."** Claiming
  *most pairs moved* ("90% movers") is a **global-shift confound** — don't. But
  *of the pairs that moved, the fraction that traced* is the legitimate
  conditional that lifts the number above the global 0.20 — **provided it clears
  the matched-strength null** (which rules out the global-shift explanation for
  that concentrated effect). Report ρ **with** its ×null ratio and the
  %held / %struct decomposition; ×null alone ≠ trace.
