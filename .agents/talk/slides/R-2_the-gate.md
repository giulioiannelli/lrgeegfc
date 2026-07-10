---
name: talk-slide-R2-the-gate
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: R-2
part: III — Results
duration: ~75 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# R-2 · The gate — which bands hold the trace

## Slide placeholder (copy into Canva)

**Main point.** At the cohort level, **β and α hold the trace**; γ_low / δ / γ_high
don't reach the cohort gate; **θ is the clean cohort-negative**. But **no band is
trace-free** — every band carries per-patient variance (γ_low even has the
*largest single-patient* trace) — so the gate certifies **cohort-level hold via
the matched-strength null**, not presence-vs-absence. Peak ≠ verdict.

**Concepts to land.**
- **The gate is cohort-level hold** (Wilcoxon vs matched-strength) — **not "has a
  trace / doesn't."** *Every band carries some per-patient variance; no band is
  trace-free.* The question is whether the hold is **cohort-wide**.
- Cohort verdicts (ρ_sym): **β +0.20, 6/10, p = 0.032 (16.6×)** and
  **α +0.10, 6/10, p = 0.024 (13.7×)** hold; **γ_low +0.08, p = 0.080** and
  δ / γ_high don't reach the gate; **θ −0.04, p = 0.784** is the clean
  cohort-negative.
- **Peak ≠ verdict — the γ_low lesson:** γ_low carries the *largest single-patient*
  trace (**Pat_05 +0.86**, vs β's max +0.54) yet **doesn't clear the cohort gate**.
  Strong individual traces ≠ cohort-level hold — certify by the **null**, not the
  peak. (Also why we never rank bands by the ×null ratio — γ_low's is *bigger*
  than β's.)

**Figures / visuals.**
- **ρ_sym "gate bloom"** (radial — the talk-native view; *lead figure*). Six band
  wedges, β at 12 o'clock; one **petal per patient**, grown from a T=0 **waterline**
  (outward = trace held, inward = reset). **Wedge glow ∝ cohort significance**
  (−log₁₀ gate p), so **only β and α bloom**; the rest hug the ring. One view carries
  the whole slide: the **petals are the per-patient spread** (every band has
  variance — no band is trace-free); the **glow is the cohort gate** (β/α hold);
  and the **longest petals in the figure — γ_low Pat_05 +0.86 and γ_high Pat_05
  +0.92 (rim chevron) — sit in *non-glowing* wedges**, which is *peak ≠ verdict*
  made visual. Path:
  `data/preprint/figures/_drafts/fig_reasoning_bloom_gate_DRAFT.pdf`
  (builder: `scripts/01_compute/figures_embedded/fig_reasoning_bloom.py --func rho_sym_gate`).
- **Per-band per-patient forest** (rigorous alternative / backup — linear axis, each
  dot own-null coloured): `data/preprint/figures/results_section1/fig_trace_a_band_forest.pdf`
  (alt / summary: `data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf`).

**References.** None — this is our result.

---

## Keep honest (content constraints, not styling)

- **No band is trace-free — there is always some per-patient variance.** The gate
  is *cohort-level hold* (Wilcoxon vs matched-strength), **not presence/absence.**
  Never say "θ has no trace" — say "θ **doesn't hold at cohort level** / is the
  clean cohort-negative."
- **Peak ≠ verdict:** γ_low has the largest single-patient trace yet fails the
  cohort gate; **never rank bands by the ×null ratio** (M-4).
- **Show the spread** (forest / violin), not one number per band — the per-patient
  variance *is* part of the result.
- ρ_sym numbers; the **6/10 counts are context, never the gate** (the Wilcoxon is);
  pair with LOO.
- α **holds the cohort gate** but is secondary to β — don't overclaim α here (its
  weaker standing / lack of an anatomical home is the band-taxonomy story, later).
- **In the gate bloom the radial value IS ρ_sym itself** (the main cross-phase
  gate, audit_150) — *not* one of the consolidation-arc trace functionals (those
  are a different, finer analysis, R-6/R-7). What certifies a band is the **glow
  (cohort gate p vs matched-strength)**, never the petal length — say this if anyone
  reads the big γ spikes as the result. Bloom uses the turbo band palette rather
  than the canonical `band_color` (aesthetic, presenter's call).
