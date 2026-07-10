---
name: talk-slide-R8-epi-coda
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: R-8
part: III — Results
duration: ~75 s (optional / cuttable)
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - epi_allcontacts_marker_2026_06_23
  - soz_diffusion_community_layout_2026_07_08
---

# R-8 · Coda (N3) — the same operator reads epileptogenic tissue

## Slide placeholder (copy into Canva)

**Main point.** The *identical* diffusion view we built for cognition doubles as a
**clinical marker**: seizure-onset contacts form a **strength-independent diffusion
community**, and the cognitive bands **split over it** — **β spares** the epileptogenic
core, **α recruits** it. Same operator, second payoff — a cognitive → clinical bridge,
for free.

**Concepts to land.**
- **Seizure tissue is a diffusion community.** Under the same heat-kernel propagator,
  seizure-onset (SOZ) contacts co-diffuse. Matched-strength-gated **co-diffusion AUC**:
  **δ 0.80, γ_low 0.74, β 0.69**; off-shaft **δ LOSO 0.72, 8/10** — a
  **strength-independent** marker (not just "where connectivity is strong").
- **A calibrated detector.** A seed-based 6-band detector reaches **mean AUC 0.81,
  precision@5 ≈ 60 % (~7× chance), 9/10** patients above chance.
- **The cognitive bands split over the clinical core.** **β spares** the SOZ (the
  cognitive trace lives *off* the epileptogenic tissue), while **α recruits** it — the
  cognitive → clinical dissociation on the *same* map.

**Figures / visuals.**
- **The bridge — β spares / α recruits the SOZ** (the cognitive → clinical panel):
  `data/preprint/figures/results_section1/fig_trace_f_soz_divergence.pdf`
- **The SOZ as a diffusion community** (force-directed heat-kernel layout):
  `data/preprint/figures/results_section1/fig_soz_diffusion_community.pdf`
- **The marker + the detector:**
  `data/reports/results_section3/fig_epi_a_relational_marker.pdf` ·
  `data/reports/results_section3/fig_epi_b_calibrated_detector.pdf`
  (compound alternatives:
  `data/preprint/figures/all_bands/fig_epi_compound_probability.pdf`,
  `…/fig_epi_soz_marker.pdf`).

**References.** None external — this is our own diffusion marker (audit_132 all-contacts
propagator; audit_89 diffusion community). LRG operator = Villegas 2023/2025 (M-3).

---

## Keep honest (content constraints, not styling)

- **Scope = clinical SOZ labels, NOT surgical outcome.** State it once: this predicts
  *clinically labelled* seizure-onset contacts, not resection success. Don't imply a
  surgical claim.
- **The clinical band is δ, the cognitive band is β — same OPERATOR, different band.**
  The strongest epi marker is **δ (AUC 0.80)**; the cognitive trace is **β**. Don't
  conflate the δ clinical read-out with the β cognitive trace — the point is *one
  operator, two independent read-outs*, and β **spares** the SOZ.
- **Strength-independent** — every AUC quoted is **matched-strength-gated**; the marker
  is not a restatement of "seizure tissue is highly connected."
- **AUC ~0.8 is a marker, not a diagnostic** — a strength-independent enrichment, not
  a clinical-grade classifier. Keep the register honest.
- **This slide is CUTTABLE.** The talk stands on N1 + N2 alone; N3 is the "same tool,
  second use" bonus. Cut it first if you run long.
