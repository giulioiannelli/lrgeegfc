---
name: talk-slide-20-epilepsy-marker
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 20
status: draft
updated: 2026-07-13
canva: page 19 · ~30% (cuttable; ⚠ 2026-07-13 added τ role = ranking-not-precision + mst@0.20 robustness)
---

# Slide 20 — Coda: the same operator reads epileptogenic tissue

1. TITLE
Coda — the same operator reads epileptogenic tissue

2. MAIN CONCEPT
- The identical diffusion view, no retuning: seizure-onset (SOZ) contacts form a strength-independent diffusion community.
- Matched-strength-gated co-diffusion AUC: δ 0.80, γ_low 0.74, β 0.69; off-shaft δ LOSO 0.72, 8/10 — strength-independent, not "where connectivity is strong".
- The cognitive bands split over the clinical core: β spares the SOZ (the cognitive trace lives off the epileptogenic tissue), α recruits it.
- One operator, two independent read-outs — a cognitive → clinical bridge for free.
- τ's role here (honest): sweeping the diffusion time *deepens* the picture — the seizure zone separates better as a co-diffusing community, so multiscale AUC beats single-scale in most patients. But the top-of-list precision (~60% @5) comes from *fusing the bands*, not from τ. So τ buys **ranking, not precision** — a marker, not a τ-boosted detector. And the marker reproduces on the sparsified `mst@0.20` backbone → it's a property of the operator, robust to how the graph is thinned.

3. ON-SLIDE TEXT
SOZ = a strength-independent diffusion community
co-diffusion AUC: δ 0.80 · γ_low 0.74 · β 0.69
β spares the SOZ · α recruits it
(SOZ labels, not surgical outcome)

4. SPEECH
And the same tool, with nothing changed, pays a second dividend. The seizure-onset contacts form their own diffusion community — strength-independent, so it isn't just "where connectivity is strong". The clinical band is δ, area under the curve 0.80; the cognitive band is β. And on the same map they dissociate: β spares the epileptogenic core, α recruits it. One operator, two independent read-outs — the cognitive lens doubles as a clinical one.

Careful: scope is clinically-labelled SOZ, not surgical outcome. δ is the clinical band, β the cognitive one — same operator, different band; β spares the SOZ. AUC ~0.8 is a marker, not a diagnostic. This slide is CUTTABLE if you run long.

5. FIGURES
- β spares / α recruits the SOZ — data/preprint/figures/results_section1/fig_trace_f_soz_divergence.pdf
- SOZ as a diffusion community — data/preprint/figures/results_section1/fig_soz_diffusion_community.pdf
- Marker + detector — data/reports/results_section3/fig_epi_a_relational_marker.pdf · fig_epi_b_calibrated_detector.pdf

6. REFERENCES
- Ours (audit_132 all-contacts propagator; audit_89 diffusion community).

7. CANVA STATUS
Page 19 · ~30% ("Coda N3"). Have: the concept + limitations text. Missing: the SOZ figures, compressed on-slide text, presenter notes. Cuttable.
