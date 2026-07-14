---
name: talk-slide-11b-the-measure
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym · mst@0.20 backbone)
slide: 11b   # S2 in the 10-slot results restructure (.agents/talk/2026-07-13_results-restructure-9slot.md); sits between pipeline (11) and taxonomy (12)
status: draft
updated: 2026-07-13
canva: page 17 (deck) · math column LOCKED to the on-slide LaTeX below; figure = reinstatement colored circular dendrograms. RE-RENDER pending decisions on (i) the extra conceptual image, (ii) whether the per-patient panel b stays here.
---

# Slide 11b (S2) — The Measure

1. TITLE
The Measure

2. MAIN CONCEPT
- **Head (the juice):** the pipeline turned each phase into a *hierarchy*; the measure asks **how similar two hierarchies are**, at a given scale — one rank number, ρ^coph(τ), built from cophenetic distances and read against a **matched-strength null**. Give only the *sense* here; the selectivity, scale-signature and full null argument come next.
- **From tree to number, three steps (the left column is the whole story):**
  - **cophenetic distance** — the dendrogram assigns every pair of contacts D^coph_ij(τ) = the merge height of their lowest common ancestor. A scale-aware fingerprint of the hierarchy, one number per pair.
  - **displacement Δ** — how a phase moved the hierarchy relative to another: Δ^task = D_task − D_ref, Δ^rest = D_post − D_ref.
  - **ρ^coph** — rank-correlate two phases' cophenetic distances (Spearman, ∈ [−1, 1]); for the trace we correlate the two *displacements* against a split rest_pre baseline (halves A, B; a different half per displacement so they share no baseline noise), symmetrized (A↔B).
- **The colour IS the number (right column).** Each contact is shaded by its per-contact cophenetic preservation to the reference phase — the per-leaf term of ρ^coph. Green = it kept its place in the hierarchy; grey = it was reshuffled. The whole-tree bar under each dendrogram is the single ρ^coph number.
- **Two axes for everything that follows** (plant lightly): ρ is read at **every scale τ** (→ scale-signature slide) and it **compares phases** (→ taxonomy slide, next).
- **The baseline is a matched-strength null.** The per-patient panel reads each ρ^coph against a strength-preserving surrogate — name it here (it is the *only* null in the talk); the "why matched-strength / what it destroys" argument is unpacked on the nulls slide.

3. ON-SLIDE TEXT
Left column (bullets + compiled LaTeX, top → bottom) — LOCKED to the Canva deck:

- recap: timeseries → |ImCoh| → diffusion → dendrogram (linkage ℒ) → merge height

  D^{\mathrm{coph}}_{ij}(\tau)=h_{\mathcal{L}[\mathcal{D}(\tau)]}\!\bigl(\mathrm{LCA}(i,j)\bigr)

- ρ = rank-correlation of phases' cophenetic D ∈ [−1, 1]
- Δ = how a phase moved the hierarchy wrt another

  \Delta^{\mathrm{task}}_{X}=D^{\mathrm{coph}}_{\mathrm{task}}-D^{\mathrm{coph}}_{X},\\
  \Delta^{\mathrm{rest}}_{X}=D^{\mathrm{coph}}_{\mathrm{post}}-D^{\mathrm{coph}}_{X}

- compare two phases: rank-correlate their cophenetic distances at scale

  \rho^{\rm coph}(\tau)=\tfrac{1}{2}\Bigl[\rho_S\bigl(\Delta^{\mathrm{task}}_{\mathrm{rsPre}A},\,\Delta^{\mathrm{rest}}_{\mathrm{rsPre}B}\bigr)+\\+(A\leftrightarrow B)\Bigr]

Small caption near panel b: "per-patient ρ^coph vs a matched-strength null (the sole null — unpacked next)".

4. SPEECH
The pipeline turned each phase into a hierarchy. The measure asks one thing: how similar are two of those hierarchies, at a given scale. Three steps, all on the left. First, the tree hands every pair of contacts a cophenetic distance — the height at which the two first merge into one branch. That's a fingerprint of the whole hierarchy, one number per pair. Second, a displacement: how far a phase moved that fingerprint away from a reference — Δ for the task, Δ for the rest afterwards. Third, the measure itself: rank-correlate the two displacements. If they agree, the task's reshaping is still there at rest. We measure both against a split rest-pre baseline, a different half for each side so they share no noise, then swap and average — that's the symmetrized rho-coph. And the colour on the right is that number made visible: every contact is shaded by how well it kept its cophenetic place — green kept it, grey lost it — and the bar under each tree is the single number. Keep two things in mind because they organize the whole rest of the talk: we read this at every scale, and it compares phases. One last thing — the per-patient bars at the bottom are read against a matched-strength null, the one null in this talk; I'll show you exactly what that destroys in a couple of slides.

Careful: give only the SENSE of the measure here — do NOT front-run the selectivity (which bands), the scale-signature (β vs α), or the full matched-strength argument (all later). ρ is SPEARMAN on cophenetic distances (rank), never Pearson, never a partition metric (ARI/NMI/fcluster-at-k). The symmetrization uses the two rest_pre halves A,B with a DIFFERENT half per displacement (bias-fix — not "averaging two runs"). Read PER-SCALE τ, never a single best-τ scalar. Matched-strength is the SOLE null (drift retired) — if panel b stays, NAME the null (its p is on the x-axis) but keep the "what it destroys" cartoon for the nulls slide. The tanglegram + 4 fates belong to the NEXT slide (S3) — do NOT pre-empt.

5. FIGURES
- **MAIN — the multiscale reinstatement across scales (user pick 2026-07-13).** Use
  `data/outputs/figures/talk/slide_measure_reinstatement_scales.png`
  (gen `scripts/01_compute/figures_embedded/new_results_sec1/fig_reinstatement_scales_measure.py`):
  a **3×4 grid** — rows = three scales (s = 1.0, 5.6, 31.9), columns = the four phases
  (rest_pre · task_learn · task_test · rest_post), every tree shaded by per-contact ρ^coph TO
  task_test (grey→green), with a similarity bar + "ρ^coph to test" under each. ONE figure that
  carries BOTH of the measure's two axes at once: **colour = the number on the tree** (compares
  phases — rest_post rides greener / closer to test than rest_pre) AND **read at every scale** (the
  three rows). SINGLE-PATIENT ILLUSTRATIVE — a measure demo, NOT the cohort result; the cohort
  ρ_sym(τ) curves + nulls + survivors are the whole of slide 13. Name the matched-strength null in
  WORDS here; it is shown quantitatively on 13.
  - Optional small companion `fig_measure_rhocoph.png` (reference / real 0.53 / matched-strength
    shuffle 0.23) if you want the null shown *visually* on 11b; otherwise it is named in speech.
- **ADDITIONAL CONCEPTUAL IMAGE — the binned correlation + honest calibration gauge.**
  `data/outputs/figures/talk/fig_measure_binned.png` (gen `scripts/07_figures/talk_fig_measure_binned.py`,
  Pat_08 β, s=5.6). RIGHT: a **smoothed density map** (log-scaled, peak-normalized) of pairs ranked
  by Δ^task vs Δ^rest — the density tilts along the diagonal (ρ^coph=0.46; strongest-moving pairs
  agree most), "the task's reshaping predicts rest's." (Replaces the point-cloud scatter and the
  binned-line curve — both read as noisy; deprecated.) LEFT: a **vertical filled calibration bar**
  cut to the reachable max — 0 at bottom, **reachable-max ≈0.53 at top** (the hierarchy's own
  reproducibility; 1.0 impossible), green fill to **0.48** (≈91% of the bar). The figure itself is
  stripped to essentials (only ρ^coph, axis labels, colorbar, and the bar's 0/0.48/0.53); the
  calibration numbers answering "isn't 0.5 low?" live in the narration — **it is NOT low, and the
  load-bearing reason is the NULL**: **0.48 is ≈8× the matched-strength null median and past its
  95th percentile**, and it sits just below the model-free split-half reliability (0.53) so 1.0 is
  unreachable. The honest
  reason 1.0 is unreachable is the model-free **split-half reliability 0.53** (two halves of the
  SAME rest recording only reproduce the hierarchy at 0.53 → 0.48 sits just below the instrument's
  own reproducibility). The "perfect-persistence ceiling" 0.62 is shown ONLY as a *lighter generous
  upper bound* (it shares the D_task term, so it's inflated — the tell is that 0.62 > 0.53; the true
  trace ceiling is lower). Cohort strip flags Pat_08 (0.48) as **top-of-range** (cohort median 0.20).
  ⚠ VERIFIED 2026-07-13 by a 3-lens adversarial panel (verdict "reframe"): do NOT lead with "≈¾ of
  the ceiling" (spin — ratio of two baseline-biased quantities); lead with the null; the COHORT trace
  claim rests on the null-based per-patient panel (10/10, Wilcoxon p=0.001, 16/16 scales), NOT this
  exemplar's magnitude. See [[rho_sym_magnitude_calibration_2026_07_13]].
- **OPTIONAL grounding glyph — the merge-height inset.** A tiny rectangular dendrogram with leaves i, j, their LCA and the height h marked → makes D^coph_ij = h(LCA) tangible next to the first formula. Build in Canva or as a small vector; low priority.
- **DECISION FLAGS:** (i) extra conceptual image = the Δ-scatter (build below); (ii) per-patient
  panel b — **RESOLVED 2026-07-13: MOVED to slide 13** (the lasting-trace punchline). 11b is now a
  pure sense-of-measure slide via `fig_measure_rhocoph.png`; the null is *named* here (the shuffle
  panel), the cohort evidence (panel b, 10/10) lives on 13. Ensure 13 / nulls / scale slides ADD
  (multiscale ρ_sym(τ) curves, band comparison, scale-signature), not repeat.

6. REFERENCES
- Cophenetic distance / UPGMA tree: standard hierarchical-clustering quantity (Sokal & Rohlf 1962, cophenetic correlation). The split-half symmetrized ρ^coph(τ) trace estimator + the matched-strength null are ours.
- Spearman rank correlation — standard.

7. CANVA STATUS
Page 17 · math column FINAL (LaTeX above locked to the deck). Figure = the reinstatement colored circular dendrograms (panel a) with per-patient panel b. OPEN (this session): (i) add the Δ-scatter conceptual image (built → `fig_measure_scatter.png`); (ii) decide whether panel b stays (recommend keep + name the matched-strength null in one line, defer the "what it destroys" cartoon to the nulls slide). ⚠ Coordination: S2 owning the ρ^coph formula means slide 11 (pipeline) should NOT also carry the ρ estimator formula (2) — keep only the D^coph definition there.
