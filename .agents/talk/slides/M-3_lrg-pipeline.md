---
name: talk-slide-M3-lrg-pipeline
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: M-3
part: II — Methods
duration: ~70 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - .agents/talk/slides/I-5_diffusion-lrg.md
---

# M-3 · The LRG pipeline — matrix to hierarchy

## Slide placeholder (copy into Canva)

**Main point.** The pipeline made concrete: build the diffusion operator, let heat
flow for time τ, read the **communication hierarchy**. The heat kernel sees the
whole web (all indirect paths, Σₙ L̂ⁿ), not one edge. We show it **recovers**
known multiscale structure on a synthetic hierarchical network, then apply it to a
**real complete-weighted brain connectome**.

**Concepts to land.**
- The chain (visual, minimal algebra):
  **Â (imcoh_abs) → L̂ = D̂ − Â → K̂(τ) = e^{−τL̂} → communication distance
  D_ij = 1/K_ij → UPGMA dendrogram at τ = 1/λ_max → cophenetic distance**
  (the merge height of the lowest common ancestor).
- The heat kernel **sees the whole web** — indirect paths, walks of every length,
  Σₙ L̂ⁿ — the I-4 idea made concrete.
- **Validation on known ground truth:** on a synthetic **hierarchical modular
  network (HMN)** the reprojected propagator lights up **nested blocks**, and
  diffusion picks out the *known* levels in turn → LRG demonstrably **recovers**
  multiscale structure.
- **The payoff — real data:** a real **complete weighted connectome** (imcoh_abs,
  ~100% dense — the Villegas *outlier* regime from I-5) yields the **same kind of
  multiscale communities** under diffusion. So the method is **not** confined to
  sparse modular graphs — it reads hierarchy out of a dense weighted correlator,
  which is *why* it's legitimate on imcoh_abs FC at all.
- **"Multiscale" = the dendrogram's own nested cuts at one τ.** No τ-sweep needed:
  at τ = 1/λ_max the propagator already encodes the whole tree; you read scales by
  sliding the **cut height** (I-5's continuous-spectrum case, shown concretely).

**Figures / visuals.**
- **The combined demo — the single M-3 slide.** 2×4 grid: **row 1 = sparse HMN,
  diffusion-τ sweep** (a new tree for every τ) vs **row 2 = real β FC, one tree at
  τ_max + sliding cut**; LRG-native bubble layout; block matrix direct-coloured by
  community.
  - dark slide: `data/outputs/figures/lrg_diffusion_zoom/lrg_multiscale_combined.mp4`
  - light slide: `…/lrg_multiscale_combined_transparent.gif`
  - static contact sheet: `…/lrg_multiscale_combined_snapshots.pdf`
- Siblings (for a deeper beat): `…/lrg_diffusion_zoom.gif` (the I-5 teaser),
  `…/lrg_dendrogram_cut.gif` (+ `…_snapshots.pdf`).
- Real-data grounding: `data/outputs/figures/network_templates/matrix_plus_network/`
  and `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf`.

**References.**
- Villegas, Gili, Caldarelli & Gabrielli 2023, *Nature Physics* — L̂, K̂, ρ̂.
  https://doi.org/10.1038/s41567-022-01866-8
- Villegas, Gabrielli, Poggialini & Gili 2025, *Phys. Rev. Research* —
  communication distance, ultrametric, UPGMA.
  https://doi.org/10.1103/PhysRevResearch.7.013065

---

## Keep honest (content constraints, not styling)

- **Do NOT show C(τ) / specific-heat curves** — flat for our continuous spectrum;
  the multiscale content is the **dendrogram at τ = 1/λ_max** (see I-5).
- **Combinatorial Laplacian only** (L̂ = D̂ − Â); τ = 1/λ_max is load-bearing,
  never `optimal_threshold`.
- **The real β connectome is weakly hierarchical** — one dominant core + a few
  peripheral communities. Present as "the same readout applies to real data," NOT
  deep nesting.
- **Don't conflate τ with a cut height.** In the combined demo the sparse row's
  tree *morphs with τ* (a family of trees); the dense row is *one tree + a sliding
  cut*. That contrast IS the teaching point (I-5).
- imcoh_abs is **second-order**; LRG yields higher-order **graph** features
  ("emergent from pairwise connectivity") — never "nonlinear."
