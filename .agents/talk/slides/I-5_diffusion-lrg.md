---
name: talk-slide-I5-diffusion-lrg
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-09
slide: I-5
part: I — Introduction
duration: ~50 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# I-5 · The turn — renormalization by diffusion (LRG)

## Slide placeholder (copy into Canva)

**Main point.** Let **heat diffuse** on the network; diffusion time **τ is a
continuous zoom**. LRG was built to read **sparse** networks — where **spectral
gaps** hand you one scale after another. Our FC is the opposite case: a
**complete** weighted graph with a **continuous** spectrum and no gaps — so we
read the hierarchy by **clustering the propagator as a relational distance**, and
the **continuum of scales lives in the dendrogram, not in τ**.

**Concepts to land.**
- Put **heat** on the network and let it **diffuse**; **τ = a continuous zoom**
  knob, local → global. Coarse-graining *by diffusion*.
- **LRG was created for sparse networks** to detect their structure: there,
  successive **spectral gaps** reveal **subsequent scales** — nested blocks light
  up one level at a time as diffusion proceeds. *(the animation shows exactly
  this, on a known-ground-truth network.)*
- **Our FC is the outlier case:** a **complete weighted graph** with a
  **near-continuous spectrum** — **no clean spectral gaps**, so no discrete
  sequence of scales.
- **So we read the hierarchy differently:** hierarchical **clustering** using the
  **propagator as a relational distance** (communication distance
  D_ij = 1/K_ij at τ = 1/λ_max). The **continuum of scales is reflected in the
  dendrogram**, not in τ.
- **A deliberate choice — we did *not* sparsify.** Sparsifying would have forced a
  clean block structure (clean gaps, discrete scales) — but we keep the
  **complete** network to **perturb it as little as possible**. Honest trade:
  minimal intervention, at the cost of no discrete spectral scales.
- Payoff of the I-1 wink: **Gabrielli** is an author of LRG.

**Figures / visuals.**
- **Animation (sparse case — spectral gaps → scales):** the diffusion-zoom on the
  synthetic HMN, nested blocks lighting up at successive τ.
  `data/outputs/figures/lrg_diffusion_zoom/lrg_diffusion_zoom.gif`
  (still: `…/lrg_diffusion_zoom_snapshots.pdf`; GIF only — sharp MP4 on request).
- The **complete-graph** payoff (dendrogram read off the propagator distance) is
  developed in **M-3** (combined demo).

**References.**
- LRG: Villegas, Gili, Caldarelli & Gabrielli 2023, *Nature Physics*
  (https://doi.org/10.1038/s41567-022-01866-8); Villegas, Gabrielli, Poggialini &
  Gili 2025, *Phys. Rev. Research* (https://doi.org/10.1103/PhysRevResearch.7.013065).
- Lineage: De Domenico & Biamonte 2016, *PRX* — spectral entropy / density-matrix
  formalism (https://doi.org/10.1103/PhysRevX.6.041062); building on Braunstein,
  Ghosh & Severini 2006, *Annals of Combinatorics* 10:291–317 — "The Laplacian of
  a graph as a density matrix" (the graph-Laplacian / heat-kernel foundation cited
  in that PRX). *DOI to verify (candidate 10.1007/s00026-006-0289-3).*

---

## Keep honest (content constraints, not styling)

- **This is the crux of why we read scales from the dendrogram, not C(τ):** with
  no spectral gaps the specific heat C(τ) / entropy curve is **flat and
  uninformative**, so the multiscale content is carried by the **dendrogram at
  τ = 1/λ_max**. Do **not** show C(τ) curves.
- **Not sparsifying is a deliberate, honest choice** (minimal perturbation of the
  complete weighted graph), and it puts us in the **Villegas *outlier* regime** —
  say so; don't present the complete-graph use as the canonical LRG application.
- **Combinatorial Laplacian only** (L̂ = D̂ − Â). τ is a diffusion time — never
  `optimal_threshold`. τ = 1/λ_max is load-bearing for the dendrogram.
- The tree is hierarchical **by construction** — the guarded "mirror" is I-6.

## Open decision (content)

- This slide now carries a real methods nuance (sparse-gaps vs continuous-
  dendrogram + the no-sparsify choice). If it feels heavy for the intro, the
  mechanics can move to **M-3** (which already shows the sparse-vs-complete demo)
  and I-5 stays the lighter "diffusion = zoom" turn.
