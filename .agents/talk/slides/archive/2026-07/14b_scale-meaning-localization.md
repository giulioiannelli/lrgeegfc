---
name: talk-slide-14b-scale-meaning-localization
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 14b   # sits between 14 (detect ≠ discriminate + resolving power + higher order) and 15 (inference)
status: draft
updated: 2026-07-14
canva: REPURPOSED 2026-07-14 (user regroup). This slide = the CONCEPT ONLY: what the scale parameter s
  physically means (diffusion time → coarse-grains topology AND space, micro → meso → macro). The old
  "and where it points" localization half is CUT — anatomical localization was judged not presentable and
  the slide-16 localization reveal was dropped (16 is now epilepsy). No forward-pointer to a localization
  slide. The old figure-2 (fine-atlas-vs-hemisphere ruler-mismatch teaser) is CUT.
---

# Slide 14b — What the scale parameter means

1. TITLE
The scale parameter isn't a knob — it coarse-grains the network, in topology and in space

2. MAIN CONCEPT
- **Head (the juice):** slide 14 showed what the scale-swept read *buys* — discrimination, scale-resolution,
  a higher-order view. Here is *what the parameter is*: the diffusion scale s = τ·λ_max coarse-grains the
  network — **not just its topology, but its space**. Small s reads many fine, local populations; large s
  reads a few coarse, global blocks. That single idea is the spine of the whole "multiscale" thesis.
- **What s means.** The propagator e^(−τL̂) lets each node's activity diffuse for a time τ; we report the
  dimensionless s = τ·λ_max. **Small s → fine, local structure; large s → coarse, global structure.**
  Sweeping s reads the *same* network at every resolution; a single-scale measure is one horizontal slice.
- **Coarse-graining, made concrete (the main figure).** One fixed implant, read at three diffusion scales:
  MICRO (~O(n) single-population clusters), MESO (~n/10 balanced ~10-node groups), MACRO (~2 blocks). As s
  grows the local populations merge into progressively coarser communities — that IS the coarse-graining.
- **Functional, not spatial (the reach inset).** The physical length the diffusion connects grows from the
  finest functional cluster (~15 mm — already well ABOVE the 3.5 mm electrode pitch, because |ImCoh| nulls
  same-shaft coupling) to the implant span (~9 cm) by s ≈ 2. So the coarse-graining is over a FUNCTIONAL
  graph: even the finest unit is delocalised. (Honest: the scale is a topological ruler; the mm mapping is
  descriptive geometry, not a clean scale↔granularity diagonal.)

3. ON-SLIDE TEXT
what is s?  —  s = τ·λ_max, the diffusion time of  e^(−τ L̂)
  small s → fine / local / many     ·     large s → coarse / global / few     ·     one network, every resolution

the key idea:  diffusion coarse-grains the network  —  micro → meso → macro
  single populations → ~10-node groups → ~2 blocks   as s grows

(inset) physical reach: finest functional cluster ≈ 15 mm (≫ 3.5 mm pitch) → implant span ≈ 9 cm — functional, not spatial

4. SPEECH
One more thing about the scale parameter, because it's the spine of everything I've called "multiscale."
What is it, really? We let each contact's activity diffuse across the network for a time — small scale reads
the fine, local structure, large scale reads the coarse, global structure — and a single number, s = τ·λ_max,
sets that resolution. Read the same implant at three scales and you see it directly: at the finest, dozens of
single-population clusters; at a middle scale, a handful of ~ten-node groups; at the coarsest, two blocks.
That progressive merging is coarse-graining, running on the connectivity graph. And one honest detail from
the little inset: this coarse-graining is *functional*, not spatial — even the finest functional cluster is
already about fifteen millimetres across, well above the electrode spacing, because imaginary coherence
throws away the near, same-shaft coupling. So "scale" here is a ruler over the network's communication
structure, and that's the lens for everything that follows.

Careful: CONCEPT slide — no new gated statistic (numbers live on 14). Do NOT reintroduce the localization /
"coarse home → β-left" bridge — anatomical localization is not presentable and there is no localization slide
(16 is epilepsy). The mm reach is DESCRIPTIVE geometry, not a hypothesis test and not a scale↔granularity
diagonal — say "functional, not spatial", don't claim a physical-resolution ladder. The micro/meso/macro
panel is on the DENSE graph (only the dense graph gives balanced 3-regime coarse-graining; the sparse
backbone is core-dominated) — tag "illustrative / mechanism", not a data claim.

5. FIGURES
- **MAIN — the "what s does" coarse-graining panel (BUILT).**
  `data/outputs/figures/talk/fig_scale_coarse_grains_space.png`
  (gen `scripts/07_figures/talk_fig_scale_coarse_grains_space.py`, transparent, grey glass shell, Canva-ready).
  One fixed implant (Pat_05, β, rest_pre) read at three diffusion scales: MICRO (~46 single-population
  clusters, avg ~2-3 nodes) · MESO (12 balanced ~10-node groups) · MACRO (2 blocks). Coloured contacts +
  strongest intra-community edges bowed inward (contacts read on top); size-rank colour so the biggest
  community keeps its hue across scales. Illustrative mechanism (dense graph — the only graph that yields a
  balanced 3-regime coarse-graining), NOT a data claim. No dendrogram, no text on the asset.
- **INSET (user places) — the spatial-reach curve.**
  `data/preprint/figures/new_results_sec1/fig_spatial_reach.png` as a small inset axis on the panel.
  ℓ(s) = physical length the diffusion connects, cohort median ± IQR: rises from the finest functional
  cluster (~15 mm ≫ the 3.5 mm pitch — |ImCoh| has no strong same-shaft edges) to the implant span (~9 cm)
  by s ≈ 2. Makes "functional, not spatial" concrete. NOTE: the curve deliberately never reaches the 3.5 mm
  pitch line — that gap is the point; for the small inset consider dropping the pitch reference line.
- No gated-statistic figure belongs here (evidence is on 14).
- CUT: the old figure-2 ruler-mismatch / fine-atlas-vs-hemisphere localization teaser (points at a dropped
  slide).

6. REFERENCES
- Scale s = τ·λ_max = the LRG diffusion time (Villegas et al. 2023, *Nat. Phys.*; Villegas et al. 2025,
  *PRR*) — mechanism cited on the pipeline slide.
- Spatial-reach geometry: `scripts/01_compute/sparsified_arc/26_diffusion_spatial_reach_mst020.py`
  (descriptive, no null).

7. CANVA STATUS
REPURPOSED (user regroup 2026-07-14). Concept-only slide: what s physically means (diffusion time, local →
global, one network at every resolution) made concrete by the micro → meso → macro coarse-graining panel,
with the spatial-reach curve as an inset showing the coarse-graining is functional (finest cluster ~15 mm ≫
pitch) not spatial. The localization "where it points" half + figure-2 are CUT; drop any "see slide 16"
pointer. Main panel BUILT (`fig_scale_coarse_grains_space.png`, grey shell). Missing on deck: place the panel
+ reach inset, terse on-slide text, presenter notes. Also action on slide 15: remove its "the where is slide
16" forward-pointers (no localization slide).
