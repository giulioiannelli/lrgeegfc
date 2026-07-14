---
name: talk-slide-14b-scale-meaning-localization
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 14b   # sits between 14 (detect ≠ discriminate + resolving power + higher order) and 15 (inference)
status: draft
updated: 2026-07-13
canva: NEW (2026-07-13 PM, user regroup). This slide = the CONCEPT: what the scale parameter s physically
  means, and the link to localization (s coarse-grains SPACE → a coarse home). The empirical resolving-power
  + higher-order EVIDENCE lives on 14.
---

# Slide 14b — What the scale parameter means — and where it points

> ⚠ **UNDER REVIEW (2026-07-14):** the localization slide this bridges to ("16 · Where the trace lives") was **DROPPED** — the β "coarse left-hemisphere home" is incoherent (mesoscale *diffusion* scale conflated with *macro* anatomy) and only marginal. Every "slide 16" / β-LEFT pointer below is **STALE**. Cut this slide, or repurpose it to keep only the *what the scale means* content (discrimination / scale-resolution) with the localization bridge removed.

1. TITLE
The scale parameter isn't a knob — it coarse-grains the network, in topology and in space

2. MAIN CONCEPT
- **Head (the juice):** slide 14 showed what the scale-swept read *buys* — discrimination, scale-resolution,
  a higher-order view. Here is *what the parameter is*, and the one idea to carry forward: the diffusion
  scale s = τ·λ_max coarse-grains the network — **not just its topology, but its space** — so the scale at
  which the trace lives also sets the *spatial resolution* at which its home appears. A coarse-scale trace
  has a **coarse home**. That is the bridge to "where does the trace live."
- **What s means.** The propagator e^(−τL̂) lets each node's activity diffuse for a time τ; we report the
  dimensionless s = τ·λ_max. **Small s → fine, local structure; large s → coarse, global structure.**
  Sweeping s reads the *same* network at every resolution, and the hierarchy is exactly that sweep stacked
  into a tree — merges at small s are fine sub-communities, merges at large s are the coarse ones. A
  single-scale measure is one horizontal slice of that tree.
- **Diffusion coarse-grains *space*, not only topology.** As s grows, nodes that co-diffuse merge — and on
  an implanted brain, co-diffusing contacts are also *spatially* extended groups. So increasing s literally
  zooms the spatial resolution out: from a single contact, to a local cluster, to a lobe, to a hemisphere.
  The scale parameter is a spatial ruler as much as a topological one.
- **The link to localization (the bridge to slide 16).** This resolves a puzzle we hit head-on next: read
  the β trace against a **fine** anatomical atlas and it looks **placeless** — no single region survives.
  But that is a *ruler mismatch*: a coarse-scale trace tested at fine spatial resolution has nowhere to
  land. Match the parcellation to the scale — read at the **hemisphere** level — and β resolves cleanly to
  the **left hemisphere**. So "where does the trace live" is a scale-matched question, and we'll answer it
  that way. *(Honest: this is a motivating lens, not a proven one-to-one scale↔granularity diagonal — β's
  hemisphere signal sits at the mesoscale, not the coarsest.)*

3. ON-SLIDE TEXT
what is s?  —  s = τ·λ_max, the diffusion time of  e^(−τ L̂)
  small s → fine / local     ·     large s → coarse / global     ·     the hierarchy = the whole sweep

the key idea:  diffusion coarse-grains SPACE, not just topology
  contact → cluster → lobe → hemisphere  as s grows  →  a coarse-scale trace has a COARSE spatial home

→ so "where does the trace live?" is a SCALE-MATCHED question:
  fine atlas → β looks placeless   ·   matched coarse (hemisphere) map → β resolves LEFT   (slide 16)

4. SPEECH
One more thing about the scale parameter, because it's the idea that unlocks the next act. What is it,
really? We let each contact's activity diffuse across the network for a time — small scale reads the fine,
local structure, large scale reads the coarse, global structure — and the hierarchy is just that whole
sweep, stacked into a tree: the fine merges near the bottom, the coarse ones near the top. Now here's the
part I want you to hold onto. This diffusion doesn't only coarse-grain the *topology* — it coarse-grains
*space*. As the scale grows, the contacts that merge are also spatially extended groups: you go from a
single contact, to a local cluster, to a lobe, to a whole hemisphere. The scale is a spatial ruler. And
that resolves a puzzle we're about to walk straight into: if you ask where the beta trace lives using a
fine anatomical atlas, it looks placeless — nothing survives. But that's the wrong ruler. A coarse-scale
trace has a coarse home. Match the map to the scale — read at the hemisphere level — and beta resolves to
the left hemisphere. So "where does the trace live" isn't a dead end; it's a scale-matched question, and
that's exactly how we'll answer it.

Careful: this is the CONCEPT slide — no new gated statistic here (the resolving-power / higher-order NUMBERS
are on 14; the localization NUMBERS are on 16). "s = τλ_max / e^(−τL̂)" is a recap from the pipeline slide —
keep it light, don't re-derive. The scale↔space link is a MOTIVATING lens, NOT a proven clean
scale↔granularity diagonal (β's hemisphere signal is at the mesoscale s≈11, not the coarsest) — say "coarse
spatial home" and flag the diagonal is not claimed. Do NOT pre-empt slide 16's β-left numbers / SOZ-
independence / low-γ→SOZ — plant only "fine → placeless, coarse → left hemisphere". Matched-strength is the
sole null. NO inference here (→15).

5. FIGURES
- **MAIN — the "what s does" schematic (BUILT).** `data/outputs/figures/talk/fig_scale_coarse_grains_space.png`
  (gen `scripts/07_figures/talk_fig_scale_coarse_grains_space.py`, transparent, Canva-ready). One fixed
  brain-like network read at three diffusion scales: small s = ≈12 local blobs, mid s = ≈6 lobes, large s =
  2 hemispheres, left = blue family / right = orange family so the eye tracks "left stays left, it just
  merges" (the endpoint that motivates β-left on slide 16). Under it a single dendrogram with three cut lines
  (`large s → 2 · mid s → 6 · small s → 12`) makes "the hierarchy IS the sweep" literal — the three panels
  are three cuts of one tree. Tagged "schematic" (illustrative mechanism, NOT a data claim; no clean
  scale↔granularity diagonal asserted). A ward hierarchy on synthetic node positions — a cartoon of the
  mechanism, honest to the slide's caveat.
- **OPTIONAL — the ruler-mismatch teaser.** A faded two-panel: fine atlas (β nothing lights) vs hemisphere
  map (β left glows), as the forward pointer to slide 16. Keep it faint / schematic so it doesn't spend the
  slide-16 reveal; or drop and let slide 16 own it.
- No gated-statistic figure belongs here (evidence is on 14 and 16).

6. REFERENCES
- Scale s = τ·λ_max = the LRG diffusion time (Villegas et al. 2023, *Nat. Phys.*; Villegas et al. 2025,
  *PRR*) — mechanism cited on the pipeline slide. The scale-matched-parcellation localization idea is ours
  (`23_localization_statistic_bakeoff_mst020`, ladder).

7. CANVA STATUS
NEW (user regroup 2026-07-13 PM). The concept + bridge slide: what s physically means (diffusion time, local
→ global, the hierarchy = the sweep) and the one idea to carry — diffusion coarse-grains SPACE, so a
coarse-scale trace has a coarse spatial home → the scale-matched localization question answered on slide 16.
The "what s does" schematic is BUILT (`fig_scale_coarse_grains_space.png`) and ready to drop. Missing on deck:
place the schematic, optional ruler-mismatch teaser, terse on-slide text, presenter notes. Conceptual breather
between the evidence (14) and the results (15 inference, 16 localization); it plants the scale↔space bridge
that slide 16 cashes.
