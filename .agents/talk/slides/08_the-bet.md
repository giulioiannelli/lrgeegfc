---
name: talk-slide-08-the-bet
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 8
status: locked
updated: 2026-07-11
canva: page 8 · ~55%
---

# Slide 8 — The bet: form mirrors form

1. TITLE
The bet — form mirrors form

2. MAIN CONCEPT
- The intuition: the brain's ability to infer relations may be written in the multiscale structure of the functional relations between its areas.
- What if the learned relational map lives in the higher-order, multiscale structure of the network — nested communities, a whole hierarchy — rather than in single edges?
- Not a brand-new idea that cognition resides in multiscale features of the brain — but untested with this lens.
- Our bet: the trace of inference lives in the multiscale topology — invisible to ordinary network metrics and to spectral clustering of the FC.
- The figure IS the bet: the functional network on the brain gives rise to a cophenetic tree of relations — network form → hierarchical form — and we wager it mirrors the task's relational hierarchy.
- Guardrail: form mirroring form is a hypothesis to test, not a proof — the nulls decide (strength · time).
- Hand-off: first we need to extract FC from these signals without artefacts.

3. ON-SLIDE TEXT
the task is a relational hierarchy — is the brain's network one too?
the bet: the inferred map lives in the multiscale topology
   ✗ not single edges (ordinary graph metrics)
   ✗ not one fixed partition (spectral clustering)
   ✓ nested communities — a whole hierarchy
predicted: a trace that is held · band-specific · multiscale
form mirrors form — a hypothesis, not proof; the nulls decide (strength · time)

4. SPEECH
Here's the bet. The brain is known to build relational maps for this kind of inference — but where does such a map live? Our wager: in the multiscale shape of the functional network — the way it organizes into nested communities, a whole hierarchy, not single connections. Picture the network on the brain giving rise to this tree of relations. It isn't a wild idea — cognition has been tied to multiscale brain features before — but nobody has looked with this lens. So we bet the trace of inference lives in the multiscale topology, where plain graph metrics and spectral clustering are blind. Form mirrors form — a hypothesis, not a proof, and the rest of the talk is us trying to break it. First, though, a methods problem: how do we even measure connectivity here without fooling ourselves?

5. FIGURES
- HERO — this figure IS the bet ("form mirrors form"): the whole machine in ONE 3D scene — the functional network on a translucent brain (contact spheres + |ImCoh| arcs) at the base, the cophenetic dendrogram rising above the head, diffusion root-streams tying each contact to its leaf; community colours link the spatial network to the tree clades so network form → hierarchical form is legible at a glance. data/outputs/figures/talk/pipeline_brain_tree_3d.png (interactive: .html). NO on-figure text (labels in Canva) → method-neutral by construction. gen `scripts/07_figures/gen_pipeline_brain_tree_3d.py` (Pat_05 β rest_pre; camera = "idea rising out of the head" 3/4 view).
  - ⚠ REUSE with slide 11: this is ALSO the planned pipeline hero. On the BET show it as the VISION/teaser ("what we wager the brain looks like"); on slide 11 show it as the step-by-step BUILD (timeseries → |ImCoh| → diffusion → dendrogram). Coordinate with the slide-11 owner so it doesn't read as a plain repeat — or give each a distinct camera. [Conductor flag — do NOT let the same still sit unchanged on both.]
  - ⚠ Background is BLACK (dramatic); the deck is LIGHT → presenter decides full-bleed dark hero vs a white-matte variant (gen has a `--qa` white-matte PNG).
- (Optional small callback) the race bracket (task's relational hierarchy) in a top corner + a "≟" glyph, to make the task↔brain mirror explicit and call back to slide 2 — data/outputs/figures/talk/ti_race_bracket.pdf.
- (Dropped 2026-07-11) external cognitive-map figures (Behrens / Park) — user: "the figure from the other papers are completely wrong… not the one we are looking for." "Form mirrors form" is carried by OUR brain network→tree, not a borrowed grid-code schematic.

6. REFERENCES
Prior art that cognition builds relational representations (cited as motivation only — NO borrowed figure on-slide):
- Park, S. A., Miller, D. S., Boorman, E. D. (2021), "Inferences on a multidimensional social hierarchy use a grid-like code", Nature Neuroscience 24(9), 1292–1301.
- (Behrens et al. 2018 "What Is a Cognitive Map?", Neuron 100(2), 490–509 — concept only; its figure was rejected for this slide.)

7. CANVA STATUS
Page 8 · ~55% (LOCKED 2026-07-11, revised — 3D hero + expanded text). Figure LOCKED: the 3D "form mirrors form" hero (pipeline_brain_tree_3d.png) — brain network at the base, cophenetic tree rising above, diffusion streams linking them, community colours as the glue. Optional small race-bracket callback + "≟". Text expanded to the 7-line block above; speech references the 3D vision. ⚠ REUSE with slide 11 (pipeline) — coordinate (teaser vs build). ⚠ black bg on a light deck — full-bleed dark or white-matte variant (presenter call). Dropped the Behrens/Park external figures (user). Remaining to 100%: presenter places the hero (+ optional bracket) in Canva; resolve the 8↔11 reuse.
