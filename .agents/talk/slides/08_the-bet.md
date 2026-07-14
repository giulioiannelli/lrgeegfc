---
name: talk-slide-08-the-bet
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 8
status: locked
updated: 2026-07-13
canva: page 8 · ~55%
---

# Slide 8 — The bet: form mirrors form

1. TITLE
The bet — form mirrors form

2. MAIN CONCEPT
- Head (the juice): the task is a relational hierarchy; we bet the resting network holds one too — and that the trace of reasoning is written not in single edges but in the *multiscale shape* of the functional network, a structured object that simple connectivity can only reduce to a yes/no.
- The reframe that runs the whole talk: simple pairwise FC gives a **scalar verdict per band** — "did the reorganization lean back toward the task?" — and several bands will say yes. That scalar is the ceiling of what one number can carry.
- The multiscale diffusion read gives a **structured object** instead: a hierarchy per phase, a similarity curve per band. So it can *characterize* the trace — **which** bands carry it (selectivity), and **at what scale** it lives (scale-shape) — questions the scalar cannot even pose.
- Our bet, stated as information-DEPTH (not exclusivity): the inferred relational map is written in the multiscale topology — nested communities, a whole hierarchy — and the multiscale lens is what lets us read its band-selectivity and its scale-signature, beyond the bare "a band moved."
- Not a brand-new idea that cognition resides in multiscale brain features — but untested with this lens on this task.
- The figure IS the bet: the functional network on the brain gives rise to a cophenetic tree of relations — network form → hierarchical form — and we wager it mirrors the task's relational hierarchy.
- Guardrail: form mirroring form is a hypothesis to test, not a proof — the null decides (matched-strength).
- Hand-off: first we need to extract FC from these signals without artefacts.

3. ON-SLIDE TEXT
the task is a relational hierarchy — is the resting network one too?

simple FC → a scalar per band: "did rest lean back toward the task?"  (several say yes)
the multiscale read → a structured object: a hierarchy per phase, a curve per band
   → WHICH bands hold the trace   (selectivity)
   → AT WHAT SCALE it lives        (scale-shape)

the bet: the inferred map is written in the multiscale shape of the network
predicted: a trace that is held · band-selective · scale-structured
form mirrors form — a hypothesis, not proof; the null decides (matched-strength)

4. SPEECH
Here's the bet. The brain builds relational maps for this kind of inference — but where does such a map live, and how would we even read it? Simple pairwise connectivity can only hand us a scalar per frequency band: did the resting network lean back toward the task — yes or no. Several bands will say yes, and that's real, but a single number is the ceiling of what it can tell us. Our wager is that the trace of reasoning is written in something richer — the multiscale shape of the functional network, the way it nests into communities, a whole hierarchy read across scales. Picture the network on the brain giving rise to this tree of relations. If the bet is right, we don't just learn that a band moved; we learn which bands carry the trace and at what scale it lives. It isn't a wild idea — cognition has been tied to multiscale brain features before — but nobody has looked with this lens on this task. Form mirrors form: the task is a relational hierarchy, and we bet the resting network holds one too. It's a hypothesis, not a proof — and the rest of the talk is us trying to break it against the null. First, though, a methods problem: how do we even measure connectivity here without fooling ourselves?

Careful (do NOT say):
- This is the BET (the vision), not a result — do not present held / band-selective / multiscale as already proven here.
- Do NOT say the trace is "invisible to" or "only visible to" plain metrics or spectral clustering. Raw FC DOES detect the β trace (p=.024). The multiscale edge is information-DEPTH — which bands, at what scale — NOT exclusive detection.
- Do NOT say "the nulls" / "strength · time" / "survives both nulls." One null: matched-strength (drift retired 2026-07-12).
- No band-specific or anatomy claims yet on this slide — no "β → OFC", no "inference β-only." Those are for the results act, and OFC survives only as the encoding anchor.

5. FIGURES
- HERO — this figure IS the bet ("form mirrors form"): the whole machine in ONE 3D scene — the functional network on a translucent brain (contact spheres + |ImCoh| arcs) at the base, the cophenetic dendrogram rising above the head, diffusion root-streams tying each contact to its leaf; community colours link the spatial network to the tree clades so network form → hierarchical form is legible at a glance. data/outputs/figures/talk/pipeline_brain_tree_3d.png (interactive: .html). NO on-figure text (labels in Canva) → method-neutral by construction. gen `scripts/07_figures/gen_pipeline_brain_tree_3d.py` (Pat_05 β rest_pre; camera = "idea rising out of the head" 3/4 view).
  - ✅ STILL VALID under the refocus: this hero encodes no exclusive-detection and no anatomy/OFC claim — it is a method-neutral vision figure (network form → hierarchical form). No rebuild needed for slide 8. Update the *spoken* framing only (information-depth, not "plain metrics are blind").
  - ⚠ REUSE with slide 11: this is ALSO the planned pipeline hero. On the BET show it as the VISION/teaser ("what we wager the brain looks like"); on slide 11 show it as the step-by-step BUILD (timeseries → |ImCoh| → diffusion → dendrogram). Coordinate with the slide-11 owner so it doesn't read as a plain repeat — or give each a distinct camera. [Conductor flag — do NOT let the same still sit unchanged on both.]
  - ⚠ Background is BLACK (dramatic); the deck is LIGHT → presenter decides full-bleed dark hero vs a white-matte variant (gen has a `--qa` white-matte PNG).
- (Optional small callback) the race bracket (task's relational hierarchy) in a top corner + a "≟" glyph, to make the task↔brain mirror explicit and call back to slide 2 — data/outputs/figures/talk/ti_race_bracket.pdf. Still valid.
- REFOCUS FLAG (figures on OTHER slides, not this one — for the conductor): retire the 3-co-equal-reads synthesis and relabel `fig_before_nulls` on slide 13; do NOT rebuild an OFC-localization fig for the β trace (localization slide DROPPED 2026-07-14; β is DELOCALIZED — use the laterality brain instead; OFC belongs only to the encoding panel on slide 15).
- (Dropped 2026-07-11) external cognitive-map figures (Behrens / Park) — user: "the figure from the other papers are completely wrong… not the one we are looking for." "Form mirrors form" is carried by OUR brain network→tree, not a borrowed grid-code schematic.

6. REFERENCES
Prior art that cognition builds relational representations (cited as motivation only — NO borrowed figure on-slide):
- Park, S. A., Miller, D. S., Boorman, E. D. (2021), "Inferences on a multidimensional social hierarchy use a grid-like code", Nature Neuroscience 24(9), 1292–1301.
- (Behrens et al. 2018 "What Is a Cognitive Map?", Neuron 100(2), 490–509 — concept only; its figure was rejected for this slide.)

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): reframed the bet from exclusive detection ("invisible to ordinary graph metrics / spectral clustering") to information-DEPTH (scalar-per-band vs structured object; which bands / at what scale); fixed the null to matched-strength only (was "strength · time"). RE-RENDER on deck.

Page 8 · ~55% (LOCKED 2026-07-11, revised 2026-07-13 to the information-depth frame). Figure LOCKED: the 3D "form mirrors form" hero (pipeline_brain_tree_3d.png) — brain network at the base, cophenetic tree rising above, diffusion streams linking them, community colours as the glue; still valid (method-neutral, no anatomy claim). Optional small race-bracket callback + "≟". On-slide text now contrasts simple FC's scalar-per-band with the multiscale structured object (selectivity + scale-shape) rather than negating ordinary metrics/spectral clustering; the null line reads "matched-strength." ⚠ REUSE with slide 11 (pipeline) — coordinate (teaser vs build). ⚠ black bg on a light deck — full-bleed dark or white-matte variant (presenter call). Dropped the Behrens/Park external figures (user). Remaining to 100%: presenter updates the text block + speech framing in Canva, places the hero (+ optional bracket), resolve the 8↔11 reuse.
