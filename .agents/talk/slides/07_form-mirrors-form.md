---
name: talk-slide-07-form-mirrors-form
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 7
status: locked
updated: 2026-07-16
canva: page 8 · ~55%
---

# Slide 7 — Form mirrors form

1. TITLE
Form Mirrors Form
SUBTITLE (McNaughton 2010, presenter-added in Canva): "Cortical hierarchies … and the extraction of knowledge from memory"

2. MAIN CONCEPT
- Head (the juice): the task is a relational hierarchy — and this bet is not a shot in the dark. The consolidation literature (slide 6) says the brain builds structure OFFLINE, at rest, and lays cortical knowledge down as a HIERARCHY (McNaughton 2010). So we do not *guess* the resting network is hierarchical — we *expect* it. Our one added move is to read that hierarchy as a **multiscale** object.
- Why HIERARCHICAL is the expectation, not our bet: cortical knowledge is extracted from memory offline and organized as a hierarchy (McNaughton 2010, "Cortical hierarchies … and the extraction of knowledge from memory") — and slide 6 already put the offline/rest window and the structured cognitive map on the table (Ellenbogen, replay, Behrens). Hierarchy is what that literature predicts the consolidated trace should look like; we inherit it, we don't invent it.
- What IS our bet: that the hierarchy is **multiscale** — nested communities read across diffusion scales — not reducible to a single edge or to one fixed partition. The multiscale diffusion lens (slide 10, introduced next in the methods) is the new instrument; consolidation, seen *through that lens*, is what we point it at.
- Form mirrors form: the task's relational hierarchy ↔ the resting network's multiscale hierarchy. The figure IS the bet — the functional network on the brain gives rise to a cophenetic tree of relations, network form → hierarchical form, and we wager that tree mirrors the task's relational structure.
- Depth, not exclusivity (keep the July-13 honesty): reading a *structured hierarchy* rather than a scalar is what later lets us say WHICH bands carry the trace and AT WHAT SCALE it lives. Plain pairwise FC does detect a change; the multiscale read *characterizes* it. Frame the edge as information-DEPTH — never "ordinary metrics are blind."
- No null on this slide (user, 2026-07-16): the matched-strength guardrail is a results-phase gate — it lives on slide 6 (motivation) and in the results act, not on the bet/vision slide where there is no claim yet to defend. Hand-off: first we extract FC from these signals without artefacts.

3. ON-SLIDE TEXT
the task is a relational hierarchy — is the resting network one too?

   « Cortical hierarchies … and the extraction of knowledge from memory »   — McNaughton 2010

• why we expect a hierarchy (slide 6): knowledge is extracted OFFLINE, at rest — the cortex lays it down as a HIERARCHY  →  the resting network should carry one too
• our step: read that hierarchy as a MULTISCALE object
      – not a single edge · not one fixed partition · the whole nested shape, across scales
• form mirrors form — the trace we look for: held · band-selective · scale-structured

4. SPEECH  (~95 s — OPENS by departing from why-it-matters: consolidation/replay → HIERARCHY → the bet)
Everything on the last slide — consolidation, replay, the cognitive map — is the brain doing its work offline, at rest. Here's the link that turns it into a bet you can test. When the cortex sorts out what it learned, it doesn't just tuck facts away — it lays them down as a *hierarchy*. That's McNaughton's point, in the line above [gesture the McNaughton subtitle]: knowledge is extracted from memory offline and laid down hierarchically.

Now look at our task — it *is* a relational hierarchy, a chain of ordered relations [gesture the small task tree, top-right]. So the two forms should meet: if consolidation lays knowledge down as a hierarchy, and the task already is one, the resting network afterwards should carry that same form. That's the wager — form mirrors form. We're not guessing it's hierarchical; the consolidation literature predicts it.

What we *add* is a single move: read that hierarchy as a *multiscale* object — the whole nested shape across scales, not a single edge and not one fixed partition. Picture the functional network on the brain giving rise to this tree of relations [gesture brain → dendrogram]. If we're right, we don't just learn that something moved at rest — we'll say which frequency bands carry the trace, and at what scale it lives, which a single number never could.

First, though, a methods problem: how do we even measure connectivity in these signals without fooling ourselves?

Careful (do NOT say):
- Keep the two halves distinct: HIERARCHY is the EXPECTATION (McNaughton 2010 + slide 6's consolidation literature); MULTISCALE is our added hypothesis/lens. Don't present "multiscale" as something the literature already established, and don't present "hierarchical" as our novel guess. The inspiration does the motivating; the multiscale read is the contribution.
- This is the BET (the vision), not a result — do not present held / band-selective / scale-structured as already proven here.
- Do NOT say the trace is "invisible to" or "only visible to" plain metrics or spectral clustering. Raw FC DOES detect the β trace (p=.024). The "not a single edge, not one partition" line describes the OBJECT's depth (what we read), NOT exclusive detection. The edge is information-DEPTH — which bands, at what scale.
- Do NOT say "the nulls" / "strength · time" / "survives both nulls." One null: matched-strength (drift retired 2026-07-12).
- No band-specific or anatomy claims yet on this slide — no "β → OFC", no "inference β-only." Those are for the results act, and OFC survives only as the encoding anchor.
- McNaughton 2010 here is the CORTICAL-HIERARCHIES / offline-extraction paper — distinct from Wilson & McNaughton 1994 (the founding replay paper cited on slide 6). Don't conflate the two McNaughton citations.

5. FIGURES
- HERO — this figure IS the bet ("form mirrors form"): the whole machine in ONE 3D scene — the functional network on a translucent brain (contact spheres + |ImCoh| arcs) at the base, the cophenetic dendrogram rising above the head, diffusion root-streams tying each contact to its leaf; community colours link the spatial network to the tree clades so network form → hierarchical form is legible at a glance. data/outputs/figures/talk/pipeline_brain_tree_3d.png (interactive: .html). NO on-figure text (labels in Canva) → method-neutral by construction. gen `scripts/07_figures/gen_pipeline_brain_tree_3d.py` (Pat_05 β rest_pre; camera = "idea rising out of the head" 3/4 view).
  - ✅ STILL VALID under the refocus: this hero encodes no exclusive-detection and no anatomy/OFC claim — it is a method-neutral vision figure (network form → hierarchical form). No rebuild needed for slide 7. Update the *spoken* framing only (information-depth, not "plain metrics are blind").
  - ⚠ REUSE with slide 11: this is ALSO the planned pipeline hero. On the BET show it as the VISION/teaser ("what we wager the brain looks like"); on slide 11 show it as the step-by-step BUILD (timeseries → |ImCoh| → diffusion → dendrogram). Coordinate with the slide-11 owner so it doesn't read as a plain repeat — or give each a distinct camera. [Conductor flag — do NOT let the same still sit unchanged on both.]
  - ⚠ Background is BLACK (dramatic); the deck is LIGHT → presenter decides full-bleed dark hero vs a white-matte variant (gen has a `--qa` white-matte PNG).
- (Optional small callback) the race bracket (task's relational hierarchy) in a top corner + a "≟" glyph, to make the task↔brain mirror explicit and call back to slide 2 — data/outputs/figures/talk/ti_race_bracket.pdf. Still valid.
- REFOCUS FLAG (figures on OTHER slides, not this one — for the conductor): retire the 3-co-equal-reads synthesis and relabel `fig_before_nulls` on slide 14; do NOT rebuild an OFC-localization fig for the β trace (localization slide DROPPED 2026-07-14; β is DELOCALIZED — use the laterality brain instead; OFC belongs only to the encoding panel on slide 16).
- (Dropped 2026-07-11) external cognitive-map figures (Behrens / Park) — user: "the figure from the other papers are completely wrong… not the one we are looking for." "Form mirrors form" is carried by OUR brain network→tree, not a borrowed grid-code schematic.

6. REFERENCES
Why "hierarchical" is the expectation (the inspiration this slide runs on):
- ★ McNaughton, B. L. (2010), "Cortical hierarchies, sleep, and the extraction of knowledge from memory", Artificial Intelligence 174(2), 205–214 — the neocortex extracts the invariant structure of experience offline and organizes it as a hierarchy; this is why a hierarchical resting trace is the *expectation*, not our guess. The on-slide subtitle is drawn from this paper. (Distinct from Wilson & McNaughton 1994, the founding replay paper cited on slide 6 — do not conflate.)
- Slide 6 already established the offline/rest window and the structured cognitive map (Ellenbogen 2007; Liu–Behrens 2019 replay; Behrens 2018). This slide builds on those; it does not re-cite them on-slide.

Prior art that cognition builds relational representations (cited as motivation only — NO borrowed figure on-slide):
- Park, S. A., Miller, D. S., Boorman, E. D. (2021), "Inferences on a multidimensional social hierarchy use a grid-like code", Nature Neuroscience 24(9), 1292–1301.
- (Behrens et al. 2018 "What Is a Cognitive Map?", Neuron 100(2), 490–509 — concept only; its figure was rejected for this slide.)

7. CANVA STATUS
⚠ 2026-07-16 SPEECH bridge + RENAME (user): slide renamed "The bet" → "Form Mirrors Form" (file `07_form-mirrors-form.md`, matches Canva title). The speech now OPENS by departing from the why-it-matters slide — consolidation / replay / the cognitive map → the LINK that they lay knowledge down as a HIERARCHY (McNaughton) → the task is itself a relational hierarchy → "form mirrors form" → our multiscale add. Figure gestures added (McNaughton subtitle · the small task tree · brain → dendrogram). Body/on-slide text unchanged (user: "mostly good"). Guardrails intact (hierarchy=expectation vs multiscale=our add; bet not result; depth not exclusive detection; no null on this slide; no anatomy).
⚠ 2026-07-16 REWEIGHTED to the inspiration (user): the slide now spends its energy on WHY we expect a hierarchy — slide 6's consolidation literature + McNaughton 2010 (cortex extracts knowledge offline, organized as a hierarchy) — and treats MULTISCALE as our one added hypothesis/lens. Dropped the "working hypothesis" hedging down to a single matched-strength guardrail. Presenter added a McNaughton 2010 subtitle in Canva; the on-slide subtitle here is a placeholder (the paper's title phrase) — match it to the exact Canva wording. The old ✗/✓ list ("✗ not single edges (ordinary graph metrics) / ✗ not one fixed partition (spectral clustering) / ✓ multiscale hierarchy") is RETAINED as a trichotomy but reframed to describe the OBJECT's depth ("not a single edge · not one fixed partition · the whole nested shape, across scales") — NOT "ordinary metrics are blind" (that exclusive-detection reading stays retired; raw FC detects the β trace). Kept the July-13 information-depth honesty. The matched-strength null line was REMOVED from this slide (user) — no claim is made here to defend; the guardrail lives on slide 6 and in the results act, not on the bet/vision slide. RE-RENDER on deck.
⚠ 2026-07-13 REFOCUSED (settled three results): reframed the bet from exclusive detection ("invisible to ordinary graph metrics / spectral clustering") to information-DEPTH (scalar-per-band vs structured object; which bands / at what scale); fixed the null to matched-strength only (was "strength · time"). RE-RENDER on deck.

Page 8 · ~55% (LOCKED 2026-07-11, revised 2026-07-13 to the information-depth frame). Figure LOCKED: the 3D "form mirrors form" hero (pipeline_brain_tree_3d.png) — brain network at the base, cophenetic tree rising above, diffusion streams linking them, community colours as the glue; still valid (method-neutral, no anatomy claim). Optional small race-bracket callback + "≟". On-slide text now contrasts simple FC's scalar-per-band with the multiscale structured object (selectivity + scale-shape) rather than negating ordinary metrics/spectral clustering; the null line reads "matched-strength." ⚠ REUSE with slide 11 (pipeline) — coordinate (teaser vs build). ⚠ black bg on a light deck — full-bleed dark or white-matte variant (presenter call). Dropped the Behrens/Park external figures (user). Remaining to 100%: presenter updates the text block + speech framing in Canva, places the hero (+ optional bracket), resolve the 8↔11 reuse.
