---
name: talk-slide-13-tanglegram-categories
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 13
status: draft
updated: 2026-07-11
canva: NEW — no Canva page yet
---

# Slide 13 — Task-induced taxonomy of neuronal populations (trace / anchor / reset / reorganized)

1. TITLE
Task-induced taxonomy of neuronal populations
(the "A 4-fold …" long form is fine too; keep it to the one title line)

2. MAIN CONCEPT
- Ground the scale first: each sEEG contact samples a ~mm-scale neuronal population, and the hierarchy is a hierarchy OF these populations. What follows is the fate of a population's PLACE in that hierarchy (its cophenetic relationships), not of its raw activity.
- Before the cohort statistics, one visual: how does a single patient's hierarchy change from phase to phase?
- A tanglegram pairs two dendrograms and links the same contacts — you see directly whether the structure held, moved, or reverted (lines parallel = order kept, crossing = reorganized).
- Four outcomes a pattern can fall into, one illustrative patient/band each:
  - trace — task reorganized it AND the change persists into rest (Pat_06, δ) — what we hunt
  - reset — task reorganized it AND it reverts in rest (Pat_10, β — our genuine β resetter; drawn with rest-post in the CENTRE so the reverted rest-pre ≈ rest-post pair sits adjacent and reads as parallel, task is the lone scrambled column)
  - anchor — the same hierarchy in every phase (Pat_08, low-γ)
  - reorganized — it keeps changing, never settling (different in every phase) (Pat_08, δ)
- This is the vocabulary for the results: we hunt traces, read against anchors, resets and reorganizations.

3. ON-SLIDE TEXT
TITLE ONLY + the figure. No bullets, no reading key, no "we hunt traces" caption. The single-panel
figure already carries the four category titles, the phase labels, and the ρ + crossing annotations,
so the visual defines the taxonomy on its own; everything else is spoken (see §4). This is the
deliberate design: one title, one figure, let the picture make the point.

(Everything below documents the choice — it is NOT text to place on the slide.)

4. SPEECH
Each of these contacts picks up a millimetre-scale population of neurons, and the hierarchy I've been showing you is a hierarchy of those populations. So the question is what happens to a population's place in that hierarchy. Before any statistics, here's what we're actually looking at. Take one patient, one band, and lay the hierarchy from the task next to the hierarchy from the rest afterwards; link the same contacts between them. Lines that stay parallel mean the order was kept; lines that cross mean it reorganized. Four things can happen. It can be reorganized by the task and hold into rest — that's a trace, what we're after. It can reorganize and then revert — a reset. It can stay put throughout — an anchor. Or it can keep changing in every phase, never settling — reorganized. These four words are the vocabulary for the rest of the talk.

Careful: these four are the ρ^coph cross-phase-similarity classes shown by the tanglegrams — trace / anchor / reset / reorganized (do NOT relabel the 4th "emergent"; emergent is a separate community-membership concept). The annotated similarity on each figure is LOCAL to the shown clade (≠ whole-tree ρ^coph) and the clades are illustrative, disclosed-selected — NOT a gated test. Trace requires BOTH task-reorganization AND persistence.

SPEECH → GRAPHICS AUDIT (does any spoken point need to be ON the figure?). Verdict: no text needs
adding — the graphic is self-defining for a network audience. Point by point:
- "parallel = retained, crossing = reorganised" (the reading key) — this is the ONE borderline
  item; it was the removed grey footnote. It is, however, already double-encoded: the ρ values are
  coloured green when held / red when reorganised, and each gap prints its crossing count. A network
  audience reads crossing-ribbons as reordering on sight. Left to speech.
- the four definitions ("reorganised at task then retained", etc.) — the category TITLE names the
  fate; the definition is spoken. Putting the descriptors back would re-clutter exactly what we
  stripped. Keep spoken.
- "each contact ≈ a ~mm neuronal population" and "we hunt traces" — framing / punchline, spoken.
If the slide ever has to stand ALONE (handout, no narration), the only thing worth adding is a tiny
reading-key glyph — two mini ribbons, one parallel labelled "retained", one crossing labelled
"reorganised" — in a corner. Not needed for the live talk.

5. FIGURES
- SINGLE-PANEL (current, import this one): taxonomy_quad_tanglegram.png (transparent, 300 DPI,
  ~4112×3091 px) + .pdf (vector) in data/outputs/figures/talk/. The four fates side by side,
  EQUAL HEIGHT, columns L→R = RESET · REORGANISED · ANCHOR · TRACE. Text minimised to exactly
  four things: category title, per-column phase labels, per-gap ρ + crossing counts, in-bead
  contact labels. Removed: the Pat/N/descriptor italic subtitle, the grey ρ_pre,post note, the
  sub-clade colour legend, and both honesty footnotes. Exemplars/clades are byte-identical to the
  four per-fate PNGs below (same patient·band·clade). RESET keeps its rest-pre | rest-post | task
  column reorder (gap 1 = pre↔post ρ 0.89 parallel/green, gap 2 = post↔task ρ 0.22 scrambled/red).
  - (generator: scripts/07_figures/talk_fig_taxonomy_quad.py — composites
    fig_cophenetic_tanglegram.build_row across four fixed exemplars; DRY.)
- Per-fate PNGs (the earlier four separate slides, still valid if you want one fate per slide):
  - trace — tanglegram_trace_Pat_06_delta.png (ρ pre→task 0.16 → task→post 0.71: reorganized then retained)
  - reset — tanglegram_reset_Pat_10_beta.png (COLUMNS REORDERED: rest-pre | rest-post | task, so gap 1 = pre↔post ρ 0.89 reads parallel/green and gap 2 = post↔task ρ 0.22 is the scrambled/red departure — the reversion is now visible, not just annotated)
  - anchor — tanglegram_anchor_Pat_08_low_gamma.png
  - reorganized — tanglegram_reorganized_Pat_08_delta.png
- Vector sources (the _DRAFT set, FINALIZE for the preprint):
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_trace_Pat_06_delta_N23_DRAFT.pdf
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_reset_Pat_10_beta_N30_DRAFT.pdf (post-in-centre reorder; regenerate via --classes reset --patients Pat_10 --bands beta --topk 1)
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_anchor_Pat_08_low_gamma_N23_DRAFT.pdf
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_reorganized_Pat_08_delta_N20_DRAFT.pdf
- (per-fate generator: scripts/01_compute/figures_embedded/fig_cophenetic_tanglegram_class.py)

6. REFERENCES
- Taxonomy / measure are ours (ρ^coph cross-phase similarity).

7. CANVA STATUS
NEW — no Canva page yet. Build right after the pipeline (ρ^coph) slide; this is the visual glance
of the four categories. Import the SINGLE-PANEL taxonomy_quad_tanglegram.png (transparent) and drop
it full-bleed — it already carries the four titles, phase labels, and ρ/crossing annotations, so no
extra on-slide text is needed beyond the slide title.
