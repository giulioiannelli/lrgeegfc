---
name: talk-slide-12-tanglegram-categories
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 12
status: draft
updated: 2026-07-11
canva: NEW — no Canva page yet
---

# Slide 12 — A taxonomy of neuronal populations (trace / anchor / reset / reorganized)

1. TITLE
A taxonomy of neuronal populations — trace, anchor, reset, reorganized

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
each sEEG contact ≈ a ~mm neuronal population
tanglegram = two dendrograms + matched-contact links
trace (held) · reset (reverted) · anchor (unchanged) · reorganized (never settles)
we hunt traces

4. SPEECH
Each of these contacts picks up a millimetre-scale population of neurons, and the hierarchy I've been showing you is a hierarchy of those populations. So the question is what happens to a population's place in that hierarchy. Before any statistics, here's what we're actually looking at. Take one patient, one band, and lay the hierarchy from the task next to the hierarchy from the rest afterwards; link the same contacts between them. Lines that stay parallel mean the order was kept; lines that cross mean it reorganized. Four things can happen. It can be reorganized by the task and hold into rest — that's a trace, what we're after. It can reorganize and then revert — a reset. It can stay put throughout — an anchor. Or it can keep changing in every phase, never settling — reorganized. These four words are the vocabulary for the rest of the talk.

Careful: these four are the ρ^coph cross-phase-similarity classes shown by the tanglegrams — trace / anchor / reset / reorganized (do NOT relabel the 4th "emergent"; emergent is a separate community-membership concept). The annotated similarity on each figure is LOCAL to the shown clade (≠ whole-tree ρ^coph) and the clades are illustrative, disclosed-selected — NOT a gated test. Trace requires BOTH task-reorganization AND persistence.

5. FIGURES
- Canva-ready PNGs (300 DPI, 2256×2857 px) in data/outputs/figures/talk/ — import these directly:
  - trace — tanglegram_trace_Pat_06_delta.png (ρ pre→task 0.16 → task→post 0.71: reorganized then retained)
  - reset — tanglegram_reset_Pat_10_beta.png (COLUMNS REORDERED: rest-pre | rest-post | task, so gap 1 = pre↔post ρ 0.89 reads parallel/green and gap 2 = post↔task ρ 0.22 is the scrambled/red departure — the reversion is now visible, not just annotated)
  - anchor — tanglegram_anchor_Pat_08_low_gamma.png
  - reorganized — tanglegram_reorganized_Pat_08_delta.png
- Vector sources (the _DRAFT set, FINALIZE for the preprint):
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_trace_Pat_06_delta_N23_DRAFT.pdf
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_reset_Pat_10_beta_N30_DRAFT.pdf (post-in-centre reorder; regenerate via --classes reset --patients Pat_10 --bands beta --topk 1)
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_anchor_Pat_08_low_gamma_N23_DRAFT.pdf
  - data/preprint/figures/_drafts/fig_cophenetic_tanglegram_reorganized_Pat_08_delta_N20_DRAFT.pdf
- (generator: scripts/01_compute/figures_embedded/fig_cophenetic_tanglegram_class.py)

6. REFERENCES
- Taxonomy / measure are ours (ρ^coph cross-phase similarity).

7. CANVA STATUS
NEW — no Canva page yet. Build right after the pipeline (ρ^coph) slide; this is the visual glance of the four categories.
