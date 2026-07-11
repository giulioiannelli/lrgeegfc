---
name: talk-slide-17-encoding-vs-inference
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 17
status: draft
updated: 2026-07-11
canva: pages 16–17 · ~20%
---

# Slide 17 — Encoding vs inference: the decomposition (setup)

1. TITLE
Encoding vs inference — the decomposition

2. MAIN CONCEPT
- The four phases let us cut the trace in two: what the brain was shown (encoding) vs what it had to reason out (inference).
- encoding e = D(learn) − D(pre); inference-specific f = D(test) − D(learn); persistence p = D(post) − D(pre).
- Flagship quantity: T_infspec = partialcorr(f, p | e) — how much the inferred change predicts what persists, with encoding partialled out.
- Sets up one question: does rest hold the pairs it saw, or the order it inferred? (Next slide answers.)

3. ON-SLIDE TEXT
encoding = pairs seen (learn − pre)
inference = order reasoned out (test − learn), encoding controlled out
persistence = what lasts into rest (post − pre)
T_infspec = partialcorr(inference, persistence | encoding)
seen pairs, or inferred order?

4. SPEECH
The β trace holds — but does the brain keep the pairs it was shown, or the order it worked out? The four-phase design separates them. Encoding is the reorganization while it sees the adjacent pairs; inference-specific is the extra reorganization in the test phase — the transitive structure it was never shown — with encoding statistically controlled out. Persistence is what survives into rest. The one quantity we care about is the partial correlation between the inference change and what persists, holding encoding constant. That "holding encoding constant" is the whole point: it isolates the inferred relations from the seen ones.

Careful: this is the SETUP — don't preview the numbers. The partial correlation is the control that licenses the word "inferred". Relational chaining motivates the multiscale read; the result is multiscale ROBUSTNESS, not a peak at a "scale of integration". No behavioural data.

5. FIGURES
- Decomposition schematic (build in Canva): rest_pre → task_learn → task_test → rest_post, with e / f / p as three arrows and T_infspec = partialcorr(f, p | e) as the caption.
- Rendered decomposition forest — data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf

6. REFERENCES
- Design is ours; cognitive-map framing from slide 2.

7. CANVA STATUS
Pages 16–17 · ~20% (two Canva pages that collapse into this one). Have: title, auto-filler bullets. Missing: the decomposition schematic, compressed on-slide text, presenter notes.
