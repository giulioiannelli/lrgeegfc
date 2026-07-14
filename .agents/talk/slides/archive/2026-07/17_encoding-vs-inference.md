---
name: talk-slide-17-encoding-vs-inference
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 17
status: draft
updated: 2026-07-13
canva: pages 16–17 · ~20%
---

# Slide 17 — Encoding vs inference: the decomposition (setup)

1. TITLE
Encoding vs inference — the decomposition

2. MAIN CONCEPT
- Head: the held reorganization is not one thing — it splits into what the brain was **shown** (encoding: the adjacent pairs) and what it **inferred** (the transitive order it was never shown), and the four-phase design lets us cut it cleanly. This slide sets up the question; the next slide shows that **both persist**, in **different cortex**.
- The four phases give three difference-of-hierarchy quantities: encoding `e = D(learn) − D(pre)` (reorganization while the pairs are on screen), inference-specific `f = D(test) − D(learn)` (the extra reorganization in the test phase, the structure never shown), persistence `p = D(post) − D(pre)` (what survives into rest).
- Two co-equal persistence questions — not one signal and one nuisance: (i) does rest keep the encoding reorganization, and (ii) does rest keep the inference-specific reorganization?
- The partial correlation `T_infspec = partialcorr(f, p | e)` isolates the inference-specific persistence from the encoding one — it holds encoding constant so we are measuring the order the brain constructed, not re-counting the pairs it saw. Encoding is measured on its own terms in parallel, not thrown away.
- One-line: does rest hold the **pairs it saw**, or the **order it inferred** — and where does each live? (Next slide answers both.)

3. ON-SLIDE TEXT
encoding = pairs shown (learn − pre)
inference = order reasoned out (test − learn), encoding held constant
persistence = what lasts into rest (post − pre)
T_infspec = partialcorr(inference, persistence | encoding)
does rest keep the pairs seen — and the order inferred?

4. SPEECH
The β trace holds — but holding what? In the task the patient did two things: they were shown adjacent pairs, and from those pairs they had to work out an order they were never shown. Encoding versus inference. The four-phase design separates them. Encoding is the reorganization that appears while the pairs are on screen. Inference-specific is the extra reorganization in the test phase — the transitive structure — with encoding held constant, so we are not just re-counting the pairs the brain saw. Then, for each, we ask the same thing: does it survive into rest? These are two co-equal questions — encoding is a real result in its own right, not a nuisance we scrub away. The partial correlation is simply what earns the word "inferred": it strips out the seen relations and leaves the ones the brain had to construct. And the answer — next slide — is that rest keeps both, but not in the same place.

Careful: this is the SETUP — no verdict numbers here. Encoding is a genuine persisting result, NOT a nuisance to partial out; the partial correlation only isolates inference-specific from encoding, it does not demote encoding. Do NOT say inference is "β-only" (it is δ/α/β; and don't preview the δ cell at all). Do NOT frame the multiscale read as belonging to inference or as "multiscale robustness / a scale of integration" — inference's multiscale value is localization, not detection, and that is a next-slide point. Matched-strength is the only null. No behavioural data.

5. FIGURES
- Decomposition schematic (build in Canva): rest_pre → task_learn → task_test → rest_post, with `e` / `f` / `p` as three arrows and `T_infspec = partialcorr(f, p | e)` as the caption. Draw encoding and inference as two co-equal branches feeding persistence — NOT encoding as a subtracted term.
- Rendered decomposition forest — data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf. FLAG (verify before use): keep only if it shows encoding and inference as co-equal decomposition and carries NO verdict numbers / NO "β-only" or "inference = multiscale" framing. If it embeds any of those, rebuild as a neutral setup schematic. PNG export is fine for the talk.

6. REFERENCES
- Design is ours; cognitive-map / transitive-inference framing carried from slide 2.

7. CANVA STATUS
⚠ 2026-07-13 REFOCUSED (settled three results): recast so encoding and inference are co-equal persisting results (encoding is no longer a nuisance to partial out), added the "both persist, different cortex" preview, and removed the "inference = multiscale robustness" implication. RE-RENDER on deck.
Pages 16–17 · ~20% (two Canva pages that collapse into this one). Have: title, auto-filler bullets. Missing: the co-equal decomposition schematic, compressed on-slide text, presenter notes.
