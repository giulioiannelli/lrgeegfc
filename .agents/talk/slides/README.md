---
name: talk-slides-index
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-09
updated: 2026-07-11
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
---

# Talk slides — per-slide expansions

## Head

One markdown file per slide of the 20-minute talk, expanding each row of
the [structure file](../../reports/2026-07-07_talk-structure-20min.md) into
a build-ready spec: **what goes on the Canva slide**, **what you say**, the
figure to drop in, the transition, and the guardrails. The structure file
stays the single source of truth for the arc and the numbers; these files
are the shooting script. Slides are built in **Canva** (no HTML deck).

## Per-slide file format (plain — standardized 2026-07-11)

Each slide file is a copy-paste content placeholder for Canva. Keep it PLAIN —
minimal markdown, no bold-soup, no nested bullets. It carries content only, never
aesthetics (layout, colours, fonts, positions stay the presenter's call in Canva).

Seven parts, in this order:

1. TITLE — the slide's working title (one line).
2. MAIN CONCEPT — what the slide must convey, as a short, simple bullet list.
3. ON-SLIDE TEXT — the text that actually appears on the slide. As compressed and
   minimal as possible (a few words / short phrases — not sentences).
4. SPEECH — what you say out loud on this slide (the spoken script).
5. FIGURES — the figure(s) shown, each as: what it is — file path (or "from <paper>").
6. REFERENCES — bibliographical references for the slide.
7. CANVA STATUS — current state on the deck "The Multiscale Shape of Neural
   Inference": page N · ~X% · have: … · missing: …. Always re-read the deck to score.

Completion rubric for part 7 — four equal quarters (~25% each): title placed ·
on-slide text final · speech/notes present · figure(s) placed.

One allowed tail: a single "Careful:" line right after SPEECH, only when the slide
has a real overclaim landmine (e.g. "held", not "tighter"). One line, omitted otherwise.
R-7 is the exception — it carries several locked guardrails, so its Careful is longer.

The reference exemplar of this format is [R-7](R-7_inference-persists-beta.md).

## Slide index (18 slides — 2026-07-16 reorder: diffusion moved AFTER the connectivity methods; 11b→12; compound closer at 18)

Part I — Introduction
- [x] [01 · Title](01_title.md)
- [x] [02 · Transitive inference — the puzzle](02_ti-puzzle.md)
- [x] [03 · The sEEG dataset](03_seeg-dataset.md)
- [x] [04 · Functional connectivity](04_functional-connectivity.md)
- [x] [05 · Multiscale organization — a higher-order feature](05_higher-order-multiscale.md)
- [x] [06 · Why it matters — from detecting a change to characterizing it](06_why-it-matters.md)
- [x] [07 · Form Mirrors Form](07_form-mirrors-form.md)

Part II — Methods
- [x] [08 · Imaginary coherence](08_imaginary-coherence.md)
- [x] [09 · Frequency bands](09_frequency-bands.md)
- [x] [10 · Revealing multiscale structure with diffusion](10_diffusion-lrg.md)
- [x] [11 · The network-analysis pipeline](11_pipeline.md)
- [x] [12 · The measure — ρ^coph](12_the-measure.md)

Part III — Results
- [x] [13 · A taxonomy of neuronal populations — trace / anchor / reset / reorganized (tanglegrams)](13_tanglegram-categories.md)
- [x] [14 · A lasting, multiscale trace](14_lasting-trace.md)
- [x] [15 · Detection and discrimination](15_detect-vs-discriminate.md) *(the multiscale read resolves what edges can't: discriminate bands + scale-resolution; θ the non-circularity control; raw = robust complementary detector)*
- [x] [16 · Encoding vs inference — held offline](16_encoding-vs-inference.md) *(both persist; inference α/β not β-only, mesoscale-emergent)*
- [x] [17 · From cognition to pathology (one operator, epilepsy)](17_one-operator-epilepsy.md) *(band dissociation + SOZ marker; β the bridge)*

Part IV — Conclusion
- [x] [18 · Take-homes · ongoing · outlook (compound closer)](18_closing-takeon-ongoing-outlook.md)

**2026-07-16 — deck reorder (this session, user-directed).** Two moves: (1) the **diffusion/LRG** slide moved from slot 7 to **slot 10** — *after* the connectivity methods (imaginary coherence 8, frequency bands 9); this fixes the bet's "everything on the last slide" back-reference (now correctly slide 6) and makes the diffusion slide's three-axis synthesis a clean backward close, while **Form Mirrors Form** shifts up to **7**. (2) **11b (the measure) → 12**, cascading the results act (taxonomy 13 · lasting-trace 14 · detection/discrimination 15 · encoding-vs-inference 16 · epilepsy 17) so the deck is a clean **1–18** with no sub-numbered slide; the compound closer is **18**. Files renamed; headers, frontmatter slugs, and ~40 cross-references updated across all slides. `Slides.md` is a stale scratch outline (left as-is). ⚠ The Canva deck pages still need re-ordering to match this sequence.

**Restructure log.** 2026-07-13 PM (mst@0.20 settled results, per
`.agents/talk/2026-07-13_results-restructure-9slot.md` + the PM checkpoint):
the post-trace act was re-compressed 9→6 result slides. Old 14 (selectivity) + 15 (scale-shape)
folded into **13** (other agent); old 16 (β held/placeless/lateralized) superseded by the new
localization slide **16** (β coarse-left, not placeless); old 17+18 merged into **15**
(encoding vs inference, dissociation now by SCALE — OFC/cingulate localization RETIRED); old 19+20
merged into **17** (epilepsy); old 21→**18**, old 22→**19**. Superseded specs live in
`archive/2026-07/`. Older I-/M-/R-/C- files retained in git history.

**2026-07-14 — localization slide DROPPED.** Slide 16 ("Where the trace lives") retired to `archive/2026-07/`: the β "coarse LEFT-hemisphere home" is incoherent (it welds the *mesoscale diffusion* scale onto *macro* anatomy — a whole hemisphere) and only marginal; there is no presentable cognitive localization (the one focal address, low-γ→SOZ, is the disease band, carried by the epilepsy slide). Slides 17→16, 18→17, 19→18 renumbered; the figure `fig_localization_spectrum` was removed. **14b** was FOLDED into slide 07 (2026-07-14) and retired to `archive/2026-07/`: a standalone scale-concept slide was under-used once localization dropped, so its micro→meso→macro coarse-graining figure is now slide 07's main visual (the reach curve → 07 Q&A backup). `15`/`16` localization pointers reconciled 2026-07-14.
