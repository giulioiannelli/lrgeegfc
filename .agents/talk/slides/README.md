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

## Slide index (18 slides — results act RE-COMPRESSED 2026-07-13 PM to the mst@0.20 settled results)

Part I — Introduction
- [x] [01 · Title](01_title.md)
- [x] [02 · Transitive inference — the puzzle](02_ti-puzzle.md)
- [x] [03 · The sEEG dataset](03_seeg-dataset.md)
- [x] [04 · Functional connectivity](04_functional-connectivity.md)
- [x] [05 · Higher-order & the forgotten multiscale](05_higher-order-multiscale.md)
- [x] [06 · Why this matters — signatures across scales](06_why-it-matters.md)
- [x] [07 · Revealing multiscale structure with diffusion](07_diffusion-lrg.md)
- [x] [08 · The bet — form mirrors form](08_the-bet.md)

Part II — Methods
- [x] [09 · Imaginary coherence](09_imaginary-coherence.md)
- [x] [10 · Frequency bands](10_frequency-bands.md)
- [x] [11 · The network-analysis pipeline](11_pipeline.md)
- [x] [11b · The measure — ρ^coph](11b_the-measure.md)

Part III — Results
- [x] [12 · A taxonomy of neuronal populations — trace / anchor / reset / reorganized (tanglegrams)](12_tanglegram-categories.md)
- [~] [13 · A lasting, multiscale trace](13_lasting-trace.md) *(other agent — folds in selectivity + scale-shape)*
- [x] [14 · Detect ≠ discriminate](14_detect-vs-discriminate.md) *(the multiscale read resolves what edges can't: discriminate bands + scale-resolution + higher-order; θ the non-circularity control; raw = robust complementary detector)*
- [⚠] [14b · What the scale parameter means](14b_scale-meaning-localization.md) *(s coarse-grains topology AND space. ⚠ UNDER REVIEW 2026-07-14 — its payoff "coarse home → β-LEFT on slide 16" points at the now-DROPPED localization slide; the "coarse spatial home" bridge needs cutting or repurposing since anatomical localization is not presentable)*
- [x] [15 · Encoding vs inference — held offline](15_encoding-vs-inference.md) *(both persist; inference α/β not β-only, mesoscale-emergent)*
- [x] [16 · One operator, two readouts](16_one-operator-epilepsy.md) *(band dissociation + SOZ marker; β the bridge)*

Part IV — Conclusion
- [x] [17 · Take-homes](17_take-homes.md)
- [x] [18 · Outlook & thanks](18_outlook-thanks.md)

**Restructure log.** 2026-07-13 PM (mst@0.20 settled results, per
`.agents/talk/2026-07-13_results-restructure-9slot.md` + the PM checkpoint):
the post-trace act was re-compressed 9→6 result slides. Old 14 (selectivity) + 15 (scale-shape)
folded into **13** (other agent); old 16 (β held/placeless/lateralized) superseded by the new
localization slide **16** (β coarse-left, not placeless); old 17+18 merged into **15**
(encoding vs inference, dissociation now by SCALE — OFC/cingulate localization RETIRED); old 19+20
merged into **17** (epilepsy); old 21→**18**, old 22→**19**. Superseded specs live in
`archive/2026-07/`. Older I-/M-/R-/C- files retained in git history.

**2026-07-14 — localization slide DROPPED.** Slide 16 ("Where the trace lives") retired to `archive/2026-07/`: the β "coarse LEFT-hemisphere home" is incoherent (it welds the *mesoscale diffusion* scale onto *macro* anatomy — a whole hemisphere) and only marginal; there is no presentable cognitive localization (the one focal address, low-γ→SOZ, is the disease band, carried by the epilepsy slide). Slides 17→16, 18→17, 19→18 renumbered; the figure `fig_localization_spectrum` was removed. **14b** (scale→localization bridge) is UNDER REVIEW — cut or repurpose. `15`/`16` localization pointers reconciled 2026-07-14.
