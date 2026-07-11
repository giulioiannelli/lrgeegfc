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

## Slide index (22 slides — restructured 2026-07-11; 14+15+16 → two, 17+18 → one, 18+19 → one)

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
- [x] [11 · The network-analysis pipeline (ρ^coph)](11_pipeline.md)

Part III — Results
- [x] [12 · A taxonomy of neuronal populations — trace / anchor / reset / reorganized (tanglegrams)](12_tanglegram-categories.md)
- [x] [13 · A lasting, multiscale trace](13_lasting-trace.md)
- [x] [14 · Who's right? — the nulls decide](14_the-nulls.md) *(merge: old 14 who's-right + 15 the-nulls)*
- [x] [15 · Who survives — it must be multiscale](15_survivors-table.md)
- [x] [16 · β → OFC — a held consolidation](16_beta-ofc-consolidation.md) *(merge: old 17 β→OFC + 18 consolidation)*
- [x] [17 · Encoding vs inference — the decomposition](17_encoding-vs-inference.md)
- [x] [18 · Inference persists (β only), and localizes — encoding→OFC, inference→cingulate](18_inference-persists.md) *(merge: old inference-persists + localization)*
- [x] [19 · The band taxonomy — β is the flagship](19_per-band-taxonomy.md)
- [x] [20 · Coda — epilepsy marker](20_epilepsy-marker.md) *(cuttable)*

Part IV — Conclusion
- [x] [21 · Take-homes](21_take-homes.md)
- [x] [22 · Outlook & thanks](22_outlook-thanks.md)

Old I-/M-/R-/C- files superseded by this numbered set (content rewritten in plain format;
git history retains the originals). 2026-07-11: 14+15+16 → the nulls (14) + survivors table (15);
17+18 → β→OFC / held-not-replayed (16); 18+19 → inference-persists + localization (18); → 22 slides.
