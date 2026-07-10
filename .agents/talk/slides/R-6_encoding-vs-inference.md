---
name: talk-slide-R6-encoding-vs-inference
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: R-6
part: III — Results
duration: ~75 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - .agents/preprint/headlines/02_encoding_vs_inference.md
---

# R-6 · The flagship (N2) setup — encoding vs inference

## Slide placeholder (copy into Canva)

**Main point.** The four-phase design lets us cut the trace into two pieces: what the
brain was **shown** (encoding — the adjacent pairs, `task_learn`) and what it had to
**reason out** (inference-specific — the transitive structure in `task_test`, with
encoding *controlled out*). Then one question: **which of the two persists into rest?**
This is the setup for the climax.

**Concepts to land.**
- **Three phase-contrasts, one clean decomposition:**
  - **e = D(learn) − D(pre)** — *encoding* reorganization (seeing the pairs)
  - **f = D(test) − D(learn)** — *inference-specific* reorganization (the reasoning,
    beyond what was seen)
  - **p = D(post) − D(pre)** — *persistence* into rest
  - The flagship quantity is **T_infspec = partialcorr(f, p | e)** — how much the
    *inference-specific* change predicts what persists, **with encoding partialled
    out**. That "with encoding out" is the whole point: it isolates the *inferred*
    relations from the *seen* ones.
- **Why the multiscale read is the right instrument for it.** Relational chaining —
  A > B, B > C ⇒ A > D — **integrates multi-step paths**; that is *motivation* for
  looking across scales (the diffusion propagator sums walks of every length), not a
  claim about a particular scale.
- **The question, stated:** does the brain merely hold the **pairs it saw**, or the
  **order it inferred**? R-7 answers.

**Figures / visuals.**
- **Decomposition schematic** (build in Canva) — the 4-phase timeline
  (rest_pre → task_learn → task_test → rest_post) with **e / f / p** drawn as the
  three arrows, and T_infspec = partialcorr(f, p | e) as the one-line caption.
- **Rendered decomposition forest** (the real per-patient e/f/p, if you want data on
  the setup slide):
  `data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf`
- *(Context, optional) encoding already anchors in OFC:*
  `data/preprint/figures/results_section2/fig_arc_c_encoding_ofc_anchor.pdf`

**References.** None new — the design is ours; the cognitive-map / inference framing
was seeded in I-2 (Behrens 2018; Wilson 2014 / Schuck 2016 for OFC).

---

## Keep honest (content constraints, not styling)

- **This slide is the SETUP, not the result.** Don't preview R-7's numbers — land the
  *decomposition* and the *question* (seen vs inferred), then pay it off in R-7.
- **The partial correlation is the control, not a flourish.** Inference-specific =
  `task_test` reorganization **with encoding partialled out** — that is exactly what
  licenses the word "inferred." State it plainly.
- **Multi-step integration is MOTIVATION, not a located scale.** The N2 result is
  **multiscale robustness** (clears the null at fine *and* meso scales), **not** a
  peak at a "scale of integration" — that peak claim was dropped 2026-07-07. Don't
  reintroduce it here.
- **No behavioural data** — "encoding vs inference" is a *design* contrast on the
  connectivity, not a brain–behaviour split. Say it once.
