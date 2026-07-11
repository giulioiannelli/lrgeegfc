---
name: talk-slide-13-lasting-trace
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 13
status: draft
updated: 2026-07-11
canva: page 13
---

# Slide 13 — A lasting trace, read three ways

TITLE
: A lasting trace — read three ways

ONE-LINE
: The task leaves a trace that outlasts it — but the three ways of reading the connectivity give three different answers, agreeing only on β. Who's right? (the null decides, next slide.)

---

## 3. ON-SLIDE TEXT  (each bullet = one speech beat)

- **the trace holds** — post-task rest keeps the task's hierarchy (10/10 patients closer to task than to pre-task rest)
- **raw edges** — reorganize at task, lean back at rest — but the *same blunt drift in every band*
- **two multiscale lenses separate the bands — and disagree** — hierarchy: β, α · subspace: β, γ_low, δ
- **all three agree on one thing** — β is the loudest reorganization by every read
- **so who's right?** → the null decides (next)

---

## 4. SPEECH  (~75 s)

Did the task leave a trace? Read through the hierarchy, yes — the resting brain afterwards doesn't go back to where it started. Its tree still looks like the task's, and across all ten patients post-task rest sits closer to the task state than pre-task rest did.

But whether you *see* that trace depends on the lens. Look at the raw connectivity — the edges themselves — and the network clearly reorganizes during the task and leans back at rest. The trouble is that this blunt drift looks the same in every band; the raw edges can't tell a band that carries the trace from one that doesn't.

So we read the same connectivity through two multiscale lenses — the cophenetic hierarchy and the Grassmann subspace — and now the bands separate. But the two lenses don't agree: the hierarchy points to beta and alpha, the subspace to beta, low-gamma and delta.

Put all three side by side and there's exactly one thing they agree on — beta is the loudest reorganization by every measure. Everywhere else they disagree. So before we run a single control: who's right? That's what the null is for.

---

## 5. FIGURES  (all built — PNGs in data/outputs/figures/talk/)

Reveal in speech order; suggested build = trace → raw → subspace → synthesis.

| # | role | file |
|---|------|------|
| 1 | **the trace holds** (cophenetic reinstatement: 3 phase-trees + 10/10 strip) | `trace_reinstatement_coph.png` |
| 2 | **raw edges** read (3 FC matrices: rest_pre · task_test · rest_post, Pat_08 β) | `fig_raw_fc_three_phase_Pat_08_beta.png` |
| 3 | **Grassmann subspace** read (per-band subspace trace-mass heatmaps) | `fig_grassmann_heatmap_bands.png` |
| 4 | **synthesis — who's right** (per band × per measure, before nulls; β loudest under all three, reads disagree) | `fig_before_nulls_three_measure.png` |

Generators (repo): `scripts/01_compute/figures_embedded/fig_talk_raw_fc_three_phase.py`, `fig_talk_before_nulls_three_measure.py`; reinstatement PDF `data/preprint/figures/results_section1/fig_trace_e_reinstatement.pdf`; Grassmann PDF `data/preprint/figures/all_bands/fig_bands_grassmann_heatmap.pdf`.

Number to verify on the deck: the reinstatement panel prints **10/10, p = 0.003**; the finalized compound `fig_trace1.pdf` reports **p = 0.001** (LOO 0.002) — use the finalized number, or say "10/10" without the exact p.

Fig 4 honesty note: bars = observed cohort effect (native value annotated), *before* the matched-strength / drift nulls that adjudicate on the next slide. The ✓ is each read's own basic cohort test, and it is deliberately not uniform — raw & coph: one-sided Wilcoxon of per-patient effect vs 0 (no surrogate); Grassmann: observed contiguous cohort-significant run > null p95. Note raw β is loudest by median yet fails its plain Wilcoxon (p = 0.053, no ✓) — the per-patient spread swamps it, which is exactly why the raw edges can't certify a band and the null is needed.

---

## 6. REFERENCES
- All three read-outs are ours (raw |ImCoh| edges · Grassmann chordal distance · cophenetic ρ^coph, both LRG distances from `e^{−τL̂}`). No external citation on this slide.

## 7. CANVA STATUS
Page 13 · figures done. Rewritten M 07-11: numbering fixed (was "Slide 11"); slide now ends on "who's right?" (folds in the old slide-14 "who's right" beat — see below); all four figures built and PNG'd. On the deck: place the four figures in build order, paste the 5 on-slide bullets, add the speech as presenter notes.

NOTE (structure): this slide now carries what slide 14 ("Who's right?") used to. Recommend deleting/merging slide 14, or repurposing it as the null-setup transition. The old slide-14 framing ("the three measures can't all be right") was also wrong for coph-vs-Grassmann (per `results_sec_1.tex` they are *complementary*, sort by kind); the genuine tension is **raw fires everywhere vs the multiscale reads carve bands**, which this slide now states correctly and the null slide resolves.
