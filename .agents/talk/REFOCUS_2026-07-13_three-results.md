---
name: talk-refocus-three-results
type: refocus-plan
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
supersedes: results-act rows (slides 13–20) of REVIEW_2026-07-13_tau-sweep-propagation.md
source_of_truth: .agents/preprint/established_results/2026-07-13_settled-three-results.md
pointers:
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md
  - .agents/talk/REVIEW_2026-07-13_tau-sweep-propagation.md
---

# Talk refocus — the three settled results

**Head.** The results have settled, and they move the talk off two things it currently
leans on: **(a) "β trace → OFC"** (now DELOCALIZED — subpoint 1.8 ✅ supersedes it as a
degenerate-graph artifact; OFC survives only as the *encoding* anchor) and **(b) the
"it must be multiscale / who's right" exclusive-detection spine** (retired: raw FC *and*
Grassmann both detect β; the multiscale edge is **information-depth**, not exclusivity).
This plan refocuses **slides 8, 13–21** around the settled three results and the
information-depth frame, keeping slide numbers/slots stable. Numbering/drift housekeeping
on 5/6/7/9/10/12 + hubs stays in the propagation-review punch-list.

## The spine (say this once; every slide inherits it)

**Thesis — information-DEPTH, not exclusive detection.** Simple pairwise FC gives a
**scalar verdict per band** ("did the reorganization lean back toward the task?") and
several bands say yes — that scalar is the ceiling of its information. The multiscale
diffusion read gives a **structured object** (a hierarchy per phase; a ρ_sym(τ) curve per
band) that characterizes what the scalar can't pose:
- **WHICH bands — selectivity.** Raw FC fires 4 bands; cophenetic fires only **α/β**;
  spectral resistance is dead in all. (Selectivity, *not* "only multiscale sees it".)
- **AT WHAT SCALE — scale-shape.** ρ_sym(τ) grades each firing band: **β scale-invariant
  (genuinely multiscale)**, **α single mesoscale peak (single-scale)**, δ fires-but-shallow.
- **WHERE — localization (WEAK on mst@0.20; de-emphasize).** β trace is delocalized;
  enc→OFC / inf→cingulate are whole-graph trends, re-derivation pending.
- **THROUGH WHICH PATHWAYS — higher-order.** The propagator sums walks of all lengths.

**The three results.**
1. Reasoning leaves a **held, genuinely multiscale trace** in β (scale-invariant) + a
   single-scale one in α. β is delocalized, held (not replayed); inter-patient spread ∝ implant laterality.
2. What persists = **encoding** (pairs shown) + **inference** (order inferred), NOT β-only;
   they sit in different cortex (OFC vs cingulate).
3. The **same operator** localizes epileptogenic tissue (δ/low-γ/β, strongest δ).

**Band-dissociation device (the unifying table, R2→R3 bridge + a take-home):**
cognition = **α/β** · epilepsy = **δ/low-γ/β** · **β = the bridge**.

**Guardrails (❌ — violating these re-creates the contradictions):**
- ❌ "only LRG/multiscale sees the trace" / "invisible to simple methods" / "multiscale-exclusive".
- ❌ "β trace → OFC" (delocalized; OFC = ENCODING anchor only).
- ❌ "inference β-only" (it is δ/α/β); ❌ the δ inference-specific effect (partial-corr artifact — cut).
- ❌ drift / "second null" / "survives both nulls"; ❌ ranking bands by ×-null ratio.
- Null = **matched-strength only**. Read **per-scale**, never best-scale.

Status tags from the settled doc: ✅ re-verified on mst@0.20 · ⚠️ whole-graph, re-derivation
pending (present honestly) · ❌ do not claim.

## Canonical numbers (source of truth; every fixer uses these)
- **Trace:** β 16/16 scales p .014→.001 (scale-invariant); α 12/16 p .024→.007 (mesoscale, peaks s≈5);
  θ 0/16, low-γ 0/16 (null); δ 8/16 patchy & high-γ 4/16 coarse (artifacts, not claimed).
- **β trace DELOCALIZED** (1.8 ✅): nothing clears BH at s=1/2.83/5.65. **Held, not replayed:** reinstatement 10/10 p=.001 (LOO .002) ⚠️.
- **Laterality (1.12 ✅):** β trace ∝ left-contact fraction ρ=+0.685 p=.029; α not lateralized (ρ=0.16 ns).
- **Encoding:** β 16/16 (meso .003); α 3/16 (τmin .042, meso .032; raw-visible: raw .032, clustering .010).
- **Inference (NOT β-only):** β 7/16 (.005), α 7/16 (.024); δ inference = CUT ❌. Value = localization not detection (raw detects β inf .010, placeless).
- **Localization (⚠️ whole-graph):** encoding→OFC q .010/.040; inference→cingulate q .030/.040; low-γ encoding→cingulate ρ_sym +0.25 q .035 (8/8).
- **Epilepsy:** AUC δ **.83** / low-γ **.82** / β **.745** / α .61; beats fake-SOZ null δ7/10, low-γ8/10, β8/10;
  β spares SOZ / α recruits (α SOZ–SOZ +0.41 p=.005); τ = ranking not precision; fused prec@5 = 60%;
  two populations (8 community AUC↑.99 + 2 right-hemi hub .49/.57); LOPO median AUC .87 (9/10>chance); robust dense→mst@0.20 (.80→.83, .74→.82, .69→.745).

---

## Deck re-map (slots stable; job of each slide changes)

**R1 — the held multiscale trace = 13·14·15·16** (was: trace / who's-right / ladder / β→OFC):
- **13 — the trace holds + the frame.** Kill "read three ways / who's right". Beats: (1) the
  hierarchy holds — 10/10 rest_post closer to task than rest_pre (1.1); (2) the reframe —
  raw FC = a scalar per band ("several leaned back"); the multiscale read = a structured
  object that asks *which bands, at what scale*. Grassmann = an honest whole-graph companion
  (β only; do NOT attribute to mst@0.20). Sets up 14 (selectivity) + 15 (scale-shape). Soften
  "drift" wording. FIG: reinstatement (trace holds) + a "scalar vs structured object" concept;
  RETIRE / relabel fig_before_nulls (drop the 3-co-equal-reads synthesis; fix the footer).
- **14 — selectivity (folds in the controls ladder).** Matched-strength = the sole null
  (drift retired). SELECTIVITY: raw FC fires 4 bands (δ/α/β/low-γ) → non-selective; cophenetic
  fires only **α/β**; spectral resistance dead in all (1.6). Honesty: raw *detects* β (p=.024) —
  the edge is **selectivity, not exclusivity**. This is old-15's controls-ladder REFRAMED as
  selectivity (not "it must be multiscale"). FIG: the controls-ladder / selectivity matrix.
- **15 — scale-shape (NEW job; τ-curves).** The sweep grades each firing band: β scale-invariant
  (16/16, .014→.001, *deep*), α single mesoscale peak (12/16, .024→.007, *single-scale*), δ
  fires-but-shallow (not claimed). "Only β is genuinely multiscale." Only the sweep shows this
  (1.7). FIG: ρ_sym(τ) curves per band (β flat, α peaked) — the multiscale payoff figure.
- **16 — what kind of trace: held, placeless, lateralized.** KILL "β→OFC". Beats: (1) HELD not
  replayed — reinstatement 10/10 p=.001 (1.10 ⚠️); (2) DELOCALIZED — no single address at any
  scale (1.8 ✅; be honest it overturns the prior OFC claim); (3) LATERALIZED carriers — spread
  ∝ left-fraction ρ=.685 p=.029 (1.12 ✅), resetters right-only, α not lateralized. FIG:
  reinstatement + laterality brain (NOT the OFC-localization fig).

**R2 — encoding vs inference = 17·18** (keep slots):
- **17 — the decomposition (setup).** Both encoding & inference persist (don't cast encoding as
  a nuisance to control out). Question: seen pairs vs inferred order — and they'll differ in
  *where*. No numbers. Light edit.
- **18 — both persist, different cortex.** Fix per settled doc: encoding β 16/16 (meso .003), α
  3/16; inference β 7/16 (.005) + α 7/16 (.024) → **NOT β-only**; CUT δ inference (2.4 ❌).
  Encoding→OFC, inference→cingulate (⚠️ whole-graph, present honestly; localization = the
  value-add, detection tied). low-γ encoding→cingulate (2.8 ⚠️). Not a duration effect (2.9).
  Scrub the +0.091/p=0.0098 scalar → per-scale; scrub the drift analogy; retitle off "β only".

**R3 — the shared operator = 19·20** (keep slots):
- **19 — the band dissociation (bridge).** Repurpose the taxonomy slide into the settled doc's
  selling point: cognition α/β · epilepsy δ/low-γ/β · **β the bridge**. α → mesoscale 12/16 .007;
  θ-silence ⇒ selectivity is real; keep low-γ→cingulate encoding (⚠️). One-glance band table;
  the transition from cognition into the clinical coda.
- **20 — the epileptogenic marker.** Update AUCs → δ .83 / low-γ .82 / β .745. Add: strength-
  independent community; beats fake-SOZ null (δ7/low-γ8/β8 of 10); β spares / α recruits; τ =
  ranking not precision; precision from band fusion (60% prec@5); two populations (community
  AUC↑.99 + 2 right-hemi hubs); LOPO median AUC .87; robust dense→mst@0.20; scope = labels not
  outcome, marker not detector.

**Frame + close:**
- **8 — the bet (light reframe).** Drop "invisible to single-edge metrics / spectral clustering"
  (exclusive-detection overclaim — raw detects β). The bet becomes: the trace is *characterized*
  by multiscale structure (its scale-shape), a structured object beyond the scalar. Keep it a
  hypothesis the nulls test (matched-strength only).
- **21 — take-homes (align).** Ensure no β→OFC and no exclusive-detection language; carry the
  information-depth frame + the band-dissociation take-home + honest inference (δ/α/β, localization).

## What is NOT in this refocus (stays in the propagation-review punch-list)
5, 6, 7, 9, 10, 12 (numbering/drift/exemplar housekeeping) + Slides.md/README/PLAN hub re-sync.

## Execution
Parallel rewrite of slots 8, 13–21 (each agent: read settled doc + this plan + its current
slide → rewrite the .md in the 7-part format to the settled framing). Then a consistency pass
(shared vocabulary; OFC only on 18; band-dissociation only on 19; no exclusive-detection anywhere).
