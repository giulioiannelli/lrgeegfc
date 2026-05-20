---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: pointer to the parallel epilepsy track (Directions A–D); explicitly NOT folded into notes_imcoh.tex
---

# Epilepsy track — pointer

**Head.** The user is running a parallel research track on the epileptic dynamics of the same cohort. The track is **not folded into `notes_imcoh.tex`** — the manuscript is about the cognitive-task trace pattern (TARR taxonomy), and the epilepsy track investigates graph-theoretic primitives at epi-vs-non-epi nodes without invoking the trace frame. This file states the connection points where the manuscript could optionally cite the parallel track in a future revision, and explicitly rules out folding the epi work into the §5.5 Bonferroni or §5.6 taxonomy machinery.

## Parallel track artefacts

### Plan

- `.agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md` — 4-direction scope plan (A: cross-phase rigidity; B: virtual resection in `ρ̂(τ)`; C: eigenmode IPR; D: `ρ̂`-leakage propagation-zone candidates).

### Scope reports

- `.agents/guides/task-persistence-investigation/2026-05-08_epi-cross-phase-rigidity.md` — Direction A (cross-phase rigidity of the epi sub-network: KC distances on induced subtrees at λ ∈ {0, 0.5, 1}; size-matched random-leaf null R=100; cross-probe variant).
- `.agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md` — Direction C (mode projection mass `m^E_k = Σ_{i ∈ E_p} v_k(i)²` binned on normalized eigenvalue; degree-corrected null via strength-stratified shuffling).

### First-pass results

- `.agents/reports/2026-05-08_direction-a-null-induced-subtree.md` — Direction A first-pass: δ-heights RESET 75% IQR, β-heights RESET 43% IQR (dissociation from §5 β trace), α-topology TRACE 55–58% IQR. **Phenomenology, not trace verdict.** The framing is "structured per-band/per-λ regularities" (per `feedback_regularity_over_bh_null.md`), not BH-FDR rejection.

### Memory pointers

- `feedback_epilepsy_not_trace_locked.md` — epilepsy investigations are NOT trace-locked; map graph-theoretic primitives, don't fold into TARR.
- `epilepsy_lrg_research_directions.md` — 4-direction summary (current).
- `epileptic_imcoh_universal.md` — n=10 verdict (2026-05-07): δ cross-probe 1.55× is known-biology *confirmation* not discovery; β V1 retired.
- Canonical writeup: `.agents/reports/2026-05-07_epileptic-n10-revisit.md`.

## Why the epilepsy track is parallel and not folded

1. **Different scientific question.** The manuscript investigates whether cognitive-task reorganization leaves a structural memory in `rsPost`. The epilepsy track investigates whether graph-theoretic primitives at epi-vs-non-epi nodes differ in ways the LRG framework can detect (cross-phase rigidity, eigenmode localization, ρ̂-leakage to propagation-zone candidates).
2. **Different baseline / null.** The §5 cohort claim is gated by within-rsPre split-half nulls. The epilepsy track uses size-matched random-leaf-subset nulls — incompatible with the trace-direction frame.
3. **Different cell unit.** §5 cells are (probe, band) for the LRG-probe family or (band, region) for the anatomy family. Epilepsy track cells are (band, λ_KC) on the **induced subtree** restricted to epi leaves `E_p`, with a size-matched null built from random leaf subsets.
4. **Trace-locked frame would be wrong.** Per `feedback_epilepsy_not_trace_locked.md`: the epilepsy investigation is a **structural** investigation of how the LRG framework reads the epi sub-network, not a hypothesis test about whether epi nodes participate in the cognitive trace. The TARR taxonomy applied to epi leaves would silently impose the cognitive-task interpretation onto a baseline graph-property question.

## Connection points where the manuscript could OPTIONALLY cite the parallel track

These are forward-look mentions, not folding the analysis into §5. Each is a one-line citation in §5.5 / §5.6 / §6.3 if the parallel track produces a publishable result before the manuscript is finalized.

### §5.5 (Hippocampus + fusiform regions involve epi contacts)

The §5.5 anatomical regions where the trace-leaf cell counts are highest (β Hippocampus 5/27, β fusiform 5/38, γ_l fusiform 13/38) overlap with the typical clinical sEEG sampling for epilepsy localization. **Optional citation**: a one-paragraph note in §5.5 stating "the contact fraction labeled epileptic in these regions is X%" once item #10 of `15_new_analyses_owed.md` is computed. This is a **descriptive disclosure**, not a hypothesis test:

> "Across the cohort, the contact pool at γ_l ctx-lh-fusiform contains X% of contacts labeled clinically epileptogenic; the contact pool at β Hippocampus contains Y%; β fusiform Z%. The trace-leaf rate within and outside these epi contacts (across the cortical contact pool) is reported in Supplementary Table N. The cohort claim of §5.3 is computed on the full contact pool and does not condition on epi labeling; the disclosure here is for downstream interpretive context only."

### §5.6 anchor class (epi-driven anatomical anchors)

Epi contacts that recur across phases as same-probe anchor clades (per item #1 of `15_new_analyses_owed.md`) overlap with the surgical-target literature on focal epileptogenic zones. **Optional citation**: a one-line note in §5.6 anchor-class prose stating "anchor clades restricted to epi-contact leaf sets are reported separately in the parallel epilepsy track (`.agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md`, Direction A)". This citation does **not** import the TARR frame; it just points to the parallel investigation.

### §6.3 outlook (forward direction)

The §6.3 outlook can mention the parallel epilepsy track as a forward direction, framed as "graph-theoretic primitives at epi-vs-non-epi nodes" not as "do epi nodes trace?":

> "A parallel investigation reads the same |ImCoh|-LRG cohort under graph-theoretic primitives restricted to clinically-labeled epileptogenic contact subsets — eigenmode localization on the leading Laplacian eigenmodes, cross-phase rigidity of the induced subtree distance, and propagation-zone candidacy via ρ̂-leakage. These probes ask whether the LRG framework resolves epi-specific structure independent of the cognitive trace pattern reported here. They are scoped under `.agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md` and reported separately."

This is a forward-look-and-defer, not a result.

## What the verification package does NOT propose

### NOT — hypergeometric epi × β-trace test

A hypergeometric test of "are β trace-leaves over-represented at epi-labeled contacts?" would import the trace frame onto the epi-vs-non-epi contrast. This is exactly what `feedback_epilepsy_not_trace_locked.md` rules out. **The verification package does NOT propose this test.**

If the parallel track wants to address the question "do epi nodes participate in the cognitive trace?", that is a separate scoping exercise and would need its own scope report under `.agents/guides/task-persistence-investigation/` framed without the TARR taxonomy.

### NOT — folding epi metadata into §5.5 anatomy Bonferroni

The §5.5 family `m=48^{(anat)}` is over (band, region) cells where regions are Desikan-Killiany cortical parcels. Adding an epi-vs-non-epi axis would expand this family, change the multiple-comparison structure, and conflate the trace anatomy question with the epi-vs-non-epi question. **The verification package does NOT propose this expansion.**

### NOT — folding epi metadata into §5.6 taxonomy class assignment

The §5.6 taxonomy assigns each leaf to one of {trace, reset, rearrange, anchor, diffuse} based on its dendrogram-clade behavior. Adding an epi-vs-non-epi covariate would either (a) require a 2 × 5 product taxonomy (10 classes) or (b) require stratified hypothesis tests per epi label. **The verification package does NOT propose this stratification** — the taxonomy is structural-descriptive and is not designed to support epi-stratified inference.

## Action — next agent / next session

1. **Read this file before referencing epi material in the manuscript.** The parallel track is forward-look context, not a manuscript section.
2. **Item #10 in `15_new_analyses_owed.md`** (epi contact fraction disclosure) is the only direct connection that could land in `notes_imcoh.tex` in the next revision. It is a descriptive disclosure, not a hypothesis test.
3. **Scope reports for Directions A and C** (`2026-05-08_epi-cross-phase-rigidity.md`, `2026-05-08_epi-eigenmode-localization.md`) are the canonical home for the parallel track's analyses. New audit scripts (`audit_54_*`, `audit_55_*`) implement them. Their outputs go to `data/audit/epi_cross_phase_rigidity/` and `data/audit/epi_eigenmode_localization/` — **not** to `data/audit/section5_v2_round3_redo/` (which is the §5 manuscript directory).
4. **The acceptance gate for the parallel track** (per the plan): at least one of the four directions must produce a signal not derivable from raw |ImCoh| cross-probe enrichment. If all four are null after this control, the conclusion is that LRG adds visual / descriptive value but no inferential content beyond raw FC for the epi question. If one or more passes, the parallel track yields a separate publication, not a section in `notes_imcoh.tex`.

## Final framing for any agent

The manuscript is about the **cognitive trace** in `rsPost`. The parallel track is about **epileptic dynamics** in the LRG framework. The two share the cohort and the FC method but **not** the scientific question. Do not fold them. Do not import TARR into the epi track. Do not import epi-stratified tests into the manuscript without scoping them as descriptive disclosures.
