---
name: writing-directive-results-section-full-draft
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: A COMPLETE scaffold ("fake") draft of the entire Results section (R1–R3, every paragraph written, with four tables and real numerics) to hand to the LaTeX writing agent as the structural model + content directive. It is a starting scaffold, NOT final prose — every \TODO must be filled from the cited CSV and every number re-verified before submission. Measure definitions belong to Methods, not Results. Changes no locked verdict.
created: 2026-06-24
companion: directives/writing_directive_2026-06-24_results-section-full-draft.tex (the draft itself — compilable LaTeX)
verified_against: data/audit/ctm_triangle/cohort_summary.csv (audit_33); data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (audit_63, the C3 referee); data/audit/grassmann_cluster_extent/cohort_summary.csv (audit_70); data/audit/spectral_distance_swap/cohort_summary.csv (audit_143); data/audit/localization_atlas/matched_strength_R1000_{include,exclude}.csv (audit_83/92); data/audit/epi_stratified/ (audit_77/85); data/audit/replay_states/ (audit_124–141); data/audit/consolidation_arc/ (audit_103/103d); data/audit/inference_localization/ (audit_110/112/113b); data/audit/epi_propagator_detector/ (audit_117)
---

# Writing directive — full Results-section scaffold draft

**Head.** The companion `.tex` is a **complete, long, drop-in draft of the whole
Results section** (R1 trace · R2 flagship/abstraction · R3 epilepsy) — every
paragraph already written in Nature-Neuroscience register, with four tables and
the real cohort numerics woven in. **Use it as the skeleton you write *into* and
expand, not as finished text.** Its job is to fix the structure, the paragraph
titles, the claims, and the evidence so the writing agent never has to re-derive
them — only polish prose, verify numbers, and lengthen.

> **This is a scaffold, not a submission.** Two hard duties before any of it
> ships: (1) **fill every `\TODO{…}`** from the CSV it names — do **not** invent;
> (2) **re-verify every number** against the `verified_against` CSVs (frontmatter).
> The draft's numbers are transcribed from the headline files / ledgers, but the
> cached CSV is the source of truth.

---

## 1. What this is, and how to use it

- The `.tex` compiles standalone (generic `article` class) so the writing agent
  can see it rendered immediately, then port it into the real journal template.
- It is the **"fake paper"** version: real structure, real claims, real numbers,
  placeholder authorship/figures. Treat it as the **directive for what the
  Results must say and in what order**, plus a first full prose pass.
- Workflow for the writing agent: (a) drop into the manuscript template;
  (b) resolve all `\TODO`; (c) re-verify numerics against the cached CSVs;
  (d) expand each paragraph to journal length (the draft is dense, not yet
  long-form); (e) produce the five figures from the figure-legends section.

## 2. Numerics status — real vs. to-fill

- **Real (transcribed, re-verify against CSV):** the held-state cohort result
  (10/10, p=0.001), the inference-specific β arc (p≈0.007), the encoding echoes
  (α +0.22 / β +0.25, p=0.014), the mesoscale nulls (τ≈2.6 p=0.014, τ≈6.8
  p=0.010), OFC localization (R=1000, BH q≈0.009–0.013; duration-survive
  q=0.050), low-γ cingulate (ρ≈0.39, 8/8, q≈0.035), the detector (AUC≈0.81,
  precision@5≈60%, 9/10), the band×probe dissociation table.
- **`\TODO` (not in scope when drafted — fill from CSV/config):** patient IDs &
  demographics, exact band edges, phase durations, the per-band statistics behind
  the verdict matrix (from `locked/VERDICT_LEDGER.md`), the detector feature list
  (from `epi_propagator_detector/README.md`), data/code availability.

## 3. Non-negotiable rules baked into the draft (keep them)

1. **Results stated in the LRG framework** (cophenetic / Grassmann); raw FC is the
   baseline, never a result. (`feedback_results_only_in_laplacian_framework`.)
2. **Matched-strength (C3) is the referee** — every cohort claim names it.
   (`feedback_matched_strength_mandatory`.)
3. **OFC is a *concentration*, not a container** — "brain-wide, yet concentrated";
   never "the trace lives in OFC". It washes out at lobe/hemisphere.
4. **Inference → cingulate is a directional hint only** (duration control demoted
   it); never an established location.
5. **No behaviour, ever** — performance data is unavailable; "abstraction" is an
   interpretation of a persistence result, never tied to inference success.
6. **β-scoping for tissue** — "core spared" is the β result; α *recruits* the core
   (R1.6 bridge). Do not generalize "healthy cortex, not the core" across bands.
7. **Held, not bursty / not "tighter"** — the trace is a sustained, non-ergodic
   state; the transient-replay reading is a clean negative; do not claim a tighter
   attractor (only a weak trend, p=0.08).
8. **Methods are not described in Results** — the framework, the two probes, the
   substrate, and the C1–C5 controls live in Methods; Results name them only as
   properties of the finding.

## 4. Structure map (paragraph titles, for quick reference)

- **R1 — A held, multiscale β-band trace of learning and reasoning**
  (R1.1 imprint rest does not erase · R1.2 written across scales, not single links ·
  R1.3 brain-wide yet concentrated in the orbitofrontal map · R1.4 β trace rides
  healthy cortex, spares the seizure core · R1.5 a state the brain dwells in, not a
  flash · R1.6 *bridge* — over the seizure core the memory and inference bands part
  ways). Table 1 = band×probe dissociation.
- **R2 — Offline abstraction of a learned structure (FLAGSHIP)**
  (R2.1 carries what was inferred, not just what was seen · R2.2 mesoscale of
  multi-step integration · R2.3 orbitofrontal anchor set by learning · R2.4 memory
  and inference ride different rhythms · R2.5 hidden low-γ cingulate memory trace ·
  R2.6 inference leans cingulate [hint] · R2.7 cingulate switches content by
  rhythm). Tables 2 (arc by band) + 3 (localization + duration).
- **R3 — Propagator-inspired markers of epileptogenic tissue**
  (R3.1 relational not contact-by-contact · R3.2 cross-band calibrated P(SOZ) ·
  R3.3 seed-free frontier · R3.4 two populations · R3.5 scope/ceiling). Table 4 =
  detector performance.

## 5. Sources of truth (read alongside)

- Outline / TOC: `headlines/RESULTS_OUTLINE.md` (same paragraph titles + flags).
- Per-headline detail + provenance: `headlines/01_trace.md`,
  `02_encoding_vs_inference.md`, `03_epileptogenic_markers.md`.
- Locked verdicts: `locked/VERDICT_LEDGER.md`, `locked/CONTROLS.md`,
  `locked/ANATOMY_LEDGER.md`. The draft must not contradict these.
- Per-band briefs: `bands/01_beta.md` … (operational source for β/α numbers).

## 6. Writing-agent checklist

- [ ] Port `.tex` into the journal template; keep paragraph titles.
- [ ] Resolve every `\TODO`; transcribe missing numerics from the cited CSVs.
- [ ] Re-verify each stated number against the `verified_against` CSVs.
- [ ] Expand each paragraph to full length (draft is dense, not long-form).
- [ ] Confirm no rule in §3 is violated anywhere in the expanded text.
- [ ] Build the five figures (legends section at the end of the `.tex`); PDF
      vector only, X-epi variants where noted, no C4 figure.
- [ ] Cross-check against `locked/` ledgers; flag any drift before submission.
