---
name: results-draft-review-handoff
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-04
pointers:
  - .agents/preprint/directives/writing_directive_2026-06-24_results-section-full-draft.tex   # THE artifact
  - .agents/preprint/directives/writing_directive_2026-06-24_results-section-full-draft.md     # its directive/rules
  - .agents/preprint/headlines/RESULTS_OUTLINE.md      # TOC (R1 biology-led)
  - .agents/preprint/headlines/01_trace.md             # N1 home (biology-led)
  - .agents/reports/2026-06-25_per-band-consistency-taxonomy-lock.md   # taxonomy LOCK
  - .agents/reports/2026-06-25_multiphase-snr-reliability.md           # Q1 (residual biology)
  - .agents/reports/2026-06-25_raw-vs-multiscale-trace.md              # Q2 (multiscale claims)
  - data/audit/localization_atlas/carrier_loo.csv                      # γ_low LOO
---

> **Head — resume token: "review the Results `.tex` paragraph by paragraph and fix."**
> The full long Results draft is
> `.agents/preprint/directives/writing_directive_2026-06-24_results-section-full-draft.tex`
> (R1–R3, 4 tables, real numerics + provenance comments). **R1 was rebuilt
> biology-led this session** (cophenetic = primary probe; Grassmann demoted to a
> β *confirmation*; per-band probe dissociation → Methods). R2 (flagship
> encoding/inference dissociation) and R3 (epi markers) are the earlier detailed
> pass and are **consistent** with the new R1. **Next: go paragraph by paragraph,
> PI reviews, agent fixes.** Honor the guardrails in §2; use the numbers in §3
> (do not re-derive); the known fixes are §4.

## 1. The artifact + workflow

- Edit the **`.tex`** in place, paragraph by paragraph (R1.1 → … → R3.5). Keep the
  provenance `% source:` comments. Register = Nature Neuroscience; clean-forward
  prose (no "previously/now" diary language — it's an article, not a changelog).
- The **`.md` companion** is the directive: it lists the 8 baked-in rules, the
  `\TODO` inventory, and the re-verify-before-submission duty. Read it once.
- `RESULTS_OUTLINE.md` (TOC) and `01_trace.md` (N1 home) already carry the
  biology-led R1; keep the `.tex` in sync with them.

## 2. Framing decisions to HONOR (locked this session)

- **Biology-led, method-as-lens.** R1.2 = the OFC/PFC cognitive-map consolidation
  (the heart). Do **not** headline the method or the probe dissociation.
- **Cophenetic is the primary probe. Grassmann is a β *confirmation* only** — the
  per-band readout × band dissociation lives in **Methods** (a table), not as an
  R1 headline. Rationale: the dissociation reads to a referee as "your two methods
  disagree"; the *convergence at β* is the sellable part.
- **OFC = concentration / over-expression above baseline, NOT "the trace lives in
  OFC"** (it washes out at lobe/hemisphere; within-OFC lean ≈0.59; hotspot not
  container). Same for γ_low→PFC (corroborating, a rung below β).
- **No behaviour, ever.** Performance data is unavailable and won't come.
  "Abstraction/consolidation" is an interpretation of a *persistence* result,
  never tied to inference success.
- **inference→cingulate = directional hint only** (fails the duration control).
  encoding→OFC is the solid localization.
- **Cohort-level over a *variable* cohort** — never "every patient shows it." The
  per-patient spread is real biology (Q1), not a reliability artifact; the SNR/
  reliability analysis is a **control that lives in Methods** + one body clause.
- **Results in the LRG (cophenetic/Grassmann) framework; raw FC = baseline only.**
- **Matched-strength (C3) is the mandatory referee** — every cohort claim names it.
- **Q2 FORBIDS:** "invisible to edge-by-edge comparison", "stronger/better-resolved
  in the hierarchy than raw" (raw magnitude is *larger*; ultrametric compresses),
  "α is multiscale-only" (α clears in raw too). ALLOWED: the band taxonomy is a
  hierarchy property (only cophenetic expresses absent-θ / present-β); at β the
  hierarchy clears matched-strength (p=.005) where raw falls just short (p=.053).

## 3. Key numbers (do NOT re-derive)

**Per-band trace (cophenetic gate · Grassmann · localization · tier):**

| band | coph gate (n>p95 / med ρ / p) | Grassmann | localization | tier |
|---|---|---|---|---|
| **β** | 7/10 / +0.22 / **.005** | **.005** (run 29) | **OFC**, LOO 10/10 top, q≤.05; sensorimotor depleted | consistent — **flagship** |
| **γ_low** | 5/10 / +0.083 / .116 | **.005** | **PFC**, LOO sig 10/10, top **8/10** (Pat_08/14 dethrone) | consistent — 2nd, a rung below β |
| **α** | 5/10 / +0.105 / **.002** | .35 (null) | none — scattered/patient-specific | patient-specific |
| **δ** | 4/10 / +0.008 / .28 | **.005** | none | patient-specific (Grassmann-only cell) |
| **γ_high** | 4/10 / ~0 / .25 | .06 (null) | parietal (net-null) | not-traced |
| **θ** | 2/10 / −0.04 / .72 | .16 (null) | MTL (net-null) | **absent** |

- **Probe dissociation:** α cophenetic-only, δ Grassmann-only, **β both = the
  convergence/cross-check**, θ/γ_high neither. Convergence at β is meaningful
  *because* the probes otherwise dissociate.
- **Grassmann baseline caveat (PARKED):** Grassmann = triangle
  `T_G = d(pre_A,task) − d(task,post)` (audit_66/70), single rest_pre; cophenetic =
  split-baseline `ρ_split` (audit_63). **Not baseline-matched** → the δ
  Grassmann-only cell is the exposed one. β convergence is safe (Grassmann run-29 +
  β passes the rigorous cophenetic independently). Fix = split-baseline Grassmann
  (compute follow-up) OR Methods caveat. Not on R1's critical path.
- **Q1 (multiphase SNR, audit_145):** heterogeneity is **NOT** pure reliability —
  a residual biological who-traces axis survives (R²=.064; who-traces ⊥ reliability;
  **Pat_15 = most reliable, lowest tracer**). So SNR is a control, not a result.
- **R2 (arc):** T_infspec·e **β p≈.007** (β-only, LOO-robust, matched-strength);
  T_learn **α +0.22 / β +0.25, p=.014** (learning leaves its own trace, β→OFC same
  hotspot q=.010); mesoscale τ≈2.6 p=.014, τ≈6.8 p=.010; encoding→OFC q=.010;
  **inference→cingulate = HINT** (q≈.42–.76, fails duration); low-γ **cingulate
  encoding** ρ≈.39, 8/8, q=.035. **α = memory-only; β = memory + inference.**
  Duration control: test 1.35–2.53× longer than learn; β ρ≈+0.25 n.s. (duration-
  clean), α ρ≈+0.53 (duration-suspect). audit_103b truncation null RETIRED (invalid).
- **R3 (epi):** relational co-diffusion community (not contact-by-contact); distant
  off-shaft SOZ marker δ heat LOSO AUC .72; **seed-based detector LOPO AUC≈.81,
  prec@5≈60% (~7× chance), 9/10** above chance; **two populations** (tight-community
  vs hub, Pat_10/15 marginal); scope = clinical SOZ labels, **not** surgical outcome
  (triage tool). ~41% of SOZ contacts are WM (a WM filter hurts).

## 4. Known fixes pending (the "what's there to fix")

- **Methods** needs two new paragraphs to match the biology-led R1: (a) the SNR/
  reliability **control**, (b) the Grassmann **baseline-asymmetry** caveat (triangle
  vs split-baseline). Not re-checked this session; Discussion also unreviewed.
- **OFC seam:** OFC now appears in R1.2 (overall β trace) *and* R2.3 (encoding
  component) — true, not contradictory, but add one smoothing sentence.
- **`01_trace.md` N1.2 / N1.2b** (spectral-PCA / raw-edge method-superiority
  sub-sections) → **demote to Methods** (still pending; the "α invisible to
  spectral" claim narrowed under Q2). PI hasn't green-lit yet.
- **`\TODO` inventory:** patient IDs/demographics, exact band edges, phase
  durations, per-band verdict-matrix supp. tables (from VERDICT_LEDGER), detector
  feature list, data/code availability.
- **Re-verify every stated number** against the cached CSVs before submission
  (the draft's own rule).

## 5. Source-of-truth (read alongside)

- Locked verdicts: `.agents/preprint/locked/{VERDICT_LEDGER,ANATOMY_LEDGER,CONTROLS}.md`
  (do **not** edit — cascaded separately by the compute session).
- Taxonomy lock + LOO: `2026-06-25_per-band-consistency-taxonomy-lock.md` +
  `data/audit/localization_atlas/carrier_loo.csv`.
- Q1/Q2: `2026-06-25_multiphase-snr-reliability.md`,
  `2026-06-25_raw-vs-multiscale-trace.md`.
- Per-headline detail: `headlines/{01_trace,02_encoding_vs_inference,03_epileptogenic_markers}.md`.
