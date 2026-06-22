---
name: lowgamma-focal-cingulate-encoding
type: report
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-06-22
headline: N2.4 (.agents/preprint/headlines/02_encoding_vs_inference.md)
pointers:
  - scripts/01_compute/audit/audit_112_within_system_trace.py
  - data/audit/inference_localization/within_system_trace_include.csv
  - data/audit/inference_localization/within_system_trace_exclude.csv
---

# A focal low-γ *memorising* trace in the cingulate, hidden by whole-brain averaging

**Head.** Low-γ has **no whole-brain consolidation trace** — averaged over the
network it looks empty. Yet a region-local test finds a **genuine, cohort-wide
memorising (encoding) trace in the cingulate** that survives the strict null and
the removal of diseased contacts. It is a real signal that cohort-wide averaging
dilutes to nothing — and it is the cleanest argument in the paper for reading the
trace **per region** rather than as one global number. This is the N2.4 backing
result; it is **encoding-only and therefore immune to the duration confound** that
complicated the inference story (encoding never touches the longer `task_test`).

## The result in plain terms

The β and α rhythms carry a consolidation trace you can see across the whole brain.
Low-γ does not — its whole-brain trace is null. The naive conclusion would be
"low-γ does nothing offline." That conclusion is wrong. When we stop averaging over
the whole network and instead ask, region by region, *"does this region's own
wiring keep a copy of the learning phase into the final rest?"*, the **cingulate**
answers yes, and only for the **memorised** content (what was shown), not the
inferred content. The signal was there all along; averaging over hundreds of
uninvolved region-pairs washed it out.

## Technical statement

Measure: the **absolute** within-system trace (`audit_112`) — for every region,
the Spearman concordance between (a) how the region's pairwise cophenetic distances
moved during learning and (b) how they moved into the final rest, with **no
per-patient demeaning** (so it tests an *absolute* trace, not a relative
concentration), tested one-sided against the matched-strength surrogate over the
region's incident pairs.

- **Cingulate, encoding, epi-included:** median concordance ρ ≈ +0.39, **8/8
  sampling patients positive** (6/8 individually significant), Wilcoxon p ≈ 0.004,
  BH q ≈ 0.035. (K = 8 implanted patients — the cingulate is the best-sampled
  system.)
- **Cingulate, encoding, epi-excluded** (diseased contacts removed): ρ ≈ +0.51,
  7/8 positive, p ≈ 0.008, q ≈ 0.070 — a **larger** effect on **fewer** contacts,
  so the trace is **not** an artifact of the epileptic core; it is power-limited,
  not epi-driven.
- **Cingulate, inference:** null (ρ ≈ −0.10, p ≈ 0.98) — the focal low-γ trace is
  **encoding-specific**. Whatever low-γ consolidates in the cingulate, it is the
  *memorised* material, not the *inferred* structure.
- Other systems' low-γ hotspots from the demeaned localizer (insula, OFC,
  sensorimotor) do **not** survive this absolute test → they were demeaning
  structure, not absolute traces.

## Why it matters

It is the existence proof for the method's core claim that a **per-pair, per-region
multiscale read-out** sees structure a global summary cannot. A reviewer who asks
"why not just take one whole-brain number?" is answered by this panel: the
whole-brain low-γ number is zero, and the cingulate trace is real.

## Critical issues (lead with the weaknesses)

- **It is encoding, not inference.** This is a memorising trace, not the headline
  inference result. Frame it as a *band-and-region completeness* finding (low-γ is
  not inert; the cingulate consolidates memorised content there), not as evidence
  for inference.
- **Power floor.** The absolute within-system Wilcoxon **cannot reach BH for
  systems sampled by ≤ 5 patients** (e.g. β-OFC encoding is large, ρ ≈ +0.56, but
  K = 5 → cannot clear). So "only the cingulate clears the absolute test" is partly
  a sampling fact — the cingulate is simply the best-sampled system (K = 8) — **not**
  proof that only the cingulate has a focal trace. Do not over-read exclusivity.
- **Single diffusion scale.** Computed at τ = 1/λ_max only (see the τ-sweep scope).
- **No behaviour.** As with all of N2, this is the neural residue of
  reorganisation, not a correlate of memory *success* (no performance data exists).

## Provenance

- Data → `data/audit/inference_localization/within_system_trace_{include,exclude}.csv`
  · `audit_112_within_system_trace.py` · 2026-06-18 (recomputed/verified 2026-06-22).
- Whole-brain low-γ null (the contrast) → `data/audit/consolidation_arc/arc_null_per_patient.csv`
  (low-γ `T_infspec_pe` p ≈ 0.12, `T_learn` p ≈ 0.12) · `audit_103 --null`.
- Headline home: N2.4 in `.agents/preprint/headlines/02_encoding_vs_inference.md`.
- Memory: [[inference_mark_localization_2026_06_18]].
