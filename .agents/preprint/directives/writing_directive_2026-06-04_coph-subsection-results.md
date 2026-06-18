---
name: writing-directive-coph-subsection-results
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Results encapsulation for the cophenetic (LRG communication-distance) per-pair subsection (ssec title "LRG-based cophenetic distance for the disentanglement of multiscale per-pair communication task-induced trace"). What the measure SHOWS at cohort level + how it clears/does-not-clear controls + per-band reading of the first figure. Measure definition belongs to Methods, NOT here. Changes no locked verdict.
created: 2026-06-04
companion: directives/writing_directive_2026-06-03_raw-copula-figure-swap.md (the raw twin; same copula render). Figures: data/preprint/figures/all_bands/fig_bands_joint_density_copula_coph.pdf (first), fig_bands_null_triangle_coph.pdf (controls).
verified_against: data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (gate p, obs_median_rho, surr p95); data/reports/imcoh_continuous_trace/controls_summary.csv (drift paired Wilcoxon + cross-probe, recomputed this session).
---

# Cophenetic per-pair subsection — results encapsulation

**Head.** Where the raw substrate reported a per-pair trace in *every* band,
the LRG cophenetic (communication-distance) substrate **disentangles** it:
only **β (primary) and α (secondary)** carry a genuine per-pair trace. θ turns
**anti**, and δ, γ_l, γ_h fall to **null**. The α/β trace **clears the
matched-strength gate** (the mandatory control) and is corroborated by the
drift and cross-probe controls; the raw substrate's signal does not survive
this scrutiny. This is the disentanglement the subsection title promises.

> **Do not describe the measure here** (cophenetic distance, the LRG flow, the
> split-baseline ρ^coph) — that is Methods. This subsection reports *what the
> measure shows* and *how it behaves under the controls*. Same copula render and
> justification as the raw twin (do not re-justify the analytic choice at
> length — one back-reference suffices).

---

## 1. The result to encapsulate (one short paragraph)

Frame it as **disentanglement / structural enrichment**, not resolution:
moving from edge weights to the communication-distance hierarchy does not
"reveal hidden bands" — it **strips the non-specific, drift-driven component**
that made the raw substrate light up everywhere, leaving the per-pair trace
concentrated where it survives a proper null: **β and α**. β is the larger,
tightest effect (the primary trace band); α is secondary but equally
gate-significant. This is consistent with the project framing that the LRG
step is structural enrichment, not a confirmation-only or resolution step.

## 2. First figure — `fig_bands_joint_density_copula_coph.pdf` (per-band comments)

Same Gaussian-copula render as the raw figure (cohort-median ρ^coph against the
matched-strength p95 floor), so the colour field already encodes the gate. Short
reading, band by band:

- **β** — tightest green diagonal; the strongest, primary per-pair trace.
- **α** — clear green diagonal, moderate; the secondary trace.
- **θ** — **anti**: red on the diagonal, green on the anti-diagonal (negative ρ^coph).
- **δ, γ_h** — blank/cream: the observed correlation sits at or below the
  matched-strength floor → no trace.
- **γ_l** — also cream **at this floor**: its correlation does not clear the
  matched-strength gate (see §3 — it clears only the weaker controls).

The one-sentence contrast to make explicit: the raw figure's diagonal in every
band collapses here to **two bands**.

## 3. Controls — what clears, what does not (verified, n=10)

| control | clears | does NOT clear | what it rules out |
|---|---|---|---|
| **matched-strength gate** (mandatory) | **α (p=0.002), β (p=0.005)** | γ_l (0.116), δ (0.278), θ (0.722), γ_h (0.246) | strength/topology-matched noise — the level reachable from electrode-strength structure alone |
| **within-session drift** | α (0.007), β (0.014), **γ_l (0.010)** | δ (0.246), θ (0.278), γ_h (0.423) | slow session drift mimicking a trace |
| **cross-probe** (ρ_cross/ρ_full) | α 0.91, β 1.00, γ_l 1.02 (≈1, undegraded) | — | same-probe / shared-shaft geometry artifact |

**Verdict = the gate = α/β.** Both clear it decisively, and both are
corroborated by the drift and cross-probe controls — so the α/β trace is robust
across the whole battery.

**The γ_l caveat — state it honestly, do NOT promote γ_l to a trace band.**
γ_l clears the drift and cross-probe controls but **fails the mandatory
matched-strength gate**; it therefore reads as null in the first (gate-floored)
figure. It is a sub-threshold band, not a third trace band. (Per
`feedback_matched_strength_mandatory`: the gate is the verdict; drift/cross-probe
are corroborating, not promoting.)

## 4. Guards

- Headline is **"mainly β, secondarily α."** No third trace band.
- **No anatomy, no per-patient localization** — both retracted/null this era
  (the trace is spatially delocalized). Keep this subsection cohort-level.
- Do not re-state the raw→coph difference as a count change or "lines turning
  gray"; it is **selectivity** — the signal concentrates onto α/β.
- One back-reference for the copula render ("rendered as in
  \FigRef{fig:rawfc_joint_density}") instead of re-justifying the analytic
  function.
- Controls battery wording must match the raw subsection's (drift null,
  cross-probe, matched-strength). The second figure is the cophenetic
  null-triangle `fig_bands_null_triangle_coph.pdf`.

## 5. Verified numbers (sources)

Figure + gate — `matched_strength_surrogate_split_baseline/cohort_summary.csv`:

| band | obs_median ρ^coph | surr p95 | gate p | first-figure |
|------|------------------|----------|--------|--------------|
| δ    | +0.008 | 0.140 | 0.278 | cream |
| θ    | −0.040 | 0.032 | 0.722 | anti |
| α    | +0.105 | 0.047 | **0.002** | green |
| β    | +0.221 | 0.076 | **0.005** | green (tightest) |
| γ_l  | +0.083 | 0.109 | 0.116 | cream |
| γ_h  | +0.000 | 0.100 | 0.246 | cream |

Drift + cross-probe — `imcoh_continuous_trace/controls_summary.csv` (cohort,
paired Wilcoxon split > drift, one-sided; cross/full ratio of median ρ^coph):
α drift 0.007 / cross 0.91 · β drift 0.014 / cross 1.00 · γ_l drift 0.010 /
cross 1.02 · δ/θ/γ_h drift ≥0.25 (fail). Note: the controls file's per-patient
ρ^coph differs slightly from the matched-strength file's (independent
half-split computations) — cite ρ magnitudes from the figure/gate file, control
pass/fail and p-values from the controls file; the qualitative verdict (α/β
gate, +γ_l drift only) is identical in both.
