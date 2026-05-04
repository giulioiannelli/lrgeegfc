---
name: 2026-04-29_e1-cohort-verdict
type: report
era: COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/plans/active/2026-04-29_eigenvector-direct-pivot-plan.md
  - .agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md
  - .agents/reports/2026-04-29_critical-post-mortem.md
  - data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv
  - data/audit/cohort_coverage_matrix/triangulation_n10_imcoh_abs.csv
---

# E1 cohort verdict — first eigenvector-direct probe lands at Outcome B

## Renormalization head

**E1 (Grassmann chordal distance between top-k eigenspaces of L̂ across phases) does not surface a cohort-positive task trace in any band at n=10 IMCOH_ABS. δ remains partition-resolution-locked at L5_k; all other bands ergodic. Theta is negatively positioned (mean Δ_E1 = +0.171, frac_pos=0.20, classified `V=negative`); low_gamma is the strongest near-miss (frac_pos=0.70, mean Δ_E1 = −0.142). The eigenvector-direct hypothesis from the pivot plan is not supported by the cleanest probe — pivot plan §7 Outcome B applies. The negative result is now stronger because it covers both tree-based (L1, L5_k, L4, L7) and eigenvector-direct probes.**

---

## TOC

1. **No band crosses the 0.8 cohort threshold at any k ∈ {2, 3, 5, 8, 13, 21}.**
2. **Theta is the only band with V_E1 = `negative` (frac_pos=0.20, mean Δ > 0).**
3. **Low_gamma is the strongest near-miss: 7/10 patients positive, most negative mean Δ_E1.**
4. **Wilcoxon q across the six bands sits at 0.93 (low_gamma at 0.51) — no signal stronger than chance.**
5. **No new triangulation predicate fires: zero `eigenvector-confirmed` bands.**
6. **The δ partition-resolution-locked verdict is unchanged.**
7. **Pivot plan acceptance gate: Outcome B (negative result strengthened).**

---

## 1. Per-band E1 results

Frac_pos at each k (number of patients positive, out of 10):

| band | k=2 | k=3 | k=5 | k=8 | k=13 | k=21 | max | mean Δ_E1 | V_E1 |
|---|---|---|---|---|---|---|---|---|---|
| theta | 1 | 1 | 1 | 2 | 1 | 1 | 0.20 | +0.171 | **negative** |
| alpha | 2 | 2 | 3 | 2 | 3 | 4 | 0.40 | +0.003 | silent |
| delta | 4 | 2 | 3 | 4 | 3 | 4 | 0.40 | −0.035 | silent |
| high_gamma | 4 | 3 | 3 | 5 | 3 | 2 | 0.50 | −0.010 | silent |
| beta | 3 | 2 | 4 | 3 | 4 | 6 | 0.60 | +0.039 | silent |
| **low_gamma** | 7 | 5 | 3 | 6 | **7** | **7** | **0.70** | **−0.142** | silent |

**Cohort threshold = 0.8 (8/10).** No band passes. Low_gamma is the closest, with consistent trace direction (mean negative) and three k-values at 0.7. Wilcoxon raw p for low_gamma at the band level = 0.084, BH-FDR q over 6 bands = 0.51 — does not reach significance.

## 2. Triangulation matrix (v1.2)

```
band         V_L1     V_L5_k      V_E1       T
delta        silent   positive    silent     partition-resolution-locked
theta        silent   silent      negative   ergodic
alpha        silent   silent      silent     ergodic (dissenter-suppressed)
beta         silent   silent      silent     ergodic
low_gamma    silent   silent      silent     ergodic
high_gamma   silent   silent      silent     ergodic (dissenter-suppressed)
```

Zero bands satisfy `eigenvector-confirmed = (V_L5_k = positive AND V_E1 = positive)`. δ keeps its partition-resolution-locked verdict (the L5_k modular swap at k=20–32 does not show as a subspace rotation at the chosen E1 grid).

## 3. What this confirms and what it doesn't

**Confirmed.**
- δ trace is real but partition-resolution-locked (consistent with 2026-04-27 §8.8).
- The cohort task-trace hypothesis at the cohort-wide level (≥ 8/10 in any band) is not supported by a clean eigenvector-direct probe.
- The methodological pivot was worth running — the same continuous-spectrum FC graphs that killed scalar-spectrum probes (S, C, Ψ, Ψ_L, τ*) also do not present a cohort-wide subspace rotation at the chosen k-grid.

**Not confirmed.**
- A finer k-grid or a different τ might still surface signal — E1 was registered at τ = 1/λ_max with a sparse Fibonacci-ish k-grid {2, 3, 5, 8, 13, 21}. A sensitivity sweep over k or τ remains an open question (E1 scope §10).
- E2 (per-node diffusion-map drift, Wasserstein on per-node embedding distances) and E3 (IPR distribution shifts) test different geometric reductions of the same eigenbasis. Whether they surface signal that E1 averages out is not yet known. By pivot plan §7 Outcome B they remain optional / diagnostic; the user can decide whether to run them or freeze.

**Anti-finding (the strongest negative).** Theta has a *positive* mean Δ_E1 across the cohort (Δ_E1(test, post) > Δ_E1(pre, post)) — task-test pushes the diffusion subspace *farther* from rest_post than pre was, not closer. This is consistent with theta's L5_k silence; theta task is not a subspace-rotating event in any direction.

## 4. Acceptance gate decision (pivot plan §7)

**Outcome B.** "If E1 confirms ergodic for α / β / γ_h, the negative result is much stronger now (it covers both tree-based and eigenvector-direct probes). Paper headline is 'we tested both tree-based and eigenvector-direct LRG probes; only δ partition-resolution-locked survives controls.'"

This is now the supported reading. δ remains the only band with a defensible task-trace claim, and its claim is narrowly scoped: a partition-label reorganisation at integer-k that does not project onto either pair-distance (L1) or subspace rotation (E1).

## 5. Recommended next steps

1. **Decide on E2 / E3.** They test orthogonal geometric reductions (per-node diffusion drift, IPR distribution shift). If the user wants the negative result tightened, run them; if Outcome B is acceptable as-is, freeze the pipeline at v1.2.
2. **Sensitivity sweep on E1 (optional).** Re-run E1 with a finer k-grid (e.g. `range(2, N//2, 2)`) and report whether low_gamma's near-miss tightens or relaxes. If a narrow ridge appears at intermediate k, the v1.2 verdict deserves reconsideration.
3. **Lock the v1.2 paper headline.** With no eigenvector-confirmed band and δ alone surviving, the paper becomes: (a) δ is real and locked; (b) every other band is honestly ergodic at n=10 IMCOH_ABS; (c) we tested both tree-based and eigenvector-direct LRG probes; (d) the LRG framework's spectrum-scalar ladder (S, C, Ψ) is non-informative for our continuous-spectrum FC outlier case (already concluded 2026-04-29).

## 6. Provenance

- E1 producer: `scripts/01_compute/audit/audit_27_e1_subspace_alignment.py` (sanity gate passed at k=3: within=0.07, shuffle=1.72, reflex=1.1e-7).
- E1 cohort CSV: `data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv` (360 rows = 10 patients × 6 bands × 6 k).
- Matrix integrator: `scripts/01_compute/audit/audit_24_cohort_coverage_matrix.py` extended with `populate_e1_subspace`, RUNG_ORDER + E1_subspace, `eigenvector-confirmed` triangulation predicate.
- Run metadata: `data/audit/cohort_coverage_matrix/audit_24_run_metadata.json`.
- Headline figure: `data/outputs/figures/section6/cohort_coverage_matrix_n10_imcoh_abs.pdf` (now 9 columns wide; E1_subspace adjacent to L7).
- Phase 0 cache: full-phase + halves LRG re-cached with eigenvectors persisted (240 + 240 NPZs); spectral identity verified by audit_26 smoke test.
