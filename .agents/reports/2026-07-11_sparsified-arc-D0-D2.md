---
name: sparsified-arc-D0-D2
title: Sparsified multiscale arc — D0 protocol + D1 diagnostics + D2 trace (exploratory)
date: 2026-07-11
status: current
era: IMCOH_ABS x COHORT_N10 (percolation-backbone multiscale)
tags: [sparsification, percolation, tau-sweep, rho_sym, band-selectivity, exploratory]
scope: recover the LRG multiscale on a percolation backbone after audit_174 showed the fully-connected propagator is degenerate
---

**Head.** On a parameter-free percolation backbone the LRG diffusion τ-sweep is
**no longer degenerate** (D1: 32–62 % of cells develop a multiscale C(τ) ladder,
vs ~0 % fully-connected), and the cross-phase trace becomes **τ-resolved**. The
clean scientific win is **β**: its trace is *scale-localized* — absent at
τ_min (s=1, p=.22) and emerging at a **mesoscale s≈14–17** (matched-strength
scale-max p=.042), robust across percolation + two fixed densities. **BUT there
is a red flag the user named explicitly: band-selectivity collapses on the
backbone** — δ, θ, β, γ_high all clear the matched-strength gate, and by raw
p-value δ/γ_high are *more* robust than β. If every band traces, that is a
**negative** result (no band-specificity = the contribution is gone). The
decisive open test is the **drift null**: does a stronger control restore
selectivity (kill δ/θ/γ_high, spare β)? All results live under
`data/sparsified_arc/` (NOT the preprint). Matched-strength stays mandatory.

---

## What was built (library-first, one canonical kernel)
- `src/lrg_eegfc/utils/fc/backbone.py::percolation_sweep` — P∞/E∞ vs θ curve.
- `src/lrg_eegfc/utils/fc/heat_multiscale.py` — `laplacian_eig`,
  `entropy_specific_heat` (Ŝ∈[0,1] auto-ranged), `specific_heat_peaks`,
  `scale_grid` (s=τλmax∈[1, s(Ŝ=.02)] cap 200), `cophenetic_at_tau/scale`,
  `rho_sym`, `rho_sym_over_scales`.
- `src/lrg_eegfc/utils/metrics/surrogate.py::matched_strength_shuffle`
  (numba, byte-copied from audit_150 → bit-identical anchor).
- Scripts `scripts/01_compute/sparsified_arc/0{1..6}_*.py`. Scope report:
  `.agents/guides/task-persistence-investigation/2026-07-11_sparsified-arc-master-scope.md`.
- Infra note: run with the env python directly
  (`/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`) —
  `conda run -n lapbrain` fails on a cross-compiler activate.d hook.

## D0 — percolation protocol (VALIDATED, with one honest caveat)
`data/sparsified_arc/percolation/` (+figs). 240 cells in 21 s.
- **0/240 tree-like, 0/240 disconnected** — always connected & cycle-rich
  (cycle-rank median ~1100, min 6). Good.
- **θ\* is a single-edge quantity** (weakest max-spanning-tree edge), so the
  surviving-edge fraction swings **2 %→70 %** across cells. This is **structured
  by patient, not noise**: Pat_14 near-tree (median 2.5 %, its vendor-replaced
  FC), Pat_15 near-dense (43 %, RH-only implant); 75 % of cells sit in a clean
  5–35 % moderate range. Verdict: **usable as primary**, with a **fixed-density
  cross-check** (d0.15/d0.20) so no claim rests on the single-edge bottleneck.
  Fig `percolation_stability.pdf`.

## D1 — entropy / specific heat (multiscale partially restored)
`data/sparsified_arc/entropy/` (+figs). 240 cells in 2 s.
- C(τ) median 1 peak most bands, **2 for γ_high**; **32–62 % of cells >1 peak**
  (fully-connected was ~0 %). Sparsification does move cells off the
  Wigner-semicircle single-scale regime (Villegas 2025). Dominant scale is
  **coarse, s≈50–180**, i.e. away from τ_min — "result at a different τ."
  Fig `entropy_summary.pdf`, `entropy_curves_{band}.pdf`.

## D2 — τ-resolved trace ρ_sym(s) (the headline + the red flag)
`data/sparsified_arc/trace_arc/` (percolation, primary) +
`trace_arc_d015/`, `trace_arc_d020/` (robustness). R=200, 16 scales.
Estimator ρ_sym; matched-strength null sparsified the same way; scale-max gate
(null maxes over its own best scale). s=1 = τ_min anchor.

**Band nature (percolation):** β single-mesoscale (p_smax=.042 @ s≈17, FAILS
s=1 p=.22); δ/γ_high "continuous-multiscale" (fire at ~all scales incl s=1);
θ single-scale (.0098); α none (.116); γ_low none (.080).

**Robustness (clears + LOO across percolation/d0.15/d0.20):**
| band | verdict | gate perc/.15/.20 | LOO perc/.15/.20 |
|---|---|---|---|
| delta | **ROBUST** | .014/.007/.002 | .027/.014/.004 |
| high_gamma | **ROBUST** | .001/.019/.019 | .002/.037/.037 |
| beta | robust-gate, LOO-borderline@perc | .042/.002/.002 | .082/.004/.004 |
| alpha | partial 2/3 (fixed only) | .116/.019/.019 | — |
| theta | partial 1/3 (perc only) | .010/.116/.461 | — |
| low_gamma | none | .080/.097/.116 | — |

**The red flag (user, 2026-07-11): "if all bands trace, that is terrible news —
no band-variability is a negative result."** On the full graph the trace was
**band-selective (α/β only)** — the manuscript's contribution. On the backbone
under matched-strength, **selectivity collapses** (δ/θ/β/γ_high all clear; δ/γ_high
strongest). This is the failure mode, not a richer result.

**Why β is still the honest anchor — mechanism, not p-value:**
1. **Scale-localization.** β alone is *absent at τ_min and emerges at a mesoscale*
   — the genuine multiscale value-add ("resolves what no single scale can").
   δ/γ_high are *uniformly elevated at every scale incl τ_min* — the signature of
   a weakened null (a global offset), and γ_high clears at τ_min on the backbone
   where the fully-connected graph did **not**.
2. **Biology.** β's per-patient structure is the established 6-trace /
   Pat_07+15-anti / Pat_14-weak tier (`rho_sym_beta_tier`), not a new pattern.
3. **Prior drift-robustness.** β is the band with established drift-robustness
   (C2 β.0137; β infspec duration-robust .007); δ/γ_high are prime drift/EMG
   artifact candidates and drift-UNTESTED here.

## Reconciliation with prior (1/K(τ_min)) results
- audit_150 full-graph gate: α .024, β .032 CLEAR (band-selective). audit_174/176:
  at τ_min the propagator≈adjacency → that gate was effectively **single-scale**.
- On the backbone β **survives and gains a mesoscale signature** (the recovery
  goal). α **does not survive percolation** (density-sensitive; clears only fixed
  density) — a genuine loss/softening for α. δ/θ/γ_high **newly appear** —
  provisionally an artifact of null-weakening on sparse/dense backbones.
- **Net:** the LRG multiscale is recovered as *non-degenerate* and β is *more*
  interesting (mesoscale). But band-selectivity is **not** preserved under
  matched-strength alone; establishing it requires the drift null.

## Open / next (in flight or pending)
- **DRIFT null (decisive).** Does it restore selectivity (spare β, kill
  δ/θ/γ_high)? Recipe scout dispatched; whole-task drift valid for T_test, INVALID
  for conditional T_infspec_pe (use MS). This is the single most important
  remaining control.
- D3 enc/inf arc (running, `enc_inf_arc/`): T_test/T_learn/T_infspec/T_infspec_pe(s);
  headline is β-only T_infspec_pe (prior full-graph ρ_sym p=.0098).
- D4 epi marker (running, `epi_arc/`): multiscale vs single-scale heat_t5, precision-first.
- β localization on the backbone (OFC over-expression before p-values) — not yet run.
- D1 networks per cell — not yet plotted.
- D5 full reconciliation report.

## Caveats
Exploratory. Matched-strength alone does not gate a cohort claim here because it
loses band-specificity; drift is required. Do NOT migrate anything to the preprint
until selectivity is re-established. Per-patient plots exist for every step
(control against silly mistakes).
