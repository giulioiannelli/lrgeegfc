---
name: sparsified-arc-master-scope
title: Sparsified multiscale arc — master scope (protocol + 3-point arc, τ-resolved)
date: 2026-07-11
status: current
era: IMCOH_ABS x COHORT_N10 (percolation-backbone multiscale)
supersedes_in_spirit: 2026-07-11_sparsification-propagator-necessity.md
---

**Head.** On the fully-connected `imcoh_abs` graph the LRG propagator at
`τ=1/λmax` is degenerate: `D=1/K ≈ D=1/A`, so every prior result was effectively
*single-scale* (audit_174/176). This scope rebuilds the whole arc on a
**parameter-free percolation backbone** (keep every edge ≥ the connectivity
bottleneck `θ*` = weakest maximum-spanning-tree edge — the largest global
threshold at which the giant component still spans all N nodes), where the graph
is **connected and cycle-rich**, and makes **τ (equivalently the dimensionless
scale `s ≡ τ·λmax`) the one and only swept parameter**. A result may now be
*single-scale* (fires at one s), *multiscale* (fires at several separated s), or
*continuously-multiscale* (fires across the whole s-band). The matched-strength
null remains mandatory for every cohort claim.

---

## 1. Notation

- `W_p` — cleaned symmetric zero-diagonal `imcoh_abs` FC for phase
  `p ∈ {A, B, task_test, rest_post}` (A,B = split-half of `rest_pre`).
- `θ` — global weight threshold. `A_t = W ⊙ 𝟙[W ≥ θ]`.
- `P∞(θ)` — node fraction in the largest connected component (LCC) of `A_t`.
- `E∞(θ)` — edge fraction *inside* the LCC / all edges.
- `θ*` — bottleneck = `max{θ : P∞(θ)=1}` = `min` over maximum-spanning-tree edges.
- `B_p = W_p ⊙ 𝟙[W_p ≥ θ*_p]` — the **percolation backbone** (per phase, its own θ*).
- `L̂ = D̂ − B` — combinatorial Laplacian; `ev,V = eigh(L̂)`; `λmax = ev[-1]`.
- `s ≡ τ·λmax` — dimensionless diffusion scale (`s=1` ⇔ `τ=τ_min=1/λmax`).
- `K(τ)=e^{−τL̂} ≥ 0` exactly (Metzler; underflow floored at `+1e-30`);
  `ρ=K/Tr K`; `D_ij(τ)=(1−δ_ij)/ρ_ij`.
- `Ŝ(τ)=−Tr[ρ ln ρ]/ln N ∈ [0,1]`; `C(τ)=−dŜ/d log₁₀τ`.
- `D^coph_p(s)` — condensed cophenetic distances, `UPGMA(D(τ=s/λmax_p))`.
- `ρ_sym(s) = ½[ Spearman(D_task−D_A, D_post−D_B)
            + Spearman(D_task−D_B, D_post−D_A) ]` on `D^coph(s)`.

## 2. Predicates / classification (per (patient,band) cell, per arm)

At each `s` on a pre-committed grid, matched-strength gate → per-cell `p(s)`.
Define the fire-set `F = {s : p(s) < 0.05 ∧ ρ_sym(s) > surrogate median}`:
- **none** — `F = ∅`.
- **single-scale** — `F` is one contiguous band spanning < ⅓ of the grid.
- **multiscale** — `F` has ≥2 separated bands.
- **continuous-multiscale** — `F` spans ≥ ⅔ of the grid.
Cohort verdict per band = Wilcoxon(obs−surr_p50, greater) at each `s`, plus a
scale-max gate (null also maxes over `s`) to control scale-selection.

## 3. τ / s range (the "S from 0 to 1" rule)

- Entropy/specific-heat **plots**: grid chosen so `Ŝ` spans ~[0.999, 0.001];
  mark `τ_min=1/λmax` and `τ_max=1/λ2` as references. `C(τ)` peaks = the
  data's characteristic scale(s) (Villegas: 1 peak = collapsed single scale;
  >1 = multiscale ladder).
- Arc **sweep**: `s ∈ [1, s_max]` where `s_max` is where `Ŝ` of phase A drops
  to a floor (0.02, "graph collapsed to one cluster"), capped `[1,200]`, 16
  log-spaced. `τ_max = 1/λ_min = 1/0 = ∞` is undefined, hence the entropy floor.

## 4. Nulls (mandatory)

- **Matched-strength** (primary): 4-cycle ±δ strength-preserving shuffle
  (`matched_strength_shuffle`, numba, bit-identical to audit_150), **sparsified
  the same way** (its own θ*) at every `s`. Upper-tail `p=mean(surr≥obs)`.
- **Drift** (secondary, whole-task windows): reused from the established drift
  harness where a conditional (enc/inf) claim needs it.
- Same per-cell RNG seed (`BASE_SEED+idx`) and phase order as audit_150 so the
  `s=1` full-graph anchor reproduces the canonical gate bit-for-bit.

## 5. Five-point critical preamble

1. **Claim.** On the percolation backbone the cross-phase trace is band-selective
   and **not reducible to a single scale** — some bands fire only at a
   characteristic mesoscale `s>1` that the raw adjacency (single scale) cannot
   resolve — and survives matched-strength.
2. **Null.** The backbone changes nothing (same picture as `s=1`), OR the trace
   is a sparsification artifact reproduced by the matched-strength surrogate.
3. **Strongest alternative the null must kill.** (a) scale-selection winner's
   curse → scale-max gate lets the null pick its own best `s`; (b) global
   strength shift masquerading as trace → matched-strength holds strength fixed;
   (c) density p-hacking → θ* is data-set, not chosen, and reported per cell;
   (d) shared-draw artifact → identical surrogate draws feed every arm, `s=1`
   reproduces audit_150.
4. **What it cannot reject.** A backbone is a modeling choice; a percolation-only
   result must be cross-checked against the moderate-density robustness sweep.
   Matched-strength cannot reject a genuine strength-independent reorganization
   that is *also* single-scale — that is a real (if less exciting) outcome.
5. **Falsifiers.** If every band that clears clears at `s=1` with a flat `ρ_sym(s)`,
   there is no multiscale content (back to single-scale). If β dies on the
   backbone, the trace lived in the weak-edge bulk.

## 6. Deliverables (drives the compute order)

D0 protocol: `P∞,E∞ vs θ` per patient/band/phase + stability report.
D1 diagnostics: `Ŝ(τ),C(τ)` and the backbone **networks** per patient/band/phase.
D2 trace arc: `ρ_sym(s)` + gates + per-cell class + β localization (over-expression
   before p-values).
D3 enc-vs-inf arc, τ-resolved, per-band nature + patient spread.
D4 epilepsy marker/detector, τ-resolved, **precision-first**, seeded + unseeded,
   vs state of the art.
D5 reconciliation report: what stands, what is now single-scale, what is retracted.

## 7. Library kernel (single source of truth)

- `lrg_eegfc.utils.fc.backbone.percolation_backbone / percolation_sweep`
- `lrg_eegfc.utils.fc.heat_multiscale` — `laplacian_eig`, `entropy_specific_heat`,
  `specific_heat_peaks`, `scale_grid`, `cophenetic_at_tau`, `rho_sym`,
  `rho_sym_over_scales`.
- `lrg_eegfc.utils.metrics.surrogate.matched_strength_shuffle` (numba,
  byte-copied from audit_150 `_swap_loop` → bit-identical anchor).
- Data loading stays script-side (`load_phase` in audit_150; dataset-specific).

## 8. Output location

All artifacts under `data/sparsified_arc/` (NOT the preprint). Each heavy cache
ships a `config.json` (cohort, bands, θ*-rule, s-grid, R, seed, git SHA).

## 9. Open directions (log as found)
- Does the backbone cophenetic **localize** β→OFC better than raw/geodesic
  (audit_171 axis on the backbone)?
- Unseeded epileptic-node detection now that the graph is thresholded (the dense
  graph may have hidden it).
- Per-band "why it fires/doesn't": local vs global, patient spread as signal.
