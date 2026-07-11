---
name: ms-relative-scale-sweep
type: scope
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-11
implements: audit_172 (fulfils the audit_159-scoped follow-up of 2026-07-07_r1-tau-multiscale-reexamination.md §6)
estimator: rho_sym (canonical, audit_150)
pointers:
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - scripts/01_compute/audit/audit_121_tau_sweep_cophenetic_trace.py
  - .agents/reports/2026-07-07_r1-tau-multiscale-reexamination.md
  - .agents/guides/02_methods/lrg-framework-guide.md
memory: tau_sensitivity_trace_2026_06_22
---

# MS-relative scale sweep of the cophenetic trace — does a band that fails at τ_min pass at its own scale?

## Head

audit_121 mapped the bare cophenetic trace across the diffusion scale
`s ≡ τ·λ_max` for all six bands but carried only the cheap trace-free
**placebo** referee — never the mandatory **matched-strength** null at any
`s ≠ 1`. It left exactly one live lead: the **α-band trace roughly doubles at
s ≈ 1.5–2** in the collapse-free window. This scope runs the referee audit_121
lacked: the canonical **ρ_sym gate (audit_150) swept over s with a
matched-strength null rebuilt at every s**, so the honest quantity is the
signal-minus-null separation `z(s)`, not raw `ρ_sym(s)`. The question:
**is there an intrinsic band ↔ scale match — does a band that fails the cohort
gate at `s = 1` (α first, then δ, low_γ) clear matched-strength at its own
coarser scale — or does the null bump identically and kill the apparent gain?**

## Notation

- `s ≡ τ · λ_max` — dimensionless diffusion scale. `s = 1 ⇔ τ = τ_min = 1/λ_max`
  (finest). `s_F ≡ λ_max/λ_gap ≈ 3` (Fiedler / interior). `s* ≈ 8–15` (the single
  C(τ) susceptibility peak = collapse onset; guide §5.2). Read below `s*`.
  Symbol chosen to avoid collision with the α frequency band (2026-07-11).
- Per matrix `W` (giant component, |ImCoh| FC): `L̂ = D̂ − W`, spectrum
  `{λ_i, v_i} = eigh(L̂)`. Diffusion density
  `ρ(s) = V diag(e^{−(s/λ_max)λ}) Vᵀ / Z`. Communication distance
  `D(s) = 1/ρ(s)` (symmetrised, zero diag, non-finite → row max).
- Cophenetic image `D_coph(s) = cophenet(UPGMA(D(s)))` (average linkage,
  scipy). Condensed upper-triangle vector.
- Trace statistic (canonical, audit_150):
  `ρ_sym(s) = ½[ρ_S(D_task(s)−D_A(s), D_post(s)−D_B(s)) +
                ρ_S(D_task(s)−D_B(s), D_post(s)−D_A(s))]`
  over the four phase cophenetic vectors `{A, B, task_test, rest_post}`
  (A,B = independent rest_pre split-halves, raw-FC halves cache).

## Predicates (what we compute per patient × band × s)

- `ρ_sym^obs(s)` — observed trace at scale s.
- `{ρ_sym^(r)(s)}_{r=1..R}` — matched-strength surrogates: each of the four
  phase matrices independently rewired by the **canonical 4-cycle ±δ
  strength-preserving swap** (audit_150 `_swap_loop`, numba, `SWAP_FACTOR=20`,
  `W_MAX=1`), eigendecomposed once, cophenetic reformed at each s.
- `z(s) = (ρ_sym^obs(s) − mean_r ρ_sym^(r)(s)) / std_r ρ_sym^(r)(s)`.
- `p_1sided(s) = mean_r[ ρ_sym^(r)(s) ≥ ρ_sym^obs(s) ]` (per-cell upper tail).
- Cohort gate `gate_p(band, s) = Wilcoxon(ρ_sym^obs − ρ_sym^surr_p50,
  alternative='greater')` over the 10 patients — **identical form to
  audit_150**, now a function of s.

## Anchor (correctness gate before any interpretation)

Seeds and swap match audit_150 exactly (`BASE_SEED=20260706 + cell_idx`,
`R=200`, same phase iteration order A→B→task→post). Therefore the `s = 1`
column MUST reproduce audit_150's `per_patient_per_band.csv` `obs_rho`,
`surr_*`, and `cohort_summary.csv` `gate_p_sym` **bit-for-bit**. Assert
`|Δ| < 1e-10` on obs and on the surrogate mean at s=1 for every cell; abort on
mismatch. This is the audit_121-style anchor: any s ≠ 1 movement is real scale
physics, not a code difference.

## Properties / expected behaviour

- `s = 1`: β, α CLEAR the gate; δ/θ/low_γ/high_γ FAIL (reproduces audit_150).
- β: `z(s)` broad and positive finest → Fiedler (scale-broad, τ-robust).
- α (the test): if the doubling is real, `z(s)` and `−log gate_p(s)` **peak at
  s ≈ 1.5–2.5** and the gate crosses 0.05 there while the surrogate stays flat.
  If it is generic coarse-graining geometry, the surrogate mean rises with the
  observed and `z(s)` stays flat/declining — no gate crossing.
- Negative-control discipline: δ/θ/high_γ are dead-null at s=1; if any "passes"
  at a coarse s the surface is manufacturing significance (the audit_103b
  lesson) → that s is past-trust, not a result.

## 5-point critical preamble (mandatory)

1. **Claim.** A band that fails the cohort trace gate at `s = 1` has a genuine,
   matched-strength-clean cophenetic trace at its own coarser scale
   (`α` at `s ≈ 1.5–2.5`); equivalently, there is a band ↔ scale match rather
   than one universal `τ_min`.
2. **Null.** Canonical matched-strength (4-cycle ±δ strength-preserving,
   `R = 200`), cophenetic reformed at each s from the surrogate's own
   eigenpairs. Cohort Wilcoxon of (obs − surr median), one-sided greater.
3. **Strongest alternative the null must kill.** A mesoscale `ρ_sym` bump is
   generic diffusion-geometry of *any* strength-heterogeneous dense graph
   (coarse-graining pulls every phase onto the slow Fiedler mode, inflating
   cross-phase concordance). A matched-strength surrogate has the *same*
   strength sequence and therefore the *same* slow-mode geometry, so it bumps
   identically → `z(s)` stays flat. This is exactly what audit_121's placebo
   could only weakly probe. The discriminator is **obs − surr**, per band, and
   **α vs the dead-null control bands** at the same s.
4. **What the null cannot do / neg-control discipline.** It cannot certify a
   reading taken past collapse onset (geometry degenerate, Spearman variance →
   0): restrict trust to `s < s*` (≈ 8–15) and flag/exclude any s with a
   zero-variance guard trip. δ/θ/high_γ must stay null at every s; if one
   "clears" at coarse s the reading is collapse, not trace — retire it.
5. **Falsification.** PASS = α clears matched-strength at `s ≈ 1.5–2.5`
   (gate p < 0.05, leave-one-patient-out robust, `z(s)` peaked, geometry not
   collapsed) while every control band stays null. FAIL = α does not clear at
   any collapse-free s, OR a negative-control band clears → the α mesoscale is
   observed-only / manufactured, the headline stays at `s = 1`. β is reported at
   `s = 1` regardless (its value is a robustness reframe, not a new gate).

## Pseudocode

```
SGRID = unique([0.5,0.7,0.85] ∪ geomspace(1.0, 12.0, 13))     # dense near 1–3
for (idx, pat, band) in cells:                                # 10×6 = 60
    Ws  = {A,B,task_test,rest_post}  (raw-FC halves + loader)
    Eig = {ph: eigh(diag(sumW)-W)}                            # once per phase
    obs[s] = rho_sym( coph_s(Eig[ph], s) for ph ), all s      # reuse eigenpairs
    rng = default_rng(BASE_SEED + idx)                        # == audit_150
    for r in 1..R:
        Es = {ph: eigh(shuffle(Ws[ph], n_swaps, rng))}        # canonical swap
        for s in SGRID: surr[s,r] = rho_sym( coph_s(Es[ph], s) for ph )
    z[s], p1[s] from obs[s], surr[s,:]
assert obs[s=1], surr[s=1] == audit_150 (bit-exact)
cohort: gate_p(band,s) = wilcoxon(obs - surr_p50, 'greater') over patients
```

Optimisation (per feedback_optimize_time_and_surface_progress): eigendecompose
each matrix **once** and reuse across all s (the D(s) = V e^{−(s/λ_max)λ} Vᵀ
trick, exactly audit_121); numba swap reused verbatim from audit_150 (RNG drawn
outside the jitted loop → bit-identical); multiprocessing Pool over the 60 cells;
live `[i/N] patient/band` progress with `flush=True`; pilot one cell and
extrapolate before the full launch.

## Visualisation (after verdict)

Per-band `z(s)` and `−log10 gate_p(s)` vs s on one axis, with s=1 / s_F / s*
guides; α highlighted; control bands greyed; collapse-onset shaded. Only if α
PASSES does a figure/supplement become relevant (else it is documented-negative
housekeeping, per the 2026-07-07 decision).

## Connection to prior tools

- Generalises **audit_150** (ρ_sym gate) from `s = 1` to a scale sweep.
- Supplies the matched-strength referee that **audit_121** (observed-only +
  placebo) explicitly deferred; same eigenpair-reconstruction engine.
- Mirrors **audit_157** (arc ρ_sym × τ × MS at 3 τ) but on the *bare* trace and
  a *dense* s-grid; reuses the same swap ensemble methodology (not the same
  cached seed — we regenerate with audit_150's seed to anchor the gate).

## Open questions (explicitly out of scope here)

- Encoding / inference components over s (`T_learn`, `T_infspec_pe`) — Part B,
  separate run; audit_157 already has a 3-point MS version favouring s ≈ 2.6.
- β→OFC per-node localisation at α's preferred s — only if α PASSES.
- Whether a certified α mesoscale changes the Methods "τ-sweep is degenerate"
  language — a writing decision downstream of the verdict, not settled here.

## Part B — encoding / inference components over a dense s-grid (audit_173)

Same construction, applied to the **decomposition** functionals instead of the
bare trace: `T_learn = ρ_sym(task_learn−pre, post−pre)` (encoding echo) and
`T_infspec_pe = partial-ρ_sym(inference, post | encoding)` (inference-specific,
the β-only headline), over the 5-phase arc
`{rest_pre_A, rest_pre_B, task_learn, task_test, rest_post}`.

- **Generalises audit_157** (ρ_sym × τ × matched-strength at 3 τ: 1.0, 2.610,
  6.813) to a **dense s-grid** — the exact analogue of audit_172 ⊃ audit_150.
  Reuses audit_157's machinery verbatim (`_sym_functionals`,
  `_surr_functionals_sym`, `_cophenet_at_taumult`) and the **cached
  seed-20260511 surrogate eigenpairs** (cophenetic reformed at each s, NO
  regeneration → cheap).
- **Anchors:** arm1 at s ∈ {1.0, 2.610, 6.813} reproduces audit_103d (ρ_split);
  sym at s = 1 reproduces audit_152 (ρ_sym fine-scale). Both asserted.
- **5-point preamble = audit_157's** (inference-specific consolidation is
  multiscale, mildly mesoscale-favouring; MS null at each s; neg-control bands
  δ/θ/γ must stay null at every s; read below collapse onset; falsify if β
  T_infspec_pe fails MS at every collapse-free s or a control clears).
- **Read (same discipline as Part A):** off the collapse-free window and the
  *shape* — a scale-tuned component **peaks** (audit_157 hinted β inference peaks
  at s ≈ 2.6); a collapse channel **climbs monotonically**. The full curve shows
  which. Output: `data/audit/arc_scale_sweep_rhosym/`.
