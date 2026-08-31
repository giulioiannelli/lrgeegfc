---
name: pipeline-contract
kind: locked-contract
era: PAPER_FINALIZATION (Wave 0, lane W0-A)
status: locked_2026-08-31
created: 2026-08-25
scope: The locked substrate every analysis lane inherits — the ImCoh magnitude transform, the sparsified backbone (a knob-integrated window, not a single fraction), the phase set, and the diffusion-scale grid. States the pre-registered rule that chose each, the evidence, and an explicit list of what is NOT settled.
supersedes: the bare "mst@0.20" choice recorded in .agents/guides/02_methods/sparsification-choice.md §4a
companion: .agents/reports/2026-08-25_w0a-substrate-contract.md
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
backbone: mst-union, knob-integrated over f in [0.07, 0.20]
pointers:
  - scripts/01_compute/paper_final/w0a_00_preregistration.md
  - src/lrg_eegfc/workflow/substrate.py
  - data/paper_final/w0a_substrate/
---

# Pipeline contract — the locked substrate

## Head

The substrate is **`imcoh_abs`** on an **mst-union backbone reported over the window `f ∈ [0.07, 0.20]`**, never at a single fraction. That window is where α and β both carry a cohort trace while θ does not, it is contiguous and 1.51 octaves wide, and it was found by sweeping the knob rather than by choosing it. Two things a reviewer will ask are answered honestly and negatively: **no parameter-free filter qualifies** (TMFG α p = 0.141, β p = 0.086; percolation α p = 0.312, β p = 0.076 — neither clears, so the free parameter cannot be eliminated, only integrated over), and **the full band-selectivity claim is narrower than the trace claim** — requiring low_γ to be null as well shrinks the window to `f ∈ [0.14, 0.20]`, 0.51 octaves, which is below the pre-registered one-octave bar. The incumbent `f = 0.20` sits at the **upper edge** of both windows, which is exactly the signature of an outcome-selected knob; the contract moves the representative fraction to the interior.

---

## 1. The contract

| choice | locked value | how it was chosen |
|---|---|---|
| FC transform | `imcoh_abs` = `⟨\|Im C\|⟩_f` | R1 selected `imcoh_sq` on split-half reliability (5/6 bands, p = 5.0e-6, margin ≈ 1.2 % relative). `imcoh_sq` then **failed the frozen R2 admissibility gate**: its admissible window is a single fraction, 0.00 octaves. β — the paper's strongest claim — is positive at 6/10 swept configs under `abs` and only 3/10 under `sq`. See §4.1; this is the one place my pre-registration was under-specified and the resolution is flagged as such. |
| backbone family | mst-union (`MST ∪ strongest f of edges`) | The only family in the panel with an octave-wide admissible window. |
| backbone density | **knob-integrated over `f ∈ [0.07, 0.20]`**, fractions `(0.07, 0.10, 0.14, 0.20)`; representative single graph `f = 0.14` | R2.1 (parameter-free) failed for every candidate → R2.3 fallback. `f = 0.14` is an interior point of both windows; the incumbent `0.20` is at the edge of both. |
| phase set | 5 phases, `task_learn` first-class | The paradigm separates learning a structure from applying it. §2. |
| scale grid | `s = τ·λ_max ∈ logspace(0, log10 180, 16)` | Incumbent grid, unchanged, so this contract stays comparable to existing artifacts. §3. |
| cohort gate | sign-flip cluster mass over the whole scale axis | The 16-scale axis is worth **n_eff = 1.51** independent tests (mean cross-scale r = +0.79), so per-scale multiplicity corrections over-charge ~11×. §5. |
| null | matched-strength on the dense FC, then sparsify | Incumbent. **Not sufficient** — see §6. |

### How to use it

Every lane calls one function. Do not re-derive the graph inline, do not hardcode a phase tuple, do not hardcode a fraction.

```python
from lrg_eegfc.workflow.substrate import (
    CANONICAL, canonical_graph, canonical_graph_ensemble,
    canonical_phase_eigs, canonical_scale_grid,
)
from lrg_eegfc.utils.fc.heat_multiscale import cross_phase_functionals_over_scales

eigs = canonical_phase_eigs("Pat_02", "beta")        # five phases, from CANONICAL_PHASES
T    = cross_phase_functionals_over_scales(eigs, canonical_scale_grid())
# -> {"T_probe", "T_encode", "T_probespec", "T_probespec_pe"}, each length 16
```

**Any number that goes in the paper must be knob-integrated**, i.e. computed at every fraction of `CANONICAL.plateau_fracs` and reported as the median with the across-fraction spread. `canonical_graph_ensemble(patient, phase, band)` yields `(frac, graph)` over exactly those fractions. A single-fraction number is a tuned number and will be treated as one.

A four-phase caller passes `phases=("A", "B", "task_test", "rest_post")` to the same function and gets `{"T_probe"}` alone — identical to the incumbent `rho_sym_over_scales` to 1e-16. A raw-FC baseline is `canonical_graph(..., dense=True)`: a visible argument, not a different code path.

---

## 2. The phase set (locked)

Five phases, because the paradigm is transitive inference and separating **learning a structure** from **applying it** is part of the science:

| role | phase | what it is |
|---|---|---|
| `baseline_a` | `A` | first contiguous half of `rest_pre` |
| `baseline_b` | `B` | second contiguous half of `rest_pre` |
| `encode` | `task_learn` | ordered premises presented |
| `probe` | `task_test` | novel non-adjacent pairs judged — answerable only by inference |
| `follow` | `rest_post` | resting state after the task |

Four cross-phase functionals, all symmetrised over the arbitrary A/B arm assignment and all **cross-baseline** (an A-referenced change is always paired with a B-referenced persistence, so no arm shares a baseline with itself):

| functional | definition | question |
|---|---|---|
| `T_probe` | `½[ρ(D_test−D_A, D_post−D_B) + ρ(D_test−D_B, D_post−D_A)]` | does the total task reorganisation persist? (the standard trace) |
| `T_encode` | `½[ρ(D_learn−D_A, D_post−D_B) + ρ(D_learn−D_B, D_post−D_A)]` | does the **encoding** reorganisation persist? |
| `T_probespec` | `½[ρ(D_test−D_learn, D_post−D_B) + ρ(D_test−D_learn, D_post−D_A)]` | does what the probe adds **on top of** encoding persist? |
| `T_probespec_pe` | `½[pr(f,p\|e) + pr(f,p2\|e2)]` | the same, with the encoding component partialled out |

`T_probe` is numerically identical to the incumbent `rho_sym`. The phase set is passed as **data** everywhere; adding or removing a phase requires no code change.

**`T_probespec_pe` is UNCALIBRATED and must not be gated.** It is a conditional (partial-correlation) statistic, and lane W0-C demonstrated that such a statistic can return a significantly positive result (p = 0.007) on a sham arc built entirely inside pre-task rest, where no consolidation can exist. Until a data-based placebo has been run for it, its cells are descriptive only.

---

## 3. The diffusion scale grid (locked)

`s = τ·λ_max ∈ logspace(0, log10 180, 16)` — from `s = 1` (`τ = 1/λ_max`, the fastest mode) to `s = 180`, 16 log-spaced points.

**Report per-scale, never best-scale.** No number produced on this substrate may be collapsed to a scale-maximum before it is reported. `s` is a *functional / topological* scale, not a physical distance: never label a scale micro / meso / macro without quoting `ℓ(s)` or a cluster count.

---

## 4. What chose this

### 4.1 The transform, and the one place the pre-registration was under-specified

R1 ranked the two transforms on split-half reliability of the FC matrix itself — a property of the estimator, independent of our hypothesis. `imcoh_sq` won: 5 of 6 bands, cohort-median Δ = −0.0078, Wilcoxon p = 5.0e-6 over 60 (patient, band) cells. The margin is ~1.2 % relative on reliabilities of 0.55–0.74. A post-hoc diagnostic ruled out the obvious confound: the win is not inherited from coherence magnitude (`Spearman(observed, band-averaged |C|)` differs by ≤ 0.005 in five of six bands).

Addendum R1-bis, written after A1 and **before any trace number**, said: run A2 under both transforms; if the verdict is transform-invariant, keep `abs` for continuity; if not, adopt `sq`. The verdict is **not** transform-invariant — but not in the direction R1-bis anticipated. Under `sq` the **admissible window collapses to a single fraction, 0.00 octaves**: α survives (7/10 shared configs, identical to `abs`) but **β does not** (6/10 under `abs` → 3/10 under `sq`; on `T_encode`, 4/10 → 1/10), and low_γ leaks more, so band-selectivity is worse too.

R1-bis's literal text would have me adopt a transform under which **no stable substrate exists at all**. That is an under-specification in my own pre-registration, and I am flagging it rather than quietly re-writing it. The resolution: **R1 is a preference ordering; R2 is a hard admissibility gate**, frozen before any number, and a (transform, backbone) pair that admits no octave-wide invariant region is not admissible whatever its reliability. `imcoh_sq` fails R2. `imcoh_abs` passes it at 1.51 octaves. So `imcoh_abs` is retained — **not because it won R1 (it lost) but because it is the only one of the two that admits a stable substrate.**

The honest consequence, which belongs in Methods: **the β trace is specific to `⟨|Im C|⟩` and does not survive `⟨(Im C)²⟩`.** That is a real fragility, not a nuisance, and it should be stated rather than discovered by a referee.

### 4.2 The backbone

The knob was swept, not chosen: 18 configurations spanning `mst-union` at 10 fractions from 0.02 to 1.00, a plain global threshold at 3 fractions (the mechanism contrast — same edge budget, no connectivity policy), the disparity filter at 3 α, TMFG, and percolation; × 6 bands × 4 functionals × 16 scales × 10 patients, with matched-strength surrogate draws **shared across configurations within a cell** so that config-to-config differences are attributable to the backbone and not to surrogate noise.

A fraction is **admissible** when both signal bands are positive and the null band is not, on the same graph:

| f | density | θ | α | β | low_γ | admissible (θ) | admissible (θ + low_γ) |
|---|---|---|---|---|---|---|---|
| 0.02 | 0.028 | 0.311 | 0.064 | **0.025** | **0.020** | | |
| 0.04 | 0.044 | 0.278 | **0.022** | 0.094 | 0.173 | | |
| 0.05 | 0.053 | 1.000 | **0.008** | 0.082 | 0.062 | | |
| 0.07 | 0.072 | 1.000 | **0.009** | **0.022** | **0.013** | ✓ | |
| 0.10 | 0.101 | 0.286 | **0.028** | **0.031** | **0.022** | ✓ | |
| 0.14 | 0.140 | 0.393 | **0.017** | **0.026** | 0.111 | ✓ | ✓ |
| 0.20 | 0.200 | 0.426 | **0.018** | **0.004** | 0.098 | ✓ | ✓ |
| 0.28 | 0.280 | 0.250 | 0.166 | **0.004** | 0.091 | | |
| 0.40 | 0.400 | 0.402 | 0.093 | **0.003** | 0.068 | | |
| 1.00 | 1.000 | 0.322 | 0.120 | **0.016** | 0.061 | | |

`T_probe` cluster p; bold = p < 0.05. **Admissible window `f ∈ [0.07, 0.20]`, contiguous, 1.51 octaves** — passes the pre-registered ≥ 1-octave bar. Adding low_γ to the required-null set shrinks it to `f ∈ [0.14, 0.20]`, **0.51 octaves — below the bar**.

Over the wide window all four fractions agree: α 4/4 (median p = 0.017, margin +0.125), β 4/4 (p = 0.024, +0.127), θ 0/4 (p = 0.409, −0.012).

**α is recoverable.** Its plateau is `f ∈ [0.04, 0.20]`, 2.32 octaves, positive at 6/10 fractions — *invariantly-positive*. β's is wider, `f ∈ [0.07, 1.00]`, 3.84 octaves, 8/10 fractions, holding all the way to the dense graph. θ is *invariantly-null* across the entire 5.64-octave grid.

**No parameter-free filter qualifies.** TMFG: α p = 0.141, β p = 0.086. Percolation: α p = 0.312, β p = 0.076. Neither clears on `T_probe` in any band, so neither sits inside the plateau and R2.1 fails. The free parameter therefore cannot be eliminated — only integrated over, which is what R2.3 and `canonical_graph_ensemble` do.

**Planarity, not density, is what kills α.** At matched density (TMFG 0.050 vs mst@0.05 0.053, ratio 0.947) the α margin falls from +0.202 to +0.043 — a 4.7× drop — and p from 0.0055 to 0.134. The two graphs keep the same number of edges and differ only in which ones, so the loss is the topological prior, not sparsity.

---

## 4b. What this substrate does and does not license

This section exists because it is the single most likely thing for a referee to attack, and because sibling lane W0-B established the facts that make the old wording indefensible.

**Matched-strength is an *independence* null, not a lag null.** W0-B showed it agrees cell-for-cell with an independent segment-lattice null (r = 0.983 across 800 cells), so it tests whether the cross-phase structure could arise from independently-drawn graphs with the same strength sequence. It says nothing about whether the structure is carried by coupling *magnitude* or by *lag*. Every verdict in §4.2 is therefore an independence verdict.

**The trace fails the lag-destroying nulls.** With the backbone edge set pinned to the observed `|ImCoh|` one — the control built to give the trace its best chance — W0-B found a cohort median margin of **−0.003 for β and +0.004 for α**, against +0.214 and +0.130 under matched-strength. No separation, at any scale, in any band. The α and β traces clear independence nulls and survive order-preserving session repartitioning (β 10/16 scales, best p = 0.0049, not BH-significant, min q = 0.186), but they do not clear a null that preserves coherence magnitude and destroys only lag.

Two rules follow, and they bind on every lane inheriting this contract:

1. **Every cohort claim needs three rungs**, reported together: matched-strength / N1b (independence), **N1 with the backbone fixed** (lag), and **N3 order-preserving** (session structure). A claim resting on matched-strength alone is an independence claim and must be labelled as one.
2. **Where N1 is not cleared — currently α, β, and all four five-phase functionals — the claim must be stated as being about *coherence structure*, not about lagged interaction, and volume-conduction immunity must not be offered as its justification.** The substrate is `|ImCoh|`, which is *constructed* to be VC-immune; but a result that does not separate from a magnitude-preserving, lag-destroying surrogate has not demonstrated that its content is the lagged part. The conservative construction is a property of the estimator, not evidence about the finding.

**No functional may ever be tested against zero.** W0-B found that none of the four five-phase functionals is zero-centred: all four, `T_probe` included, are significantly positive at 16/16 scales for α and β on a sham arc carved from a single resting recording with temporal order destroyed. This does not void surrogate-referenced results — the surrogate inherits the same construction — but it does mean the only admissible statistic is a **margin against a per-configuration surrogate**. Everything in this contract is computed that way (§5); anything downstream must be too.

---

## 5. The cohort gate used here

One **sign-flip cluster-mass test** (Maris-Oostenveld; whole-patient-profile flips, which preserves the along-axis correlation) per (config, band, functional) over the entire 16-scale margin profile. This collapses the swept axis to a single test, so it does not inherit the per-scale multiplicity problem at all.

That choice is forced by a measurement, not a preference: on these margin matrices the 16-scale axis is worth **n_eff = 1.51** independent tests (mean cross-scale correlation +0.79), so treating the scales as 16 independent tests over-charges by ~11×. Sibling lane W0-C reached the same conclusion independently on a 28-point sweep (n_eff 1.2–1.9, mean r +0.70 to +0.93).

Everything is gated on the **margin** `obs − surr_p50`, never a raw observed value, because per-patient observed values co-vary with the height of that patient's own null.

Per-scale cleared-counts are reported as **secondary** and under two families side by side (raw α = 0.05, and whole-grid BH). Verdict agreement with the primary gate is 0.729 and 0.750 respectively over 432 cells — the family matters, which is precisely why the plateau is defined on the family-free test.

**Reconciliation note for integration:** this lane's gate is the axis-cluster test, which corresponds to W0-C's `axis_cluster_gate`. W0-C's locked default is whole-grid BH across the (band × scale) family, with the cluster test offered as the alternative that respects axis smoothness. The gate was deliberately **not** swapped mid-sweep, so every row here is internally comparable. Where this contract and W0-C's ledger disagree on a per-cell verdict, the difference is the multiplicity family and not the data.

---

## 6. What is NOT settled by this contract

1. **The trace does not clear a lag-destroying null, and that is now the largest open item in the project.** Matched-strength is an independence null (§4b); W0-B's N1 with the backbone pinned gives a cohort median margin of −0.003 (β) and +0.004 (α). Every verdict here is an **independence** verdict about **coherence structure**. Whether any of it is carried by lagged interaction is unresolved, and §4b's wording rules are mandatory until it is.
2. **The full band-selectivity claim is knob-dependent.** low_γ is null only at `f ≥ 0.14`, giving 0.51 octaves — below the pre-registered bar. Any claim of the form "the hierarchy rejects low_γ" must be reported with that window and that width, or not made.
3. **The free parameter is integrated over, not eliminated.** No parameter-free filter passed. A referee can still ask why mst-union; the answer is the surface in §4.2, not a principle.
4. **`T_probespec_pe` is uncalibrated** (§2) and cannot support a claim until a data-based placebo exists.
5. **`T_probespec` (inference-specific) is invariantly null for α across the whole grid** (0/10 fractions), and for low_γ and high_γ too. Whatever the inference-specific story becomes, it is not a knob-robust cross-phase trace in those bands.
6. **δ and high_γ are knob-dependent** and no verdict is offered for them here.
7. **The transform comparison rests on two candidates.** `⟨Im C⟩` (signed) and other coherence measures were not run.
8. **PMFG was evaluated observed-only**; it is O(N³) planarity testing and does not sit inside a 100-surrogate null, so the planar arm establishes agreement of the readout, not that PMFG clears a null.
9. **n = 10.** The cohort signed-rank p-floor is 1/1024, and the cluster test inherits the same ceiling. Nothing here has the headroom for a heavily corrected family.

---

## 7. Provenance

- Pre-registration (frozen before any number, with the 5-point critical preamble): `scripts/01_compute/paper_final/w0a_00_preregistration.md`
- Transform head-to-head: `scripts/01_compute/paper_final/w0a_01_transform_headtohead.py` → `data/paper_final/w0a_substrate/a1_transform/`
- Split-half cache, both transforms: `w0a_01b_cache_halves_both_transforms.py` → `IMCOH_HALVES_CACHE`. Reproduces the legacy `imcoh_halves_fc` cache to ≤ 3.0e-8.
- Stability sweep: `w0a_02_sparsification_stability.py` → `data/paper_final/w0a_substrate/a2_stability/{abs,sq}/`. R = 100 surrogates, 18 configs (abs) / 10 (sq).
- Surface analysis: `w0a_03_plateau_analysis.py` (10 000 sign-flip permutations).
- Parameter-free filters: `w0a_04_parameter_free_filters.py` → `data/paper_final/w0a_substrate/a3_parameter_free/`
- Equivalence checks, on real FC: `w0a_verify_functionals.py`, `w0a_verify_canonical_graph.py`.
- Library: `src/lrg_eegfc/workflow/substrate.py`, `utils/fc/heat_multiscale.py`, `utils/fc/split_half.py`, `utils/fc/backbone.py`, `utils/metrics/graph_descriptors.py`, `config/paths.py`.
