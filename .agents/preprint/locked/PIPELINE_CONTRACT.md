---
name: pipeline-contract
kind: locked-contract
era: PAPER_FINALIZATION (Wave 0, lane W0-A)
status: DRAFT — pending the A2 stability surface
created: 2026-08-25
scope: The locked substrate every analysis lane inherits — the ImCoh magnitude transform, the sparsified backbone (or the knob-integrated rule that replaces a single fraction), the phase set, and the diffusion-scale grid. States the pre-registered rule that chose each, the evidence, and an explicit list of what is NOT settled.
supersedes: the bare "mst@0.20" choice recorded in .agents/guides/02_methods/sparsification-choice.md §4a
companion: .agents/reports/2026-08-25_w0a-substrate-contract.md
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
pointers:
  - scripts/01_compute/paper_final/w0a_00_preregistration.md
  - src/lrg_eegfc/workflow/substrate.py
  - data/paper_final/w0a_substrate/
---

# Pipeline contract — the locked substrate

## Head

PLACEHOLDER-HEAD

---

## 1. The contract

PLACEHOLDER-CONTRACT-TABLE

### How to use it

Every lane calls one function. Do not re-derive the graph inline; do not hardcode a phase tuple; do not hardcode a fraction.

```python
from lrg_eegfc.workflow.substrate import (
    CANONICAL, canonical_graph, canonical_phase_eigs, canonical_scale_grid,
)
from lrg_eegfc.utils.fc.heat_multiscale import cross_phase_functionals_over_scales

eigs = canonical_phase_eigs("Pat_02", "beta")          # five phases, from CANONICAL_PHASES
T    = cross_phase_functionals_over_scales(eigs, canonical_scale_grid())
# -> {"T_probe", "T_encode", "T_probespec", "T_probespec_pe"}, each length 16
```

A four-phase caller passes `phases=("A", "B", "task_test", "rest_post")` to the same function and gets `{"T_probe"}` alone — the standard trace, identical to the incumbent `rho_sym_over_scales` to 1e-16. A raw-FC baseline is `canonical_graph(..., dense=True)`: a visible argument, not a different code path.

---

## 2. The phase set (locked)

Five phases, because the paradigm is transitive inference and separating **learning a structure** from **applying it** is part of the science, not an optional extra:

| role | phase | what it is |
|---|---|---|
| `baseline_a` | `A` | first contiguous half of `rest_pre` |
| `baseline_b` | `B` | second contiguous half of `rest_pre` |
| `encode` | `task_learn` | ordered premises presented (A>B, B>C, …) |
| `probe` | `task_test` | novel non-adjacent pairs judged — answerable only by inference |
| `follow` | `rest_post` | resting state after the task |

The four cross-phase functionals, all symmetrised over the arbitrary A/B arm assignment and all **cross-baseline** (an A-referenced change is always paired with a B-referenced persistence, so no arm shares a baseline with itself):

| functional | definition | question |
|---|---|---|
| `T_probe` | `½[ρ(D_test−D_A, D_post−D_B) + ρ(D_test−D_B, D_post−D_A)]` | does the total task reorganisation persist? (the standard trace) |
| `T_encode` | `½[ρ(D_learn−D_A, D_post−D_B) + ρ(D_learn−D_B, D_post−D_A)]` | does the **encoding** reorganisation persist? |
| `T_probespec` | `½[ρ(D_test−D_learn, D_post−D_B) + ρ(D_test−D_learn, D_post−D_A)]` | does what the probe adds **on top of** encoding persist? |
| `T_probespec_pe` | `½[pr(f,p\|e) + pr(f,p2\|e2)]` | the same, with the encoding component partialled out |

`T_probe` is numerically identical to the incumbent `rho_sym`. The phase set is passed as **data** everywhere; adding or removing a phase requires no code change.

---

## 3. The diffusion scale grid (locked)

`s = τ·λ_max ∈ logspace(0, log10 180, 16)`, i.e. `s` from 1 (τ = 1/λ_max, the fastest mode) to 180, 16 log-spaced points. The incumbent grid, unchanged, so this contract is comparable to every artifact built on it.

**Report per-scale, never best-scale.** No number produced on this substrate may be collapsed to a scale-maximum before it is reported. `s` is a *functional / topological* scale, not a physical distance: never label a scale micro / meso / macro without quoting `ℓ(s)` or a cluster count.

---

## 4. What chose this, in one paragraph per choice

PLACEHOLDER-RATIONALE

---

## 5. What is NOT settled by this contract

PLACEHOLDER-OPEN

---

## 6. Provenance

PLACEHOLDER-PROVENANCE
