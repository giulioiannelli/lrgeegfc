---
name: methods-displacement-taxonomy
era: IMCOH_ABS_COHORT_N10
status: draft
kind: methods
scope: per-pair cross-phase displacement portrait decomposing the validated ρ_split^coph into an anchor/reset/trace/random/rest-drift per-band character — DESCRIPTIVE enrichment of the locked "no trace" verdicts, NOT a new probe or a verdict change
---

# Cross-phase displacement taxonomy — decomposing ρ_split into anchor / reset / trace / random

**Head + limitations (read first).** This is a **descriptive decomposition** of
the already-validated per-pair measure `ρ_split^coph` (audit_63,
matched-strength resilient) — not a new probe, not a verdict change. It splits
the binary "trace / no-trace" into a per-(patient, band) reading on two axes:
*did the dendrogram reorganize* (displacement magnitude `m_task`, `m_rest` in
baseline-noise units) and *did the change persist* (`ρ_split`, the validated
quantity). **Unverified caveats, stated up front per project rule:**

1. The persistence axis (`ρ_split`) is matched-strength-validated; the
   **magnitude axis (`m_task`) is referenced only to the baseline split-half
   noise floor `‖D_coph(preA) − D_coph(preB)‖`, NOT to a matched-strength
   null.** "Task reorganized the tree" is therefore baseline-referenced, not
   strength-controlled.
2. Category labels are **soft-threshold** (the `m ≈ 1` cutoff wobbles at the
   margins — e.g. Pat_02 β flips trace↔weak-trace on `m_rest = 0.98`). The
   robust content is the two continuous numbers, not the bin.
3. Cohort "character" is a **plurality, not unanimity** (θ = 5/10 random).
4. It is a **whole-tree per-patient reading, not a per-pair partition.** The
   trace is a broad, low-amplitude, sign-coherent drift across the *bulk* of
   pairs (diagnostic 2026-06-03: `ρ` among the bottom-50 % movers ≈ `ρ_full`,
   while among the top-10 % movers it is often lower or negative). There is no
   clean subset of "trace pairs."

It therefore does **not** enter `locked/VERDICT_LEDGER.md` or the per-band
briefs as fact, and changes no verdict. The per-band trace verdicts stand
exactly as locked. What this adds is a *mechanism* for the existing verdicts —
especially a positive identity for the two locked "no trace" bands. Promotion
path in §5.

---

## 1. What it decomposes

`ρ_split^coph` is the Spearman rank correlation of a scatter whose every dot is
**one pair** `(i, j)`, split-baseline (audit_63):

```
x = Δ_task(i,j) = D_coph_task(i,j)  − D_coph_preA(i,j)
y = Δ_rest(i,j) = D_coph_post(i,j)  − D_coph_preB(i,j)
ρ_split        = Spearman(x, y)              # the validated tilt of the cloud
```

The dots that *make* the validated number already carry the information; the
taxonomy reads it instead of only its tilt. Pipeline is bit-exact with audit_63
(`load_phase_fc` + `lrg_ultrametric_condensed`, τ = 1/λ_max imported verbatim);
`ρ_split` reproduces the locked `matched_strength_surrogate_split_baseline/
per_patient_per_band.csv` to < 1e-3 for every cell (hard integrity gate in the
script).

## 2. The three quantities (per patient, per band)

```
m_task = median|Δ_task| / median|D_coph(preA) − D_coph(preB)|   # reorg vs baseline noise
m_rest = median|Δ_rest| / median|D_coph(preA) − D_coph(preB)|   # rest displacement
ρ_split                                                          # persistence (validated)
```

Descriptive reading (soft, not a gate):

| reading | condition | plain meaning |
|---|---|---|
| **anchor** | `m_task < 1` and `m_rest < 1` | task barely moved the tree |
| **reset** | `m_task ≥ 1`, `m_rest < 1` | task moved it, rsPost returned to baseline |
| **rest-drift** | `m_task < 1`, `m_rest ≥ 1` | task quiet, rsPost drifted on its own |
| **trace** | both `≥ 1`, `ρ_split ≥ 0.15` | reorganized and persisted |
| **random** | both `≥ 1`, `ρ_split < 0.15` | reorganized but rsPost went elsewhere |
| **weak-trace** | low magnitude but `ρ_split ≥ 0.15` | faint, rank-coherent drift |

## 3. Per-band cohort character (DRAFT — descriptive)

Per-patient readings tallied across the n = 10 cohort (data:
`data/audit/pair_displacement_taxonomy/per_patient_<band>.csv`):

| band | persist (trace+weak) | character | enriches locked verdict |
|---|---|---|---|
| **β** | **7/10** | **TRACE** (4 trace + 3 weak; dissent = 1 random, 1 anchor, 1 rest-drift) | confirms `strong trace` |
| **θ** | 1/10 | **RANDOM** (5 random + 2 anchor) | explains `no trace` → transient-scattered, not restorative |
| **γ_h** | 4/10 | **BIMODAL** (3 coherent trace vs 2 anti + 2 reset + 2 random) | explains `no trace` → cohort cancellation of large opposing effects |
| α | 4/10 | mixed, reset-leaning (3 reset) | `strong trace, only D_coph` carried by 2 trace + 2 weak |
| δ | 4/10 | mixed (3 rest-drift, 2 reset; 1 massive reorganizer Pat_06 m=13) | `weak trace, only Grassmann` |
| γ_l | 5/10 | split trace/random (4 random; 1 massive reorganizer Pat_05 m=18) | `strong trace, only Grassmann` |

**The two clean cohort characters** are β (TRACE) and θ (RANDOM). θ is the
informative one: in 7/10 patients task *did* reorganize the θ tree
(`m_task > 1`), but the change did not persist (`ρ ≈ 0`) and rsPost moved to an
*unrelated* configuration (5 random) rather than back to baseline (only 1
reset). So the locked θ "no trace" has a positive identity — **task-induced θ
reorganization is real but transient and scattered, not homeostatic**. γ_h's
"no trace" is the opposite failure: large per-patient effects pointing every
direction (3 strong coherent traces, e.g. Pat_05 ρ = 0.92, vs Pat_15 ρ = −0.47
and clean resets Pat_10/Pat_13 with `m_task` 7/4 snapping back), cancelling to a
cohort null.

## 4. β-outlier resolution (the payoff)

The three β patients that fail `ρ_split^coph` get a mechanism each, consistent
with the cross-probe picture:

| patient | m_task | ρ_split | reading | cross-probe note |
|---|---|---|---|---|
| Pat_10 | 1.58 | −0.09 | **random** | reorganized but scattered → caught by Grassmann subspace, missed per-pair |
| Pat_14 | 0.93 | −0.05 | **anchor** | task barely moved the tree |
| Pat_15 | 0.47 | +0.08 | **rest-drift** | task quiet (R-hemi-only implant), rsPost drifts on its own |

This dovetails with `methods_grassmann_cluster_extent.md`: Pat_10's β trace
lives in the leading-subspace rotation, not in per-pair merge-height ranks —
exactly the "random at the pair level, coherent in the subspace" signature.

## 5. Relationship to the canonical taxonomy + promotion path

- **trace / anchor / reset** map directly onto the project's canonical
  `trace / anchor / reset / emergent` taxonomy
  (`.agents/guides/01_project/terminology.md`). **random** and **rest-drift**
  are pair-field extensions of the continuous decomposition (rsPost reorganizes
  orthogonally to task / moves only in rest). **emergent** is a module-level
  (subtree) phenomenon, orthogonal to this pair-level reading.
- **To promote θ = RANDOM and γ_h = BIMODAL into the briefs** (`bands/04_theta.md`,
  `bands/05_gammah.md`, `bands/00_cohort.md`) as enrichments of the locked
  "no trace" verdict, the **magnitude axis must be validated under
  matched-strength** — i.e. show that `m_task` exceeds what a strength-preserving
  surrogate produces, so "task reorganized the tree" is strength-controlled, not
  baseline-referenced. Then add dated revision entries (the verdicts themselves
  do not change — only their narrative gains the mechanism). Until that null is
  run, this stays **draft** and outside the locked artifacts.

## 6. Provenance

- Scope report: `.agents/guides/task-persistence-investigation/2026-06-01_cross-phase-displacement-taxonomy.md`
- Script: `scripts/01_compute/audit/audit_76_pair_displacement_taxonomy.py` (bit-exact with audit_63; integrity-gated)
- Data: `data/audit/pair_displacement_taxonomy/` — `figures/<band>_pair_clouds.pdf`, `per_patient_<band>.csv`
- Companion report: `.agents/reports/2026-06-01_trace-concordance-vs-blind-fc.md`
- Source measure: `data/audit/matched_strength_surrogate_split_baseline/` (audit_63, `ρ_split^coph`)
