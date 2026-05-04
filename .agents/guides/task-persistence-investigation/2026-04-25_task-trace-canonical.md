---
name: task-trace-canonical
type: scope
era: COHORT_N10
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - .agents/reports/2026-04-24_h1-h4-vi-results.md
  - src/lrg_eegfc/utils/metrics/tree.py
  - src/lrg_eegfc/utils/metrics/tree_distance.py
  - src/lrg_eegfc/utils/metrics/hypothesis.py
---

# Task-trace canonical reformalization (n=10, IMCOH_ABS)

**Single canonical statement of the task-trace question — node-role
regimes (Persistence / Trace / Reset / Rearrange) defined at the
leaf-set level, claim resolution at `(band × scale)` cohort level
with the bar `≥ 8/10` patients sign-correct, and the central claim
explicitly stated as a *residual* assertion (rest_post differs from
rest_pre **in the direction of** task_test) rather than a global
"rest_post becomes task-like" claim. Supersedes the fragmented
`2026-04-24_multiscale-task-trace.md` framing now that Pat_14 is
restored. Reuses MRL primitives (`J_min = 0.9`, `k_min = 3`, native
heights) so cross-checks with MRL/CBR are straightforward.**

---

## 1. Notation

Indices and sets:

| Symbol | Domain | Meaning |
|---|---|---|
| `p ∈ P` | `\|P\| = 10` | patient — Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15 |
| `b ∈ B` | `\|B\| = 6` | band — δ, θ, α, β, γ_l, γ_h |
| `φ ∈ Φ` | `\|Φ\| = 4` | phase — `rest_pre`, `task_learn`, `task_test`, `rest_post` |
| `L_p` | `\|L_p\| = N_p` | leaf-index set for patient `p` (113 for Pat_10; raw N elsewhere) |
| `Z^{p,b,φ}` | UPGMA linkage matrix | one per `(p, b, φ)` under `imcoh_abs` |
| `T^{p,b,φ}` | binary tree | from `Z^{p,b,φ}`; `N_p` leaves, `N_p − 1` internal nodes |
| `Internal(T)` | set of internal nodes | excludes leaves and root for candidacy |
| `leaves(v)` | `2^{L_p}` | leaf-set of internal node `v` |
| `size(v)` | `ℕ` | `\|leaves(v)\|` |
| `h_rel(v)` | `(0, 1]` | fractional native merge height — only scale variable |

Match primitives (identical to MRL):

| Symbol | Type | Definition |
|---|---|---|
| `J(A, B)` | `[0, 1]` | `\|A ∩ B\| / \|A ∪ B\|`; `0` if both empty |
| `J*(S, T)` | `[0, 1]` | `max_{u ∈ Internal(T)} J(leaves(u), S)` |
| `J_min` | `(0, 1]` | match threshold; **default `0.9`** (locked, matches MRL Q1) |
| `match(S, T; J_min)` | predicate | `J*(S, T) ≥ J_min` |
| `k_min` | `ℕ_{≥ 2}` | subtree-size floor; **default `3`** (locked, matches MRL) |

**No new symbols, no new τ, no new size floor.** Cross-checks with MRL
require identical primitives.

## 2. Operational regimes (P / T / R / RA)

For a candidate subtree `v` with leafset `S = leaves(v)` and
`size(v) ≥ k_min`, three Boolean predicates are evaluated against the
patient's three other LRG trees at the same band:

```
present_pre(v;  J_min) ⇔ match(S, T^{p,b,rest_pre};  J_min)
present_task(v; J_min) ⇔ match(S, T^{p,b,task_test}; J_min)
present_post(v; J_min) ⇔ match(S, T^{p,b,rest_post}; J_min)
```

The four regimes partition candidate space:

| Regime | Predicate combination | Meaning |
|---|---|---|
| **P (persistence)** | `present_pre ∧ present_task ∧ present_post` | anchored: same module across all phases |
| **T (trace)** ★ | `¬present_pre ∧ present_task ∧ present_post` | task-induced module that persists into rest_post |
| **R (reset)** | `¬present_pre ∧ present_task ∧ ¬present_post` | task-induced module that dissolves back |
| **RA (rearrange)** | otherwise | ergodic / no consistent role |

★ **T is the central scientific claim.**

The candidate set itself is **task-anchored** (matches MRL):

```
V_task(p, b) := { v ∈ Internal(T^{p,b,task_learn}) ∪ Internal(T^{p,b,task_test})
                  : k_min ≤ size(v) ≤ N_p − 1 }
```

V_task is a multiset (a leafset present as an internal node in both
task trees contributes twice). Rationale: a leafset that emerges in
both task phases is more "task-anchored" than one in only one. This
matches the MRL pooling decision.

### 2.1 Electrode-level lift

Electrode-level statements (e.g. "node X participates in T") are
**derived**, never primary:

```
electrode_in_T(ℓ; b, h_rel-window) :=
    ∃ v ∈ V_task(p, b) : ℓ ∈ leaves(v)  ∧  h_rel(v) ∈ window
                       ∧ regime(v) = T
```

Aggregation across patients is **count of patients in which `ℓ`
participates in T at `(b, window)`**, not a continuous score. The
sensor-space heatmap visual (Stage 4 of the audit-and-recovery plan)
renders this aggregate.

## 3. Claim levels

The hypothesis question collapses three *resolution levels* that prior
reports occasionally conflated. Reformalised here:

| Level | Object | Operationalization | Cohort-claim threshold |
|---|---|---|---|
| **Subtree** | one `(p, b, v)` row | `regime(v) = T` | n/a (raw observation) |
| **Cell** | one `(b, h_rel-bin)` | fraction of patients with ≥ 1 T subtree in cell | **`Π(b, bin) ≥ 0.8` (8/10 patients)** |
| **Band** | one `b` | non-empty set of cohort-wide cells in band | descriptive band ranking, no scalar gate |

The `(band × scale)` map (Stage 2 of audit-and-recovery) is the
**cell-level** claim. The "Pure cohort-wide multiscale band-specific
trace" claim of the project lives at the **band level**, sustained by
the cell-level pattern.

**No FDR over bands × hypotheses.** The `2026-04-24_post-mortem-scalar-session.md`
established that scalar gates with `q < 0.05` FDR `m = 6` are
incompatible with `n = 10` effect sizes (`r_rb ≈ 0.6`). The cell-level
unanimity bar `≥ 8/10` is the operative gate; cluster-permutation in
the `k`-axis is reported as a secondary, descriptive overlay (matches
the `cluster_stats` helper already in the library).

## 4. Residual vs global — the H2-RAW failure

The historical confusion: a casual reading of the trace claim is
"rest_post becomes task-like". The data refute this — H2-RAW
(direct correlation of upper-triangle ultrametric distances) and
H2-FROB (Frobenius distance ratio) **fail at every band** (0/6 q < 0.05).
In several bands rest_post is *farther* from task than from rest_pre
(`r > 1` in H2-FROB).

The correct framing is **residual**:

> The pattern of rest_post − rest_pre alignment with task − rest_pre
> alignment, restricted to **specific subtrees / specific scales /
> specific bands**.

Equivalently in MRL/CBR terms: the trace observable counts subtrees
that newly form during task and persist into rest_post. The subtree
count can be cohort-wide and band-specific while the *global*
geometry of rest_post stays anatomy-dominated (and therefore close to
rest_pre). The two are not contradictory once the residual nature of
the claim is named.

## 5. n=10 cohort — what changes vs n=9

Pat_14 was excluded from COHORT_N9 (locked 2026-04-22) because its
`task_test` was corrupt at vendor import. Vendor replaced 2026-04-25;
the timeseries is now valid. Cohort returns to n=10 (Pat_14 back).

**Implications for prior n=9 results:**

1. The cluster-permutation results in `2026-04-24_multiscale-task-trace.md`
   (δ k=23–31 p=0.014, α k=2–4 conditional-H p=0.050) are **n=9**.
   They must be re-evaluated at n=10 before being cited. If Pat_14
   sign-aligns, the bar tightens to `≥ 8/10` (~78% same as 7/9). If
   Pat_14 cuts against, signal may weaken — report honestly.
2. `H2c` and `H2d` were universal at all 6 bands (n=9). Pat_14
   inclusion is unlikely to flip them globally but per-cell unanimity
   counts must be recomputed.
3. `H2a` was demoted to supplementary at n=9 (2/708 unanimous cells).
   Pat_14 at n=10 is unlikely to rescue it; H2a stays supplementary.
4. `KC / MC / wRF` n=5-only artefacts are unaffected; they remain
   `stale-numeric` and must not be cited as current.
5. The legacy `mean_contrast_heatmap_4panel.pdf` (n=5) stays archived.

## 6. Properties

| # | Statement |
|---|---|
| 6.1 | The four regimes P/T/R/RA partition `V_task(p, b)` exactly: `(present_pre, present_task, present_post) ∈ {0,1}³`, with `present_task = 1` always (candidates are drawn from task trees). 4 of the 8 cells = P/T/R/RA; the other 4 are vacuous (`present_task = 0`). |
| 6.2 | `T` is **conditional / residual**, not global. By construction `T` is silent on overall rest_post ↔ rest_pre similarity. |
| 6.3 | Range of cell-level statistic: `Π(b, bin) ∈ [0, 1]`. Cohort-wide cell ⇔ `Π ≥ 0.8`. |
| 6.4 | Threshold-monotonicity: `T` is **non-monotone** in `J_min` (matches MRL property 5.4). The conjunction `¬present_pre ∧ present_post` flips opposing predicates. Sensitivity sweep `J_min ∈ {0.85, 0.9, 0.95}` is informative. |
| 6.5 | Independence of bands: every analysis is per-band; no cross-band pooling at the cohort-claim level. |
| 6.6 | Anchored entirely in existing primitives — `tree_internal_nodes`, `jaccard_leafsets`, `dmax_from_Z`, `h_log_grid` (plot-time only), `cluster_stats` (descriptive overlay only). |

### 6.7 What this reformalization does NOT measure

- **Global tree similarity** between phases — failed direction
  (H2-RAW / H2-FROB / H2a-topo, see `2026-04-24_multiscale-task-trace.md` §7c).
  Out of scope by design.
- **Causal direction** of the trace (is rest_post "more like task"
  because of task, or because of session drift?) — H2e split-half
  drift floor at n=9 says drift `≈ 0`, but a residual H2e at n=10
  with Pat_14 included is needed before any causal language.
- **Statistical-significance gate** at FDR `q < 0.05 m = 6` — closed
  loop, see post-mortem (archived).
- **Disruption** (rest_pre modules that fail in rest_post) — symmetric
  observable, dropped per the same Q4 directive that retired the MRL
  `M̄_−` mirror.
- **task_test-specific** vs `task_learn`-specific provenance — pooled
  candidate set conflates them; we cannot, at this n, separate
  learning carry-over from test-specific module formation.

## 7. Caveats & failure modes

### 7.1 Same-probe bias at coarse `h_rel`

Same-probe MSC bias is mostly mitigated under `imcoh_abs` (≥ 50%
reduction per `probe-bias-guide.md`), not eliminated. At
`h_rel ≳ 0.7`, merges follow probe geometry. Mitigation: the
`(band × h_rel)` map hatches the `h_rel ≥ 0.7` strip; the cohort-wide
cells under load-bearing claim must lie at fine/mid scales
(`h_rel ∈ [0.05, 0.5]`).

### 7.2 Anatomy-dominated rest invariance

`T_rest_pre` and `T_rest_post` share anatomy — `J*(S, T_rest_pre) ≈
J*(S, T_rest_post)` for most leafsets `S`. The conjunction
`¬present_pre ∧ present_post` is therefore a **small-numerator**
observable. Cell-level `Π ≥ 0.8` is hard to fake.

### 7.3 Pat_03 outlier (1024 Hz)

`h_rel` normalisation handles scale comparability. Sensitivity
analysis: report `Π_with_Pat_03` and `Π_without_Pat_03` per cell;
flag any cell where the difference exceeds `0.1` — Pat_03 is
load-bearing in that cell and the claim is fragile.

### 7.4 Pat_14 newly restored

Pat_14 has been excluded since 2026-04-22 (n=9 lock). Spot-check
required: open one Pat_14 dendrogram per band; verify `dmax`,
`h_rel` distribution, and `tree_internal_nodes` count are
comparable to other patients before trusting Pat_14 contributions.
Flag in the audit-and-recovery handoff.

### 7.5 k-artefact at extremes

VI(k) per-cut analyses get eaten at `k ≈ N` (singleton-dominated)
and at `k ≈ O(1)` (giant-cluster-dominated). The user's stated
complaint. Reformalization here uses **native heights** for
candidate enumeration (not integer-k cuts), eliminating the VI(k)
extremes problem at the *subtree-level* claim. The k-axis re-enters
only as a **descriptive overlay** in the `(band × k)` map (Stage 2
of audit-and-recovery), where a separate companion diagnostic
(`scripts/01_compute/diagnostics/diag_k_artefact.py`) reports
`n_eff(k)`, singleton-fraction, max-cluster-fraction so a reader can
visually discount extremes. **No headline-figure masking** — that
was the user's explicit choice.

### 7.6 Threshold sensitivity (`J_min`)

Primary `J_min = 0.9` matches MRL. Sensitivity sweep
`{0.85, 0.9, 0.95}` reported as supplementary panels. Cells where
the band ranking flips across `J_min` are unreliable and flagged.

## 8. Pseudocode

```
Input:
  P, B                              # cohort, bands
  Φ = {rest_pre, task_learn, task_test, rest_post}
  J_min ∈ {0.85, 0.9, 0.95}         # primary 0.9
  k_min = 3
  H = (h_0, …, h_K), K = 12         # log-spaced viz bin axis (h_log_grid)
  fc_method = "imcoh_abs"
  cohort_threshold = 0.8            # 8/10

# === Per-patient pass ===
For each (p, b) ∈ P × B:
  Z_pre   ← load_lrg_result(p, "rest_pre",   b, fc_method).linkage
  Z_TL    ← load_lrg_result(p, "task_learn", b, fc_method).linkage
  Z_TT    ← load_lrg_result(p, "task_test",  b, fc_method).linkage
  Z_post  ← load_lrg_result(p, "rest_post",  b, fc_method).linkage

  V_pre   ← [u ∈ tree_internal_nodes(Z_pre)  : size(u) ≥ k_min]
  V_post  ← [u ∈ tree_internal_nodes(Z_post) : size(u) ≥ k_min]

  V_task ← []
  for (Z, src) ∈ {(Z_TL, "TL"), (Z_TT, "TT")}:
    for v ∈ tree_internal_nodes(Z):
      if k_min ≤ size(v) ≤ N_p − 1:
        v.source ← src
        V_task.append(v)

  for v ∈ V_task:
    S      ← leaves(v)
    J_pre  ← max(jaccard_leafsets(S, leaves(u)) for u in V_pre)
    J_post ← max(jaccard_leafsets(S, leaves(u)) for u in V_post)
    for J_min ∈ {0.85, 0.9, 0.95}:
      pres_pre  ← (J_pre  ≥ J_min)
      pres_post ← (J_post ≥ J_min)
      regime[J_min] ← classify(pres_pre, True, pres_post)
        # P: T T T,  T: F T T,  R: F T F,  RA: T T F
    emit row: (p, b, v.source, v.h_rel, v.size,
               J_pre, J_post,
               regime[0.85], regime[0.90], regime[0.95])

# === Cohort aggregation (plot-time) ===
For J_min ∈ {0.85, 0.9, 0.95}:
  For (b, bin_k) ∈ B × {1, …, K}:
    R       ← rows where row.b = b ∧ h_{k−1} ≤ row.h_rel < h_k
    P_e     ← unique row.p in R                           # eligible patients
    Π_T     ← (1/|P_e|) · |{p ∈ P_e : ∃ row ∈ R with row.p=p ∧ regime=T}|
    Π_P, Π_R, Π_RA computed analogously
    Cohort-wide cell flag: Π_T ≥ cohort_threshold

Output:
  data/audit/trace_modules/regime_per_node.csv
  data/audit/trace_modules/cohort_pi_by_regime.csv
  (figures produced by Stage 4 of audit-and-recovery)
```

## 9. Connection to prior tools

| Prior measure | Relation |
|---|---|
| H2c (`h2c_ultrametric_drift`) | continuous drift, **complementary** — the directional Spearman ρ on cophenetic-distance shifts at the **patient** level. Should pass when `T` is non-empty in the band; sanity check, not a substitute. |
| H2d (`h2d_coactivation_persistence`) | **subset** — the size-2-and-integer-k restriction of `T`. Captures "pair-level" trace; this reformalization captures arbitrary-size subtree trace. |
| H1-topo / H2a-topo (`h2_topology_directed`) | unconditional bipartition overlap. `T` adds the conditional `¬present_pre`; H1-topo is a coarser sibling, *not* a substitute. |
| KC / MC / wRF (`tree_distance`) | scalar tree distances integrating over scales. **Failed direction** for cohort claim at n=9 (post-mortem). Stays scalar, this reformalization is multi-resolution. |
| Partition-multiscale cluster-perm (`h2_partition_multiscale`) | k-resolved partition contrast. The δ k=23–31 (n=9, `Δ_VI` cluster-perm `p = 0.014`) and α k=2–4 (n=9, conditional H `p = 0.050`) cells should manifest as cohort-wide cells in the `(band × h_rel)` map at compatible scales (δ moderate `h_rel ≈ 0.2–0.4`; α near-root `h_rel ≈ 0.5–0.8`). Cross-check after Stage 2 of recovery. |
| MRL (`2026-04-25_module-retention-landscape.md`) | **identical primitives**, focused observable: `retained = T`. MRL reports `Π(b, bin)` for the `T` regime only, with `M̄_step / M̄_smooth` as continuous companions. This reformalization is the **superset** that names all four regimes; MRL is the "T-only" projection of it. |
| CBR umbrella (`2026-04-25_cbr-investigation.md`) | sibling family of single-module-matching variants. Each variant differs in anchor (rest_post vs task) and similarity function (Jaccard vs containment vs exact equality). All speak the same regime vocabulary; this canonical doc names the regimes once, CBR variants instantiate them with different operators. |

## 10. Implementation plan

This canonical doc is **definitional**. Two downstream consumers in
the audit-and-recovery plan:

1. **Stage 4 of audit-and-recovery** (`scripts/01_compute/audit/audit_15_trace_modules.py`):
   produces `data/audit/trace_modules/regime_per_node.csv` and
   `cohort_pi_by_regime.csv` per the pseudocode above. Reuses
   `tree_internal_nodes`, `jaccard_leafsets`, `dmax_from_Z`,
   `h_log_grid`. No new helpers.
2. **Stage 4 figure** (`scripts/01_compute/figures_embedded/fig_trace_modules.py`):
   cohort summary PDF + per-patient mini-dendrograms (the visual
   proof figure, mirrored on MRL Figure 2).

**Library promotion (deferred):** if a third consumer (e.g. CBR
cohesion variant) needs the regime classifier, promote
`(classify_regime, regime_membership)` into a new
`lrg_eegfc.utils.metrics.tree_persistence` module per the
coding-rules promotion-on-second-use principle.

## 11. Open questions

| Q | Note |
|---|---|
| Q1 — Should the cell `bin_k` axis be in `h_rel` (matches MRL) or in `k` (matches partition-multiscale)? | The audit-and-recovery Stage 2 figure uses `k` to mirror the n=5 4-panel layout; this canonical doc and Stage 4 use `h_rel` to match MRL/CBR. Both views are produced. Cross-pollination via the connection-to-prior-tools table in §9. |
| Q2 — What happens if `\|P_e(b, bin)\| < 8`? | Underpowered for the ≥ 8/10 claim. Cell flagged with grey-dotted border; not counted toward "cohort-wide" cells. Same convention as MRL §4. |
| Q3 — Permutation null at the regime level? | Deferred per MRL Q5. Descriptive only in v1. |
| Q4 — Should `task_test`-anchored variant be reported separately from `task_learn`-pooled? | Default pooled (matches MRL Q3). Sensitivity panel: `task_test`-only candidate set (`V_task` restricted to `T_TT` internals) — confirms pooling does not invent T modules. |

---

## Locked decisions

| # | Decision | Date |
|---|---|---|
| L1 | `J_min = 0.9` primary; sensitivity sweep `{0.85, 0.9, 0.95}` | 2026-04-25 |
| L2 | `k_min = 3` (matches MRL) | 2026-04-25 |
| L3 | Native heights at per-patient level; log-spaced bins for cohort viz only | 2026-04-25 |
| L4 | Pooled task candidate set `V_task = V(T_TL) ∪ V(T_TT)` | 2026-04-25 |
| L5 | Cohort-wide threshold `Π ≥ 8/10 = 0.8` | 2026-04-25 |
| L6 | No FDR / no scalar gate. Cluster-permutation reported as descriptive overlay only. | 2026-04-25 |
| L7 | T regime is the central claim; P/R/RA reported as auxiliary fields, not headline | 2026-04-25 |

## Pass criterion (descriptive)

The trace claim is sustained **at the band level** if, for some
`b ∈ B`, the `(band × h_rel)` `Π_T` map shows:

1. **Vertical persistence** — at least two adjacent `h_rel` bins with
   `Π_T ≥ 0.8` (multiscale, not a single isolated cell).
2. **Lying at non-probe-bias scale** — bins with mid-edge
   `h_rel < 0.7`.
3. **Robust across `J_min ∈ {0.85, 0.9, 0.95}`** — qualitative
   band ranking preserved.
4. **Robust to Pat_03 in/out** — `\|Π_with − Π_without\| < 0.1` per
   cell at `Π_T ≥ 0.5`.

If 1–4 hold for at least one band, the band-specific multiscale
trace claim has descriptive evidence at n=10. The headline
interpretation reads:

> *"Task-induced modules at fractional dendrogram heights `h_rel ∈ [a, b]`
> in band `X` persist into `rest_post` in `≥ 8/10` patients (Jaccard
> match floor `J_min = 0.9`). Equivalent re-statement at integer-k
> resolution: cohort-wide unanimity in the `(band × k)` `(band × h_rel)`
> map (Stage 2 of audit-and-recovery)."*
