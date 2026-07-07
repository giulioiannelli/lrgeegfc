---
name: module-retention-landscape
type: scope
era: COHORT_N9
status: superseded
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_mrl-vs-cbr-reconciliation.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - .agents/reports/2026-04-24_h1-h4-vi-results.md
  - .agents/plans/active/2026-04-25_surface-multiscale-trace.md
  - src/lrg_eegfc/utils/metrics/tree.py
  - src/lrg_eegfc/utils/metrics/tree_distance.py
---

> **Superseded 2026-04-25** by Cohesion-CBR for the per-patient
> discrete claim and by the (yet-to-be-written) continuous-trace-matrix
> scope for the cohort claim. See
> [`2026-04-25_mrl-vs-cbr-reconciliation.md`](2026-04-25_mrl-vs-cbr-reconciliation.md)
> for why MRL's hard-threshold Jaccard predicates filtered out the
> limbo zone where most real signals live, conflated cohesion with
> size match, and produced a scalar field that hid per-leaf identity.
> The TAM analysis showed the (J_pre, J_post) cloud sits on the
> diagonal at high Jaccard with the trace zone empty — the discrete
> trace claim does not hold under this operationalisation. The
> continuous H2c signal remains the genuine cohort-wide finding;
> Cohesion-CBR (`audit_12`, soft affinities + per-leaf colour) is the
> right tool for visualising per-patient discrete events. This scope
> kept for historical traceability.

# Module-Retention Landscape (MRL)

**Per-internal-node descriptive observable on the pooled task-phase trees
`T_task_learn ∪ T_task_test`. For every node `v` we report (a) a strict
0/1 indicator `retained(v) ⇔ leafset(v) is absent from rest_pre AND
present in rest_post` at Jaccard match floor `J_min = 0.9` (step), and
(b) the continuous Jaccard improvement
`ΔJ(v) := J*(v, T_rest_post) − J*(v, T_rest_pre)` (smooth). Both
observables live at the node's native fractional merge height
`h_rel(v) = h(v) / dmax(T)`. No synthetic query grid, no tolerance
window. Cohort views are histograms over a log-spaced binning axis
applied at plot-time only.**

---

## 1. Notation

Indices and sets:

| Symbol | Domain | Meaning |
|---|---|---|
| `p ∈ P` | `|P| = 9` | patient (`Pat_02, 03, 05, 06, 07, 08, 10, 13, 15`) |
| `b ∈ B` | `|B| = 6` | band (`δ, θ, α, β, γ_l, γ_h`) |
| `φ ∈ Φ` | `|Φ| = 4` | phase (`rest_pre, task_learn, task_test, rest_post`) |
| `L_p` | `|L_p| = N_p` | leaf-index set for patient `p` (113 for Pat_10; raw N for the rest) |

LRG output objects (per `(p, b, φ)`, all on `imcoh_abs`):

- `D^{p,b,φ} ∈ ℝ_{≥0}^{N_p × N_p}` — ultrametric distance matrix.
- `Z^{p,b,φ} ∈ ℝ^{(N_p − 1) × 4}` — UPGMA linkage matrix (scipy convention).
- `T^{p,b,φ}` — binary tree from `Z^{p,b,φ}`; `N_p` leaves, `N_p − 1`
  internal nodes.

For each internal node `v ∈ Internal(T)`:

- `leaves(v) ⊆ L_p` — leafset.
- `size(v) := |leaves(v)| ∈ {2, …, N_p}`.
- `h(v) := Z[i, 2]` — native merge height at the linkage row that produced `v`.
- `dmax(T) := max_v h(v) = Z[N_p − 2, 2]` — root height.
- `h_rel(v) := h(v) / dmax(T) ∈ (0, 1]` — fractional native height.
  **The only "scale" variable in this spec.** No synthetic query grid.

Functions:

| Symbol | Type | Definition |
|---|---|---|
| `J(A, B)` | `2^{L_p} × 2^{L_p} → [0, 1]` | `|A ∩ B| / |A ∪ B|`; `0` if both empty |
| `J*(S, T)` | `2^{L_p} × Trees(L_p) → [0, 1]` | `max_{u ∈ Internal(T)} J(leaves(u), S)` |
| `J_min` | `(0, 1]` | match threshold; **default `0.9`**; sensitivity `{0.85, 0.9, 0.95}` |
| `match(S, T; J_min)` | predicate | `J*(S, T) ≥ J_min` |
| `k_min` | `ℕ_{≥ 2}` | subtree-size floor; **default `3`** |

Cohort aggregators (plot-time only):

| Symbol | Type | Definition |
|---|---|---|
| `H = (h_0, …, h_K)` | log-spaced bin edges in `(0, 1]` | `h_0 = 0.05`, `h_K = 1`, `K = 12` |
| `bin_k` | `[h_{k−1}, h_k]` | k-th bin |
| `M̄_step(b, bin_k)` | `[0, 1]` | cohort mean of `retained(v)` over (p, v) in cell |
| `M̄_smooth(b, bin_k)` | `[0, 1]` | cohort mean of `max(ΔJ(v), 0)` over (p, v) in cell |
| `Π(b, bin_k)` | `[0, 1]` | proportion of patients with `≥ 1` retained `v` in cell (replaces `W` to avoid collision with Kendall's `W` in H4) |

**Symbols intentionally NOT reused**: `τ` (LRG diffusion time
`e^{-τL}`), `δ` (delta band), `ξ`/`Ξ`/`ε_ξ` (would imply a synthetic
query grid — there is none). `W` reserved for Kendall.

---

## 2. Predicates

### 2.1 Jaccard

```
J(A, B) := |A ∩ B| / |A ∪ B|     if A ∪ B ≠ ∅
        := 0                       if A = B = ∅
```

Symmetric, `J(A, A) = 1`, `J(A, B) = 0 ⇔ A ∩ B = ∅`. Implemented in
`lrg_eegfc.utils.metrics.tree.jaccard_leafsets`.

### 2.2 Best-match against a tree

```
J*(S, T) := max_{u ∈ Internal(T)} J(leaves(u), S)
```

Convention: the root is in `Internal(T)`, hence `J*(L_p, T) = 1` for
every `T` over `L_p`. Root is excluded from candidate sets at v-side
via `size(v) ≤ N_p − 1`; the root cannot be a "task-induced module"
because it has the same leafset as all other roots.

### 2.3 Match threshold

```
match(S, T; J_min) ⇔ J*(S, T) ≥ J_min
```

`J_min = 0.9` (locked). Interpretation: a match permits at most ~10%
leafset disagreement — for a subtree of size 20 that's at most 1
leaf swapped in and 1 leaf swapped out (`J = 18/22 ≈ 0.82` would
fail, `J = 19/21 ≈ 0.905` would pass). This is strict by design:
the user's directive is *"a group of leaves should really be the
same group of leaves, otherwise how can you justify it"*.

### 2.4 Pooled task-phase candidate set

```
V_task(p, b) := { v ∈ Internal(T^{p,b,task_learn})
                ∪ Internal(T^{p,b,task_test})
                  : k_min ≤ size(v) ≤ N_p − 1 }
```

Each `v ∈ V_task` carries its native `h_rel(v)` from its source tree
plus a `source(v) ∈ {TL, TT}` tag. The pooled set is a **multiset**:
if the same leafset `S` appears as an internal node in both `T_TL`
and `T_TT` (with possibly different `h_rel`), both occurrences are
kept. This is intentional — a leafset that emerges as an internal
node in BOTH task trees is more "task-anchored" than one appearing
in only one. A leafset-deduplicated variant is offered as a
sensitivity panel (§6.7).

### 2.5 retained — step indicator

For `v ∈ V_task` with `S = leaves(v)`,

```
task_induced(v; J_min) ⇔ J*(S, T^{p,b,rest_pre})  < J_min
persisting(v;   J_min) ⇔ J*(S, T^{p,b,rest_post}) ≥ J_min
retained(v;     J_min) ⇔ task_induced(v) ∧ persisting(v)
```

`retained ∈ {0, 1}`. Note: there is **no condition on the other task
phase**. A subtree from `T_TT` is allowed to be present in `T_TL`
(and vice versa); the trace claim is "any task → rest_post is OK"
per the user directive.

### 2.6 ΔJ — smooth (continuous) variant

```
J_pre(v)  := J*(leaves(v), T^{p,b,rest_pre})
J_post(v) := J*(leaves(v), T^{p,b,rest_post})
ΔJ(v)     := J_post(v) − J_pre(v)   ∈ [−1, +1]
```

Positive `ΔJ` = leafset is more closely matched in `rest_post` than
in `rest_pre`. Continuous in node positions, no `J_min`. Cohort
aggregation uses `max(ΔJ, 0)` so the cohort field stays in `[0, 1]`
and treats disruption (`ΔJ < 0`) as zero — disruption side is out of
scope per the user directive (Q4 lock).

---

## 3. Per-patient observables

For every `(p, b)`, enumerate `V_task(p, b)`. For each `v ∈ V_task`
record one row:

```
(p, b, source(v), h_rel(v), size(v),
 J_pre(v), J_post(v),
 retained(v; J_min), ΔJ(v))
```

This is the patient-level data product. **No grid, no tolerance** —
each subtree is evaluated at its own native height. Output:
long-format CSV
`data/reports/imcoh_mrl/mrl_per_node_tau{85,90,95}.csv`.

### 3.1 Quantity per (p, b)

`|V_task(p, b)|` is the count of pooled task internal nodes
satisfying `k_min ≤ size ≤ N_p − 1`. With `N_p ≈ 115` and `k_min = 3`:

```
|V_task| ≤ 2 · (N_p − 1 − N_kmin)   ≈ 2 · 110 = 220 candidates per (p, b)
```

Roughly half are size-2 leaf pairs that get filtered by `k_min ≥ 3`,
leaving `~150–180` candidates per (p, b). Cohort total:
`9 patients × 6 bands × ~170 candidates ≈ 9 × 10³` rows. Trivial.

---

## 4. Cohort aggregation (plot-time)

The cohort heatmap requires binning along `h_rel` for visualization.
The bins are a *display* choice; they do not enter the per-patient
data product.

Define the log-spaced bin axis

```
H = (h_0, h_1, …, h_K),    h_k = 0.05 · (1 / 0.05)^{k/K},   K = 12,
bin_k = [h_{k−1}, h_k),       k ∈ {1, …, K}.
```

Computed via `lrg_eegfc.utils.metrics.tree.h_log_grid([1.0], n=K+1)`.
The argument `[1.0]` is a placeholder; we want the grid in `h_rel`
units, not absolute distance. (See `h_log_grid` docstring for the
existing project convention.)

Per `(b, bin_k)`, define the eligible-row sets:

```
R(b, bin_k)   = { (p, v) : v ∈ V_task(p, b),  h_rel(v) ∈ bin_k }
R_p(b, bin_k) = { v      : v ∈ V_task(p, b),  h_rel(v) ∈ bin_k }      ∀ p
P_e(b, bin_k) = { p ∈ P  : |R_p(b, bin_k)| ≥ 1 }                       (eligible patients)
```

Cohort fields:

```
M̄_step(b, bin_k; J_min)  := (1 / |R|) · Σ_{(p, v) ∈ R} retained(v; J_min)
M̄_smooth(b, bin_k)       := (1 / |R|) · Σ_{(p, v) ∈ R} max(ΔJ(v), 0)
Π(b, bin_k; J_min)        := (1 / |P_e|) · |{ p ∈ P_e : Σ_{v ∈ R_p} retained(v; J_min) ≥ 1 }|
```

`M̄_step, M̄_smooth, Π ∈ [0, 1]`. A cell is flagged
**underpowered for the cohort claim** when `|P_e(b, bin_k)| < 7` —
the ≥ 7/9 cohort threshold cannot be evaluated.

### 4.1 Per-patient quantisation note

Because `M̄_step` is a count divided by `|R|` (pooled across patients),
it is not visibly quantised at typical cell-counts (`|R| ≈ 50–150`).
The patient-level analogue is more discrete:

```
m_p_step(b, bin_k) := (1 / |R_p|) · Σ_{v ∈ R_p} retained(v),
   range = {0, 1/|R_p|, …, 1}
```

If `|R_p| ≤ 5`, `m_p` is coarsely quantised — visualization should
annotate the per-patient cell-fill where `|R_p| < 5` so a reader can
discount low-resolution cells.

---

## 5. Properties

| # | Statement |
|---|---|
| 5.1 | `M̄_step, M̄_smooth, Π ∈ [0, 1]` |
| 5.2 | Native-height design: per-patient evaluation has no synthetic grid, no tolerance, no parameter `ε`. |
| 5.3 | Band-row-independence: bands are computed independently; no cross-band pooling. |
| 5.4 | Monotonicity in `J_min` is **opposing** between the two component predicates: `task_induced(v; J_min) ⇔ J*(S, T_pre) < J_min` is monotone **increasing** in `J_min` (stricter ⇒ more nodes flagged task_induced); `persisting(v; J_min) ⇔ J*(S, T_post) ≥ J_min` is monotone **decreasing** in `J_min` (stricter ⇒ fewer nodes flagged persisting). The conjunction `retained = task_induced ∧ persisting` is therefore **non-monotone** in `J_min` in general — a node can be retained at one threshold and not at a stricter or looser one. Implication: the {0.85, 0.90, 0.95} sweep is informative (not redundant); `J_min = 0.9` is the locked primary, the other two are robustness-only. Watch for cells where the band ranking flips across `J_min` — those are unreliable. |
| 5.5 | Asymmetry: `M̄_+` only; mirror `M̄_−` (rest_pre modules failing into rest_post) is dropped per Q4 directive. |
| 5.6 | NOT a metric. `retained` is a directed retention indicator; not a tree distance. |
| 5.7 | NOT statistical. No FDR, no q-values. Permutation null deferred to validation phase per Q5 directive. |
| 5.8 | Anchored in primitives — `tree_internal_nodes`, `jaccard_leafsets`, `dmax_from_Z`, `h_log_grid` (plot-time only). All in `lrg_eegfc.utils.metrics.tree`. |

### 5.9 What MRL does NOT measure

- **Global tree similarity** between phases (failed direction H2-RAW / H2-FROB / H2a-topo, see `2026-04-24_multiscale-task-trace.md` §7c).
- **Pair-level co-clustering** at integer cuts — that is H2d. MRL is the subtree-level lift; H2d is the `size = 2` ∧ integer-`k` restriction.
- **Direction of continuous distance shift** — that is H2c.
- **Disruption** of pre-existing modules — `M̄_−` mirror, dropped per Q4.
- **Statistical significance** at the (b, bin_k) level — null model deferred per Q5.
- **task_test-specific provenance** — by Q3 directive, presence in either task phase counts; we cannot distinguish learning carry-over from test-specific module formation at this n.

---

## 6. Caveats & failure modes

### 6.1 Same-probe bias at coarse `h_rel`

Same-probe MSC bias mostly mitigated under `imcoh_abs` but not eliminated
(`.agents/guides/02_methods/probe-bias-guide.md`). At `h_rel ≳ 0.7`, merges
follow probe geometry — anatomy dominates, task signal is small. Match
probabilities are then driven by anatomy in all three trees simultaneously,
which mostly *reduces* `task_induced` (rest_pre also matches), keeping
`M̄_step` low — but spurious `M̄_smooth` spikes are possible if a probe
configuration drifts artifactually.

**Mitigation.** Hatched overlay on bins with mid-edge `h_rel ≳ 0.7`.
Confirm load-bearing claims at fine/mid scales (`h_rel ∈ [0.05, 0.5]`).

### 6.2 Anatomy-dominated rest invariance

`T_rest_pre` and `T_rest_post` share anatomy — they are typically very
similar. `J*(S, T_rest_pre) ≈ J*(S, T_rest_post)` for most leafsets `S`.
The conjunction `task_induced ∧ persisting` is therefore necessarily a
*small-numerator* observable. MRL is conservative under anatomy-shared
baselines; high `M̄_step` is hard to fake.

### 6.3 Strict `J_min = 0.9` ⇒ sparse step landscape

At `J_min = 0.9` many `(b, bin_k)` cells will have `M̄_step = 0`. This
is by design — strict justifiability (Q1 directive). The smooth
co-primary `M̄_smooth` captures the gradient where `M̄_step` is zero.
Together they answer "where do task-induced modules persist exactly
(step)" and "where does post tend toward task-like membership at all
(smooth)".

### 6.4 Pooling task_learn and task_test

`V_task = V(T_TL) ∪ V(T_TT)` is a multiset. A leafset shared between
both task trees contributes twice (once per source tree). This is
intentional — repeated presence is more task-anchored. Sensitivity:
re-compute on leafset-deduplicated `V_task'` (per (p, b), unique
leafsets only; `h_rel` taken as the mean across sources where it
appears) and report side-by-side. If `M̄_step(dedup) ≈ M̄_step(pooled)`
within ± 0.05 cohort-wide, dedup is safe to ignore.

### 6.5 Pat_03 outlier

1024 Hz vs 2048 Hz cohort. `h_rel` normalisation protects scale
comparability. Report `M̄_with_Pat_03` and `M̄_without_Pat_03` per
(b, bin_k); flag any cell where `|M̄_with − M̄_without| > 0.1` —
Pat_03 dominates that cell, claim is fragile.

### 6.6 Root and singleton exclusion

`size(v) ≤ N_p − 1` excludes the root (would always match every
tree, `J*(L_p, ·) = 1`, hence `task_induced` always false at root).
`size(v) ≥ k_min = 3` excludes leaf pairs (high chance match by
coincidence at low `J_min`; with `J_min = 0.9`, a size-2 cluster
must share both leaves with some `T_rest_post` size-2 cluster to
count). Sensitivity sweep `k_min ∈ {2, 3, 5}` reports stability.

### 6.7 Ties in heights

Rare (`< 0.5%` empirically on `imcoh_abs` UPGMA). When two internal
nodes share an exact merge height, both contribute to the same bin
along `h_rel`. No bias — bin-counts simply increment by 2 instead
of 1.

### 6.8 What `M̄ > 0` cannot distinguish

- (i) Genuinely band-and-scale-specific modules forming in task and
  surviving — the desired interpretation.
- (ii) Consistent measurement drift / epoch-length effect that
  produces a "task-novel" subtree in rest_post by accident.

H2e (`2026-04-24_multiscale-task-trace.md` §7a) showed within-session
drift `≈ 0` for the H2c continuous direction; an MRL analogue (split
each phase, recompute `M̄_step`/`M̄_smooth` cohort-wide on the half-data
trees, compare to full-data MRL) is the right reviewer-defense
companion. Deferred to validation phase per Q5.

### 6.9 Identifiability under empty `R_p`

If `|R_p(b, bin_k)| = 0` for some patient (no eligible task subtree
falls in that bin), `p` does not contribute to the cell. The
`|P_e|` map is reported alongside `M̄`; cells with `|P_e| < 7` are
flagged underpowered for the ≥ 7/9 cohort claim.

---

## 7. Pseudocode

```
Input:
  P, B, Φ                           # cohort, bands, phases
  J_min ∈ {0.85, 0.9, 0.95}         # primary 0.9; sweep three values
  k_min = 3                         # subtree size floor
  H = (h_0, …, h_K), K = 12          # log-spaced viz bin axis
  fc_method = "imcoh_abs"

# === Per-patient pass (no grid, no tolerance) ===
For each (p, b) ∈ P × B:
  Z_TL   ← load_lrg_result(p, "task_learn",  b, fc_method).linkage
  Z_TT   ← load_lrg_result(p, "task_test",   b, fc_method).linkage
  Z_pre  ← load_lrg_result(p, "rest_pre",    b, fc_method).linkage
  Z_post ← load_lrg_result(p, "rest_post",   b, fc_method).linkage

  V_pre  ← [u ∈ tree_internal_nodes(Z_pre)  : size(u) ≥ k_min]
  V_post ← [u ∈ tree_internal_nodes(Z_post) : size(u) ≥ k_min]

  V_task ← []
  for (Z, src) ∈ {(Z_TL, "TL"), (Z_TT, "TT")}:
    for v ∈ tree_internal_nodes(Z):
      if k_min ≤ size(v) ≤ N_p − 1:
        v.source ← src
        V_task.append(v)

  for v ∈ V_task:
    S ← leaves(v)
    J_pre  ← max(jaccard_leafsets(S, leaves(u)) for u in V_pre)
    J_post ← max(jaccard_leafsets(S, leaves(u)) for u in V_post)
    for J_min ∈ {0.85, 0.9, 0.95}:
      retained[J_min] ← (J_pre < J_min) ∧ (J_post ≥ J_min)
    ΔJ ← J_post − J_pre
    emit row: (p, b, v.source, v.h_rel, v.size, J_pre, J_post,
               retained[0.85], retained[0.90], retained[0.95], ΔJ)

# === Cohort aggregation (plot-time only) ===
For each J_min ∈ {0.85, 0.9, 0.95}:
  For each (b, bin_k) ∈ B × {1, …, K}:
    R     ← rows where row.b = b, h_{k−1} ≤ row.h_rel < h_k
    P_e   ← unique row.p in R
    M̄_step[J_min, b, bin_k]   ← mean(row.retained[J_min] for row in R)
    M̄_smooth[b, bin_k]         ← mean(max(row.ΔJ, 0) for row in R)
    Π[J_min, b, bin_k]         ← (#patients in P_e with at least one
                                   retained[J_min] = 1 row in this cell)
                                   / |P_e|
    (Note: M̄_smooth does not depend on J_min.)

Output:
  data/reports/imcoh_mrl/mrl_per_node_tau{85,90,95}.csv  (long format)
  data/reports/imcoh_mrl/mrl_cohort_tau{85,90,95}.csv     (cohort fields)
  data/reports/imcoh_mrl/figures/mrl_landscape.{pdf,png}
```

### 7.1 Complexity

Per `(p, b)`: enumerating internal nodes is `O(N_p)`; per node `v`,
computing `J*` against `V_pre` and `V_post` is `O(N_p · N_p)` set ops
on hashed leafsets. Total per `(p, b)`: `O(N_p^3)` worst case.
With `N_p ≈ 115`: `~ 1.5 × 10^6` ops. Cohort: 54 (p, b) cells →
`~ 8 × 10^7` ops, well under a second in vectorised numpy.

### 7.2 Optimisation notes

- Cache leafsets as `frozenset[int]` or `uint8[N_p]` bitmasks for
  fast intersection.
- Early-exit `J*` when any candidate exceeds `J_min` (turns
  `O(N_p^2)` into expected `O(N_p)` per node).
- Pre-prune candidate `u` by size constraint:
  `J(leaves(v), leaves(u)) ≥ J_min` requires
  `size(u) / size(v) ∈ [J_min, 1/J_min]` — for `J_min = 0.9` this
  is `[0.9, 1.111]`, eliminating most candidates.

---

## 8. Visualization

The visualization has **two figures**, in service of two questions:

- **"Where (band × scale) does the trace live cohort-wide?"** — Figure 1.
- **"What does the trace *look like* on the actual data?"** — Figure 2.

The heatmap (Figure 1) is the cohort observable; the dendrogram
reprojection (Figure 2) is the visual proof on the network's
hierarchical organization. The MRL self-selects which subtrees to
draw on Figure 2 — the "retained" leafsets at the cell of interest.

### 8.1 Figure 1 — H3-style cohort consensus heatmap (single panel)

Layout deliberately mirrors `H3_imcoh_abs.pdf`:

- One panel, 6 rows (`δ, θ, α, β, γ_l, γ_h` from top), x-axis
  log-spaced `h_rel` bins (default K=12, range `[0.05, 1)`).
- Cell value: `Π(b, bin_k; J_min = 0.9)` — fraction of patients with
  at least one retained subtree in the cell.
- Sequential white→orange→red colormap matching H3 exactly:
  - `Π < 0.5` → white-to-pale
  - `Π ∈ [0.5, 7/9)` → orange ("majority")
  - `Π ∈ [7/9, 1)` → red ("≥ 7/9, cohort-wide")
  - `Π = 1` → dark red ("unanimous")
- Colorbar with the same threshold ticks (`majority`, `5/9`, `7/9`,
  `unanimous`).
- Hatched overlay on bins with mid-edge `> 0.7` ("probe-bias suspect").
- Cell border (red, 1.6 pt) wherever `Π ≥ 7/9` to make the cohort-wide
  cells unmistakable.
- Cell border (grey dotted, 0.8 pt) wherever `n_pat < 7` (underpowered).
- No `fig.suptitle`; companion `.md` carries the caption.

A single panel forces a single read: "where in (band, h_rel) does the
cohort agree that retention occurs". A horizontal red ridge in row α
between `h_rel ∈ [0.1, 0.2]` would say *"fine-scale α modules trace
cohort-wide"*; a blank β row would say *"no β trace"*.

### 8.2 Figure 2 — Dendrogram reprojection (the visual proof)

Driven by Figure 1: pick the band × `h_rel` bin with the strongest
cohort signal (highest `Π` among the cohort-wide cells, at non-probe-bias
`h_rel`). For 3 representative patients (priority: highest
patient-level `M_step` in that cell):

- One row of 4 dendrograms per patient: `rest_pre | task_learn |
  task_test | rest_post`, sharing the leaf set `L_p` but with
  phase-specific orderings.
- For each patient, identify the *strongest retained module* at that
  scale: among `v ∈ V_task(p, b)` with `h_rel(v)` in the target bin
  AND `retained(v; J_min = 0.9) = 1`, take the one with the highest
  `J*(leaves(v), T_post)` (most strongly persistent). Call it `v*_p`,
  with leafset `S*_p ⊆ L_p`.
- **Colour the leaves in `S*_p`** with a distinct, per-patient colour
  (consistent across the 4 phase panels of that patient). All other
  leaves rendered in neutral grey.
- **Mark the matched subtree height** in each phase: find the internal
  node `u_φ` of `T^{p,b,φ}` with the highest `J(leaves(u_φ), S*_p)`.
  Draw a thin horizontal bar in the patient's colour at `h_rel(u_φ)`,
  spanning the leaves of `u_φ`. Annotate the bar's `J` to the right
  ("J = 0.94"). The reader sees the matched subtree at-a-glance.
- Y-axis log-scale, dendrogram limits per the `tmin/tmax` rule
  (`merge_heights[0]*0.8`, `merge_heights[-1]*1.05`).

The reader narrative:
- `rest_pre` panel: the patient's coloured leaves are scattered across
  the dendrogram, the matching internal node has a low `J` (e.g. 0.4)
  and is small / shallow.
- `task_learn` panel: the same leaves coalesce into a single coherent
  branch at intermediate height (`J` jumps to ~0.85+).
- `task_test` panel: same. (Confirms task-anchoring.)
- `rest_post` panel: the leaves are AGAIN coherent, with `J ≥ 0.9` to
  the same module — the trace.

Three patient rows × 4 phase columns = 12 dendrogram panels per
figure. Per-band figures live at
`data/reports/imcoh_mrl/figures/mrl_reproject_<band>_h<bin>.pdf`.

### 8.3 Reading rules

- **Figure 1 ridge in band-row** → multiscale cohort agreement on
  retention in that band.
- **Figure 1 horizontal contrast across rows** → band-specificity.
- **Figure 2 — leaf coherence across the 4 phases** = the trace, on the
  actual hierarchy. If the coloured leaves are coherent in pre too,
  the patient is NOT a positive contributor — and `retained(v) = 0`
  would have prevented us from selecting them in the first place; if
  coherent in pre is what you see, debug.
- **Discrepancy between Figure 1 and Figure 2** → diagnostic; expect
  that high-`Π` cells *must* support 3+ visually-clean reprojections
  in Figure 2. If they don't, the heatmap is misleading us.

### 8.4 Supplementary panels

- `J_min ∈ {0.85, 0.95}` cohort heatmaps (single-row strip of three
  H3-style panels) — robustness.
- Per-patient strip: 9-row × K-bin `m_p_step(b, bin)` per band,
  6 panels (one per band) — cohort-wideness disaggregated.
- `Pat_03 in/out` sensitivity: two H3-style heatmaps side-by-side at
  `J_min = 0.9`.
- Smooth co-primary as a separate H3-style panel using
  "fraction of patients with mean `ΔJ > 0` in the bin" — captures the
  gradient when the strict step is sparse.

---

## 9. Connection to prior tools

| Prior measure | Relation to MRL |
|---|---|
| H2d (`h2d_coactivation_persistence.py`) | MRL is the subtree-level lift from pairs (`size = 2`) and integer cuts (`k`) to subtrees of arbitrary size at native fractional heights. H2d ≅ MRL with `size(v) = 2` and an integer-`k` partition cut. |
| H1-topo / H2a-topo (bipartition overlap, `h2_topology_directed.py`) | MRL adds the conditional "absent in rest_pre" predicate. H1-topo is unconditional bipartition overlap; MRL is *residual* / *conditional* on baseline absence. |
| H2c (`h2c_ultrametric_drift.py`) | H2c operates on continuous distance shifts (Spearman ρ). MRL operates on discrete subtree identity. Complementary observables on the same residual structure. |
| KC / MC / wRF (`utils.metrics.tree_distance`) | These are scalar tree distances integrating over scales and assuming a fixed leaf labeling. MRL keeps the scale axis explicit (vector-valued in `h_rel`) and uses leafset matching directly — does not aggregate to a scalar. The integration-to-scalar direction is the failed scalar-gate effort (post-mortem, archived). |
| partition-multiscale cluster-perm (`h2_partition_multiscale.py`) | Cluster-perm on partition divergence, k-resolved. MRL is also h-resolved but on subtree identity. Complementary. The δ k=23–31 (`Δ_VI p=0.014`) and α k=2–4 (`Δ_H p=0.050`) cluster-perm hits should appear in `M̄` at compatible `h_rel` ranges — δ at moderate `h_rel ≈ 0.2–0.4` (k=23–31 of an N≈115-leaf tree), α at near-root `h_rel ≈ 0.5–0.8` (k=2–4). Sanity-check pass after first compute. |
| CBR umbrella (`task-persistence-investigation/2026-04-25_cbr-investigation.md`, README) | MRL is the *task-anchored Jaccard* variant of the broader CBR (cluster-birth-retention) family the README umbrella points to. Sibling variants (containment, consensus-subtree, rest_post-anchored legacy) live in separate scope reports under that umbrella; primitives are shared. |

---

## 10. Implementation plan

### 10.1 Script

`scripts/01_compute/hypothesis_tests/mrl_landscape.py` — new, ~200 lines.

CLI:
```
python scripts/01_compute/hypothesis_tests/mrl_landscape.py \
    [--patients Pat_02,Pat_03,...] \
    [--bands delta,theta,alpha,beta,low_gamma,high_gamma] \
    [--j-min-sweep 0.85,0.9,0.95] \
    [--k-min 3] \
    [--n-bins 12] \
    [--h-min 0.05] \
    [--fc-method imcoh_abs] \
    [-v]
```

Imports: `load_lrg_result`, `tree_internal_nodes`, `jaccard_leafsets`,
`dmax_from_Z`, `h_log_grid` from existing library. No new helpers in
v1 (per `coding-rules.md`, promotion to library waits for second consumer).

Outputs:
- `data/reports/imcoh_mrl/mrl_per_node.csv` — long format, all rows
  (p, b, source, h_rel, size, J_pre, J_post, retained_85, retained_90,
  retained_95, ΔJ).
- `data/reports/imcoh_mrl/mrl_cohort_tau{85,90,95}.csv` — cohort fields
  per (b, bin_k).
- `data/reports/imcoh_mrl/mrl_cohort_smooth.csv` — cohort `M̄_smooth`
  per (b, bin_k); separate file because no `J_min`.
- Sidecar `data/reports/imcoh_mrl/mrl_landscape.md` — short
  renormalised summary + reading rules + cell counts where `Π ≥ 0.78`.

### 10.2 Figure generator

`scripts/01_compute/figures_embedded/fig_mrl_landscape.py` — new, ~150 lines.

Reads cohort CSVs, produces:
- `data/reports/imcoh_mrl/figures/mrl_landscape.{pdf,png}` — headline
  3-panel at `J_min = 0.9`.
- `mrl_landscape_supplementary.{pdf,png}` — `J_min ∈ {0.85, 0.95}`
  variants.
- `mrl_per_patient_strips.{pdf,png}` — 6-panel per-patient strips.
- `mrl_pat03_sensitivity.{pdf,png}` — Pat_03 in/out paired heatmaps.
- `mrl_dedup.{pdf,png}` — leafset-deduplicated sensitivity.

### 10.3 Library promotion (deferred)

Per `coding-rules.md`: a helper used by ≥ 2 callers gets promoted.
If a sibling CBR variant (under `task-persistence-investigation/`)
needs `J*`, `retained`, or `ΔJ`, promote
`(jaccard_star, retained_predicate, delta_jaccard)` into a new
`lrg_eegfc.utils.metrics.tree_persistence` module in the same commit
as the second consumer.

### 10.4 Tests

`tests/test_mrl.py` covering:

| Test | Expectation |
|---|---|
| Identity: `T_pre = T_TL = T_TT = T_post`, toy tree | `M̄_step = 0` for all (b, bin); every node matches both pre and post; `task_induced` is false everywhere |
| Saturation: `T_pre = star`, `T_TT = T_post = balanced` | Every internal node `v ∈ T_TT` of `size > 1` has `J*(·, T_pre) = 0` (only root matches the leafset of `v`, and root J ≠ 1 unless `v` is root); `J*(·, T_post) = 1`; `M̄_step = 1` for `J_min ≤ 1`. Verify on toy. |
| ΔJ symmetry under (pre↔post) swap | ΔJ negates; `max(ΔJ, 0)` swaps with `max(−ΔJ, 0)`; `M̄_smooth(swap) = M̄_smooth(disrupt)`. |
| Root exclusion | `M̄` at the bin containing `h_rel = 1` is 0 (no eligible node). |
| Bin-edge inclusion | `[h_{k−1}, h_k)` half-open; node at `h_rel = h_k` falls in `bin_{k+1}`. |

---

## 11. Locked decisions (closed open questions)

| Q | Decision | Date |
|---|---|---|
| Q1 — Default `J_min` | `J_min = 0.9` primary; sweep `{0.85, 0.9, 0.95}` for sensitivity | 2026-04-25 |
| Q2 — Height-window scheme | None at per-patient level (native heights only); log-spaced bins for cohort viz only | 2026-04-25 |
| Q3 — Task-induced predicate | Pooled task search across `V(T_TL) ∪ V(T_TT)`; condition is `¬match(rest_pre)` only (no condition on the other task phase) | 2026-04-25 |
| Q4 — Mirror field `M̄_−` | Drop. Only `M̄_+` (positive trace into rest_post) reported. | 2026-04-25 |
| Q5 — Permutation null | Defer to validation phase. Descriptive-only in v1. | 2026-04-25 |
| Q6 — Smooth alternative `ΔJ` | Co-primary with step `M̄_step`. Three-panel headline (step + smooth + Π). | 2026-04-25 |

## 12. Pass criterion (descriptive)

The MRL passes the surface-evidence test if the headline three-panel
figure shows:

1. A **non-uniform landscape** in at least one of (`M̄_step`, `M̄_smooth`)
   — some `(b, bin_k)` cells materially above zero, others near zero.
   Uniformity in either direction is a null result.
2. **Band-row contrast** consistent with H2d band typology
   (α strongest, θ weakest at fine-mid `h_rel`).
3. **Vertical persistence** — multiple bins active in the same band-row
   (the nested-modules claim).
4. **Cohort-wideness** — at least one cell per "active" band has
   `Π ≥ 7/9 = 0.78`.
5. **Robustness across `J_min ∈ {0.85, 0.9, 0.95}`** — qualitative
   band ranking preserved.
6. **Pat_03 sensitivity** — removing Pat_03 changes `M̄_step` by `< 0.1`
   in every cell with `M̄_step ≥ 0.1`.

If 1–4 hold (5–6 are sensitivity controls), the MRL is a usable
descriptive observable for the synthesis figure of the paper. The
load-bearing interpretation reads:

> *"Task-induced modules at fractional dendrogram heights `h_rel ∈ [a, b]`
> in band `X` persist into `rest_post` in ≥ 7/9 patients, where
> 'task-induced' = absent from rest_pre and 'persists' = present in
> rest_post at Jaccard ≥ 0.9."*

Band-specific, multiscale, cohort-wide, no scalar gate.
