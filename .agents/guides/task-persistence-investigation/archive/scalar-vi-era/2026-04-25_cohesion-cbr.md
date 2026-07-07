---
name: cohesion-cbr
type: scope
era: COHORT_N9
status: current
created: 2026-04-25
updated: 2026-04-25
implementation_note: |
  Sanity check on the diagnostic cells (§5 of this scope) revealed that
  tightness `T(S, T_φ) = |S| / |LCA(S, T_φ)|` is too strict a similarity
  score on its own — a small number of stragglers can pull the LCA up to
  the root and crash `T_tt` to ~0.4 even when the bulk of the leafset is
  cohesive. The primary classifier in audit_12 therefore uses **Jaccard**
  scores `J_rpre, J_tl, J_tt` with `J_task := mean(J_tl, J_tt)` and
  corner-distance soft affinities on the `(J_rpre, J_task)` unit square.
  Tightness is retained as a diagnostic sidecar in
  `anchor_classification.csv` for follow-up. The 4-corner geometry,
  per-leaf affinity averaging, and confidence-graded alpha rendering are
  unchanged from §2.5–2.6. See §10 for the implementation notes.
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-classifier-failure-analysis.md
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Cohesion-CBR (variant E)

**Replaces both Jaccard and containment with *tightness* `T(S, T_φ)
:= |S| / |leaves(LCA(S, T_φ))|` — the directional cohesion score that
asks the visual question literally ("do these leaves bunch into one
subtree of T_φ?"). Drops the four hard-threshold predicates in favour
of four soft affinities so no anchor falls into a limbo zone and no
leaf defaults to rearrange because its ancestors went unclassified.
Anchored on rest_post for direct comparability with the legacy CBR
figures; per-leaf colour is the dominant affinity and saturation is the
margin (top1 − top2).**

## 1. Notation

Inherits the cohort/band/phase notation of the parent scope reports.
Recap:

- `P, B, Φ` — `COHORT_N9` patients, 6 bands, 4 phases.
- `T^{p,b,φ}` — binary tree induced by the LRG linkage on
  `imcoh_abs`, with `N_p` leaves and `N_p − 1` internal nodes.
- For each internal node `v`: `leaves(v) ⊆ L_p`, `size(v)`, `h(v)`,
  `h_rel(v) = h(v) / dmax(T^{p,b,φ})`.

Fixed parameters (kept in lockstep with all other CBR variants for
direct comparability):

- `s_min := 5`,  `s_max := 60`.
- No `τ_match`, `τ_dissim` thresholds — affinities are continuous.

## 2. Definitions

### 2.1 Least-common-ancestor leafset and tightness

For a leafset `S ⊆ L_p` and a tree `T_φ`, define the **least common
ancestor** of `S` in `T_φ` as the smallest internal node containing
`S`:

    LCA(S, T_φ) := argmin_{v ∈ Internal(T_φ), S ⊆ leaves(v)} size(v).

If `S` is a singleton, `LCA(S, T_φ)` is the leaf itself's parent
(size ≥ 2). If `S` spans the entire tree, `LCA(S, T_φ)` is the root.
The LCA exists for every non-empty `S ⊆ L_p`.

The **tightness** of `S` in `T_φ` is

    T(S, T_φ) := |S| / size(LCA(S, T_φ)).

`T : 2^{L_p} \\ ∅ × Trees(L_p) → (0, 1]`. Properties:

- `T = 1` ⇔ `S` is the leafset of an internal node of `T_φ` (S is *exactly*
  a subtree).
- `T → 0` as `|LCA(S, T_φ)| → N_p` (S spans the whole tree).
- *Asymmetric in the (S, T_φ) sense* — `T(S, T_φ)` is not equal to
  `T(leaves(v), T_S)` in general.
- *Monotonically non-increasing in `S`*: enlarging `S` can only enlarge
  or preserve the LCA, never shrink it.
- Cheap to compute: `O(|S| · |Internal(T_φ)|)` worst case using the
  set-membership scan, or `O(N · |S|)` using ancestor walking.
  Implementation reuses `tree_internal_nodes` from
  `lrg_eegfc.utils.metrics.tree`.

### 2.2 Anchor universe

This variant anchors on `rest_post` for direct comparability with the
legacy CBR figures (`audit_08`). The anchor universe is

    A(p, b) := {η ∈ Internal(T^{p,b,rest_post}) : s_min ≤ size(η) ≤ s_max}.

A task-test-anchored sister variant is straightforward by symmetry but
is *not* in this scope.

### 2.3 Per-anchor tightness vector

For each anchor `η ∈ A(p, b)` with leafset `S = leaves(η)`:

    T_rpre(η)  := T(S, T^{p,b,rest_pre})
    T_tl(η)    := T(S, T^{p,b,task_learn})
    T_tt(η)    := T(S, T^{p,b,task_test})              (≡ NaN if task_test missing)
    T_rpost(η) := T(S, T^{p,b,rest_post})              (= 1 by construction; reported for symmetry)

We collapse the two task scores into a single "task tightness" via
the symmetric-task minimum,

    T_task(η) := min(T_tl(η), T_tt(η))                 (or T_tl alone if T_tt is NaN).

Using min over the two task phases is the conservative choice: an
anchor only counts as "cohesive in task" if it is cohesive in both
task phases. This is the cohesion analog of the legacy CBR's
`min(J_tl, J_tt) ≥ τ_match` requirement.

### 2.4 Soft pattern affinities

For each anchor `η ∈ A(p, b)`, define four affinities in [0, 1]:

    a_trace(η)     := max(0, T_task(η)  − T_rpre(η))
    a_persist(η)   := min(T_rpre(η), T_task(η), T_rpost(η))
    a_reset(η)     := max(0, T_rpre(η)  − T_task(η))
    a_rearrange(η) := max(0, 1 − max(T_rpre(η), T_task(η)))

Reading:

- `a_trace ↑` when the anchor's leaves bunch in task and not in rsPre
  ("task-induced module that persists into rsPost" — `T_rpost = 1` by
  the rsPost anchor).
- `a_persist ↑` when the anchor's leaves bunch in *every* phase
  (anatomy / always-on connectivity).
- `a_reset ↑` when the leaves bunch in rsPre but not in task (the task
  disrupted a pre-existing module that returned in rsPost — hence the
  rsPost anchor).
- `a_rearrange ↑` when the leaves don't bunch in any non-rsPost phase
  (rsPost-unique structure).

The four affinities are *not* mutually exclusive and do *not* sum to 1.
That is intentional — they are independent axes; an anchor with
`(T_rpre = 0.3, T_task = 0.9, T_rpost = 1)` has high `a_trace = 0.6`,
moderate `a_persist = 0.3`, zero `a_reset`, and zero `a_rearrange`,
which correctly captures "predominantly trace-leaning with some weak
persist component".

### 2.5 Per-leaf affinity vector

For each leaf `ℓ ∈ L_p` and each pattern `p ∈ {trace, persist, reset,
rearrange}`,

    A_p(ℓ; p, b) := mean_{η ∈ A(p, b), ℓ ∈ leaves(η)} a_p(η).

Mean (uniform weighting over containing anchors) is the default;
size-weighting and h_rel-weighting are open-question variants (§9).

Edge case: a leaf may belong to no anchor of size ∈ [s_min, s_max]
(typical for terminal leaves whose smallest containing internal node
has size 2-4). For such leaves all four affinities are 0 and we
explicitly mark them `unclassified`.

### 2.6 Per-leaf colour and saturation

For each leaf `ℓ` with at least one defined affinity:

    dominant(ℓ) := argmax_p A_p(ℓ).
    confidence(ℓ) := A_dominant(ℓ) − A_secondbest(ℓ)        ∈ [0, 1].

Colour rule for the figure:

- Hue = `PATTERN_COLOR[dominant(ℓ)]` (same palette as audit_08:
  trace red, persist black, reset gray, rearrange blue).
- Alpha = `α_min + (α_max − α_min) · confidence(ℓ)` with
  `α_min = 0.30`, `α_max = 1.00`. Faint = ambiguous; saturated =
  clear.
- For unclassified leaves (no eligible ancestor): alpha = 0.10,
  hue = neutral light gray.

This is the only figure-level change relative to audit_08; the
dendrogram layout, codebar geometry, and 4-phase row remain identical.

## 3. Properties

| property | value |
|:---|:---|
| Range of `T` | `(0, 1]` |
| Symmetric in `(S, T_φ)` | **No** (directional). |
| Captures cohesion vs. mere containment | **Yes** — fixes failure mode #2 of legacy. |
| Threshold-free | **Yes** — fixes failure mode #1 (limbo zone). |
| Every leaf gets a colour decision | **Yes** — fixes failure mode #3 (default-to-blue) modulo the small-tree edge case. |
| Detects fragmentation | **No** by construction — same limit as legacy / TA-CBR / Containment. A future *multi-LCA* variant could; out of scope here. |
| Identifiable failure mode | If two patterns score equally high (large `a_trace` and large `a_persist`), the leaf is faint and reads as "ambiguous" rather than incorrectly forced into one bucket. |
| Complexity | `O(|A| · N²)` per (patient, band) for the LCA scans. For `N ≈ 120` and `|A| ≈ 80`: ~1.2 M operations per cell, ~70 M cohort-wide; trivial. |

**Negative property — what Cohesion-CBR cannot detect:**

- Anchor leafsets whose pieces are split across *multiple* sibling
  subtrees of T_φ. Their LCA is forced upward; tightness drops; the
  anchor reads as "scattered". This is the same blind spot as every
  single-match CBR variant — the natural fix is a *multi-LCA*
  generalisation that reports the K best disjoint sub-LCAs and their
  combined tightness; deferred.
- The procedure is anchored on rest_post; subtrees unique to task that
  fragment in rsPost will be missed (same as legacy). A
  task-anchored Cohesion-CBR is straightforward and a follow-up.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| LCA gives no information about *where* the missing leaves are. Two anchors with the same `T_rpre = 0.5` could be 50% of a same-LCA subtree (single contiguous block plus extras) or two halves of a split LCA. Both read identically. | Acknowledged. The figure's value is the *aggregate* per-leaf affinity; ambiguous individual anchors get washed out by the per-leaf average. |
| Pat_14 has no task_test → `T_tt` undefined. Use `T_task = T_tl` alone for that patient. | Reported in figure title and tally CSV. |
| Pat_03 1024 Hz outlier → no special handling required at the cohesion level (the anchor is in rsPost regardless of fs); flagged for the eventual cohort census. | Inherit existing flag. |
| Probe-bias coarse-scale modules can dominate `a_persist` regardless of task. | Inherit `s_max = 60` window. Diagnose by comparing `a_persist` distribution on the probe-debiased FC pipeline (when available). |
| Affinity definitions are convex-but-not-canonical: alternatives exist (`a_trace = T_rpost · (1 − T_rpre)`, sigmoidal contrasts, etc.). | This scope locks in the simplest defensible forms. Sensitivity to the choice is part of §9. |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)
       s_min=5, s_max=60

I_rpre  ← internal_nodes(Z_rpre)
I_tl    ← internal_nodes(Z_tl)
I_tt    ← internal_nodes(Z_tt)        (None if missing)
I_rpost ← internal_nodes(Z_rpost)

FUNCTION lca_size(S, I):
    best ← N
    FOR each v in I:
        IF S ⊆ leaves(v) AND size(v) < best:
            best ← size(v)
    RETURN best

FUNCTION tightness(S, I):
    RETURN |S| / lca_size(S, I)

per_anchor ← []
FOR each η in I_rpost:
    IF size(η) < s_min OR size(η) > s_max: CONTINUE
    S ← leaves(η)
    T_rpre  ← tightness(S, I_rpre)
    T_tl    ← tightness(S, I_tl)
    T_tt    ← tightness(S, I_tt) IF I_tt ELSE NaN
    T_rpost ← 1.0
    T_task  ← min(T_tl, T_tt) IF NOT NaN ELSE T_tl
    a_trace     ← max(0, T_task  − T_rpre)
    a_persist   ← min(T_rpre, T_task, T_rpost)
    a_reset     ← max(0, T_rpre  − T_task)
    a_rearrange ← max(0, 1 − max(T_rpre, T_task))
    APPEND (η, T_rpre, T_tl, T_tt, T_task, a_*) to per_anchor

leaf_aff[1..N][p ∈ patterns] ← 0
leaf_count[1..N] ← 0
FOR (η, …, a_*) in per_anchor:
    FOR ℓ in leaves(η):
        FOR p in patterns:
            leaf_aff[ℓ][p] += a_p
        leaf_count[ℓ] += 1

FOR ℓ ∈ {1..N}:
    IF leaf_count[ℓ] == 0:
        leaf_label[ℓ] ← "unclassified"
        leaf_alpha[ℓ] ← 0.10
        CONTINUE
    FOR p in patterns:
        leaf_aff[ℓ][p] /= leaf_count[ℓ]
    sorted ← sorted(leaf_aff[ℓ], descending)
    leaf_label[ℓ] ← sorted[0].pattern
    confidence ← sorted[0].value − sorted[1].value
    leaf_alpha[ℓ] ← 0.30 + 0.70 · confidence

OUTPUT per_anchor, leaf_aff, leaf_label, leaf_alpha
```

## 6. Visualization spec

One PDF per (patient, band) at
`data/audit/per_patient_hierarchy_cohesion/{Pat_XX}/{band}_hierarchy-crossphase.pdf`.

Layout — *identical* to audit_08 / variant A / variant B / variant D
so all five measures are directly comparable side-by-side:

- 4 phase dendrograms in a row (rest_pre, task_learn, task_test,
  rest_post).
- Each panel: dendrogram + codebar strip beneath.
- Codebar: one rectangle per leaf, hue = `PATTERN_COLOR[dominant(ℓ)]`,
  alpha = `leaf_alpha[ℓ]`. **Faint colour ⇔ ambiguous classification;
  saturated ⇔ clear classification.** This is the only meaningful
  layout difference vs audit_08.
- Branches: an internal-node link gets the colour of its dominant
  category iff every leaf below it has the same dominant category
  *and* mean alpha ≥ 0.5; else background gray.
- Leaf stems: overdrawn in `(hue, alpha)` of the leaf itself.
- Legend title: `"{Pat_XX} — {band} — Cohesion-CBR (rsPost-anchored,
  size ∈ [5, 60]; α = top1 − top2)"`.

Auxiliary CSVs:
- `leaf_assignment.csv` schema:
  `patient, band, leaf_id, dominant, confidence,
   A_trace, A_persist, A_reset, A_rearrange`.
- `anchor_classification.csv` schema:
  `patient, band, anchor_row, size, h_rel, T_rpre, T_tl, T_tt, T_task,
   a_trace, a_persist, a_reset, a_rearrange`.

## 7. Connection to prior tools

| existing tool | how Cohesion-CBR relates |
|:---|:---|
| Legacy CBR (rsPost-anchored Jaccard) | Same anchor; replaces J with T; replaces strict-classify with soft affinities. Same per-leaf priority spirit, but every leaf is now classified. |
| TA-CBR (variant A) | Different anchor (task vs. rsPost). The two are independent axes; could be combined in a future "dual-anchor cohesion" measure. |
| Containment-CBR (variant B) | Cohesion is the score function containment *should* have been. Containment over-permissivity is exactly the failure tightness fixes. |
| CON (variant D) | CON is the limit case `T = 1` (exact equality). Cohesion is its continuous relaxation. |
| MRL | MRL uses Jaccard via `match_τ`. A cohesion-MRL is conceivable but a separate scope. |
| H2c / H2d | Untouched. Per-pair / per-scale statistics; complementary. |

## 8. Implementation plan

- **Library reuse**: `tree_internal_nodes` from
  `src/lrg_eegfc/utils/metrics/tree.py`. The `lca_size` /
  `tightness` helpers are inlined for the first run; second caller
  triggers promotion to that module per `coding-rules.md`.
- **Script**: `scripts/01_compute/audit/audit_12_cohesion_cbr.py`.
  Reuse `plot_dendrogram_with_strip` from `audit_08` with a small
  patch to accept per-leaf alpha.
- **Outputs**: as listed in §6.
- **Pat_14**: `T_tt` undefined → use `T_task = T_tl`. Already handled
  by the conditional in §2.3.

## 9. Open questions

- **Anchor weighting** in §2.5. Current default: uniform mean over
  containing anchors. Size-weighted mean (`weight = 1/size(η)`) would
  pull harder on small specific anchors; h_rel-weighted (low-h_rel
  more) would push toward fine-scale modules. Defer choice until first
  visual review.
- **`J_task` collapse rule** (current implementation). Default: `mean(J_tl,
  J_tt)`. `min` is conservative but penalises asymmetric task signals
  (e.g. Pat_05 γ_l where `J_tl=0.85` but `J_tt=0.80` reads slightly
  asymmetric). Decision after first cohort census.
- **Multi-LCA generalisation**. To detect fragmented anchors, replace
  `LCA` with the *K* best disjoint sub-LCAs and report total cohesion
  `Σ |sub_S_k| / Σ size(LCA_k)`. Separate scope.
- **Task-anchored Cohesion-CBR**. Trivial mirror; if Pat_14-excluding
  is acceptable, run alongside the rsPost-anchored variant for
  sensitivity. Currently deferred to keep the diff vs. legacy minimal.
- **Confidence threshold on alpha**. Should leaves with `confidence <
  0.05` (essentially undecidable) be marked separately? Current rule
  fades them via alpha; a categorical "ambiguous" colour is an
  alternative.

## 10. Implementation epilogue

The scope's original primary score was *tightness*. First-pass
implementation (commit ad65; see git history) used `T(S, T_φ)` for the
classifier with `T_task := min(T_tl, T_tt)` and contrast-based affinities
`a_trace := max(0, T_task − T_rpre)` etc. The result was bad: only 9
total trace leaves cohort-wide, dominated by spurious persist /
unclassified. Two failure modes:

1. **`min` over task phases** is too conservative. Pat_05 γ_l size-42
   has `T_tl = 0.69`, `T_tt = 0.37` (the 8 outlier leaves blow up the
   task_test LCA), so `min = 0.37` reads "weak task" though Jaccard
   says clean trace. `mean(T_tl, T_tt) = 0.53` is more honest.
2. **Tightness is too strict at the LCA**. A handful of straggler
   leaves can drag the LCA all the way to the root and crash `T` to
   `|S|/N`. Real-world FC modules have 1–3 noisy boundary leaves.
   Jaccard tolerates this via the union denominator; tightness does
   not.

The corrected implementation in `audit_12_cohesion_cbr.py` therefore
uses **Jaccard** for the classifier (familiar, validated, tolerant of
boundary noise) and **corner-distance soft affinities** for the
classification step. Tightness is retained in the per-anchor CSV as a
diagnostic sidecar. The figure layout, palette, codebar geometry,
confidence-graded alpha, and per-leaf projection are unchanged from
§2.5–2.6.

A future scope (`2026-XX-XX_multi-lca-cohesion.md`) can address the
fragmentation blind spot of single-match similarity scores via a
multi-LCA tightness generalisation; that's the principled fix for
"one anchor split across multiple coherent task subtrees".
