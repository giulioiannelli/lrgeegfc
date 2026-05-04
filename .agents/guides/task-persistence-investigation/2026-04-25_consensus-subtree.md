---
name: consensus-subtree
type: scope
era: COHORT_N9
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_containment-cbr.md
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Consensus-subtree (CON)

**A leafset `L ⊆ L_p` is *consensus task-trace* iff `L` is the exact
leafset of some internal node in BOTH `T^{p,b,task_test}` AND
`T^{p,b,rest_post}` AND NOT the leafset of any internal node in
`T^{p,b,rest_pre}`. The strictest CBR variant — no thresholds, no
similarity scores, exact subtree identity — designed to surface only
the most undeniable task-induced persistent modules.**

## 1. Notation

Inherits §1 of `2026-04-25_task-anchored-cbr.md`. Recap of the
load-bearing objects:

- For each `(p, b, φ)`: a binary tree `T^{p,b,φ}` with `N_p` leaves
  and `N_p − 1` internal nodes.
- For each internal node `v ∈ T`: the leafset `leaves(v) ⊆ L_p`.
- The *leafset family* of a tree `T` is the multiset of all `leaves(v)`
  for `v ∈ Internal(T)`. By construction it has `N_p − 1` elements,
  one per internal node, all distinct (binary trees have unique
  internal-node leafsets).

Define for each `(p, b, φ)`:

    F^{p,b,φ} := { leaves(v) : v ∈ Internal(T^{p,b,φ}) } ⊆ 2^{L_p}.

`|F^{p,b,φ}| = N_p − 1`.

Size constraint (kept identical to legacy / TA-CBR for comparability):
`s_min := 5`, `s_max := 60`. Restrict each leafset family to the size
window:

    F^{p,b,φ}_w := { L ∈ F^{p,b,φ} : s_min ≤ |L| ≤ s_max }.

## 2. Definitions

### 2.1 Consensus predicate

A leafset `L ⊆ L_p` is `CON-trace` for `(p, b)` iff

    CON-trace(L; p, b) ⇔  L ∈ F^{p,b,task_test}_w
                       ∧  L ∈ F^{p,b,rest_post}_w
                       ∧  L ∉ F^{p,b,rest_pre}.

(The `∉ F^{p,b,rest_pre}` condition is *not* size-windowed: even if
a near-identical rsPre subtree exists at any size, it disqualifies.)

This is exact set-equality at the leafset level. There is no
threshold, no max, no similarity score.

### 2.2 Sibling consensus predicates

The other three CBR pattern names extend by analogy:

    CON-persist(L; p, b)   ⇔ L ∈ F^{p,b,task_test}_w ∧ L ∈ F^{p,b,rest_post}_w ∧ L ∈ F^{p,b,rest_pre}.
    CON-reset(L; p, b)     ⇔ L ∈ F^{p,b,rest_pre} ∧ L ∉ F^{p,b,task_test} ∧ L ∉ F^{p,b,rest_post}.
    CON-rearrange(L; p, b) ⇔ L ∈ F^{p,b,rest_post}_w ∧ L ∉ F^{p,b,task_test} ∧ L ∉ F^{p,b,rest_pre}.

These four are mutually exclusive. Many leafsets fall into NONE
(e.g., a leafset present only in `task_test` and `task_learn` but in
none of the resting phases). We do not track unclassified leafsets.

### 2.3 Per-leaf projection

For each leaf `ℓ ∈ L_p`:

- Collect all `L ∈ F^{p,b,rest_post}_w` with `ℓ ∈ L` and
  `pattern(L) ≠ None`.
- Apply the priority chain trace > persist > reset > rearrange.
- Default `rearrange` if no consensus-classified leafset contains `ℓ`.

Same projection rule as legacy / TA-CBR / Containment-CBR.

### 2.4 Optional ε-relaxation

Exact equality is an extreme criterion. We **do not** include a
thresholded fallback in this scope (the whole point of CON is to be
strict). If exact CON-trace returns zero across the cohort, a follow-up
scope `2026-XX-XX_consensus-subtree-relaxed.md` can introduce
`L ≃_ε L'` ⇔ `|L ∆ L'| ≤ ε`, but until then `ε = 0`.

## 3. Properties

| property | value |
|:---|:---|
| Threshold-free | **Yes** — no `τ_match`, `τ_dissim`. |
| Symmetric in `(task_test, rest_post)` | **Yes**. |
| Asymmetric in `(rest_pre, post-task)` | **Yes** — rsPre absence is required. |
| Sensitive to single-leaf changes | **Extremely**. A leafset of 9 leaves in task that loses one leaf in rsPost (becoming an 8-leaf set) is **not** CON-trace, even if the 9-leaf rsPre set is dissimilar. This is the price of strictness. |
| Detects fragmentation | **No** — same limitation as legacy CBR. |
| Detects partial persistence | **No** by construction. |
| Cohort-level expectation | We expect *very few* CON-trace leafsets — likely 0–5 per (patient, band). The measure is a high-bar existence proof, not a census. |
| Complexity | `O(\|F^{p,b,task_test}_w\| · log \|F^{p,b,rest_pre}\| + \|F^{p,b,task_test}_w\| · log \|F^{p,b,rest_post}_w\|)` using frozenset hashing. Trivial. |

**Negative property — what CON cannot detect:**
- Anything where the set membership is even one leaf off.
- Any task subtree whose persistent rsPost manifestation is a
  refinement / coarsening of itself.
- Any "trace" that lives at sizes outside `[s_min, s_max]`.

CON is the *most rigorous* of the four scopes (legacy / TA / Containment
/ CON) but also the *least sensitive*. By design.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| Likely empty result. CON might find zero CON-trace leafsets across the entire cohort. | This is the answer if it happens. The scope's *negative* result is meaningful: it tells us that if we want a non-empty descriptive census we must accept some similarity tolerance (legacy / TA / containment) or some leaf-level slack (relaxed CON). |
| Silent failure if `s_max = 60` is too tight. | We additionally report the unconstrained variant `s_max = N_p − 1` in a sibling CSV — purely as a sensitivity check; the load-bearing column is the size-windowed one. |
| Hash collisions on frozenset of leaf integers. | None — Python's `frozenset` hashing on a fixed-domain integer set is exact. |
| Pat_14 lacks `task_test`. | Cell `n/a`. Same handling as TA-CBR. |
| Probe-bias: a probe-defined module that happens to be a small subtree could land in CON-persist, polluting the persist count without affecting CON-trace. | CON-persist count is reported but not load-bearing. CON-trace cannot be probe-driven (rsPre absence requirement). |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)
       s_min=5, s_max=60

F_rpre  ← { leaves(v) : v in internal_nodes(Z_rpre)  }
F_tl    ← { leaves(v) : v in internal_nodes(Z_tl)    }
F_tt    ← { leaves(v) : v in internal_nodes(Z_tt)    }
F_rpost ← { leaves(v) : v in internal_nodes(Z_rpost) }

(use Python frozensets so set membership is O(1) hash lookup)

F_tt_w   ← { L ∈ F_tt    : s_min ≤ |L| ≤ s_max }
F_rpost_w ← { L ∈ F_rpost : s_min ≤ |L| ≤ s_max }

per_leafset ← []
FOR each L in F_tt_w:
    in_rpost ← (L ∈ F_rpost_w)
    in_rpre  ← (L ∈ F_rpre)            (no size window — strict)
    label    ← classify_CON(in_rpost, in_rpre)
    APPEND (L, "task-anchored", label) to per_leafset

FOR each L in F_rpost_w:
    in_tt    ← (L ∈ F_tt_w)
    in_rpre  ← (L ∈ F_rpre)
    label    ← classify_CON(in_tt, in_rpre)
    APPEND (L, "rpost-anchored", label) to per_leafset
    (note: every CON-trace L will appear once with each anchor tag — dedupe at the per-leafset CSV)

leaf_label[1..N] ← "rearrange"
FOR ℓ ∈ {1..N}:
    cands ← { label : (L, _, label) ∈ per_leafset, label ≠ None, ℓ ∈ L,
                       L ∈ F_rpost_w }    (priority chain on rsPost-supported labels)
    APPLY priority trace > persist > reset > rearrange.

OUTPUT per_leafset, leaf_label
```

`classify_CON(in_rpost, in_rpre)` returns `trace` (in_rpost ∧ ¬in_rpre),
`persist` (in_rpost ∧ in_rpre), `reset` (¬in_rpost ∧ in_rpre), or None
when neither `in_rpost` nor `in_rpre` is true (those are pure-task or
pure-rearrange leafsets that we don't surface).

## 6. Visualization spec

One PDF per (patient, band):
`data/audit/per_patient_hierarchy_consensus/{Pat_XX}/{band}_hierarchy-crossphase.pdf`.

**Identical layout** to TA-CBR / Containment-CBR / audit_08 so all four
measures are comparable side by side. The key visual difference:
under CON, only leaves whose containing rsPost subtree's leafset
literally appears in `F^{p,b,task_test}` get the `trace` colour. We
expect *much sparser* red strips than under any of A / B / legacy.

CSVs:
- `leaf_assignment.csv` (same schema).
- `consensus_leafsets.csv` columns:
  `patient, band, leafset_size, leafset_hash, h_rel_in_taskT, h_rel_in_rsPost, h_rel_in_rsPre, label`.
  `h_rel_in_*` is the merge height of the corresponding internal node
  (NaN if the leafset is absent in that phase). `leafset_hash` is the
  string hex of `hash(frozenset(leaves))` truncated to 16 chars — for
  cross-(patient, band) collation only.

A separate **cohort census table** (markdown) at
`data/audit/per_patient_hierarchy_consensus/CENSUS.md`:

| patient | band | N CON-trace | N CON-persist | N CON-reset | N CON-rearrange |
|---|---|---|---|---|---|

with one row per (patient, band) and a totals row. Empty cells render
as "0".

## 7. Connection to prior tools

| existing tool | how CON relates |
|:---|:---|
| Legacy CBR | Strictest possible upper bound. Every CON-trace is also legacy-trace (Jaccard would be 1.0); the converse generally fails (legacy admits J ≥ 0.75). |
| TA-CBR | Same anchoring logic available. CON's task-anchored predicate is `L ∈ F_tt_w ∧ L ∈ F_rpost_w ∧ L ∉ F_rpre` — strictly stronger than TA-trace at any threshold ≥ 0.75. |
| Containment-CBR | CON does not use containment; the two are not nested. A leafset can be CON-trace (exact match in rsPost) AND also have C_rpre = 1 (it's contained in some giant rsPre subtree) — those two facts coexist without contradiction. |
| MRL | MRL's `match_τ` predicate at `τ = 1` is *almost* CON-trace, except MRL anchors only on `task_test` and uses `J*` rather than exact equality. With `τ = 1` and no other rounding, the predicates coincide on the task-anchored side. |
| H2c / H2d | Untouched. CON is per-leafset; H2c/H2d are per-pair. |

## 8. Implementation plan

- **Library reuse**: `tree_internal_nodes`. Frozenset hashing is built-in
  to Python; no new helper needed.
- **Script**: `scripts/01_compute/audit/audit_11_consensus_subtree.py`.
  Single pass over four trees per (patient, band); no nested Jaccard
  loop, so this is the cheapest of the four scopes computationally.
- **Run order**: only after A and B (TA-CBR + Containment-CBR) are on
  screen and judged, per user instruction (2026-04-25). Implementation
  task tracked separately (#18).
- **Outputs**: as listed in §6.

## 9. Open questions

- **Should we implement the relaxed `ε`-version eagerly?** Defer until
  we see whether CON-trace is empty across the cohort. If empty for
  ≥ 50 % of (patient, band) cells, the relaxed version becomes the
  load-bearing measure and gets its own scope.
- **Should the "absent in rsPre" requirement be size-windowed too?**
  Currently no — any rsPre internal node with this leafset (at any
  size) disqualifies. A size-window on the rsPre check would relax
  the requirement; defer.
- **Cross-(patient, band) leafset-hash collation**: a leafset
  `{3, 7, 12, 19, 22}` that is CON-trace in (Pat_05, β) and also
  CON-trace in (Pat_05, γ_l) is interesting. Cross-band consistency of
  CON-trace leafsets is an obvious follow-up; not in this scope.
