---
name: task-anchored-cbr
type: scope
era: COHORT_N9
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - scripts/archive/2026-04_failed-scalar-session/gallery_cbr_modules.py
  - scripts/01_compute/audit/audit_08_per_patient_hierarchy.py
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Task-anchored Cluster-Birth-Retention (TA-CBR)

**A descriptive per-(patient, band) classifier that anchors on every
internal node of the `task_test` dendrogram (rather than the `rest_post`
one) and asks whether its leafset (a) is *not* a coherent subtree in
`rest_pre` and (b) *is* a coherent subtree in `rest_post`. Same Jaccard
machinery as the legacy CBR; the only change is the anchor — flips the
asymmetry so task subtrees that fragment into multiple smaller rsPost
subtrees become detectable.**

## 1. Notation

Inherits the cohort/band/phase notation of the MRL scope report
(§1 of `2026-04-25_module-retention-landscape.md`). Recap:

- `P, B, Φ` — patients, bands, phases (`COHORT_N9`).
- `Z^{p,b,φ}` — `scipy` linkage matrix from
  `workflow.lrg.compute_lrg_analysis(imcoh_abs)`.
- `T^{p,b,φ}` — induced binary tree, `N_p` leaves, `N_p − 1` internal nodes.
- For each internal node `v ∈ T^{p,b,φ}`: `leaves(v) ⊆ L_p`, `size(v)`,
  `h(v)`, `h_rel(v) = h(v) / dmax(T^{p,b,φ})`.

For the legacy (rsPost-anchored) CBR we will write
`Internal_rsPost(p, b) := Internal(T^{p,b,rest_post})`. For TA-CBR the
anchor set is

    Internal_taskT(p, b) := Internal(T^{p,b,task_test}).

Fixed thresholds (inherited verbatim from the legacy classifier so
results are directly comparable):

- `τ_match := 0.75`
- `τ_dissim := 0.50`
- `s_min := 5`,  `s_max := 60`

## 2. Definitions

### 2.1 Best-match Jaccard (unchanged)

For target leafset `S ⊆ L_p` and reference tree `T`,

    J*(S, T) := max_{v ∈ Internal(T)} |leaves(v) ∩ S| / |leaves(v) ∪ S|.

`J* : 2^{L_p} × Trees(L_p) → [0, 1]`. Symmetric in `S`-vs-`leaves(v)`.
Implementation: `jaccard_leafsets` from `lrg_eegfc.utils.metrics.tree`.

### 2.2 TA-CBR per-anchor predicates

Given an anchor `η ∈ Internal_taskT(p, b)` with `s_min ≤ size(η) ≤ s_max`,
define

    J_rpre(η)  := J*(leaves(η), T^{p,b,rest_pre})
    J_rpost(η) := J*(leaves(η), T^{p,b,rest_post})
    J_tl(η)    := J*(leaves(η), T^{p,b,task_learn))      (corroboration only)

The TA-CBR pattern of `η` is the unique label assigned by the priority
chain (mutually exclusive in this order):

    TA-trace(η)     ⇔ J_rpre(η)  ≤ τ_dissim ∧ J_rpost(η) ≥ τ_match
    TA-persist(η)   ⇔ J_rpre(η)  ≥ τ_match  ∧ J_rpost(η) ≥ τ_match
    TA-reset(η)     ⇔ J_rpre(η)  ≥ τ_match  ∧ J_rpost(η) ≤ τ_dissim
    TA-rearrange(η) ⇔ J_rpre(η)  ≤ τ_dissim ∧ J_rpost(η) ≤ τ_dissim

Anchors that satisfy none of the four (e.g. `J_rpre = 0.6`, `J_rpost = 0.6`)
are labelled `None` and excluded from the per-anchor census.

`J_tl(η)` is reported as a corroboration column but does **not** enter the
classifier — anchoring already lives on `task_test`, so requiring
`J_tl(η) ≥ τ_match` would over-restrict (the two task phases differ in
duration and content; there is no a-priori reason their LRG dendrograms
must contain near-identical subtrees).

### 2.3 Per-leaf projection

For each leaf `ℓ ∈ L_p`, scan all anchors `η ∈ Internal_taskT(p, b)` with
`s_min ≤ size(η) ≤ s_max` and `ℓ ∈ leaves(η)`. Collect the multiset of
TA-CBR labels. Assign `ℓ` the highest-priority pattern present
(trace > persist > reset > rearrange), defaulting to `rearrange` if no
classified anchor contains `ℓ`. Identical projection rule to
`audit_08_per_patient_hierarchy.py`, but the underlying anchor set is
`Internal_taskT(p, b)` rather than `Internal_rsPost(p, b)`.

## 3. Properties

| property | value |
|:---|:---|
| Range of each `J_*` | `[0, 1]` |
| TA-trace is symmetric in `(rest_pre, rest_post)` | **No** — pre and post enter via different predicates (low / high). |
| TA-trace is symmetric in `(task_learn, task_test)` | **No** — anchored on `task_test` only; `task_learn` is corroboration. |
| Detects task subtrees that fragment in rsPost | **Yes** if at least one rsPost subtree retains ≥ 75 % overlap with the task anchor — but it still requires a single rsPost match to clear `τ_match`. Multi-fragment rsPost re-aggregation is **not** detected. |
| Detects rsPost subtrees that have no task counterpart | **No**, by construction — the procedure only iterates over `Internal_taskT`. |
| Identifiability | A task anchor `η` may share label `TA-trace` with a different anchor `η'` if both their leafsets are well-matched in rsPost. Patterns are **per-anchor**, not per-leaf-set; correlated anchors are an expected feature, not a bug. |
| Complexity | O(\|Internal_taskT\| · \|Internal_rsPost\| + \|Internal_taskT\| · \|Internal_rsPre\|) Jaccard evaluations per (patient, band). For `N ≈ 120`: ~28k Jaccards × 60 cells = ~1.7 million per cohort. |

**Negative property — what TA-CBR cannot detect:** rsPost-only modules
that have no task-test counterpart but did not exist in rsPre either
(those would only show as `TA-trace` if a *different* task anchor happens
to land on the same leafset; otherwise invisible). Such modules are
reachable by complementary rsPost-anchored CBR — the two procedures
together cover both anchor directions.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| Anchor change does not address Jaccard's containment penalty (failure 2 from the parent menu). | Pair with measure B (containment CBR) when judging — they address orthogonal failures. |
| Best-match `J*` still uses a single reference subtree. A task anchor whose leaves are split between two adjacent rsPost subtrees is invisible. | Acknowledged; addressed only by D (consensus-subtree) and a future top-K extension. |
| Pat_14 lacks `task_test` → no anchors → empty census. Mark cell `n/a` in tally. | Out-of-scope deliverables flag. |
| Probe-bias (CLAUDE.md invariant 5): subtrees at coarse `h_rel` may be probe-defined rather than functional. The current `s_max = 60` cap already excludes the very coarsest. Scope inherits this; do **not** raise `s_max` without a probe-debias pass. | Keep `s_max = 60`. Report per-anchor `h_rel(η)` in the gallery so the reader can spot near-root anchors. |
| Symmetric Jaccard ties when `|leaves(η)| = |leaves(m)|`. With `N_p ≈ 120` and integer set sizes, exact ties are rare but possible. | First-encountered tie wins (deterministic by node id). Document in the per-anchor CSV. |
| Asymmetry in (rsPre, rsPost) means TA-trace count differs from the count of "rsPost subtrees not in rsPre and in task" — they are different predicates. | Always report which anchor was used (rsPost vs. task) in figure titles and CSV headers. |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)
       τ_match=0.75, τ_dissim=0.50, s_min=5, s_max=60

I_rpre   ← internal_nodes(Z_rpre)
I_tl     ← internal_nodes(Z_tl)
I_tt     ← internal_nodes(Z_tt)
I_rpost  ← internal_nodes(Z_rpost)

per_anchor ← []
FOR each η in I_tt:
    IF size(η) < s_min OR size(η) > s_max: CONTINUE
    J_rpre  ← max over m in I_rpre  of jaccard(leaves(η), leaves(m))
    J_rpost ← max over m in I_rpost of jaccard(leaves(η), leaves(m))
    J_tl    ← max over m in I_tl    of jaccard(leaves(η), leaves(m))
    label   ← classify_TACBR(J_rpre, J_rpost, τ_match, τ_dissim)
    APPEND (η, J_rpre, J_rpost, J_tl, label) to per_anchor

leaf_label[1..N] ← "rearrange"  (default)
FOR ℓ ∈ {1..N}:
    cands ← {label : (η, …, label) ∈ per_anchor, label ≠ None, ℓ ∈ leaves(η)}
    IF "trace"   ∈ cands: leaf_label[ℓ] ← "trace"
    ELIF "persist" ∈ cands: leaf_label[ℓ] ← "persist"
    ELIF "reset"   ∈ cands: leaf_label[ℓ] ← "reset"
    ELIF "rearrange" ∈ cands: leaf_label[ℓ] ← "rearrange"

OUTPUT per_anchor, leaf_label
```

`classify_TACBR` is the priority chain from §2.2. Identical four-pattern
schema as the legacy classifier so the colours and codebar from
`audit_08_per_patient_hierarchy.py` remain meaningful.

## 6. Visualization spec

One PDF per (patient, band) at
`data/audit/per_patient_hierarchy_taskanchored/{Pat_XX}/{band}_hierarchy-crossphase.pdf`.

Layout (identical to audit_08 so the two procedures are directly
visually comparable):

- 4 panels in a row, one per phase in canonical order
  (`rest_pre, task_learn, task_test, rest_post`).
- Each panel: dendrogram + codebar strip beneath.
- **Branches**: an internal-node link gets the colour of its category iff
  every leaf below it shares the same category; else `BACKGROUND_GRAY`.
- **Leaf stems**: overdrawn in the leaf's own category colour
  (so isolated classified leaves remain visible against gray parent links).
- **Codebar**: one rectangle per leaf (10 × 1 units), colour = leaf
  category. Leaf order = scipy `dendrogram(...)['leaves']`.
- **Colour palette** (user-pinned 2026-04-25):
  - trace = `#d62728` (red)
  - persist = `#000000` (black)
  - reset = `#7f7f7f` (gray)
  - rearrange = `#1f77b4` (blue)
- Legend at the figure top with title
  `"{Pat_XX} — {band} — TA-CBR (task_test-anchored, size ∈ [5, 60])"`.
- No `fig.suptitle` (CLAUDE.md never-list).

**Reading rule**: same RED columns clustering tightly under one phase's
dendrogram and scattered under another reveals the same task-induced /
rsPost-persisted subtrees the legacy CBR was designed to find — but now
detected even when the rsPost manifestation is the smaller of the two
matched leafsets.

Auxiliary CSV at the directory root:
`leaf_assignment.csv` with columns
`patient, band, leaf_id, category` (same schema as legacy).

Per-anchor CSV (new, since the anchor universe changes):
`anchor_classification.csv` with columns
`patient, band, anchor_row, size, h_rel, J_rpre, J_rpost, J_tl, label`.
One row per anchor `η` of `T^{p,b,task_test}` with `size ∈ [s_min, s_max]`.

## 7. Connection to prior tools

| existing tool | how TA-CBR relates |
|:---|:---|
| Legacy CBR (`gallery_cbr_modules.py`, `audit_08`) | Same Jaccard, same thresholds, same priority chain. Different anchor (`task_test` vs. `rest_post`). The two are **complementary**, not nested: a leafset can be `TA-trace` without being `legacy-trace` and vice versa. |
| Module-Retention Landscape (`2026-04-25_module-retention-landscape.md`) | MRL is the **scalar field aggregate** of TA-CBR: it counts, per `(p, b, h_rel)`, the fraction of task-anchored subtrees with `J_rpre ≤ 0.5 ∧ J_rpost ≥ 0.75`. TA-CBR is the **per-anchor** view of the same predicate, plus the persist/reset/rearrange siblings and a per-leaf projection for visual inspection. |
| H2c (continuous ultrametric drift) | H2c is global pairwise; TA-CBR is local subtree-level. H2c says "directions agree in aggregate"; TA-CBR points at *which* subtrees agree. Complementary. |
| H2d (block coactivation persistence) | H2d is per-pair-per-scale; TA-CBR is per-subtree. H2d aggregates over leaf pairs; TA-CBR aggregates over leafsets. Complementary. |
| Containment CBR (B) | Drops Jaccard for asymmetric containment. Can be combined with anchoring (B itself can be anchored on either tree). The combined "task-anchored containment" is a strict superset of TA-CBR with `τ_dissim, τ_match` interpreted on containment scores. |

## 8. Implementation plan

- **Library reuse**: `jaccard_leafsets`, `tree_internal_nodes` from
  `src/lrg_eegfc/utils/metrics/tree.py`. No new helpers required.
  Per-leaf projection lifted from
  `scripts/01_compute/audit/audit_08_per_patient_hierarchy.py:leaf_categories`.
- **Script**: `scripts/01_compute/audit/audit_09_task_anchored_cbr.py`.
  Mirrors the structure of `audit_08`. Only differences: (i) anchor scan
  iterates `I_tt` instead of `I_rpost`; (ii) classifier input is
  `(J_rpre, J_rpost)` instead of `(J_rpre, J_tl, J_tt)`; (iii) emits the
  `anchor_classification.csv`.
- **Outputs**:
  - `data/audit/per_patient_hierarchy_taskanchored/{Pat_XX}/{band}_hierarchy-crossphase.pdf`
  - `data/audit/per_patient_hierarchy_taskanchored/leaf_assignment.csv`
  - `data/audit/per_patient_hierarchy_taskanchored/anchor_classification.csv`
- **Pat_14** has no `task_test` LRG cache → cell `n/a`; emit a `missing.txt`
  listing skipped (patient, band) pairs.
- **No promotion to library** in this pass; if a third caller appears the
  helper goes to `src/lrg_eegfc/utils/metrics/tree.py`.

## 9. Open questions

- **Should `J_tl` corroboration become a hard requirement?** Currently
  reported, not classified on. If we want stricter TA-trace ("anchored on
  task_test, also coherent in task_learn, persists in rsPost"), the
  classifier becomes 3-Jaccard like the legacy. Decision deferred until
  A vs. legacy comparison is on screen.
- **Per-anchor weighting in the cohort census**: should anchor-row
  contributions be uniform, size-weighted, or `pi`-weighted (where `pi`
  is the prominence used in `gallery_cbr_modules.py`)? Default for the
  first run: uniform (each anchor counts once). Revisit after sight.
- **Threshold sweep**: the menu (failure 4) flagged hard thresholds as a
  problem. This scope keeps `(τ_match=0.75, τ_dissim=0.50)` from the
  legacy classifier so A is directly comparable. A threshold-sweep
  variant `TA-CBR(τ_match, τ_dissim)` is a separate scope (deferred).
