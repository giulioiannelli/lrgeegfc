---
name: containment-cbr
type: scope
era: COHORT_N9
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Containment-CBR

**Replaces the symmetric Jaccard in the CBR classifier with a directional
containment score `C(S, T) := max_v |S ∩ leaves(v)| / |S|` that asks
"how much of the anchor leafset is contained in *some* subtree of the
reference tree", regardless of how much extra material that subtree
carries. Restores sensitivity to small task subtrees that are clean
subsets of larger rsPost subtrees (and vice versa).**

## 1. Notation

Inherits §1 of `2026-04-25_task-anchored-cbr.md` and the MRL scope.
Recap:

- `P, B, Φ, T^{p,b,φ}, leaves(v), size(v), h_rel(v)` as before.
- The anchor tree is parameterised: `T_anchor ∈ {T^{p,b,task_test},
  T^{p,b,rest_post}}`. Both anchorings are tabulated; the containment
  rewrite is independent of the anchor.

Fixed thresholds (kept identical to legacy CBR for direct comparability;
re-interpretation differs because the score is now directional):

- `τ_match := 0.75` (containment now means: ≥ 75 % of the anchor's
  leafset lives inside *some* reference subtree).
- `τ_dissim := 0.50` (≤ 50 % of the anchor's leafset lives inside any
  single reference subtree → "spread out / no coherent counterpart").
- `s_min := 5`, `s_max := 60`.

## 2. Definitions

### 2.1 Directional containment

For a target leafset `S ⊆ L_p` and a tree `T` over the same leaf universe,

    C(S, T) := max_{v ∈ Internal(T)} |S ∩ leaves(v)| / |S|.

`C : 2^{L_p} × Trees(L_p) → [0, 1]`. Asymmetric (does *not* equal
`C(leaves(v), T_S)` in general).

Range and limits:

- `C(S, T) = 1` ⇔ `S ⊆ leaves(v)` for some `v ∈ Internal(T)`. In
  particular, the root subtree always satisfies this trivially, so we
  must restrict `v` to `size(v) ≤ s_max + κ` (we use `κ = N_p / 2` —
  i.e. forbid the matching subtree from spanning more than half the
  tree, otherwise the root trivially gives `C = 1`).
- `C(S, T) = 0` ⇔ no internal node of `T` (subject to the size cap
  above) contains a single leaf of `S`. With our `κ`, the leaves of `S`
  must collectively avoid the entire upper half of the tree's
  internal-node coverage — extremely rare.

### 2.2 Containment-CBR predicates

For an anchor `η` (in either `Internal_taskT(p, b)` or
`Internal_rsPost(p, b)` — both versions are computed and emitted):

    C_rpre(η)  := C(leaves(η), T^{p,b,rest_pre},  capped by κ)
    C_rpost(η) := C(leaves(η), T^{p,b,rest_post}, capped by κ)
    C_tl(η)    := C(leaves(η), T^{p,b,task_learn}, capped by κ)
    C_tt(η)    := C(leaves(η), T^{p,b,task_test}, capped by κ)

For the **rsPost-anchored** variant (analog of the legacy CBR; only the
score is changed):

    C-trace_rsPost(η)     ⇔ C_rpre(η)  ≤ τ_dissim ∧ min(C_tl, C_tt)(η) ≥ τ_match
    C-persist_rsPost(η)   ⇔ min(C_rpre, C_tl, C_tt)(η) ≥ τ_match
    C-reset_rsPost(η)     ⇔ C_rpre(η)  ≥ τ_match  ∧ min(C_tl, C_tt)(η) ≤ τ_dissim
    C-rearrange_rsPost(η) ⇔ max(C_rpre, C_tl, C_tt)(η) ≤ τ_dissim

For the **task-anchored** variant (analog of TA-CBR):

    C-trace_taskT(η)     ⇔ C_rpre(η)  ≤ τ_dissim ∧ C_rpost(η) ≥ τ_match
    C-persist_taskT(η)   ⇔ C_rpre(η)  ≥ τ_match  ∧ C_rpost(η) ≥ τ_match
    C-reset_taskT(η)     ⇔ C_rpre(η)  ≥ τ_match  ∧ C_rpost(η) ≤ τ_dissim
    C-rearrange_taskT(η) ⇔ C_rpre(η)  ≤ τ_dissim ∧ C_rpost(η) ≤ τ_dissim

`None` if no chain branch fires. Per-leaf projection identical to
legacy / TA-CBR (priority chain, default `rearrange`).

### 2.3 Why the size cap κ is necessary

Without a cap, the root of `T` always gives `C(S, T) = 1` because
`leaves(root) = L_p ⊇ S` for any `S ⊆ L_p`. Then every anchor would
trivially satisfy `C_rpre = C_rpost = 1`, and only `persist` would
ever fire.

The cap `size(v) ≤ κ` prevents any single ridiculously large reference
subtree from absorbing the anchor and is the natural directional
analogue of the implicit cap that Jaccard imposes via the union-set
denominator. Default `κ = ⌊N_p / 2⌋`. Reported in the per-anchor CSV
so a reader can see how often the cap binds.

## 3. Properties

| property | value |
|:---|:---|
| Range | `[0, 1]` |
| Symmetric | **No** — `C(S, T) ≠ C(leaves(v), T_S)` in general. |
| Monotonicity in subtree size | `C(S, T)` weakly increases as we consider larger reference subtrees (subject to the cap). |
| Detects small-anchor-in-larger-reference cases | **Yes**. This is the primary motivation. |
| Detects fragmentation (anchor split across multiple reference subtrees) | **No** by construction — `C` looks at one reference subtree at a time. Pair with multi-match extensions in a future scope. |
| Cap binds frequently? | Empirical question. We log how often the unconstrained max is achieved at `size(v) > κ` so we can revisit. |
| Complexity | Same as Jaccard: `O(\|I_anchor\| · Σ_φ \|I_φ\|)` per (patient, band). |

**Negative property — what containment-CBR cannot detect:**
- It cannot tell whether the *unmatched* part of a reference subtree is
  meaningful. A 50-leaf reference subtree containing a 5-leaf anchor
  scores `C = 1`, but the other 45 leaves are unaccounted-for.
- It cannot detect fragmentation (see above).
- It cannot distinguish between "anchor sits inside a tiny refinement
  of a reference subtree" (clean containment) and "anchor sits inside
  a sprawling junk subtree" (over-containment) — a refinement of `C`
  with a denominator on the *union* would address this but is a
  separate scope.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| Cap κ is a parameter, not a derivation. | Default `⌊N_p/2⌋`. Report `cap_bound` boolean per anchor. If many anchors saturate the cap, raise the issue in Open Questions. |
| Task-anchored containment-trace will count anchors that are subsets of the rsPost root subgraph but not in any compact rsPost cluster. | Mitigated by the cap. If suspicious anchors appear, inspect their `J_rpost` (Jaccard) value side-by-side — should be low if Jaccard says "no cluster", containment is then recording "the leaves are spread out across small rsPost siblings". |
| Containment is not a metric (no triangle inequality, no symmetry). All inferences are descriptive. | Acknowledged. Falsifiability tests live at the visualisation / cohort census level, not at a hypothesis-test level. |
| Threshold reuse from legacy is an analogy, not an equivalence. `J ≥ 0.75` and `C ≥ 0.75` are different statements. | Report both numerically; do not equate the two cell counts across measures. |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)
       τ_match=0.75, τ_dissim=0.50, s_min=5, s_max=60
       κ = floor(N / 2)

FUNCTION containment(S, Z, κ):
    best ← 0
    FOR each internal-node v in Z:
        IF size(v) > κ: CONTINUE
        cur ← |S ∩ leaves(v)| / |S|
        IF cur > best: best ← cur
    RETURN best

FOR variant in {rsPost-anchored, task_test-anchored}:
    I_anchor ← internal_nodes(Z_rpost) if rsPost-anchored else internal_nodes(Z_tt)

    per_anchor ← []
    FOR each η in I_anchor:
        IF size(η) < s_min OR size(η) > s_max: CONTINUE
        C_rpre  ← containment(leaves(η), Z_rpre,  κ)
        C_rpost ← containment(leaves(η), Z_rpost, κ)
        C_tl    ← containment(leaves(η), Z_tl,    κ)
        C_tt    ← containment(leaves(η), Z_tt,    κ)
        label   ← classify_C(variant, C_rpre, C_rpost, C_tl, C_tt, τ_match, τ_dissim)
        APPEND (η, C_rpre, C_rpost, C_tl, C_tt, label, variant) to per_anchor

    PROJECT per_anchor → leaf_label[1..N] using priority chain
    EMIT figure + leaf_assignment_<variant>.csv + anchor_classification_<variant>.csv
```

`classify_C` uses §2.2's predicate sets selected by `variant`.

## 6. Visualization spec

Two PDF families per (patient, band), one per anchor variant:

- `data/audit/per_patient_hierarchy_containment/rpost_anchor/{Pat_XX}/{band}_hierarchy-crossphase.pdf`
- `data/audit/per_patient_hierarchy_containment/task_anchor/{Pat_XX}/{band}_hierarchy-crossphase.pdf`

Layout, palette, codebar, leaf-stem overdraw rule, and reading rule are
**identical** to the TA-CBR / audit_08 spec — directly comparable side
by side. The legend title carries the variant tag, e.g.
`"{Pat_XX} — {band} — Containment-CBR (rsPost-anchored, κ = ⌊N/2⌋)"`.

CSVs:
- `leaf_assignment.csv` per anchor variant (same schema as TA-CBR).
- `anchor_classification.csv` columns:
  `patient, band, anchor_row, size, h_rel, C_rpre, C_rpost, C_tl, C_tt, cap_bound_rpre, cap_bound_rpost, cap_bound_tl, cap_bound_tt, label, variant`.

## 7. Connection to prior tools

| existing tool | how Containment-CBR relates |
|:---|:---|
| Legacy CBR | Same anchor (rsPost) and same priority chain in the rsPost-anchored variant; only the score function changes (`J → C`). Side-by-side comparison shows whether the legacy "no trace" verdicts in Pat_08/10/13/15 reflect genuine absence or Jaccard's containment penalty. |
| TA-CBR | Same anchor (task_test) for the task-anchored variant; same diagnostic. The two scopes (A and B) are orthogonal axes — anchor and similarity — and the four combinations are all reported. |
| MRL | MRL uses Jaccard via `J*`. A containment variant of MRL is conceivable but out-of-scope here; if Containment-CBR wins, MRL gets a follow-up scope `2026-XX-XX_module-retention-landscape-containment.md`. |
| H2c / H2d | Untouched. Containment is a per-subtree replacement for Jaccard, not a per-pair statistic. |

## 8. Implementation plan

- **Library reuse**: only `tree_internal_nodes`. The `containment(...)`
  helper is new (~6 lines). Per
  `.agents/guides/04_rules/coding-rules.md`, it ships inline in the
  audit script for the first use; second caller triggers promotion to
  `src/lrg_eegfc/utils/metrics/tree.py:containment_score`.
- **Script**: `scripts/01_compute/audit/audit_10_containment_cbr.py`.
  Both anchor variants in one run; loop over `variant ∈ {"rpost",
  "task_test"}`.
- **Outputs**: as listed in §6.
- **Pat_14**: skip task-anchored variant; rsPost-anchored variant is
  fine (rsPost cache present).
- **Cap κ**: hardcoded `floor(N / 2)` in this scope; if cap-bound rate
  > 50 %, report and revisit.

## 9. Open questions

- **Should the cap κ depend on `s_max`?** A natural alternative is
  `κ := max(s_max, ⌊N_p / 2⌋)`. Defer until first results land.
- **Is the legacy `min(C_tl, C_tt) ≥ τ_match` requirement on the
  rsPost-anchored variant too strict?** With containment scores, the
  conjunction is weaker than with Jaccard (because `C` is more
  permissive). Decision deferred.
- **Should we report the containment score's complement
  `\|leaves(v) \\ S\| / \|leaves(v)\|`?** It would let us distinguish
  "clean containment" (anchor exhausts the reference subtree) from
  "loose containment" (anchor is a small minority of the reference).
  Defer; if the figures don't disambiguate, add it.
