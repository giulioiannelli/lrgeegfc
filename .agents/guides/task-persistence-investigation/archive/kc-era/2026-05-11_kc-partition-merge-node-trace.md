---
name: kc-partition-merge-node-trace
type: scope_report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-11
updated: 2026-05-11
pointers:
  - .agents/guides/02_methods/lrg-framework-guide.md
  - scripts/01_compute/figures_embedded/fig_preprint_dendrograms_kc_trace.py
---

# KC partition merge-node trace

## Head

Per-merge-node decomposition of the KC λ-mixed trace sum. Each leaf
pair `(i,j)` is credited to **exactly its MRCA node** in the reference
phase's tree (not to every ancestor subtree). The aggregate over all
merge nodes is a strict partition of `∑_{(i,j)} s(i,j)`. Sharpens the
existing `dendrograms_kc_trace` figure by removing upward smearing of
deep-pair scores into ancestor subtrees.

## Notation

- `Z_X` — linkage matrix of phase `X ∈ {rPre, tLearn, tTest, rPost}`.
- `n` — number of leaves (constant across phases for one patient).
- `m_X(i,j)` — MRCA-edge-count for leaf pair `(i,j)` in `X` (KC `m`
  vector, integer).
- `M_X(i,j)` — MRCA merge height (KC `M` vector, float).
- `v_X^λ(i,j) = (1−λ)·m_X(i,j)/m_max + λ·M_X(i,j)/M_max` — KC
  λ-mixed pair vector, normalized using pooled maxima across all 4
  phases for the patient (Kendall & Colijn 2016, `normalize=True`).
- `s(i,j) = ½(|v_pre^λ − v_test^λ| + |v_pre^λ − v_post^λ|) −
  |v_test^λ − v_post^λ|` — per-pair trace score, tree-invariant.
- For internal node `N` of `Z_X` with children `L_left(N)`, `L_right(N)`:
  pairs *first-merging at N* = `L_left(N) × L_right(N)`.

## Predicate (per-merge-node trace)

```
trace_X(N) =  mean over (i,j) ∈ L_left(N) × L_right(N)  of  s(i,j)
```

A pair `(i,j)` belongs to exactly one merge node — its MRCA in `Z_X`.

## Properties

- **Partition** — `∪_N (L_left(N) × L_right(N)) = {(i,j) : i < j}` and
  the union is disjoint. So `∑_N |L_left(N)|·|L_right(N)|·trace_X(N)
  = ∑_{(i,j)} s(i,j)`. Per-merge-node mean is therefore a faithful
  decomposition of the full per-pair KC trace sum into local
  contributions; per-subtree mean (existing kc_trace figure) is an
  overlapping aggregation where each pair contributes to every
  ancestor subtree.
- **Range** — same as `s(i,j)` itself, `≈ [-1, +1]` after
  normalization.
- **Reference dependence** — `trace_X(N)` is computed on phase `X`'s
  tree topology; the same pair appears under different merge nodes
  in different phases.
- **Symmetric coloring** — render with `RdBu_r` divergent at
  symmetric `vlim = ±max|trace_X(N)|` pooled across phases (one
  vlim per page).

## Caveats

- **Cross-product size bias** — small `|L_left|·|L_right|` (near
  the leaves) averages few pairs and is noisier. Near the root, the
  cross-product is huge and averages many pairs (low variance, close
  to the global per-pair mean). The visualization will naturally
  show more saturation near the leaves and more "neutral" coloring
  near the root.
- **Not a quantification of `n_trace`** — this is a *visualization
  primitive*, not a hypothesis-test. Counts/statistics over
  "trace-strong merge nodes" require a within-baseline null
  comparison (analogous to the `audit_47` machinery).
- **λ > 0 mixes heights** — at λ = 0.5 the normalized `M`
  contribution adds branch-length sensitivity. The partition property
  still holds; only the per-pair `s` definition changes.
- **rPre interpretation** — `trace_pre(N) > 0` means: pairs first-
  merging at `N` in rPre meet at *different* depths in tTest and
  rPost (but tTest and rPost agree with each other). I.e., `N` is a
  rPre-specific merge "undone" by task and post.

## Pseudocode

```python
def leaves_under_each_node(Z):
    n = Z.shape[0] + 1
    out = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a, b = int(Z[i, 0]), int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out

def per_merge_node_score(Z, s_mat):
    n = Z.shape[0] + 1
    leaves = leaves_under_each_node(Z)
    out = {}
    for i in range(n - 1):
        L = list(leaves[int(Z[i, 0])])
        R = list(leaves[int(Z[i, 1])])
        out[n + i] = float(s_mat[np.ix_(L, R)].mean())
    return out
```

## Visualization

- Same layout as `dendrograms_kc_trace`: 4 panels (rPre, tLearn,
  tTest, rPost), shared y, probe-shaft barcode below each.
- Each U-shape (merge node) colored by `trace_X(N)`,  RdBu_r,
  symmetric vlim pooled across phases.
- `above_threshold_color` = neutral gray (root-trivial node).
- Bottom barcode: probe (shaft) identity, same color per leaf in all
  panels.

## Connection to prior tools

- **`dendrograms_kc_trace`** (existing): per-subtree mean (each pair
  counted in every ancestor → smearing).
- **`dendrograms_subtree_jaccard_trace`** (existing): set-overlap of
  clade leaf-sets — orthogonal axis (clade structure, not pair
  topology).
- **`dendrograms_matched_trace`** (existing): exact clade equality —
  binary partition view of `dendrograms_subtree_jaccard_trace`.
- **This (new)**: KC partition. Sharpest pair-decomposition view
  of the existing per-pair KC trace.

## Open questions

- Does the partition view visually localize the trace more clearly
  than per-subtree mean on the cohort? Comparison page-by-page on
  the same patient × band TBD.
- Cross-product size bias near the leaves — should we weight by
  `|L_left|·|L_right|` (= total contribution to the KC sum) instead
  of taking the mean? Could surface a *contribution* view rather than
  a *normalized-rate* view.

## Output

- Script: `scripts/01_compute/figures_embedded/fig_preprint_dendrograms_kc_partition_trace.py`
- Folder: `data/reports/preprint/figure1_beta_trace/dendrograms_kc_partition_trace/`
- Filename: `kc_partition_trace_lambda{lam}_{band}_all_patients.pdf`
- λ supported: 0, 0.5 (via `--lam`).
