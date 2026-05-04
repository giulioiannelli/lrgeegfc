---
name: trace-modules
type: scope
era: COHORT_N10
status: draft
created: 2026-04-25
updated: 2026-04-26
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
  - data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf
  - src/lrg_eegfc/utils/metrics/tree.py
  - src/lrg_eegfc/utils/metrics/tree_distance.py
---

# Trace-modules — module identification at significant `(band, k)` cells

**For each `(band, k_bin)` cell that the headline figure
(`task_trace_band_k_n10_imcoh_abs.pdf`) flags as cohort-wide (≥ 8/10)
*or* cluster-permutation significant, enumerate the actual leafsets
that drive the trace per patient. Apply the canonical T-regime
predicate (`¬present_pre ∧ present_task ∧ present_post` at
`J_min = 0.9`, matches MRL) and emit one CSV row per (patient, band,
k, T-subtree). Then render a cohort summary PDF and per-patient
mini-dendrograms (3 stacked: pre / task_test / post) so the
`(band × k)` cell is grounded in concrete electrode subsets, not
abstract numbers.**

---

## 1. Notation

Reuses the canonical reformalization
(`2026-04-25_task-trace-canonical.md`) verbatim. Critical symbols:

| Symbol | Meaning |
|---|---|
| `T^{p,b,φ}` | LRG tree (linkage matrix) per patient × band × phase |
| `Internal(T)` | non-leaf, non-root nodes of `T` |
| `leaves(v)` | leaf-set of internal node `v` |
| `J(A, B) = |A∩B|/|A∪B|` | symmetric Jaccard |
| `J*(S, T) = max_u J(leaves(u), S)` | best-match Jaccard `S` against tree |
| `J_min = 0.9` | match threshold (locked, matches MRL) |
| `k_min = 3` | subtree-size floor (locked, matches MRL) |
| `K_sig(b)` | set of integer-k values flagged sig in the headline for band `b` |

`K_sig(b)` is derived from the headline's per-(band, k) sign-flip
cluster-permutation results — the cluster runs whose mass exceeds the
95% null threshold (sidecar lists thresholds per panel).

## 2. Predicates

For each candidate subtree `v ∈ V_task(p, b)` evaluated against the
patient's `Z_pre` and `Z_post`:

```
present_pre(v)  ⇔ J*(leaves(v), T_pre)  ≥ J_min
present_post(v) ⇔ J*(leaves(v), T_post) ≥ J_min
T(v)            ⇔ ¬present_pre ∧ present_post     (task ⇒ post, novel)
```

Note: `present_task` is implicit — `v` is enumerated from
`Internal(T_task_learn) ∪ Internal(T_task_test)`, so `v` is by
construction "present in some task tree". The T regime is therefore
fully specified by the `pre`/`post` predicates.

The canonical doc names four regimes (P/T/R/RA); this scope report
focuses on **T only** (the central scientific claim). P/R/RA are
emitted as auxiliary fields in the same CSV row so future
cross-checks have the data without recomputation.

## 3. Properties

| # | Statement |
|---|---|
| 3.1 | Range — `T(v) ∈ {0, 1}` per row. |
| 3.2 | Independence of bands — per-band processing; no cross-band pooling. |
| 3.3 | Reuses MRL primitives — `tree_internal_nodes`, `jaccard_leafsets`, `dmax_from_Z` from `lrg_eegfc.utils.metrics.tree`. No new helpers. |
| 3.4 | Cell selection is **driven by the headline figure** (Stage 2b). Significant (band, k_bin) cells are extracted from the cluster-perm output stored in the headline's sidecar. |
| 3.5 | NOT a statistical test — descriptive enumeration only. Cohort-wide claim per cell already gated by the headline (≥ 8/10 unanimity OR cluster-perm sig). |

### 3.6 What this measure does NOT show

- **Direction of pair-level coactivation drift** — that is H2d.
- **Continuous distance-shift orientation** — that is H2c.
- **Statistical significance of T-subtrees individually** — deferred;
  permutation null is in MRL Q5 and inherited.
- **Disruption** (rest_pre modules failing in rest_post) — symmetric;
  out of scope per canonical doc §2.

## 4. Caveats & failure modes

### 4.1 Probe-bias at coarse h_rel

If a sig cell sits at `k ≤ 5` (giant-dominated), the matched subtree
likely follows probe geometry. Mitigation: per-patient mini-dendrogram
draws the subtree's leaf indices; visual comparison to the implant
geometry (`channel_labels.csv` per patient) lets the reader spot
probe-anchored "trace" subtrees and discount them.

### 4.2 Pat_03 (1024 Hz) and Pat_14 (newly restored)

Pat_03 included with marker; Pat_14 spot-check required (per canonical
doc §7.4). The cohort summary panel reports per-patient T-module
counts so a reader can flag Pat-specific dominance.

### 4.3 J_min sensitivity

`J_min = 0.9` primary. The MRL Q1 sweep `{0.85, 0.9, 0.95}` lives in
the MRL pipeline; this scope report does not re-sweep. If a future
reviewer asks, the analysis script accepts `--j-min` and can be
re-run.

### 4.4 Pooled vs task_test-only candidate set

Default pooled `V_task = V(T_TL) ∪ V(T_TT)` (matches canonical L4 and
MRL Q3). Sensitivity panel: re-run with `V_task = V(T_TT)` only and
report side-by-side counts. Cross-cell flips between the two settings
flag pools-only artefacts.

### 4.5 k vs h_rel mismatch with MRL/CBR

MRL uses native `h_rel`; the headline figure uses integer-k. Stage 4
uses **integer-k cuts** (matches the cell of interest from the
headline) but reports `h_rel(v)` per row so a reader can map a
T-subtree to the MRL native-height view if needed.

## 5. Pseudocode

```
Input:
  K_sig: dict[band → list[(k_lo, k_hi)]]  # sig clusters per band from
                                          # headline cluster-perm
  J_min = 0.9
  k_min = 3
  P = COHORT_N10
  fc_method = "imcoh_abs"

For each (b, (k_lo, k_hi)) ∈ K_sig.items():
  k_rep = round(geometric_mean(k_lo, k_hi))     # representative k for the cell
  for p in P:
    Z_pre  ← load_lrg_result(p, "rest_pre",   b, fc_method).linkage
    Z_TL   ← load_lrg_result(p, "task_learn", b, fc_method).linkage
    Z_TT   ← load_lrg_result(p, "task_test",  b, fc_method).linkage
    Z_post ← load_lrg_result(p, "rest_post",  b, fc_method).linkage

    V_pre  ← [u ∈ tree_internal_nodes(Z_pre)  : size(u) ≥ k_min]
    V_post ← [u ∈ tree_internal_nodes(Z_post) : size(u) ≥ k_min]

    V_task ← []
    for (Z, src) in {(Z_TL, "TL"), (Z_TT, "TT")}:
      V_task.extend([(v, src) for v in tree_internal_nodes(Z)
                              if k_min ≤ size(v) ≤ N_p − 1])

    for (v, src) in V_task:
      S = leaves(v)
      J_pre  = max(J(S, leaves(u)) for u in V_pre)
      J_post = max(J(S, leaves(u)) for u in V_post)
      if J_pre < J_min and J_post ≥ J_min:    # T regime
        emit row: (p, b, k_rep, k_lo, k_hi, src,
                   list(S) sorted, S_named (channel labels),
                   v.h_rel, v.size, J_pre, J_post)

Output:
  data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv

  + per (b, (k_lo, k_hi)) PDF mini-dendrograms for each patient with at
    least one T row in the cell.
```

### 5.1 Complexity

Per `(p, b)` enumerate `~110` task-tree internal nodes; per node
compute `J*` against `~110` rest_pre and `~110` rest_post candidates:
`O(N²)` set ops per node ⇒ `~10⁴` ops per (p, b). Cohort: ~6×10⁵
ops, sub-second per (band, k_bin).

## 6. Visualization

### 6.1 Cohort summary PDF

`data/outputs/figures/section6/trace_modules_summary_n10_imcoh_abs.pdf`

- **Top:** stacked bar per band — count of T-modules pooled across
  patients × bin (pre-/post-Jaccard quadrant). Hatched bars =
  cluster-perm sig cells.
- **Bottom:** sensor-space co-occurrence heatmap (channel × channel)
  of leaf overlap in T modules across all bands. Diagonal hotspots =
  channels that consistently land in T modules together.

### 6.2 Per-patient mini-dendrograms

`data/audit/trace_modules/Pat_NN_band_<band>_k<k>.pdf` (capped at
~40 PDFs: 10 patients × 2 bands × 2 k-bins).

- 3 vertically stacked dendrograms per file: rest_pre, task_test,
  rest_post (sharing leaf set `L_p`).
- For each T module emitted in this (band, k) cell, color the leaves
  consistently across the three phases. Non-T leaves grey.
- Dendrogram y-limits: `tmin = merge_heights[0]*0.8`,
  `tmax = merge_heights[-1]*1.05` (per CLAUDE.md always rule).
- Horizontal bar at the matched internal node's `h_rel(u)` in each
  phase, annotated with the J value to that T module.
- Render the dendrogram as full vector (do not call
  `set_rasterized(True)`).

### 6.3 Reading rules

- **Coloured leaves coherent in pre too** → wrong: T(v) requires
  `J_pre < 0.9`. If you see coherence in pre, the J* matched a
  different rest_pre subtree by coincidence; the row is a false
  positive (debug).
- **Coloured leaves form a tight cluster in task_test and rest_post**
  → expected (T regime).
- **Same colour leaf-set across multiple patients in the cohort
  summary** → cohort-wide module candidate.

## 7. Connection to prior tools

| Prior tool | Relation |
|---|---|
| Headline figure (Stage 2b) | drives cell selection (`K_sig(b)`); this stage grounds each cell in concrete leafsets |
| MRL (`mrl_landscape.py`) | identical primitives (`J_min`, `k_min`, `jaccard_leafsets`); MRL operates at native heights and reports `Π_T(b, bin_h_rel)`; this stage operates at integer-k cuts and reports leafset rows. The two views must agree on existence of T modules, may disagree on counts/coverage |
| CBR variants | sibling family with different anchor / similarity operators; this stage uses the canonical (rest_pre) anchor + Jaccard; CBR variants may surface different subtrees |
| H2d (block-pair persistence) | scalar precursor — `T` at `size = 2` |

## 8. Implementation plan

### 8.1 Analysis driver

`scripts/01_compute/audit/audit_15_trace_modules.py` — new, ~150 lines.

CLI:
```
python scripts/01_compute/audit/audit_15_trace_modules.py \
    [--bands delta,alpha] \
    [--k-bins '23-31,2-4'] \
    [--j-min 0.9] \
    [--k-min 3] \
    [-v]
```

Reads: LRG imcoh_abs caches per (p, b, φ).
Writes: `data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv`.

Library reuse: `tree_internal_nodes`, `jaccard_leafsets`,
`dmax_from_Z` from `lrg_eegfc.utils.metrics.tree`. **No new helpers**
(per coding-rules promotion-on-second-use).

### 8.2 Figure generator

`scripts/01_compute/figures_embedded/fig_trace_modules.py` — new, ~250 lines.

Reads `trace_subtrees_n10_imcoh_abs.csv` + LRG caches (for the
mini-dendrograms only). Writes:
- `data/outputs/figures/section6/trace_modules_summary_n10_imcoh_abs.pdf`
- `data/audit/trace_modules/Pat_NN_band_<band>_k<k>.pdf` (capped ~40)

PDF only, rasterised dendrograms.

### 8.3 Defaults for the v1 run

Driven by the n=9 era cluster-perm hits (still our best estimate
pending re-evaluation at n=10):
- δ k=23–31
- α k=2–4

Stage 4 runs on these by default; if the n=10 cluster-perm output
flags additional sig cells, append them via the `--k-bins` flag in a
follow-up commit.

## 9. Open questions

| Q | Note |
|---|---|
| Q1 — Channel-label resolution | `S_named` requires loading `channel_labels.csv` per patient. Format: comma-separated string of channel names sorted by leaf index. |
| Q2 — Per-patient PDF cap | Default 40 (10 pat × 2 bands × 2 k-bins). If a future cluster-perm result identifies a 3rd band, the cap rises proportionally. |
| Q3 — Cross-MRL validation | Sanity check: any cohort-wide (≥ 8/10) T-cell here should have a corresponding `Π_T ≥ 0.8` cell in MRL at compatible `h_rel`. Mismatch → diagnostic. |
| Q4 — `task_test`-only sensitivity panel | Default off in v1; can be re-run with `--task-only` if needed. |

## 10. Pass criterion (descriptive)

Stage 4 passes if:

1. The CSV contains ≥ 1 T-row in each significant (band, k_bin) cell
   from the headline figure.
2. At least one (band, k_bin) cell has T-rows from ≥ 8/10 patients —
   matches the headline cohort-wide claim with concrete leafsets.
3. The per-patient mini-dendrograms visibly show the T leaves clustered
   in task_test and rest_post but scattered in rest_pre (the visual
   proof of the trace).
4. The cohort summary co-occurrence heatmap surfaces ≥ 1
   off-diagonal hotspot — evidence that specific channel pairs
   recur across patients in T modules.

If 1–4 hold, the band × k headline figure is grounded in identifiable
hierarchical modules, and the user's "we lose the *which*" complaint
is closed.
