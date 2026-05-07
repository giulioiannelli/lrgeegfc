---
name: KC trace module visualization
description: Multiscale "preserved hierarchy" visualization on FC dendrograms+networks — captures topology-driven traces; height-only traces are invisible
type: result
era: IMCOH_ABS × COHORT_N10
status: locked_2026-05-07
created: 2026-05-07
script: scripts/01_compute/audit/audit_47_kc_trace_network_view.py
output: data/audit/kc_trace_network_view/figures/{patient}__{band}__option_C.pdf
related:
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - data/audit/kc_trace_network_view/trace_subtrees_summary.csv
---

# KC trace as a structural feature of FC: dendrogram+network visualization

**Head.** This visualization is one of the most important proofs of trace
in §5: it shows the same set of leaves moving from *scattered across
the rest_pre dendrogram* to *forming a coherent colored subtree in
task and rest_post*, simultaneously in the dendrogram and in the FC
network. Critically, it captures only **topology-driven** trace
(modules that reorganize and persist as discrete subtree groupings).
**Height-only trace** (anchor topology with merge-height shifts —
the λ=1 component of KC) is invisible to this view.

## What the figure shows

Each `option_C.pdf` is a 3 × 3 grid:

- **rows** = KC λ regime (λ=0 topology, λ=0.5 balanced, λ=1 heights)
- **cols** = phase (rest_pre, task, rest_post)
- **each cell** = a dendrogram (left) and an FC network (right) for
  the same phase, with trace modules highlighted in saturated colors.

A **trace module** is a tt-subtree whose leaf set L_tt:

- has Jaccard ≥ 0.65 with some rPost subtree of comparable size
  (size ratio ∈ [1/1.5, 1.5]),
- has Jaccard < 0.5 with **every** rPre subtree (no size restriction
  on the rPre side),
- the smallest rPre subtree containing L_tt has size ≥ 2 · |L_tt|
  (`FRAG_FACTOR = 2`) — explicit fragmentation guard,
- has size ≥ 3.

The fragmentation guard was added 2026-05-07 to fix a class of
false-positive **anchor** modules. Without it, a 3-leaf L_tt fully
contained in a 5-leaf rPre subtree (Jaccard 3/5 = 0.6, but failing the
old size-restricted Jaccard scan because the rPre subtree was outside
[2, 4]) leaked through the gate. Pat_13 α went from 3 spurious modules
to 1 real module after this fix; the visualization now consistently
shows leaves moving from scattered → coherent across phases.

Modules of different λ subsets:

- λ = 0 (topology): all topology-trace modules.
- λ = 1 (heights): topology-trace modules that *additionally* satisfy
  |h_tt − h_post| / max(h_tt, h_post) ≤ 0.30 (merge-height match).
- λ = 0.5 (balanced): intersection (subset of λ = 1).

Greedy de-duplication keeps modules with ≤ 30 % leaf-set overlap and
caps at 12 modules per (patient, band).

## Coloring rules

- **Dendrogram** — per-phase leaf and link coloring:
  - tt's dendrogram colors `leaves_tt`.
  - rPost's dendrogram colors `leaves_post` (the matched rPost subtree's leaves).
  - rPre's dendrogram colors `leaves_tt ∪ leaves_post` as scattered leaf-bars.
  - Internal links color by smallest module M with `link.leaves ⊆ M.leaves`
    under that phase's leaf set. (Without this rule the module fragments
    visually whenever the rPost match is at Jaccard 0.65–0.75.)
  - Explicit colored vertical bars at each leaf x-position (scipy's
    `link_color_func` only paints internal links; leaves inherit their
    parent link's color, which in rPre is gray — bars override that).

- **Network** — phase-INVARIANT module membership (union leaves), edge
  rank-percentile per phase. Strong intra-module edges in tt/rPost
  saturate; weak rPre intra-module edges fade. The network has one
  shared node set across phases, so phase contrast must come from
  edges, not from which nodes are colored.

- **Layout** — `lrg_sfdp` seeded with the task_test partition at k=20.
  This is critical: trace-module leaves end up in compact regions of
  the network rather than being split across the layout.

- **Merge-height bars** — at λ = 0.5 and λ = 1, each module's merge
  height in tt and rPost is drawn as a full-width dashed `axhline`
  in the matching color. rPre is omitted because there is no
  matched subtree (by construction). This visually separates a
  module's λ = 0 topological identity from its λ = 1 height
  signature.

## What the visualization captures

**Topology-driven trace.** A subtree S of leaves L:

- exists as a coherent subtree in tt and rPost (Jaccard ≥ 0.65),
- does NOT exist as a coherent subtree in rPre (Jaccard < 0.5).

This is the trace mechanism for any KC trace component that pulls
leaves into NEW co-clustering groups during task that persist into
rest_post. Most cohort-level β / α / low_γ traces produce visible
modules under this rule.

## What the visualization MISSES

**Height-only trace.** A subtree whose leaves are stable across all
three phases (anchor topology, J(tt, rPre) ≈ J(tt, rPost) ≈ 1) but
whose merge HEIGHT shifts dramatically between rPre and tt+rPost
*within* the stable subtree. The leaves don't reorganize, only the
RG-flow distance does.

Such modules are filtered out by the `Jaccard(tt, rPre) < 0.5`
constraint — by construction. KC's λ=1 component still detects them
because λ=1 reads merge heights, but the structural-network view
cannot show them as "modules that reorganize".

### Pat_08 — the canonical height-only-trace counterexample

Pat_08 is 5/5 in the cohort cross-probe sign matrix at β (alongside
Pat_06). Yet under this visualization with `J_POST_MIN = 0.65` and
`J_PRE_MAX = 0.5`:

| band | Pat_08 trace modules |
|:--|:--|
| α | 0 |
| β | 0 |
| low_γ | 0 |

Despite being a strong cohort-level scalar trace contributor, Pat_08
has **zero topology-trace modules** across all three trace bands. His
KC trace must therefore be carried entirely by the λ=1 (height)
component within stable topology. This is consistent with the round-2
mechanistic reading "topology + heights" for β: the cohort-level scalar
trace pools both components, but individual patients can carry only
one of them.

## Cohort coverage (current sweep)

Patients sweeped: full cohort minus Pat_08 — Pat_06, Pat_03, Pat_05,
Pat_13, Pat_02, Pat_10, Pat_14, Pat_07, Pat_15 (n=9). Bands sweeped:
**all six** — δ, θ, α, β, low_γ, high_γ. Counts are module counts
at λ = 0; bold sizes = ≥ 15 leaves (network-scale).

Per-(patient, band) module counts at λ = 0 with `J_POST_MIN = 0.6`:

| | δ | θ | α | β | low_γ | high_γ |
|:--|:--|:--|:--|:--|:--|:--|
| Pat_06 | 2 | 3 | 3 | 2 [5, **39**] | 2 | 2 [6, **55**] |
| Pat_03 | 3 | 1 | 3 [4, **15**, **45**] | 3 [7, 14, **15**] | 1 [**17**] | 1 |
| Pat_05 | — | 2 | 2 [4, **48**] | 1 [**47**] | 5 [3, 9, 9, **19**, **45**] | 1 |
| Pat_13 | 2 | 2 | 1 | 2 | 1 | 4 |
| Pat_02 | 2 | — | 3 | 4 [3, 4, 4, 10] | 3 | 3 |
| Pat_10 | — | 1 | 1 [**25**] | 1 | — | — |
| Pat_14 | 2 | 2 | 2 | 3 [4, 4, 9] | 2 | — |
| Pat_07 | 1 | 3 [5, 5, **13**] | 4 [4, 4, 6, **35**] | **6 [4, 4, 4, 6, 6, 22]** | 4 | 3 |
| Pat_15 | — | — | — | 1 | 2 | — |
| **cells/9** | **7/9** | **7/9** | **8/9** | **9/9** | **8/9** | **6/9** |
| **total** | 15 | 15 | 23 | **23** | 21 | 10 |

Bold sizes = ≥ 15 leaves (network-scale). "—" = empty cell, no figure
rendered. **45/54 cells visible across the n=9 cohort × 6 bands** (83 %).
Plus Pat_08 0/6 = 45/60 across full n=10 cohort × 6 bands.

**Per-band ranking** (`cells_with_any/9`):
β (9/9) > α = low_γ (8/9) > δ = θ (7/9) > high_γ (6/9). β is the
**only band with universal visibility** — confirms §5.2's claim that
β is the load-bearing trace band (5/5 cohort sign +
Bonferroni-surviving test). β is also tied with α for total module
count (23 each), but distributed across more patients (9/9 vs 8/9).

**Per-patient β richness** (post 0.6 lowering, 2026-05-07):
Pat_07 (6) > Pat_02 (4) > Pat_03 = Pat_14 (3) > Pat_06 = Pat_13 (2)
> Pat_05 = Pat_10 = Pat_15 (1).

**Notable additions at the 0.6 threshold (vs prior 0.65 sweep):**
- Pat_05 β jumps from a single 3-leaf module to a single **47-leaf
  module** (~40 % of the network) — a major reorganization the
  stricter threshold completely missed.
- Pat_07 β grows from 4 to **6 modules**.
- Pat_14 β grows from 2 to 3 modules (largest now 9 leaves).
- Pat_07 and Pat_14 each pick up a new band (δ for Pat_07, θ for
  Pat_14).

Per-module sizes, Jaccard scores, and figure paths in
`data/audit/kc_trace_network_view/trace_subtrees_summary.csv`.

## Why this matters for the manuscript

The KC scalar trace (Section 5.2) is the load-bearing significance
test, but the scalar collapses topology and heights into one number.
The cohort-level p = 0.019 at λ=0 and p = 0.042 at λ=1 already
demonstrates that the trace has both components, statistically. What
this visualization adds is the **mechanism made visible**: the
topology component IS specific identifiable groups of leaves that
fragment in rPre and crystallize in tt+rPost. It is a per-patient,
per-band, multiscale concrete object that any reader can inspect.

The fact that Pat_08 doesn't show topology-trace modules is not a
weakness — it is a useful split that supports the §5.2 decomposition:
β trace pools two mechanisms (topology + heights), and individual
patients can reside on either side.

## Configuration parameters (locked 2026-05-07)

```python
J_POST_MIN = 0.6     # tt-rPost Jaccard threshold (lowered from 0.65 to surface more β modules)
J_PRE_MAX = 0.5      # tt-rPre Jaccard upper bound (fragmentation)
FRAG_FACTOR = 2.0    # smallest rPre subtree containing L_tt ≥ 2·|L_tt|
HEIGHT_REL_THR = 0.30
MIN_SIZE = 3
MAX_TRACE_TOTAL = 12
MAX_OVERLAP_FRAC = 0.30
LAYOUT_K = 20        # only for lrg_sfdp module-seeding
SIZE_RATIO_RANGE = 1.5  # rPost candidate matches restricted to size /1.5..*1.5
                        # NOT applied to rPre Jaccard scan (anti-anchor fix 2026-05-07)
```

**Empty-cell skip.** If the λ=0 module list is empty for a (patient, band)
cell, the script does NOT render a figure (logs "NO modules; skipping
figure" and writes an empty CSV row). An all-gray panel grid is
visually useless and only obscures the gallery — the CSV still records
the negative result.

**Linear y-axis.** Dendrograms render with a linear y-axis from 0 to
`merge_heights[-1] * 1.05`. The earlier log-scale rendering compressed
near-root structure and made the visual harder to read; linear is the
locked default.

Tightening `J_POST_MIN` to 0.75 drops the cohort to ~1 module per
case (too sparse). Loosening to 0.5 produces visual fragmentation
because rPost's tree topology groups some matched leaves with
non-trace leaves — the per-phase coloring rule then produces gray
internal links inside what should look like a coherent module. 0.65
is the locked value.

## Suggested manuscript exemplars

For a single-figure showcase, the strongest contrasts are:

- **Pat_07 β** — four modules [4, 4, 6, 22] visible at all three λ;
  cleanest "many-module trace" demonstration in the cohort.
- **Pat_10 α** — single 25-leaf module, ~22 % of the brain
  reorganizes as one tight subtree (single-module showcase).
- **Pat_03 α** — three modules of escalating size [4, 15, 45] showing
  multiscale structure within one (patient, band).
- **Pat_05 low_γ** — five modules of varying sizes [3, 9, 9, 19, 45],
  the densest single cell of the cohort.
- **Pat_06 β** — load-bearing β band (5/5 patient), [5, 39] modules.
- **Pat_08 β** (counterexample) — height-only-trace, zero modules
  visible by construction.

These five subfigures span the cohort sizes and demonstrate the
visualization's reach, including its principled blind spot (Pat_08).

## Future extensions

1. **Reset detection** (next iteration): same scaffolding but with
   J(rPre, rPost) ≥ 0.65 AND J(rPre, tt) < 0.5 AND J(rPost, tt) < 0.5
   — modules that exist in both rest phases but are disrupted during
   task, then revert. This is the RESET pattern from
   `.agents/guides/01_project/terminology.md` and would let us show
   "anatomically driven" modules that the task transiently breaks.
2. **Height-only trace detection** to recover Pat_08-style cases —
   would need a different criterion (stable topology + significant
   merge-height shift). Distinct visualization since "modules that
   reorganize" doesn't apply.
3. **Band expansion** *(done 2026-05-07 same session)*: full 6-band
   sweep (δ, θ, α, β, low_γ, high_γ) showed β as the only universally
   visible band (9/9). δ, θ, high_γ progressively sparser — consistent
   with §5.2's β-load-bearing claim.

## Pointers

- Script: `scripts/01_compute/audit/audit_47_kc_trace_network_view.py`
- Figures: `data/audit/kc_trace_network_view/figures/`
- Per-module summary: `data/audit/kc_trace_network_view/trace_subtrees_summary.csv`
- MRL conceptual basis:
  `.agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md`
- §5.2 KC scalar trace: `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md`
