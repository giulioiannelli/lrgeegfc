---
name: network_templates / metastable_sankey_mpl
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-26
updated: 2026-05-26
pointers:
  - metastable_sankey_mpl.py
  - metastable_sankey.md
  - README.md
  - src/lrg_eegfc/visuals/metastable.py
---

# `metastable_sankey_mpl` — static matplotlib alluvial PDF

Vector-PDF twin of ``metastable_sankey`` (Plotly).  Same data prep
pipeline (``compute_clustering_across_tau`` + ``compute_sankey_flows``);
different renderer.  Cubic-bezier ribbon polygons + stacked cluster
bars, drawn entirely with ``matplotlib.patches`` — no Plotly, no
``kaleido``, no rasterisation.

Use this when you need a Sankey **for a paper figure**.  Use the
Plotly variant when you need interactive hover lists (cluster member
contacts, exact transition counts).

## What this figure shows

- **Columns** — one per τ at ``x = i · x_step``.  Left = small τ
  (fine), right = large τ (coarse).
- **Stacked bars** — clusters within each column.  Height ∝ cluster
  size in node-count units; fill = qualitative palette by cluster ID
  (matches the Plotly variant exactly so the two renders are
  visually equivalent).
- **Ribbons** — cubic-bezier polygons between consecutive columns,
  width = transition count, colour = target-bar colour at
  ``ribbon_alpha=0.55``.
- **Column headers** — ``"τ = X (k clusters)"`` text above the top
  bar, at fraction ``column_label_pad`` of total height above ``y=0``.
- **Bar labels** — optional ``"C{cluster} ({n})"`` text to the right
  of each bar.  Toggle off via ``--no-bar-labels`` for clean prints
  or when bars get too small.

## Why polygon ribbons (not nx / Plotly export)

- ``matplotlib.sankey.Sankey`` is for **energy-flow Sankey**
  (single source, branching outputs).  It does **not** support
  multi-column alluvial.  Wrong tool for this job.
- ``plotly.write_image`` requires ``kaleido``, which is unreliable
  in the lapbrain env and produces a rasterised PDF anyway.
- ``pysankey`` / ``mpl-sankey`` are third-party packages with their
  own quirks and are not in the lapbrain env.

The cubic-bezier ribbon is ~30 lines of matplotlib; we ship it in
``lrg_eegfc.visuals.metastable.make_metastable_sankey_mpl``.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of the cohort | `Pat_05` |
| `band` | one of `BRAIN_BANDS_NAMES` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | magnitude FC only | `imcoh_abs` |
| `tau_pairs` | space-separated `τ:n_clusters` pairs | `0.1:15 0.5:8 1.0:5 2.0:3 5.0:2` |
| `figsize` | inches `(W H)` | `(12, 6)` |
| `node_width` | cluster-bar width (x-axis units) | `0.06` |
| `ribbon_alpha` | ribbon opacity | `0.55` |
| `no_bar_labels` | suppress `C{c} ({n})` text | `False` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/metastable_sankey_mpl/` |

## Style choices

- **Rule compliance**: ``use_lrg_style()`` at the top, no
  ``fig.suptitle``, no ``set_rasterized(True)``, no PNG sibling.
  Axis-level title carries ``patient · band · phase · fc_method``;
  file name carries the τ-schedule.
- **Bar height in *node count* units, not normalised**: the absolute
  height of the rightmost column equals the absolute height of the
  leftmost.  This is the property an alluvial relies on for the
  "merging" reading; do not normalise to [0, 1] across columns.
- **Gap between bars within a column** is set by ``node_gap_frac``
  (default 1.2 % of total height — roughly 1 node-row at N=80).
  Decrease for chunkier bars, increase for airier columns.
- **Ribbon ordering**: within each source bar, ribbons are stacked
  top-to-bottom in **target-bar order** at the same column.  This is
  a non-crossing heuristic that works for monotone τ-schedules
  (cluster counts decreasing).  It is NOT a global anti-crossing
  optimiser; very tangled schedules may still cross visually.
- **Cubic-bezier control points** sit at the column midpoint
  vertically aligned with the source/target ends.  This is the
  Plotly default; gives a symmetric S-curve.
- **Default ``ribbon_alpha=0.55``** (vs 0.6 Plotly default) — the
  static print needs slightly more legibility margin since you
  can't disambiguate via interaction.

## How to call

```bash
# Default — 5-column schedule, Pat_05 / β / rest_pre / imcoh_abs
python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py

# Wider figure
python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py \
    --figsize 16 7

# 7-column schedule (matches gen_sankey_grid.py)
python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py \
    --tau-pairs 0.1:15 0.3:12 0.5:8 0.8:6 1.0:5 2.0:3 5.0:2

# Clean print — no bar labels, chunkier bars
python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py \
    --no-bar-labels --node-width 0.10
```

## Anti-patterns

- ❌ Calling ``matplotlib.sankey.Sankey`` instead of
  ``make_metastable_sankey_mpl`` — wrong tool (energy-flow Sankey,
  not alluvial).
- ❌ Exporting Plotly Sankey to PDF via ``kaleido`` — rasterises the
  ribbons (PDF is just a wrapped PNG); use this template instead.
- ❌ Normalising bar heights across columns — kills the "merging"
  reading.  Keep node-count units throughout.
- ❌ Passing a non-monotone ``tau_pairs`` schedule — the alluvial
  will render but ribbons will cross heavily and the LRG hierarchical-
  merging interpretation breaks.
- ❌ Setting ``rasterize=True`` on the ribbon polygons — would
  defeat the whole point of choosing matplotlib over Plotly's PDF
  export.  Vector PDF is fine: a 5-column figure produces 5–20 KB.

## Library entry points

```python
from lrg_eegfc.visuals.metastable import (
    compute_clustering_across_tau,
    compute_sankey_flows,
    make_metastable_sankey_mpl,    # static-PDF renderer
    SankeyFlowData,
    DEFAULT_TAU_NCLUST_PAIRS,
)
```

See also: ``metastable_sankey.md`` for the interactive HTML twin.
