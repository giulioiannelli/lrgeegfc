---
name: network_templates / metastable_sankey
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-26
updated: 2026-05-26
pointers:
  - metastable_sankey.py
  - metastable_sankey_mpl.md
  - README.md
  - src/lrg_eegfc/visuals/metastable.py
  - ipynb/06_presentation_figures/05_metastable_sankey.ipynb
  - scripts/07_figures/gen_sankey_grid.py
---

# `metastable_sankey` — interactive Plotly Sankey of cluster evolution

A single-panel **alluvial (Sankey) HTML** showing how the LRG
community partition evolves across a τ-schedule.  One column per τ,
one stacked bar per cluster within a column, one ribbon per
``(c_t → c_{t+1})`` transition with width ∝ number of nodes flowing.

This is the canonical "metastable nodes" visualisation in the LRG
framework — the original Pat_02 / β / rsPre reference is
``data/outputs/figures/presentation_figures/Pat_02/fig5_metastable_sankey_beta.html``.

For a static PDF alternative, see ``metastable_sankey_mpl.md``.

## What this figure shows

- **Columns** — one per τ in the schedule (left = small τ = fine
  structure; right = large τ = coarse structure).
- **Stacked bars** — clusters within a column.  Bar height ∝ number
  of nodes in the cluster.  Colour from the qualitative palette
  (Set3 + Pastel + Dark2) by cluster ID.
- **Ribbons** — width ∝ number of nodes flowing
  ``c_t → c_{t+1}``.  Colour = target-bar colour, α = 0.6.
- **Hover** — each bar shows the contact-label list and the cluster
  size; each ribbon shows the transition count.

The schedule shape encodes the multiscale lens: cluster counts
**must decrease** column-to-column so the Sankey "merges" rather than
"splits".  This is a visualisation contract, not a property of the
data — ``compute_clustering_across_tau`` cuts each τ-specific
dendrogram at a fixed cluster count (``maxclust`` criterion), so the
column shapes are externally imposed.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of the cohort | `Pat_05` |
| `band` | one of `BRAIN_BANDS_NAMES` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | magnitude FC only | `imcoh_abs` |
| `tau_pairs` | space-separated `τ:n_clusters` pairs | `0.1:15 0.5:8 1.0:5 2.0:3 5.0:2` |
| `width` | Plotly figure width (px) | `1400` |
| `height` | Plotly figure height (px) | `700` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/metastable_sankey/` |

## Style choices

- **Magnitude FC only.**  ``imcoh`` (signed) is rejected by LRG — the
  Laplacian needs non-negative weights.  Use ``imcoh_abs``,
  ``imcoh_sq``, ``msc``, or ``corr``.
- **τ-schedule = `DEFAULT_TAU_NCLUST_PAIRS`**: 5 columns at
  ``τ = 0.1, 0.5, 1.0, 2.0, 5.0`` with target ``k = 15, 8, 5, 3, 2``.
  Five columns is the upper limit for readable column labels; expand
  to seven (``gen_sankey_grid.py`` recipe) if you need finer
  resolution and accept smaller fonts.
- **Cluster-count monotonicity**.  The schedule must be monotone
  non-increasing in ``k``.  If you reverse it (fine-on-right) the
  Sankey will still render, but the visual reads "splitting" instead
  of "merging" and the LRG semantics flip.
- **Palette = Plotly Set3 + Pastel + Dark2** (36 distinct hues).
  Cluster IDs are stable across τ within a single Plotly figure but
  **not** across (patient, band, phase) calls — don't try to identify
  "cluster 3" across two separate Sankey HTMLs.
- **Hover truncates at 15 nodes** per cluster; the rest collapse to
  "+N more".  Adjust ``hover_max_nodes=`` in
  ``compute_sankey_flows`` if you need the full list.
- **Output is HTML only** for this template.  Static PDF + raster
  alternatives belong in ``metastable_sankey_mpl``.

## How to call

```bash
# Default — 5-column schedule, Pat_05 / β / rest_pre / imcoh_abs
python .agents/guides/05_plotting/network_templates/metastable_sankey.py

# 7-column schedule (matches scripts/07_figures/gen_sankey_grid.py)
python .agents/guides/05_plotting/network_templates/metastable_sankey.py \
    --tau-pairs 0.1:15 0.3:12 0.5:8 0.8:6 1.0:5 2.0:3 5.0:2

# Different cell
python .agents/guides/05_plotting/network_templates/metastable_sankey.py \
    --patient Pat_02 --band alpha --phase task_test --fc-method msc
```

## Anti-patterns

- ❌ Passing a non-monotone τ-schedule (`0.1:5 0.5:8 1.0:3`) — the
  alluvial will render but the LRG hierarchical-merging reading is
  silently broken.
- ❌ Using signed FC (`imcoh`, raw `corr` with negatives) — the
  Laplacian step will raise `ValueError`.  Always use a magnitude
  variant.
- ❌ Calling ``create_sankey_diagram`` directly in new code — that
  function is the back-compat shim around
  ``compute_sankey_flows + make_metastable_sankey_plotly``; the
  two-step path is the canonical entry point.
- ❌ Saving Plotly Sankey as a PDF via `kaleido` — ``kaleido`` is
  optional and frequently broken in the lapbrain env.  For a vector
  PDF, use ``metastable_sankey_mpl`` (pure matplotlib, fully vector,
  no extra deps).
- ❌ Reading "cluster 3" across two separate HTMLs — cluster IDs are
  τ-and-cell-specific; the palette index is just a colour, not a
  cross-figure identifier.

## Library entry points

```python
from lrg_eegfc.visuals.metastable import (
    compute_clustering_across_tau,    # τ-sweep partitions
    compute_sankey_flows,             # backend-agnostic data prep
    make_metastable_sankey_plotly,    # Plotly renderer
    SankeyFlowData,                   # dataclass
    DEFAULT_TAU_NCLUST_PAIRS,         # 5-column schedule
)
```

See also: ``metastable_sankey_mpl.md`` for the static-PDF twin of
this template.
