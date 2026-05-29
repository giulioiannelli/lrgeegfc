---
name: network_templates / circular_dendrogram_network
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-26
updated: 2026-05-26
pointers:
  - circular_dendrogram_network.py
  - chord_highlighted_edges.md
  - single_network.md
  - lrg_seeded.md
  - README.md
  - src/lrg_eegfc/visuals/network_templates.py
  - scripts/02_preprint/preprint_07_beta_rho_split_figure_test2.py
---

# `circular_dendrogram_network` — chord network + LRG dendrogram overlay

A single-axis network drawing built on graph-tool's
``get_hierarchy_control_points`` (curvy bezier bundles through a
depth-2 cut of the LRG dendrogram), with the **full** circular
dendrogram drawn as a matplotlib overlay on top.

This is the "backbone" pass of the panel-(c) recipe in
``preprint_07_beta_rho_split_figure_test2.py``, lifted to a standalone
template.  Use it when you want the LRG hierarchy itself to be the
visual subject and the FC structure as context underneath.

## What this figure shows

- **Outer ring**: every electrode contact placed at an angle given by
  ``leaves_list(Z)`` — same-cluster contacts are angularly adjacent
  by construction.  Node colour = ``shaft_colors(probes)`` (one
  ``tab20`` hue per sEEG shaft).
- **Curvy edges**: every upper-triangle pair drawn as a hierarchy-
  bundled bezier.  Same-shaft pairs in the shaft's colour, cross-probe
  in mid-gray; width and alpha scale as ``t = (|w|/|w|_max)**edge_gamma``
  (default γ=2.0, zero-floor width/alpha — the bottom ~half of edges
  becomes invisible so only the strong-edge skeleton reads as FC
  context against the dendrogram polyline overlay).
- **Dendrogram overlay**: every merge in ``lrg.linkage_matrix`` drawn
  as a radial+arc polyline pair — radial lines from each child's
  centroid inward to the merge radius, arc joining the children at
  that radius.  Merges deeper in the tree (higher ``h``) sit closer to
  the centre.

Toggle the overlay off (``--no-dendrogram``) to expose the chord
alone; toggle labels off (``--no-labels``) for a clean ring.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient`       | one of `PATIENTS_*` | `Pat_05` |
| `band`          | one of `BRAIN_BANDS` | `beta` |
| `phase`         | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method`     | magnitude FC only | `imcoh_abs` |
| `k_clusters`    | depth-2 LRG cut for chord bundling | `7` |
| `beta`          | bezier bundle tightness (0=chord, 1=fully bundled) | `0.92` |
| `render_px`     | PNG raster size (px square) | `2200` |
| `fit_view`      | graph-tool ``fit_view`` | `0.92` |
| `edge_gamma`    | γ in ``t = (|w|/|w|_max)**γ`` | `2.0` |
| `no_dendrogram` | hide overlay | `False` |
| `no_labels`     | hide radial contact labels | `False` |
| `vertex_size`   | leaf marker size (device px) | `24.0` |
| `out_dir`       | output directory | `data/outputs/figures/network_templates/circular_dendrogram_network/` |

## Style choices

- **Magnitude FC only.**  The edge recipe assumes ``|w|`` (no
  diverging colormap).  Use ``imcoh_abs``, ``imcoh_sq``, or ``msc``.
  Signed FC (``corr``, ``imcoh``) needs a separate sign-aware variant.
- **Default `k_clusters=7`** balances readability against bundle
  thickness.  Lower (4–6) → fewer, thicker bundles; higher (10–15)
  → finer fan-out.  ``preprint_07`` panel (c) uses 7.
- **`beta=0.92`** chosen empirically.  Drop to 0.7 for chord-like
  edges that follow shorter chords through the centre; push to 0.98
  to make bundles almost graph-of-paths.
- **`edge_gamma=2.0` with zero-floor `width_range=(0, 3.5)`,
  `alpha_range=(0, 0.85)`**.  Suppresses the bottom ~half of edges
  into invisibility so the dendrogram dominates while the strongest
  ~30% of FC edges remain as supporting bundles.  Earlier default
  (γ=0.55, non-zero floors) produced a dark spaghetti through the
  chord centre that drowned the dendrogram; γ=3 over-shoots and
  loses the FC context entirely.  Push to 3+ when the dendrogram
  is the only thing you want; drop to 1.0 for a denser FC backdrop.
- **`render_px=2200`** is the right balance between file size
  (~250 kB at this default) and overlay legibility.  Don't go below
  1800 — radial labels and the dendrogram polyline shrink to
  illegible at small fonts.
- **No `fig.suptitle`.**  The ax-level title carries
  ``K=…, β=…, dendro on/off``; the file name carries everything else.

## How to call

```bash
# Default — chord + full circular-dendrogram overlay
python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \
    --patient Pat_05 --band beta --phase rest_pre

# Chord-only baseline (no dendrogram overlay)
python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \
    --patient Pat_05 --band beta --phase rest_pre --no-dendrogram

# Tighter bundling, more clusters
python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \
    --patient Pat_05 --band beta --phase rest_pre --beta 0.97 --k-clusters 12

# Same chord but with the strong-edge ink restricted to the top of
# the distribution (γ=1.5 amplifies the contrast)
python .agents/guides/05_plotting/network_templates/circular_dendrogram_network.py \
    --patient Pat_05 --band beta --phase rest_pre --edge-gamma 1.5
```

## Coordinate-space caveat (read once, never read again)

graph-tool renders to a PNG with origin top-left and y growing DOWN;
matplotlib data coords have origin bottom-left and y growing UP.  The
chord render writes a PNG; the matplotlib overlay (dendrogram + labels)
is drawn in the PNG's pixel coordinates, with the y-axis flipped
implicitly via ``_chord_px_transform``.  All overlay helpers
(``draw_circular_dendrogram_overlay``, ``draw_radial_leaf_labels``)
consume the ``render_meta`` dict returned by the renderer and apply
the transform themselves — you should never compute pixel positions
by hand inside template wrappers.

## Anti-patterns

- ❌ Calling this on signed FC (``corr``, ``imcoh``) — the edge recipe
  collapses negative weights to ``|w|`` and you lose the sign.  Build
  a sign-aware variant if you need it; don't paper over with
  ``|imcoh|`` here.
- ❌ Setting ``--render-px`` below 1800 — overlay text and dendrogram
  polylines become illegible.  If you need a small chord, drop
  ``--fit-view`` (e.g. 0.78) instead.
- ❌ Sweeping ``k_clusters`` and ``beta`` independently across panels
  in a row figure — the visual identity of the chord depends on both;
  pick once and lock them across the row.
- ❌ Drawing the dendrogram overlay with default colour on a dark
  background — the dendrogram is ``(0.12, 0.12, 0.12, 0.70)`` by
  default (near-black, 70% opaque), invisible on dark plots.  Pass
  ``dendrogram_kwargs={"color": (1, 1, 1, 0.8)}`` to the
  ``plot_chord_with_dendrogram`` call site for dark backgrounds.

## Library entry points

```python
from lrg_eegfc.visuals.network_templates import (
    plot_chord_with_dendrogram,           # composite (this template)
    render_hierarchy_chord,               # graph-tool PNG render
    draw_circular_dendrogram_overlay,     # matplotlib polyline overlay
    draw_radial_leaf_labels,              # radial label ring
    build_chord_depth2_layout,            # layout primitive
    HIERARCHY_DEFAULT_K_CLUSTERS,         # = 7
    HIERARCHY_DEFAULT_BETA,               # = 0.92
)
```

See also: ``chord_highlighted_edges.md`` for the highlight-overlay
twin of this template.
