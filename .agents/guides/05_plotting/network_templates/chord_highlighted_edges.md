---
name: network_templates / chord_highlighted_edges
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-26
updated: 2026-05-26
pointers:
  - chord_highlighted_edges.py
  - circular_dendrogram_network.md
  - single_network.md
  - lrg_seeded.md
  - README.md
  - src/lrg_eegfc/visuals/network_templates.py
  - scripts/02_preprint/preprint_07_beta_rho_split_figure_test2.py
---

# `chord_highlighted_edges` — chord with selective edge highlights

A single-axis chord network where the **subject** is a small named
subset of edges drawn in vivid colour, and every other pair acts as a
faint grey "null backbone" underneath.  Generalisation of the
trace/anti overlay in panel (c) of
``preprint_07_beta_rho_split_figure_test2.py`` to any caller-supplied
edge subset(s).

## What this figure shows

- **Outer ring**: same as ``circular_dendrogram_network`` — leaves at
  ``leaves_list(Z)`` angles, shaft-coloured nodes.
- **Background**: every non-highlighted pair drawn as a faint grey
  bezier (``HIERARCHY_BG_RGB`` = `(0.62, 0.62, 0.62)` at
  ``HIERARCHY_BG_ALPHA`` = `0.085`, ``HIERARCHY_BG_PEN`` = `0.42`).
  Optional cross-probe-only filter (``cross_probe_only=True``) drops
  same-shaft pairs from the background; they're 2–8× stronger and
  otherwise dominate the grey wash visually.
- **Highlights**: one or more user-named colour classes, each drawn on
  top of the background.  Each class is a dict:

  ```python
  {
      "pairs":      np.ndarray,         # 1D flat-triu indices OR (M, 2) (i,j) array
      "rgb":        (r, g, b),          # highlight colour
      "magnitudes": np.ndarray | None,  # optional per-pair |w|; scales width/alpha
      "width_range":  (min, max),       # default (1.4, 6.4)
      "alpha_range":  (min, max),       # default (0.45, 0.95)
      "gamma":        float,            # γ in t = (m/m_max)**γ; default 0.55
      "draw_order_offset": float,       # higher = on top; default 1000*k_class
      "label":      str,                # appears in the title
  }
  ```

The default demo wires up two classes: top-K cross-probe |w| in
**forest green** + top-K same-probe |w| in **brick red** (palette
anchored to PALETTE_COPH from ``preprint_11_beta_anatomy_brain`` so
paired figures share a hue family).

For the trace/anti task-persistence use case, build your own
highlights list — typically signed ``Δ_task · Δ_rest`` magnitudes, one
class for trace (``sigma=+1``), one for anti (``sigma=-1``) — and call
``plot_chord_with_highlight`` directly from your script.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` / `band` / `phase` / `fc_method` | one panel | `Pat_05 / beta / rest_pre / imcoh_abs` |
| `highlight_mode` | `topk` is wired up here; for custom sets call the library directly | `topk` |
| `k_trace`        | primary-class size (cross-probe in demo)  | `120` |
| `k_anti`         | secondary-class size (same-probe in demo) | `40` |
| `cross_probe_bg` | cross-probe-only background grey wash     | `False` |
| `k_clusters`     | depth-2 LRG cut for chord bundling        | `7` |
| `beta`           | bezier bundle tightness                   | `0.92` |
| `render_px`      | PNG raster size                           | `2200` |
| `fit_view`       | graph-tool ``fit_view``                   | `0.92` |
| `show_dendrogram`| overlay the LRG dendrogram too (hybrid)   | `False` |
| `no_labels`      | hide radial contact labels                | `False` |
| `vertex_size`    | leaf marker size (device px)              | `24.0` |
| `out_dir`        | output directory                          | `data/outputs/figures/network_templates/chord_highlighted_edges/` |

## Style choices

- **The background should LOSE visually.**  Default ``alpha=0.085``
  and ``pen=0.42`` are chosen so the grey backbone reads as "there is
  a graph here" without competing with the highlights.  If the
  background is the visual subject, you want
  ``circular_dendrogram_network``, not this template.
- **Per-class draw order** = ``draw_order_offset + magnitude``.  By
  default class 1 sits in the 1000s, class 2 in the 2000s — last
  class wins on overlap.  Override per-class if you want the green
  on top of the red instead of the reverse.
- **Per-class γ scales WITHIN the class.**  γ=0.55 is a mild
  compression that keeps the weak edges of a class still visible.
  Push to γ=1.0–1.5 to bias hard toward the strongest members of the
  class only.
- **Two-class palette** anchored to ``preprint_11`` / ``preprint_07``:
  - primary: ``#1a7c3e`` (forest green)
  - secondary: ``#c0392b`` (brick red)
  - background: ``(0.62, 0.62, 0.62)`` at α=0.085.
- **`render_px=2200`** — same as ``circular_dendrogram_network``.
- **`show_dendrogram=False` default.**  This template's focus IS the
  highlight; the dendrogram overlay competes with the highlight ink.
  Turn on (``--show-dendrogram``) only when you actually want the
  hybrid "highlight on top of dendrogram on top of chord" look.
- **No `fig.suptitle`.**  The ax-level title carries the class labels
  ("top-120 cross-probe |w| + top-40 same-probe |w|") so the reader
  knows what's coloured.

## How to call

```bash
# Default demo — top-120 cross-probe (green) + top-40 same-probe (red)
python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \
    --patient Pat_05 --band beta --phase rest_pre

# Single-class highlight: top-200 |w| only (no red/green split)
python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \
    --patient Pat_05 --band beta --phase rest_pre --k-trace 200 --k-anti 0

# Cross-probe-only background so the same-shaft grey clump
# doesn't dominate the wash
python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \
    --patient Pat_05 --band beta --phase rest_pre --cross-probe-bg

# Hybrid with the circular-dendrogram overlay
python .agents/guides/05_plotting/network_templates/chord_highlighted_edges.py \
    --patient Pat_05 --band beta --phase rest_pre --show-dendrogram
```

## Calling from your own script (the trace/anti use case)

```python
from lrg_eegfc.visuals.network_templates import (
    plot_chord_with_highlight, _load_inputs,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

A, probes = _load_inputs("Pat_05", "beta", "rest_pre", "imcoh_abs")
lrg = load_lrg_result("Pat_05", "rest_pre", "beta", "imcoh_abs")

# your per-pair signed magnitudes (sigma=+1 trace, sigma=-1 anti)
# trace_pairs = (M_t, 2) array of (i, j); anti_pairs = (M_a, 2)
fig, ax = plt.subplots(figsize=(7.2, 7.2))
plot_chord_with_highlight(
    ax, A, probes, lrg,
    highlights=[
        {"pairs": trace_pairs, "rgb": (0.10, 0.49, 0.24),
         "magnitudes": trace_mag, "label": f"trace ({len(trace_pairs)})",
         "draw_order_offset": 2000.0},
        {"pairs": anti_pairs,  "rgb": (0.75, 0.23, 0.17),
         "magnitudes": anti_mag,  "label": f"anti ({len(anti_pairs)})",
         "draw_order_offset": 1000.0},
    ],
    cross_probe_only=True,   # background = cross-probe only
)
```

## Coordinate-space caveat

Same as ``circular_dendrogram_network`` — graph-tool's PNG is
y-down, matplotlib's data axis is y-up.  The pixel transform is
applied internally; never hand-compute pixel positions in a wrapper.

## Anti-patterns

- ❌ Passing more than 4 highlight classes — the eye stops tracking
  distinct categorical colours past three or four.  If you need more,
  rethink whether the figure should be a small-multiples grid instead.
- ❌ Making the background loud (``HIERARCHY_BG_ALPHA > 0.2``) and
  then adding highlights — both layers fight for attention and the
  reader can't tell what the figure is about.  Keep the background as
  context, not content.
- ❌ Drawing thousands of highlights — past ~300 the highlight class
  becomes its own hairball.  If your selection has that many pairs,
  it's not really a "highlight"; rebuild the figure as a
  ``circular_dendrogram_network`` chord with a magnitude-driven recipe.
- ❌ Loading ``corr`` or ``imcoh`` (signed FC) into a single-class
  highlight without distinguishing positive vs negative — split into
  two highlight classes (``rgb`` per sign) before passing.
- ❌ Forgetting ``draw_order_offset`` when stacking three or more
  classes — last-class-wins relies on the offset; without it the
  classes ravel by per-pair magnitude only.

## Library entry points

```python
from lrg_eegfc.visuals.network_templates import (
    plot_chord_with_highlight,            # composite (this template)
    render_hierarchy_chord,               # graph-tool PNG render
    draw_circular_dendrogram_overlay,     # optional overlay
    draw_radial_leaf_labels,              # radial label ring
    HIERARCHY_BG_RGB,                     # = (0.62, 0.62, 0.62)
    HIERARCHY_BG_ALPHA,                   # = 0.085
    HIERARCHY_BG_PEN,                     # = 0.42
)
```

See also: ``circular_dendrogram_network.md`` for the dendrogram-overlay
twin of this template.
