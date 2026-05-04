---
name: 05_plotting / fc_templates / README
type: figure_template_index
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-04-29
updated: 2026-04-29
pointers:
  - ../README.md
  - ../colormaps-and-styles.md
  - ../colorbars.md
  - lrgsglib/src/lrgsglib/plotlib/colorbars.py
  - src/lrg_eegfc/visuals/correlation.py
---

# FC adjacency-matrix templates

Canonical, copy-paste-ready scripts for **functional-connectivity
adjacency-matrix figures**. Use these instead of writing FC heatmap
code from scratch — every script produces a PDF that already obeys
the project's plotting rules (full vector, no `suptitle`, no
watermark by default, `imshow_colorbar_caxdivider`).

## When to read this folder

Read here whenever you are about to plot:

- a single FC matrix (one patient × one band × one phase),
- a row of FC matrices (e.g. one patient × one band × all 4 phases),
- a mosaic of FC matrices (per-band × per-phase, per-patient × per-phase, …).

If the request is "plot the adjacency matrix", "show the FC matrix",
"render the connectivity matrix", or "make a heatmap of `imcoh_abs`",
it lands here.

## Templates in this folder

| Template | Layout | Status | Script |
|---|---|---|---|
| `single_adjacency` | 1 axis (one patient × one band × one phase) | iteration in progress | [`single_adjacency.py`](single_adjacency.py) |
| `row_per_phase` | 1 row × 4 phases (one patient × one band × {RPre, TL, TT, RPost}) | iteration in progress | [`row_per_phase.py`](row_per_phase.py) |
| `row_per_band` | 1 row × 6 bands (one patient × one phase × {δ, θ, α, β, γ$_l$, γ$_h$}) | iteration in progress | [`row_per_band.py`](row_per_band.py) |
| `mosaic_band_phase` | 4 phases × 6 bands (one patient, horizontal), per-col scale | round_3_iterating | [`mosaic_band_phase.py`](mosaic_band_phase.py) |
| `mosaic_patient_phase` | n patients × 4 phases (one band), per-row scale | round_3_iterating | [`mosaic_patient_phase.py`](mosaic_patient_phase.py) |
| `mosaic_patient_band` | n patients × 6 bands (one phase), per-col scale | round_3_iterating | [`mosaic_patient_band.py`](mosaic_patient_band.py) |

Each template ships with a sibling `<name>.md` style sheet pinning
its canonical choices (colormap, vmin/vmax convention, tick policy,
file-name pattern, watermark opt-in). Update the `.md` whenever a
template's style is changed, and re-emit the example PDF.

## Where the example PDFs live

Generated PDFs go to
`data/outputs/figures/fc_templates/<template_name>/`. They are
re-generated on demand from the corresponding script in this folder.
Do **not** rely on the cached PDF being current — re-run the script
when the template is touched.

## Style guarantees (true for every template here)

- PDF only, **full vector — never `set_rasterized(True)`**.
  See [`../output-and-rasterization.md`](../output-and-rasterization.md).
- Magnitude FC (`imcoh_abs`, MSC, `corr ≥ 0`) → `cmap="magma"`,
  `vmin=0`, `vmax=row max` (or panel max for a single axis).
- `interpolation="nearest"` on every `imshow` — FC cells are
  discrete; antialiasing blurs probe-block boundaries.
- Diagonal forced to `0` for visualization (FC on the diagonal is
  not informative; leaving raw self-coherence dominates the scale).
- `imshow_colorbar_caxdivider` for every colorbar; never raw
  `fig.colorbar`.
- No `fig.suptitle`. Any context goes into a sidecar `<file>.md`
  caption only when explicitly requested.
- **No watermark / provenance footer by default.** The file name
  carries `<template_name> · <patient/cohort> · <band> ·
  <phase(s)> · <fc_method>`. Each template exposes a
  `watermark=False` kwarg; set it to `True` only when the figure
  will be detached from its file name.
- Axis labels: math indices `$i$` (x), `$j$` (y) — **never** the
  word "contact" / "channel".
- Axis ticks hidden by default. Probe-boundary ticks are an opt-in
  feature, controlled by a kwarg (see `single_adjacency.md`
  §"probe-boundary ticks").

## Routing pointer (for `/figure` and `/plotguide`)

When a figure request mentions "FC", "ImCoh", "MSC",
"correlation matrix", "adjacency matrix", "heatmap of connectivity",
"connectome", or any specific `fc_method` (`corr`, `msc`, `imcoh`,
`imcoh_abs`, `imcoh_sq`):

1. Open this folder's matching template `.md`.
2. Copy the script as a starting point.
3. Modify only the **inputs block** (patient, band, phase,
   `fc_method`, output path) — leave the styling block untouched
   unless we are intentionally iterating on the template itself.
