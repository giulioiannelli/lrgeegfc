---
name: fc_templates / mosaic_band_phase
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_3_iterating
created: 2026-05-01
updated: 2026-05-01
pointers:
  - mosaic_band_phase.py
  - row_per_band.md
  - row_per_phase.md
  - ../colorbars.md
---

# `mosaic_band_phase` — one patient · 4 phases × 6 bands (horizontal)

A 4 × 6 mosaic of FC adjacency matrices for a single patient.  Rows
iterate over phases (rest$_{\mathrm{pre}}$, task$_{\mathrm{learn}}$,
task$_{\mathrm{test}}$, rest$_{\mathrm{post}}$); columns iterate over
bands (δ, θ, α, β, γ$_l$, γ$_h$).  Horizontal layout
(wider than tall) — the natural reading direction is phase progression
along rows, spectral structure across columns.

## What this figure shows

The full band × phase decomposition for one patient at a glance.
Default colorbar mode is **per-col**: each band column has its own
scale, because γ amplitudes are much smaller than α/β under
`imcoh_abs`.  Each column's horizontal cb at the bottom is labelled
with that band's glyph (`<|imcoh_{ij}|>_δ` etc.).  Set
`--colorbar-mode shared` to compare absolute amplitudes across bands
on a single global scale (γ columns will look black).

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_4PHASE` | `Pat_05` |
| `fc_method` | `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic / index / chnames` | **`chnames`** |
| `log_scale` | logarithmic colour norm (`--no-log-scale` for linear) | **`True`** |
| `colorbar_mode` | `shared / per_row / per_col` | **`per_col`** |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/mosaic_band_phase/` |

## Style choices

- **Layout:** `add_gridspec(5, 6)` — 4 phase rows + 1 horizontal-cb
  row at the bottom × 6 band cols.
- **Per-col scale:** `vmin=0`, col-scoped `vmax = max over the 4
  phase panels of that band`.  Keeps γ columns readable.
- **Colorbars:** one per column, horizontal, in the dedicated cb row
  at the bottom.  **Every** cb carries its band-specific glyph
  (`<|imcoh_{ij}|>_δ`, … , `<|imcoh_{ij}|>_γ_h`) — the project rule
  here is "label every cb you draw" so a reader knows what scale
  they're looking at.  Rendered via `plt.colorbar(im, cax=...)` on
  manually placed axes (the project rule against `fig.colorbar` is
  about implicit `ax=` shrinkage; explicit `cax=` is the canonical
  pattern for grids).
- **Tick gating:** y-ticks on leftmost col only, x-ticks on **last
  panel row** (rest_post) only.  Panel-position is passed
  explicitly to `plot_fc_adjacency` (gridspec auto-detection would
  treat the cb row as "last", suppressing rest_post's x-ticks).
- **Row titles:** phase TeX glyphs (rotated, bold) outside the
  leftmost panel.  **Col titles:** band TeX glyphs on the top row.
- **Spacing:** `wspace=0.05`, `hspace=0.08` — panels packed tight,
  cbs visually attached to the bottom row.
- **Interpolation, vector PDF, no suptitle, watermark off** — same
  guarantees as `single_adjacency`.
- **File name:**
  `<Pat_NN>_<fc_method>_<tick_labels>_<scale>_<colorbar_mode>.pdf`.

## How to call

```bash
python .agents/guides/05_plotting/fc_templates/mosaic_band_phase.py \
    --patient Pat_05 --fc-method imcoh_abs
# add --colorbar-mode shared to compare amplitudes across bands
# add --no-log-scale for a linear colour scale
```

## Anti-patterns

- ❌ Default shared linear scale — γ rows go black.
- ❌ Per-cell colorbars (24 of them) — visual chaos.
- ❌ `fig.colorbar(im, ax=axes)` style for the cb column — that
  steals from existing axes.  Use `add_subplot(gs[:, -1])` +
  `plt.colorbar(im, cax=...)`.
