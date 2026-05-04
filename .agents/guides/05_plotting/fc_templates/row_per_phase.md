---
name: fc_templates / row_per_phase
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_2_iterating
created: 2026-04-29
updated: 2026-04-29
pointers:
  - row_per_phase.py
  - single_adjacency.md
  - ../colorbars.md
---

# `row_per_phase` — one patient × one band × four phases

A 1 × 4 row of FC adjacency matrices.  Same patient, same band, same
`fc_method`; one panel per phase (rest$_{\mathrm{pre}}$,
task$_{\mathrm{learn}}$, task$_{\mathrm{test}}$,
rest$_{\mathrm{post}}$).

## What this figure shows

Phases share a band, so the matrices are directly comparable at the
same colour scale.  The default is therefore **shared scale + one
colorbar at the right**.  Differences between panels read as
phase-induced reorganisation.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_4PHASE` | `Pat_06` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `fc_method` | `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic / index / chnames` | **`chnames`** |
| `log_scale` | logarithmic colour norm (`--no-log-scale` for linear) | **`True`** |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/row_per_phase/` |

> The default-on choices (`chnames` + log) are the canonical
> publication style for this template — visible probe geography
> + visible γ structure that linear-scale crushes.

## Style choices

- **Layout:** `subplots(1, 4, sharey=True,
  gridspec_kw={"width_ratios": [1, 1, 1, 1.06]})` — y-tick labels
  appear only on the leftmost panel.  The right-most column is pre-
  padded by `_ROW_CB_FRAC` (≈ 6%) so that the `imshow_colorbar_caxdivider`
  cax + pad eat from the *extra* width — the imshow stays the same
  size as the others.
- **Colour scale:** shared across the row;
  `vmin=0`, `vmax = max over all four matrices`.
- **Colorbar:** single, attached to the right-most panel via
  `imshow_colorbar_caxdivider` (the project's canonical helper —
  never `fig.colorbar` on multi-axis grids).  Label =
  `<|imcoh_{ij}|>_{<band>}` (or method analogue).
- **Titles:** TeX phase glyphs only; no patient / band / method
  noise on the panel.  File name carries the rest.
- **Interpolation, vector PDF, no suptitle, watermark off** — same
  guarantees as `single_adjacency`.
- **File name:**
  `<Pat_NN>_<band>_<fc_method>_<tick_labels>_<scale>.pdf`.

## How to call

```python
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row
fig, axes = plot_fc_adjacency_row(
    matrices=[M_RPre, M_TL, M_TT, M_RPost],
    titles=[r"rest$_{\mathrm{pre}}$", r"task$_{\mathrm{learn}}$",
            r"task$_{\mathrm{test}}$", r"rest$_{\mathrm{post}}$"],
    fc_method="imcoh_abs", band="beta",
)
```

```bash
python .agents/guides/05_plotting/fc_templates/row_per_phase.py \
    --patient Pat_06 --band beta --fc-method imcoh_abs \
    --tick-labels generic
# add --log-scale to compress dynamic range
# add --watermark for a slide-deck footer
```

## Anti-patterns

- ❌ Per-panel colorbars when phases share a band — wastes width
  and breaks the visual symmetry.
- ❌ `set_rasterized(True)` (vector only).
- ❌ Repeating `<patient>` / `<fc_method>` in the panel titles.
