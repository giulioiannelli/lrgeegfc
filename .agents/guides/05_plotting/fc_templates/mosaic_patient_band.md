---
name: fc_templates / mosaic_patient_band
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_3_iterating
created: 2026-05-01
updated: 2026-05-01
pointers:
  - mosaic_patient_band.py
  - mosaic_band_phase.md
  - row_per_band.md
---

# `mosaic_patient_band` — patients × bands at a fixed phase

An *n* × 6 mosaic of FC adjacency matrices.  Rows iterate over a
patient subset, columns iterate over the six bands; phase and
`fc_method` are fixed.

## What this figure shows

Cohort × spectral cube at one phase.  Comparing rows at a fixed
column shows how a single band varies across patients; comparing
columns within a row shows how connectivity reorganises across
spectral bins for one patient.

Default colorbar mode is **per-col**: each band has its own scale,
because γ amplitudes are typically much smaller than α/β under
`imcoh_abs`.  Cb axes sit at the **bottom** of the gridspec, one
horizontal cb per column.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patients` | comma-separated list | `Pat_05,Pat_06,Pat_08,Pat_13` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic / index / chnames` | **`chnames`** |
| `log_scale` | logarithmic colour norm (`--no-log-scale` for linear) | **`True`** |
| `colorbar_mode` | `shared / per_row / per_col` | **`per_col`** |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/mosaic_patient_band/` |

## Style choices

- **Layout:** `add_gridspec(n+1, 6)` — `n` patient rows × 6 band
  cols + 1 horizontal cb row at the bottom.
- **Per-col scale:** `vmin=0`, col-scoped `vmax = max over the n
  patient panels of that band`.  Keeps γ columns readable.
- **Colorbars:** one per column, horizontal, in the dedicated cb
  row at the bottom.  **Every** cb carries its band-specific glyph
  (`<|imcoh_{ij}|>_δ`, … , `<|imcoh_{ij}|>_γ_h`).
- **Tick gating:** y-ticks on leftmost col only; x-ticks on the
  **last patient row** (the bottom row of *panels*, not of the
  gridspec — the gridspec's last row hosts the cbs).  Panel-position
  is passed explicitly to `plot_fc_adjacency` to keep the bottom
  patient's x-ticks visible.
- **Row titles:** patient ids rotated bold left of the leftmost
  panel.  **Col titles:** band TeX glyphs on the top row.
- **Interpolation, vector PDF, no suptitle, watermark off** — same
  guarantees as `single_adjacency`.
- **File name:**
  `<n_or_pat>_<phase>_<fc_method>_<tick_labels>_<scale>_<colorbar_mode>.pdf`.

## How to call

```bash
python .agents/guides/05_plotting/fc_templates/mosaic_patient_band.py \
    --patients Pat_05,Pat_06,Pat_08,Pat_13 --phase rest_pre
# add --colorbar-mode per_row to give each patient its own scale instead
# add --no-log-scale for a linear colour scale
```

## Anti-patterns

- ❌ Default shared linear scale — γ columns go black.
- ❌ Cb at the right with `per_col` mode — n cbs can't all sit on
  the right column-aligned with their data; horizontal bottom cbs
  are the only correct placement.
- ❌ X-tick labels duplicated on every row — only the bottom row
  carries them.
