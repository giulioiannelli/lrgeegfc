---
name: fc_templates / mosaic_patient_phase
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_3_iterating
created: 2026-05-01
updated: 2026-05-01
pointers:
  - mosaic_patient_phase.py
  - mosaic_band_phase.md
  - row_per_phase.md
---

# `mosaic_patient_phase` — patients × phases at a fixed band

An *n* × 4 mosaic of FC adjacency matrices.  Rows iterate over a
patient subset, columns iterate over the four phases; band and
`fc_method` are fixed.

## What this figure shows

Cohort-level snapshot of how a single band's connectivity changes
across phases.  Each row is one patient; comparing the four panels
in a row exposes phase-induced reorganisation under the same
spectral bin.  Comparing rows shows cohort-level heterogeneity at
that band.

Default colorbar mode is **per-row**: each patient has its own
scale, since cohort amplitudes differ substantially.  Use
`--colorbar-mode shared` to call out absolute-amplitude differences
between patients.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patients` | comma-separated list | `Pat_05,Pat_06,Pat_08,Pat_13` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `fc_method` | `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic / index / chnames` | **`chnames`** |
| `log_scale` | logarithmic colour norm (`--no-log-scale` for linear) | **`True`** |
| `colorbar_mode` | `shared / per_row / per_col` | **`per_row`** |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/mosaic_patient_phase/` |

## Style choices

- **Layout:** `add_gridspec(n, 5)` — `n` patient rows × 4 panel
  cols + 1 cb col.  Patient rows have *different* probe layouts;
  `chnames` mode reads the per-row probe labels for the leftmost
  panel of each row.
- **Per-row scale:** `vmin=0`, row-scoped `vmax = max over the 4
  phase panels of that patient`.
- **Colorbars:** one per row, dedicated cb col at the right.  Every
  row's cb carries the band-specific label `<|imcoh_{ij}|>_<band>`
  (the row band is fixed, so each cb shows the same glyph — kept on
  every cb so a reader can confirm the scale at any row).
- **Tick gating:** y-ticks (probe letters) on leftmost col only;
  x-ticks on bottom row only — the bottom row reflects the *last*
  patient in the subset.
- **Row titles:** patient ids (`Pat_NN`) rotated bold to the left.
  **Col titles:** phase TeX glyphs on the top row.
- **Interpolation, vector PDF, no suptitle, watermark off** — same
  guarantees as `single_adjacency`.
- **File name:**
  `<n_or_pat>_<band>_<fc_method>_<tick_labels>_<scale>_<colorbar_mode>.pdf`
  where `<n_or_pat>` is `nN` for multi-patient figures, otherwise
  the single patient id.

## How to call

```bash
python .agents/guides/05_plotting/fc_templates/mosaic_patient_phase.py \
    --patients Pat_05,Pat_06,Pat_08,Pat_13 --band beta
# add --colorbar-mode shared to compare patients on a single global scale
# add --no-log-scale for a linear colour scale
```

## Anti-patterns

- ❌ A subset that always re-uses Pat_06 (or Pat_02). Cycle subsets;
  see `feedback_figure_variety.md`.
- ❌ Default shared scale across cohort — masks per-patient
  amplitude differences and crushes weak patients.
- ❌ X-tick labels using one patient's probes for ALL rows —
  meaningless for the others.  We render only the bottom row's
  x-ticks (= last patient's probes); read y-ticks per row.
