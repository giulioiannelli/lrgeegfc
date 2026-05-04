---
name: fc_templates / row_per_band
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_2_iterating
created: 2026-04-29
updated: 2026-04-29
pointers:
  - row_per_band.py
  - single_adjacency.md
  - ../colorbars.md
---

# `row_per_band` — one patient × one phase × six bands

A 1 × 6 row of FC adjacency matrices.  Same patient, same phase,
same `fc_method`; one panel per band (δ, θ, α, β, γ$_l$, γ$_h$).

## What this figure shows

Bands have **different natural amplitudes** (γ tends to be much
weaker than α/β under `imcoh_abs`).  A shared linear scale would
crush γ panels to black, so the default is **per-panel scale** with
**per-panel colorbars**.  Differences between panels read as
band-specific functional structure.

For a synoptic, scale-aware comparison, use:

- `--log-scale` — keeps per-panel scales but flattens the dynamic
  range; reveals γ structure without visually losing α/β contrast.
- `--shared-scale` — single global `vmax`, single right-side
  colorbar.  Useful when you want to call out the absolute amplitude
  drop in γ.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_4PHASE` | `Pat_06` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic / index / chnames` | **`chnames`** |
| `log_scale` | logarithmic colour norm (`--no-log-scale` for linear) | **`True`** |
| `shared_scale` | one global `vmax` + single right-side colorbar | `False` |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/row_per_band/` |

> The default-on choices (`chnames` + log) are the canonical
> publication style for this template — visible probe geography
> + visible γ structure that linear-scale crushes.

## Style choices

- **Layout:** `subplots(1, 6, sharey=True, gridspec_kw={"width_ratios": ...})`.
- **Colour scale (default):** per-panel (`vmin=0`, `vmax=panel max`).
- **Colorbars (always `imshow_colorbar_caxdivider`):** divider-based,
  compact, attached to the host axis.  *Never* `fig.colorbar` (that
  loosens spacing and breaks the project rule on multi-axis grids).
- **Per-panel mode (default):** one cb per panel.  Every panel
  surrenders the same fraction (`size="4%"`), so the row stays
  uniform.  Only the **right-most** cb carries the label
  `<|imcoh_{ij}|>_f` (generic frequency-averaged glyph; the row spans
  bands).
- **Shared mode (`--shared-scale`):** one cb on the right-most panel.
  The right-most panel's gridspec column is pre-padded by `_ROW_CB_FRAC`
  (≈ 6%) so that the cax + pad eat from the *extra* width and the
  imshow part of the last panel matches the others.
- **Titles:** TeX band glyphs from `BRAIN_BAND_TEX_DICT`.
- **Interpolation, vector PDF, no suptitle, watermark off** — same
  guarantees as `single_adjacency`.
- **File name:**
  `<Pat_NN>_<phase>_<fc_method>_<tick_labels>_<scale>_<share_tag>.pdf`
  where `share_tag ∈ {perpanel, shared}`.

## How to call

```python
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row

fig, axes = plot_fc_adjacency_row(
    matrices=[M_b for b in BRAIN_BANDS_NAMES],
    titles=[BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
    fc_method="imcoh_abs", bands=BRAIN_BANDS_NAMES,
    shared_scale=False,
)
```

```bash
python .agents/guides/05_plotting/fc_templates/row_per_band.py \
    --patient Pat_06 --phase rest_pre --fc-method imcoh_abs \
    --tick-labels generic
# add --log-scale to compress dynamic range
# add --shared-scale for one global colorbar instead of per-panel
```

## Anti-patterns

- ❌ Default shared linear scale — γ panels will appear blank.
- ❌ `set_rasterized(True)` (vector only).
- ❌ Repeating `<patient>` / `<fc_method>` in the panel titles.
