---
name: network_templates / mosaic_band_phase
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - mosaic_band_phase.py
  - row_per_phase.md
  - row_per_band.md
  - README.md
---

# `mosaic_band_phase` — one patient · 4 phases × 6 bands

A 4 × 6 mosaic of FC network drawings.  Rows iterate over phases,
columns over bands; patient and `fc_method` fixed.  The network
analogue of `fc_templates/mosaic_band_phase`.

## What this figure shows

The full band × phase decomposition of one patient's FC graph at a
glance.  Reading down a column reveals phase-induced reorganisation
within one band; reading along a row reveals spectral structure at
one phase.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `layout` | layout algorithm (LRG-seeded excluded) | **`spring`** |
| `node_color` | `shaft / uniform` | `shaft` |
| `coloring` | `probe / signed / weight` (auto-default) | `auto` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/mosaic_band_phase/` |

## Style choices

- **Independent layouts per panel.**  Each (phase, band) cell gets
  its own layout computation.  This gives the most honest per-cell
  view but loses cross-cell node correspondence.  For figures where
  positional consistency matters more than per-cell accuracy, use
  `lrg_seeded` with task-test community labels.
- **Per-panel γ-power normalisation** — each panel's edge alpha/width
  scaled to its own `|w|_max`.  γ$_l$/γ$_h$ panels won't collapse
  to black, but absolute amplitude is not comparable between cells.
- **Panel size** 2.4 × 2.4; **figure size** 14.4 × 9.6 (4 × 6).
- **Row labels** (left side, rotated): phase TeX glyphs.
  **Column titles** (top): band TeX glyphs.  No `fig.suptitle`.
- **Spacing:** `wspace=0.05`, `hspace=0.10` — packed tight.
- **File name:**
  `<Pat_NN>_<fc_method>_<layout>_<node_color>_<coloring>.pdf`.

## How to call

```bash
python .agents/guides/05_plotting/network_templates/mosaic_band_phase.py \
    --patient Pat_05 --layout kk
# Pat_06 with sfdp:
python .agents/guides/05_plotting/network_templates/mosaic_band_phase.py \
    --patient Pat_06 --layout sfdp
```

## Anti-patterns

- ❌ A single global layout shared across all 24 panels — masks
  reorganisation; nodes don't move so cells look identical.
- ❌ Drawing all 24 panels with the same `|w|_max` — γ panels go
  black.
- ❌ `fig.suptitle("...")` — band/phase glyphs in row/col labels are
  enough.
