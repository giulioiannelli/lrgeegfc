---
name: network_templates / mosaic_layout_compare
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - mosaic_layout_compare.py
  - single_network.md
  - lrg_seeded.md
  - README.md
---

# `mosaic_layout_compare` — every layout side-by-side, one case

A 3 × 4 grid of panels showing the same FC matrix rendered with every
layout in `DEFAULT_GALLERY` (`spring`, `kk`, `spectral`,
`laplacian_pca`, `community_grouped`, `backbone_guided`,
`circular_by_shaft`, `sfdp`, `arf`, `lrg_kk`, `lrg_sfdp`).  LRG-seeded
panels are dropped silently when no LRG cache exists.

## What this figure shows

Visual layout sensitivity for one (patient, band, phase, fc_method).
Modules that emerge in *most* layouts are real; modules that emerge in
only one are layout-driven artifacts.  Use this figure when you don't
yet know which layout to pick — the "test all" comparison.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `layouts` | comma-separated subset (default: all 11) | `None` (= full gallery) |
| `node_color` | `shaft / community / uniform` | `shaft` |
| `coloring` | `probe / signed / weight` (auto-default) | `auto` |
| `n_communities` | LRG cut for LRG-seeded panels | `10` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/mosaic_layout_compare/` |

## Style choices

- **Same FC matrix in every panel.**  The only thing varying is the
  layout function.
- **Same edge γ-recipe and node-color scheme** in every panel.  Only
  positions change.
- **3 × 4 grid** by default (n=11 layouts; one cell empty if all 11
  drawn, fewer if LRG layouts are dropped).
- **Panel size** 3.0 × 3.0; **figure size** 12 × 9.
- **Per-panel title** is the layout name (`spring`, `kk`, …).  No
  `fig.suptitle`.
- **File name:**
  `<Pat_NN>_<band>_<phase>_<fc_method>_<node_color>_<coloring>.pdf`.

## How to call

```bash
# Full gallery (default):
python .agents/guides/05_plotting/network_templates/mosaic_layout_compare.py \
    --patient Pat_05 --band beta --phase rest_pre

# Restrict to a subset:
python .agents/guides/05_plotting/network_templates/mosaic_layout_compare.py \
    --patient Pat_06 --band alpha --phase task_test \
    --layouts kk,spectral,sfdp,lrg_sfdp

# With community coloring (highlights LRG modules):
python .agents/guides/05_plotting/network_templates/mosaic_layout_compare.py \
    --patient Pat_05 --band beta --phase rest_pre \
    --node-color community
```

## Anti-patterns

- ❌ Using this template for paper figures — it's a *diagnostic*, not
  a publication artefact.  Once you've picked a layout, use
  `single_network` or one of the row/mosaic templates.
- ❌ Drawing only `spring` and `kk` and concluding both layouts agree —
  they share the same hairball failure mode on near-uniform weights.
  The gallery must include `spectral`, `laplacian_pca`, and
  `lrg_sfdp` to be informative.
- ❌ Calling this the "comparison figure" without labelling each
  panel with its layout name — readers won't know which is which.
