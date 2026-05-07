---
name: network_templates / single_network
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - single_network.py
  - README.md
  - mosaic_layout_compare.md
  - lrg_seeded.md
  - src/lrg_eegfc/visuals/network_templates.py
---

# `single_network` — one patient × one band × one phase

A single-axis network drawing of the FC graph at a chosen
`(patient, band, phase, fc_method)`.  The atomic template; everything
else (rows, mosaics, layout galleries) composes from this.

## What this figure shows

Edges of the full weighted FC graph rendered in 2-D with γ-power
scaled width and alpha (γ=6).  Same-probe pairs are red, cross-probe
pairs gray; this lets the reader spot probe-trivial modules at a
glance.  Nodes colored by electrode shaft by default.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `layout` | `spring / kk / spectral / laplacian_pca / community_grouped / backbone_guided / circular_by_shaft / sfdp / arf / lrg_kk / lrg_sfdp` | **`spring`** |
| `node_color` | `shaft / community / uniform` | `shaft` |
| `coloring` | `probe / signed / weight` (auto-default by `fc_method`) | `auto` |
| `n_communities` | LRG cut for `lrg_*` layouts and `node_color=community` | `10` |
| `node_size` | scatter marker size | `22` |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/single_network/` |

## Style choices

- **Edge recipe (canonical):** `t = (|w|/|w|_max)^γ`, γ=6,
  width=`(0.15, 4.0)`, alpha=`(0.03, 0.9)`.  Sort ascending →
  strongest edges drawn last.  `LineCollection` (single artist),
  **never `set_rasterized(True)`** — vector PDFs handle ~6.5k thin
  lines.
- **Probe coloring (default for magnitude FC):** same-probe red
  `(0.8, 0.2, 0.2)`, cross-probe gray `(0.35, 0.35, 0.35)`.  Reveals
  same-probe sEEG bias; if the apparent module structure is mostly
  red, it's geometry not connectivity.
- **Signed coloring (auto for `corr`/`imcoh`):** diverging `RdBu_r`
  on signed weight, `TwoSlopeNorm` centered at 0; symmetric
  `(-|w|_max, +|w|_max)`.  Width/alpha still drive off `|w|`.
- **Node coloring:** `tab20` cycled over unique probes (default) or
  over LRG community ids (`--node-color community`).  Nodes drawn at
  zorder=5, above edges (zorder=1), with white edgecolor for contrast.
- **No `fig.suptitle`, no watermark by default.** Provenance lives in
  the file name.
- **PDF only, full vector.**
- **File name:**
  `<Pat_NN>_<band>_<phase>_<fc_method>_<layout>_<node_color>_<coloring>.pdf`.

## How to call

```bash
# Default — Pat_05 β rest_pre imcoh_abs, spring layout, shaft colors
python .agents/guides/05_plotting/network_templates/single_network.py

# Different layout
python .agents/guides/05_plotting/network_templates/single_network.py \
    --patient Pat_06 --band alpha --phase task_learn --layout kk

# LRG-seeded layout + community coloring (the cluster figure variant)
python .agents/guides/05_plotting/network_templates/single_network.py \
    --patient Pat_05 --band beta --phase rest_post \
    --layout lrg_sfdp --node-color community --n-communities 10

# Signed FC (corr) — auto-switches to RdBu_r diverging coloring
python .agents/guides/05_plotting/network_templates/single_network.py \
    --patient Pat_02 --band theta --phase rest_pre --fc-method corr
```

## Anti-patterns

- ❌ `nx.draw(G, node_color="lightblue", width=widths)` — no
  zorder, no probe coloring, no γ-power; produces a gray hairball.
- ❌ `set_rasterized(True)` on the LineCollection — violates the
  no-rasterization rule.
- ❌ `fig.suptitle(...)` — context goes in the file name.
- ❌ Threshold the matrix before laying out (e.g.
  `A[A < q90] = 0`) — drops the layout's structural information.
  Use `--layout backbone_guided` if you want a sparse-layout-on-full-draw,
  or pre-process via `disparity_backbone(A)` and pass the result.
- ❌ Plain `--layout spring` on raw ImCoh without acknowledging the
  hairball — for paper figures prefer `kk`, `lrg_sfdp`, or
  `backbone_guided`.  Spring is the *default* for compatibility, not
  the recommended publication choice.
