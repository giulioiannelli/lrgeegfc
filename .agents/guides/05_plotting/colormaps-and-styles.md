---
name: 05_plotting / colormaps-and-styles
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Colormaps, palettes, line styles, fonts

## Colormaps (data → colormap)

| data class | colormap | range convention |
|---|---|---|
| FC magnitude (`imcoh_abs`, MSC, `corr` ≥ 0) | `magma` (default) | `vmin=0`, `vmax=data.max()` per row |
| Distance (≥ 0) | `viridis` | per-row vmin=0, vmax=row max |
| Signed difference | `RdBu_r` | symmetric `±99th percentile of |D|` |
| Categorical (partition labels) | `tab10` (≤ 10 cats) / `tab20` (≤ 20) | nominal |
| Phase / angle | `twilight` (cyclic) | `vmin=-π`, `vmax=π` |
| Continuous "task ↔ rest" axis | custom diverging from `lrgsglib.plotlib.colormaps` | centred at 0.5 |

Source the project's named colormaps from
`lrgsglib.plotlib.colormaps` rather than redefining inline.

## Phase-pair colour code (project-wide)

| pair | colour |
|---|---|
| RPre → TT | `tab:red` |
| RPre → RPost | `tab:blue` |
| TT → RPost | `tab:green` |
| TL → TT (within-task) | `tab:purple` |
| RPre → TL | `tab:olive` |
| TL → RPost | `tab:cyan` |
| Pat_03 marker | `tab:orange`, `marker="^"` |

These colours are already used in `audit_25/26/28` — keep them
consistent across new figures so a glance at any plot identifies the
phase pair.

## Line styles

- **Solid (`-`)**: observed quantities, primary data.
- **Dashed (`--`)**: identity / reference / null / fit lines.
- **Dotted (`:`)**: thresholds (e.g. `n_+ = 8/10`, `q = 0.05`).

## Fonts

| element | size |
|---|---|
| panel title (column header) | 10 |
| axis label | 8–9 |
| tick label | 6–7 |
| in-axis annotation | 6–7 |
| legend / footer | 6–8 |
| `fig.text(provenance)` | 6, color=gray (opt-in only; off by default) |

Smaller for high-density grids (1 row × 6 columns); slightly larger
for single-axis figures.

## Per-figure style sheets

For repeated figure families, put the canonical style for that
family under `per-figure-style-sheets/`:

- `heatmap-grid.md` — colorbar placement, vmin/vmax conventions,
  rasterisation flag.
- `scatter-with-identity.md` — identity dashed line, marker per
  cohort group, axis aspect.
- `swarm-or-strip.md` — jitter width, median bar, sign-count
  annotation.
- `network-circular.md` — node size, edge alpha by weight,
  threshold quantile.
- `4xN-phase-grid.md` — 4-phase distance heatmap convention from
  `audit_26`.

Add new sheets here whenever a new figure family is invented and
likely to be reused.
