---
name: 05_plotting / colormaps-and-styles
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-05-10
---

# Colormaps, palettes, line styles, fonts

## Colormaps (data → colormap)

**Project default:** `inferno`. Set in
[`src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle`](../../../src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle)
as `image.cmap: inferno` — every `imshow` / `pcolormesh` that does not
explicitly pass `cmap=` picks it up automatically once
`use_lrg_style()` has run. Inferno is perceptually uniform,
colour-blind friendly, and runs dark→bright sequentially (low values
dark purple, high values bright yellow).

| data class | colormap | range convention |
|---|---|---|
| **default (any sequential ≥ 0 data)** | **`inferno` (rcParam, no override needed)** | `vmin=0`, `vmax` data-driven |
| FC magnitude (`imcoh_abs`, MSC, `corr` ≥ 0) | `inferno` (default) | `vmin=0`, `vmax=data.max()` per row |
| Distance / `D̂(τ)` (≥ 0) | `inferno` (default) — small=dark, large=bright | per-page vmin/vmax (1st/99th pct) |
| Signed difference | `RdBu_r` (explicit override) | symmetric `±99th percentile of |D|` |
| Categorical (partition labels) | `tab10` (≤ 10 cats) / `tab20` (≤ 20) | nominal |
| Phase / angle | `twilight` (cyclic) | `vmin=-π`, `vmax=π` |
| Continuous "task ↔ rest" axis | custom diverging from `lrgsglib.plotlib.colormaps` | centred at 0.5 |

Override per-call with `cmap="..."` only when the data semantics
disagree with the dark→bright sequential default — e.g. `cmap="RdBu_r"`
for signed Δ matrices, `cmap="viridis_r"` if you specifically want
small distances to be bright instead of dark. Source the project's
named colormaps from `lrgsglib.plotlib.colormaps` rather than
redefining inline.

### Probe-outline overlay color

`draw_probe_outlines(ax, channel_labels, …)` defaults to **`color="cyan"`**
(see `lrg_eegfc.visuals.draw_probe_outlines`). Cyan is the cmap-complement
of inferno's bright yellow end and contrasts strongly against the dark
purple low end — visible on both regimes a same-probe diagonal block can
sit on (high FC = bright, small distance = dark). The previous
`#FFD600` (yellow) default blended with inferno's bright high-FC cells.
Override per-call only when overlaying on a non-default cmap (e.g. on
`viridis_r`, prefer `color="black"`).

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

## Project mplstyle (canonical, activate at top of every figure script)

The font / tick / output defaults above are encoded in a single
canonical matplotlib style sheet — **the project's only mplstyle**.

- **File:** [`src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle`](../../../src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle)
- **Loader:** `from lrg_eegfc.visuals.styles import use_lrg_style`
- **Usage:** call `use_lrg_style()` once at the top of any figure
  script (after imports, before the first `subplots`).

```python
from lrg_eegfc.visuals.styles import use_lrg_style
use_lrg_style()
```

What it sets (non-exhaustive — see the file for the full list):

| group | keys |
|---|---|
| fonts | `font.size=9`, `axes.titlesize=10`, `axes.labelsize=9`, `xtick.labelsize=7`, `ytick.labelsize=7`, `legend.fontsize=7` |
| ticks | `xtick.major.size=3`, `xtick.minor.size=1.5`, widths 0.6 / 0.4 |
| axes | `axes.linewidth=0.6`, `axes.titlepad=4`, `axes.labelpad=3` |
| images | `image.cmap=inferno` (default for any imshow / pcolormesh) |
| output | `pdf.fonttype=42`, `ps.fonttype=42`, `svg.fonttype=none`, `savefig.dpi=300` |

`pdf.fonttype=42` (TrueType embed) is the load-bearing entry — it
keeps PDF text editable in Illustrator / Inkscape rather than
shape-traced.

For a one-off override (e.g. a poster figure with larger ticks), use
`with rc_context({...}):` inside the script — do **not** edit the
mplstyle file for one-off needs. The mplstyle is the cohort-wide
default; per-figure deviations stay local.

### LogNorm colorbar tick policy

`_apply_factored_sci_format` (in
`src/lrg_eegfc/visuals/fc_templates.py`) is the canonical end-of-cb
call for any LogNorm colorbar created via
`imshow_colorbar_caxdivider`:

```python
from lrg_eegfc.visuals.fc_templates import _apply_factored_sci_format
_, _, clb = imshow_colorbar_caxdivider(im, ax, size="4%", pad=0.06)
clb.set_label(r"$\hat{D}(\tau) = 1/\hat{\rho}(\tau)$")
_apply_factored_sci_format(clb, axis_orientation="vertical")
```

Behaviour:

- **< 1.5 decade span** — replaces tick labels with bare mantissas
  (`1, 3, 9`) and emits one `× 10ⁿ` above the cb.
- **≥ 1.5 decade span** — keeps matplotlib's clean `10ⁿ` major-tick
  labels but forces the minor formatter to `NullFormatter` so
  inline `2 × 10ⁿ` mantissa labels never appear.

Inline `2 × 10⁻²` style minor-tick labels on a LogNorm cb are
**banned** — they crowd the cb and dwarf neighbouring panels'
labels.

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
