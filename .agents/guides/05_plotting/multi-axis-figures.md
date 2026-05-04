---
name: 05_plotting / multi-axis-figures
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Multi-axis figures — layout, shared scales, figure-level decoration

## Default layout

```python
fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(per_col_w * ncols + 0.6, per_row_h * nrows + 0.4),
    squeeze=False,                 # always 2-D axes array
    sharex="col",                  # if same x meaning per column
    sharey="row",                  # if same y meaning per row
)
```

`squeeze=False` is mandatory — never branch on
`if isinstance(axes, np.ndarray)`. Always index `axes[r, c]`.

## Per-column / per-row width budgets

- Heatmap panels: `per_col_w ≈ 2.0`, `per_row_h ≈ 2.0`.
- Scatter / swarm panels: `per_col_w ≈ 2.2`, `per_row_h ≈ 2.6`.
- Network / circular layout: `per_col_w ≈ 2.0`, `per_row_h ≈ 2.2`,
  `aspect="equal"`.

The extra `+ 0.6` / `+ 0.4` accounts for shared y-labels and shared
top titles. Adjust empirically.

## Shared axis labels (figure-level)

Per-column titles, per-row y-labels — both go on the leftmost / topmost
axis of the row / column, never `suptitle`:

```python
for c, band in enumerate(bands):
    axes[0, c].set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
for r, distance in enumerate(distances):
    axes[r, 0].set_ylabel(f"$d_{distance}$", fontsize=9)
```

Bottom x-label:

```python
for c in range(ncols):
    axes[-1, c].set_xlabel("gap (s)", fontsize=8)
```

## Provenance footer (watermark) — opt-in only

A grey provenance string in the bottom-right corner is **off by
default**. The file name already carries patient / band / phase /
fc_method, so a watermark just duplicates it. Add a watermark only
when the figure will be detached from its file name (slide deck,
screenshot, e-mail attachment).

When you do want one, use the library helper:

```python
from lrg_eegfc.visuals.layout import add_provenance_footer
add_provenance_footer(fig, f"audit_25 P3 {patient}")
```

Or expose it through a `watermark=False` kwarg on the plotting
function (see `fc_templates/single_adjacency.py` for the canonical
shape).

## Tight layout vs. constrained layout

- Use `fig.tight_layout()` as default.
- For figure-level legends or colorbars, use the `rect` argument:
  `fig.tight_layout(rect=[left, bottom, right, top])` to reserve
  margin strips. Examples in [`legends.md`](legends.md).
- Avoid `constrained_layout=True` — it conflicts with the
  `imshow_colorbar_caxdivider` divider mechanism.

## Anti-patterns

- ❌ `fig.suptitle(...)` on publication figures (already in
  `04_rules/never-always-list.md`).
- ❌ Per-axis `set_title` repeating cohort or patient ID — that's
  redundant and clutters the grid. Put it once in the file name
  (default) or in an opt-in watermark.
- ❌ A grey provenance watermark on by default — opt-in only.
- ❌ `set_rasterized(True)` — vector is the only allowed mode.
- ❌ Mixing free `subplots` and `gridspec` — pick one. `gridspec`
  is for non-uniform grids (different panel widths), `subplots` is
  for uniform grids.
- ❌ Per-axis legends when content is shared. See
  [`legends.md`](legends.md).
