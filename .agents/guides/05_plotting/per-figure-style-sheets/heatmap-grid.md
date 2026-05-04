---
name: per-figure-style-sheets / heatmap-grid
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
---

# Heatmap-grid (multi-axis `imshow` panel)

## Use when

Showing matrices side-by-side: FC matrices across phases × bands,
distance matrices across distances × bands, etc.

## Recipe

```python
from lrgsglib.plotlib import imshow_colorbar_caxdivider

fig, axes = plt.subplots(nrows, ncols,
                         figsize=(2 * ncols + 1.0, 2 * nrows + 0.4),
                         squeeze=False)

for r in range(nrows):
    last_im = None
    # determine row vmax for shared scale across the row
    row_max = max(M.max() for M in row_data[r])
    for c in range(ncols):
        ax = axes[r, c]
        im = ax.imshow(row_data[r][c], cmap="magma",
                       vmin=0, vmax=row_max, aspect="equal")
        # never set_rasterized — keep imshow vector
        last_im = im
        ax.set_xticks([]); ax.set_yticks([])
        if r == 0:
            ax.set_title(col_label[c], fontsize=10)
        if c == 0:
            ax.set_ylabel(row_label[r], fontsize=9)
    imshow_colorbar_caxdivider(last_im, axes[r, -1],
                               size="4%", pad=0.06)

fig.tight_layout()
fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
```

## Key choices

- **One colorbar per row** at the rightmost axis — not per panel,
  not figure-level (because each row may have a different scale).
- **Shared row scale** — pick `vmax = max(M.max() for M in row)` so
  panels in a row are visually comparable.
- **Full vector**. Never `im.set_rasterized(True)` — vector is the
  default and the only allowed mode for FC matrices. See
  [`../output-and-rasterization.md`](../output-and-rasterization.md).
- **No legend** for heatmaps, ever; the colorbar is the key.
- **No watermark / provenance footer by default**. If you want one
  for a slide, add it explicitly via
  `lrg_eegfc.visuals.layout.add_provenance_footer(fig, label)`.

## Variants

- For divergent data (e.g. `A − B`), use `cmap="RdBu_r"` and
  `vmin = -vmax = -np.percentile(np.abs(D), 99)`. See
  [`../colorbars.md`](../colorbars.md) §divergent.
- For 4×4 small distance matrices with cell annotations, write the
  numeric value with `ax.text(j, i, f"{v:.2f}", ha="center",
  va="center", fontsize=5, color="white" if v > row_max*0.5 else
  "black")`. See `audit_26_fc_phase_geometry.py` page 1.
