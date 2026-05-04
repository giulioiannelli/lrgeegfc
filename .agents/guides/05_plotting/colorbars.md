---
name: 05_plotting / colorbars
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Colorbars — `imshow_colorbar_caxdivider` is the default

## The rule

Every `imshow` / `pcolormesh` / `matshow` figure must use
`imshow_colorbar_caxdivider` from `lrgsglib.plotlib` (canonical) —
**not** raw `fig.colorbar(im, ax=ax)`. The helper guarantees:

- Properly spaced colorbar that does not squeeze the parent axis.
- Consistent size relative to the axis (default `5%`).
- Consistent padding (default `0.05`).
- Per-row sharing on multi-axis grids (place the colorbar on the
  rightmost axis of each row, only).

## Canonical import

```python
from lrgsglib.plotlib import imshow_colorbar_caxdivider
```

There is also a duplicate at
`src/lrg_eegfc/visuals/correlation.py` — **do not import that one**.
The lrgsglib version is the canonical one (more featured: handles
3D axes, multi-orientation, returns the divider). The duplicate
should be removed; import path migration is an open coding task.

## Recipe — single-axis

```python
fig, ax = plt.subplots(figsize=(4, 4))
im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1, aspect="equal")
# never call set_rasterized(True) — keep imshow vector
imshow_colorbar_caxdivider(im, ax, size="5%", pad=0.05)
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
fig.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close(fig)
```

## Recipe — per-row colorbar on a multi-axis grid

When rows show different distances / metrics with different scales,
attach **one colorbar per row** at the rightmost axis:

```python
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(2 * ncols + 1, 2 * nrows + 0.4),
                         squeeze=False)
for r, scale in enumerate(row_scales):
    last_im = None
    for c, datum in enumerate(row_data[r]):
        ax = axes[r, c]
        im = ax.imshow(datum, cmap="viridis",
                       vmin=scale.vmin, vmax=scale.vmax,
                       aspect="equal")
        last_im = im
        ax.set_xticks([]); ax.set_yticks([])
    # one colorbar per row, on the rightmost axis
    if last_im is not None:
        imshow_colorbar_caxdivider(last_im, axes[r, -1],
                                   size="4%", pad=0.06)
```

## Recipe — shared colorbar across all panels

When every panel shares the same scale (`vmin`/`vmax` identical),
use one figure-level colorbar at the right or bottom:

```python
fig, axes = plt.subplots(nrows, ncols, ..., squeeze=False)
for r, c in product(range(nrows), range(ncols)):
    im = axes[r, c].imshow(..., vmin=vmin, vmax=vmax)
# Reserve right strip for figure-level colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.tight_layout(rect=[0, 0, 0.91, 1])
```

`imshow_colorbar_caxdivider` is the default for per-axis or per-row
colorbars; figure-level colorbars are the exception (reserved for
the case "every panel is on the same scale").

## Sequential vs divergent colormaps

| data | colormap | rationale |
|---|---|---|
| FC magnitude (`imcoh_abs`, MSC, `corr` ≥ 0) | `magma` (default) or `viridis` | sequential, monotone-perceived |
| Differences (`A − B`) | `RdBu_r` | divergent, centred at zero |
| Categorical (partition labels) | `tab10` / `tab20` | qualitative |
| Phase (angle) | `twilight` (cyclic) | wraps cleanly at ±π |

For divergent maps, **always centre `vmin/vmax` at zero**:

```python
vmax = float(np.percentile(np.abs(D), 99)) or 1e-9
ax.imshow(D, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
```

## Anti-patterns

- ❌ `fig.colorbar(im, ax=ax)` — squeezes the parent axis, inconsistent
  sizing across panels.
- ❌ One colorbar per panel on a multi-panel figure with identical
  scales — visual noise. Use one shared colorbar.
- ❌ Magma/viridis on signed differences — the eye reads zero as
  whatever colour is in the middle of the colormap, which is rarely
  white. Use `RdBu_r`.
- ❌ Default `vmin/vmax` (auto-scaled per panel) on a multi-panel
  figure where comparison is the point. Always set shared scales.
