---
name: 05_plotting / legends
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Legends — figure-level, never axis-level for shared content

## The rule

If two or more axes in a figure share the same legend entries, the
legend lives on the **figure**, not inside any axis. Acceptable
positions, in order of preference:

1. **Bottom centre, horizontal** — for landscape (wide) figures.
2. **Top centre, horizontal** — when the figure has tight bottom
   labels (e.g. x-tick rotation already steals space below).
3. **Right centre, vertical** — only when the figure is tall and the
   right margin is empty.

Set `ncol` so the legend uses the figure's full width without
overflowing — typically `ncol = len(handles)` for ≤ 6 entries,
otherwise wrap to two rows.

## Why

Legend inside one axis breaks panel symmetry: that axis becomes
visually denser, the others "feel" empty. Tight-layout then over-pads
to make room. The whole figure looks lopsided.

## Recipe — the canonical pattern

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 6, figsize=(15, 2.6), sharey=True)
# … plotting on axes …
handles = [
    plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
               markeredgecolor="k", markersize=6, label="in-pool"),
    plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
               markeredgecolor="k", markersize=7, label="Pat_03"),
    plt.Line2D([], [], color="k", lw=0.6, ls="--", label="identity"),
]
fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.04),
    ncol=len(handles),
    frameon=False,
    fontsize=8,
)
fig.tight_layout(rect=[0, 0.04, 1, 1])  # leave bottom strip for legend
fig.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close(fig)
```

Key flags:

- `loc="lower center"` + `bbox_to_anchor=(0.5, -0.04)` — anchor on
  the figure, not an axis.
- `ncol=len(handles)` — horizontal layout, single row, balanced.
- `frameon=False` — no boxed border.
- `fig.tight_layout(rect=[0, 0.04, 1, 1])` — reserve the bottom strip
  so the legend doesn't collide with x-tick labels.

## Top-centre variant

```python
fig.legend(handles=handles, loc="upper center",
           bbox_to_anchor=(0.5, 1.02), ncol=len(handles),
           frameon=False, fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.96])
```

## Right-centre vertical variant (tall figures only)

```python
fig.legend(handles=handles, loc="center right",
           bbox_to_anchor=(1.02, 0.5), ncol=1,
           frameon=False, fontsize=8)
fig.tight_layout(rect=[0, 0, 0.94, 1])
```

## Anti-patterns

- ❌ `ax.legend(loc="upper right")` for a legend that describes other
  axes too. Always promote to `fig.legend`.
- ❌ A legend in one panel of a multi-panel figure. Symmetry breaker.
- ❌ Multi-row figure-legend at the side. Use bottom or top instead.
- ❌ `bbox_to_anchor` referring to axis coordinates instead of figure
  coordinates. Use figure coordinates: `(0.5, -0.04)` is half-way
  across the figure, slightly below it.

## Existing examples to copy from

- `src/lrg_eegfc/visuals/cross_condition.py` already uses
  `fig.legend(... bbox_to_anchor=(0.5, -0.01) ... ncol=…)` —
  reference implementation. Imitate.

## Mini-helper (not yet written; see open task)

The same five-line `fig.legend(...)` block is repeated across
scripts. Promote into `lrg_eegfc.visuals.layout.figure_legend`:

```python
def figure_legend(
    fig,
    handles,
    *,
    where: str = "bottom",   # "bottom" | "top" | "right"
    fontsize: int = 8,
    pad: float = 0.04,
) -> None:
    """Place a figure-level legend in the canonical position.

    Equivalent to fig.legend(...) with the project's defaults.
    Auto-sets ncol from len(handles) and reserves space in
    fig.tight_layout via the rect kwarg returned by the caller.
    """
    cfg = {
        "bottom": dict(loc="lower center", bbox_to_anchor=(0.5, -pad)),
        "top": dict(loc="upper center", bbox_to_anchor=(0.5, 1 + pad)),
        "right": dict(loc="center right", bbox_to_anchor=(1 + pad, 0.5)),
    }[where]
    ncol = len(handles) if where != "right" else 1
    fig.legend(handles=handles, ncol=ncol, frameon=False,
               fontsize=fontsize, **cfg)
```

When this helper exists, replace the raw `fig.legend(...)` calls
across `visuals/` and `scripts/` and remove the duplicates.
