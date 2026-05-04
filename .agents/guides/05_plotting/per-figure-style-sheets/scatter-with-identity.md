---
name: per-figure-style-sheets / scatter-with-identity
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
---

# Scatter with identity (per-band cohort scatter)

## Use when

Per-band per-patient scatter: x = quantity in phase pair A, y = same
quantity in phase pair B. Identity dashed line. Used for persistence
verdicts (`final_verdict.pdf` page 1, `Td_dS_vs_dP_scatter.pdf`,
`Td_dS_vs_dF_scatter.pdf`).

## Recipe

```python
fig, axes = plt.subplots(1, n_bands,
                         figsize=(2.2 * n_bands + 0.6, 2.6),
                         squeeze=False)
for c, band in enumerate(bands):
    ax = axes[0, c]
    sub = data[data.band == band].dropna(subset=["x", "y"])
    in_pool_excl = sub[sub.patient != "Pat_03"]
    out_p = sub[sub.patient == "Pat_03"]

    lo = float(min(sub.x.min(), sub.y.min())) * 1.05
    hi = float(max(sub.x.max(), sub.y.max())) * 1.05
    ax.plot([lo, hi], [lo, hi], color="k", lw=0.6, ls="--", alpha=0.6)
    ax.scatter(in_pool_excl.x, in_pool_excl.y, c="tab:blue",
               edgecolor="k", linewidth=0.4, s=32, zorder=3)
    if not out_p.empty:
        ax.scatter(out_p.x, out_p.y, c="tab:orange", marker="^",
                   edgecolor="k", linewidth=0.4, s=42, zorder=4)
    # patient labels (all patients)
    for _, r in sub.iterrows():
        ax.text(r.x, r.y, r.patient.replace("Pat_", ""),
                fontsize=4.5, va="bottom", ha="left", alpha=0.7)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
    ax.set_xlabel("x label", fontsize=8)
    if c == 0:
        ax.set_ylabel("y label", fontsize=9)
    ax.tick_params(labelsize=6)

# FIGURE-LEVEL legend (NEVER per-axis for shared content)
handles = [
    plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
               markeredgecolor="k", markersize=6, label="in-pool (n=9)"),
    plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
               markeredgecolor="k", markersize=7, label="Pat_03"),
    plt.Line2D([], [], color="k", lw=0.6, ls="--", label="identity"),
]
fig.legend(handles=handles, loc="lower center",
           bbox_to_anchor=(0.5, -0.04),
           ncol=len(handles), frameon=False, fontsize=8)
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
# watermark off by default; opt-in via:
#   from lrg_eegfc.visuals.layout import add_provenance_footer
#   add_provenance_footer(fig, "<script_id>")
```

## Key choices

- **Identity line** dashed black, slightly transparent, drawn first
  (lowest z-order) so points sit on top.
- **Cohort distinction**: in-pool blue circles, Pat_03 orange triangle.
- **Equal aspect** (`set_aspect("equal")`) — points are deviations
  from the identity, so visual deviation must be meaningful.
- **Patient labels** on each point in `fontsize=4.5` — readers like
  to identify outliers.
- **No `n_below` annotation inside a single panel** if the figure is
  the only one in the file. If the same scatter is repeated for
  multiple distances on different pages, add the count box per panel.
- **Figure-level legend** at the bottom centre, horizontal — see
  [`../legends.md`](../legends.md).

## Variants

- For **two-axis-line crosshairs** (e.g. zero-zero quadrants) add:
  ```python
  ax.axhline(0, color="gray", lw=0.4, alpha=0.6)
  ax.axvline(0, color="gray", lw=0.4, alpha=0.6)
  ```
  before the data.

- For **3-distance overlay** (one panel showing 3 phase-pair scatter
  series), use the project phase-pair colour code from
  [`../colormaps-and-styles.md`](../colormaps-and-styles.md).
