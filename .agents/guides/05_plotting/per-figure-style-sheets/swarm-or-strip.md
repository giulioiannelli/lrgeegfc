---
name: per-figure-style-sheets / swarm-or-strip
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
---

# Swarm / strip plot (per-band cohort dispersion)

## Use when

Showing the per-patient distribution of a scalar (e.g. T_d) across
several discrete groups (distances or bands). Median bar overlay.

## Recipe

```python
rng = np.random.default_rng(42)
fig, axes = plt.subplots(1, n_bands,
                         figsize=(2.4 * n_bands + 0.6, 2.7),
                         squeeze=False, sharey=True)
for c, band in enumerate(bands):
    ax = axes[0, c]
    sub = data[data.band == band].dropna(subset=cols)
    in_pool = sub                          # n=10
    in_pool_excl = sub[sub.patient != "Pat_03"]
    out_p = sub[sub.patient == "Pat_03"]
    for x, col in enumerate(cols):
        ys_blue = in_pool_excl[col].to_numpy()
        ys_full = in_pool[col].to_numpy()  # for stats
        xs = x + rng.uniform(-0.12, 0.12, size=len(ys_blue))
        ax.scatter(xs, ys_blue, c="tab:blue",
                   edgecolor="k", linewidth=0.3, s=24, zorder=2)
        if not out_p.empty:
            ax.scatter([x], [out_p[col].iloc[0]],
                       c="tab:orange", marker="^",
                       edgecolor="k", linewidth=0.3, s=36, zorder=4)
        # median bar (n=10)
        ax.plot([x - 0.25, x + 0.25],
                [np.median(ys_full)] * 2,
                color="k", lw=1.2, zorder=3)
    ax.axhline(0, color="gray", lw=0.4, alpha=0.6)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, fontsize=7)
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
    if c == 0:
        ax.set_ylabel("metric", fontsize=8)
    ax.tick_params(labelsize=6)

# figure-level legend
handles = [
    plt.Line2D([], [], marker="o", linestyle="", color="tab:blue",
               markeredgecolor="k", markersize=6, label="in-pool"),
    plt.Line2D([], [], marker="^", linestyle="", color="tab:orange",
               markeredgecolor="k", markersize=7, label="Pat_03"),
    plt.Line2D([], [], color="k", lw=1.2, label="median (n=10)"),
]
fig.legend(handles=handles, loc="lower center",
           bbox_to_anchor=(0.5, -0.04),
           ncol=len(handles), frameon=False, fontsize=8)
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(out_path, dpi=200, bbox_inches="tight"); plt.close(fig)
```

## Key choices

- **Pat_03 separated visually** but **statistics use n=10**: jitter
  blue dots over the in-pool 9, plot Pat_03 as orange triangle, but
  compute median and sign-count over all 10.
- **Jitter width** 0.12 inside a column (column spacing 1.0) — tight
  enough that columns don't overlap, loose enough that overlapping
  points are separable.
- **Median bar** `lw=1.2` over both subsamples — a single black
  segment at the median across all 10.
- **Sign-count annotation** (e.g. `n_neg/n_tot`) only when useful as
  the figure's headline read. Place at the top of the column with
  `bbox` around the text for legibility.

## Sign-count annotation pattern

```python
n_neg = int((ys_full < 0).sum()); n_tot = len(ys_full)
ax.text(x, ax.get_ylim()[1], f"{n_neg}/{n_tot}",
        fontsize=6, ha="center", va="top",
        bbox=dict(boxstyle="round,pad=0.18",
                  facecolor="white",
                  edgecolor="gray", linewidth=0.3))
```

## Anti-patterns

- ❌ `seaborn.swarmplot` — pulls in another plotting style; we use
  matplotlib jitter directly.
- ❌ Box-and-whiskers — too dense for n=10; show points + median.
- ❌ Sign-count over only n=9 (excluding Pat_03). The cohort is n=10;
  see [`per-figure-style-sheets/cohort-conventions.md`](cohort-conventions.md).
