---
name: per-figure-style-sheets / cohort-conventions
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
---

# Cohort visual conventions (Pat_03, in-pool, n=10)

## Pat_03 in figures

- **Always include** Pat_03 in the cohort statistics (n=10).
  Z-scores and dimensionless ratios are unaffected by the 1024 Hz
  sampling rate.
- **Mark Pat_03 distinctly**: orange triangle (`marker="^"`,
  `c="tab:orange"`, edgecolor `k`, slightly larger `s=42` than
  in-pool circles `s=32`).
- The other 9 patients are plain blue circles.
- Patient labels (last two digits, `fontsize=4.5`, slight transparency)
  on every point if the figure is per-band scatter.

## Statistic = n=10, visual = blue+orange

Cohort medians, sign counts, Wilcoxon tests, BH-FDR — all over
**n=10**. The visual split (9 blue + 1 orange) is purely a marker
distinction, not a sample-size split.

## Phase-pair colour code

(repeated from `colormaps-and-styles.md` for quick reference)

| pair | colour |
|---|---|
| RPre → TT | `tab:red` |
| RPre → RPost | `tab:blue` |
| TT → RPost | `tab:green` |
| TL → TT (within-task) | `tab:purple` |
| RPre → TL | `tab:olive` |
| TL → RPost | `tab:cyan` |

## Identity-line on cohort scatter

Always include an identity reference line, dashed black, low alpha:

```python
ax.plot([lo, hi], [lo, hi], color="k", lw=0.6, ls="--", alpha=0.6)
```

so deviations from the diagonal are immediately readable.

## Cohort threshold annotations

If the figure shows a count vs. threshold (e.g. `n_+ ≥ 8/10`), add
a dashed horizontal line at the threshold:

```python
ax.axhline(8, color="k", lw=0.6, ls="--", alpha=0.7)
```

## Provenance footer (watermark) — opt-in only

A grey footer is **off by default**. The file name already carries
the script id and cohort tag. Add a footer only when the figure will
travel detached from its file name (slide deck, screenshot).

```python
from lrg_eegfc.visuals.layout import add_provenance_footer
# opt-in:
add_provenance_footer(fig, "<script_id> P<n> cohort")
```
