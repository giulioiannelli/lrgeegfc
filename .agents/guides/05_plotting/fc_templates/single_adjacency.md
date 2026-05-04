---
name: fc_templates / single_adjacency
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-04-29
updated: 2026-04-29
pointers:
  - single_adjacency.py
  - ../colormaps-and-styles.md
  - ../colorbars.md
  - ../output-and-rasterization.md
---

# `single_adjacency` — one FC matrix, one panel

The minimal canonical figure for a functional-connectivity adjacency
matrix: one patient × one band × one phase × one `fc_method`.

> **Status:** Round 1. We are still iterating with the user. Once
> pinned, this becomes the reference for every downstream FC
> adjacency layout (`row_per_phase`, mosaics).

## What this figure shows

A square `imshow` of the symmetric (N × N) FC matrix, diagonal
zeroed for visualization. Colorbar on the right via
`imshow_colorbar_caxdivider`. **No suptitle, no watermark by
default** — patient / band / phase / fc_method live in the file
name. The watermark is opt-in via `watermark=True`.

## Inputs

| arg | meaning | default in template |
|---|---|---|
| `patient` | one of `PATIENTS_4PHASE` | `Pat_02` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | one of `corr / msc / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `tick_labels` | `generic` / `index` / `chnames` | `generic` |
| `log_scale` | logarithmic colour norm (`LogNorm`) | `False` |
| `watermark` | opt-in grey provenance footer | `False` |
| `out_dir` | output directory | `data/outputs/figures/fc_templates/single_adjacency/` |

### `tick_labels` modes

| mode | what you get |
|---|---|
| `generic` | No ticks, no axis labels. Pure matrix + colorbar. The default — clean for publication and downstream stacking. |
| `index` | Matplotlib's auto numeric ticks. Useful when you want to reference particular contact indices in a caption. |
| `chnames` | One major tick per probe family at the family midpoint, labelled with the probe prefix (`A`, `B`, `C`, …). Light minor ticks at probe boundaries. Requires `channel_labels` (parsed from `data/raw/.../channel_labels.csv`). Best for clinically anchored slides. |

## Style choices (current iteration)

- **Figure size:** `figsize=(3.4, 3.0)` square-ish, accommodates a
  vertical colorbar without squeezing the panel.
- **Colormap:** `"magma"`, `vmin=0`, `vmax=np.nanmax(M)`. Diagonal
  forced to 0 before plotting.
- **Aspect:** `aspect="equal"`. Adjacency is index-symmetric.
- **Interpolation:** `interpolation="nearest"` — FC cells are
  discrete `(i, j)` entries; the matplotlib default
  (`"antialiased"`) blurs probe-block boundaries.
- **Ticks:** controlled by `tick_labels` (see table above).
  Default `generic` hides ticks AND axis labels entirely (pure
  matrix + colorbar).  `index` adds numeric ticks with `$i$`/`$j$`
  axis labels.  `chnames` puts probe labels at family midpoints
  with light minor ticks at boundaries and ``"probe"`` axis labels.
- **Axis labels:** **never** the word "contact" / "channel" — the
  matrix is `M_{ij}` mathematically.
- **Colorbar:** `imshow_colorbar_caxdivider(im, ax,
  size="4%", pad=0.05)`; label =
  `r"$\langle |\mathrm{ImCoh}|\rangle_f$"` for `imcoh_abs`,
  `"MSC"` for MSC, `r"$|\mathrm{corr}|$"` for `corr`.
  Label fontsize 8.
- **No suptitle.** Caption in sidecar `.md` only when explicitly asked.
- **Watermark / provenance footer:** **off by default**
  (`watermark=False`). Pass `watermark=True` only when you really
  want a tiny grey provenance string in the bottom-right; the file
  name already carries patient / band / phase / fc_method.
- **Output:** PDF only, full vector. **Do NOT call
  `im.set_rasterized(True)` — FC adjacency pixels stay vector.**
- **File name:**
  `<Pat_NN>_<band>_<phase>_<fc_method>_<tick_labels>.pdf`
  (the parent directory `single_adjacency/` already labels the
  template; do not repeat that prefix in the file name).

## Anti-patterns (do not do these)

- ❌ `fig.suptitle(...)` — context belongs in the file name (and an
  optional sidecar `.md`).
- ❌ `fig.colorbar(...)` directly — always
  `imshow_colorbar_caxdivider`.
- ❌ `im.set_rasterized(True)` — FC pixels are vector; rasterising
  loses the crispness for no real file-size win.
- ❌ Watermark / provenance footer on by default — it duplicates
  the file name. Use `watermark=True` only when the figure will be
  detached from its file name (slide deck, screenshot).
- ❌ Word labels like "contact" / "channel" — use math `$i$`, `$j$`.
- ❌ Saving PNG alongside the PDF.
- ❌ Hardcoding `data/...` paths — use
  `lrg_eegfc.config.paths.FIGURES_ROOT`.
- ❌ Recomputing the FC matrix — load via `load_fc_matrix(...)`.

## How to call

The library function is the canonical entry point — keep the script
thin.

```python
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix

M = load_fc_matrix(patient="Pat_02", phase="rest_pre",
                   band="beta", fc_method="imcoh_abs")
fig, ax, im = plot_fc_adjacency(
    M, fc_method="imcoh_abs", tick_labels="generic",
)
fig.savefig("…/single_adjacency.pdf", bbox_inches="tight")
```

For the script (defaults loaded from the canonical era):

```bash
conda activate lapbrain
python .agents/guides/05_plotting/fc_templates/single_adjacency.py \
    --patient Pat_02 --band beta --phase rest_pre \
    --fc-method imcoh_abs --tick-labels generic
# tick-labels ∈ {generic, index, chnames}
# add --watermark to opt in to the bottom-right grey provenance string
```

## Open questions for next iteration

These are the dials we may still tune with the user:

1. **Tick policy** — fully hidden vs. probe-boundary minor ticks.
2. **Colorbar label units** — bare `|ImCoh|`, `<|ImCoh|>_f`, or
   "FC strength"?
3. **Diagonal handling** — zero (current) vs. NaN-mask-with-light-grey.
4. **vmax convention for single panels** — `np.nanmax(M)`,
   the 99th percentile, or a band-pinned global max?
5. **Probe-group rectangles overlay** (light frames around contacts
   sharing the same probe), opt-in?
