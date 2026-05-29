---
name: 05_plotting / README
type: guide_index
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/04_rules/coding-rules.md
  - .agents/guides/04_rules/never-always-list.md
  - lrgsglib/src/lrgsglib/plotlib/
  - src/lrg_eegfc/visuals/
  - src/lrg_eegfc/cli/plot.py
---

# Plotting guide — preferences, helpers, and the "read this first" entry

**Read this folder before producing any figure.** It is the single
source of truth for plot style, layout, colorbars, legends, fonts, and
file output in the `lrgeegfc` repo.

## TL;DR — eight rules

1. **Shared legends → figure-level, not axis-level.** Use
   `fig.legend(loc="lower center", bbox_to_anchor=(0.5, -0.02),
   ncol=…)` (or upper center) for any legend that applies to more than
   one axis. Never let a shared legend live inside one panel — it
   breaks panel symmetry and creates dead white-space. See
   [`legends.md`](legends.md).
2. **Colorbars → `imshow_colorbar_caxdivider`.** This helper is the
   default for any *single-imshow* axis. It guarantees properly spaced
   colorbars that don't squeeze the parent axis. See
   [`colorbars.md`](colorbars.md). **Scope:** the helper attaches the
   colorbar to one axis via `make_axes_locatable`; it **cannot** serve
   a colorbar shared across multiple columns in the same row (it will
   misalign or resize the wrong axis). For row-shared cbars keep the
   explicit `make_axes_locatable` / `fig.add_axes([...])` pattern.
   Locked 2026-05-28.
3. **Multi-axis layout → figure-level decoration.** Titles, legends,
   colorbars, and shared axis labels go on the *figure*, not on a
   single axis. Use `sharex` / `sharey` and per-row / per-column
   colorbars when scales differ. See
   [`multi-axis-figures.md`](multi-axis-figures.md).
4. **Library-first.** Before writing a custom plot helper, check
   `lrgsglib.plotlib` and `lrg_eegfc.visuals.*`. If a helper is used
   ≥ 2 times, promote it. See [`library-helpers.md`](library-helpers.md).
5. **PDF only, full vector — never rasterise.** No PNG siblings.
   Do **not** call `im.set_rasterized(True)` on FC matrices, audit
   panels, scatter plots, or anything else. Vector is sharper at any
   zoom and the file size for typical FC matrices is small (tens of
   KB). The previous "rasterise heavy artists" rule is withdrawn.
   See [`output-and-rasterization.md`](output-and-rasterization.md).
6. **No `fig.suptitle` on publication figures.** Context goes in the
   file name (and optionally a companion `.md` caption), not on the
   figure. (Already in `04_rules/never-always-list.md`.)
7. **No watermark / provenance footer by default.** The grey
   bottom-right string is **opt-in only** — pass `watermark=True`
   (or call `add_provenance_footer(fig, ...)`) when the figure will
   be detached from its file name. The file name itself is the
   provenance.
8. **Activate the project mplstyle.** Call
   `from lrg_eegfc.visuals.styles import use_lrg_style; use_lrg_style()`
   at the top of every figure script (after imports, before the first
   `subplots`). It pins font sizes, tick widths, and `pdf.fonttype=42`
   (TrueType-embedded PDFs). Single source:
   [`src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle`](../../../src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle).
   Do not redefine these defaults inline; do not edit the mplstyle for
   one-off needs (use `with rc_context({...}):` instead). For LogNorm
   colorbars, end with `_apply_factored_sci_format(clb, ...)` — see
   [`colormaps-and-styles.md`](colormaps-and-styles.md#lognorm-colorbar-tick-policy).

## File map

- [`legends.md`](legends.md) — figure-level legend placement and
  rationale, with copy-paste recipes.
- [`colorbars.md`](colorbars.md) — `imshow_colorbar_caxdivider` usage,
  per-row colorbars, shared colour scales, divergent vs sequential.
- [`multi-axis-figures.md`](multi-axis-figures.md) — grid layouts,
  shared axes, figure-level titles, axis-label sharing patterns.
- [`colormaps-and-styles.md`](colormaps-and-styles.md) — colormap
  preferences, qualitative palettes, line-style conventions, font
  sizes per figure size.
- [`library-helpers.md`](library-helpers.md) — index of reusable plot
  helpers in `lrgsglib.plotlib` and `lrg_eegfc.visuals` with one-liner
  descriptions and where to import each.
- [`output-and-rasterization.md`](output-and-rasterization.md) —
  PDF-only policy, **never rasterise**, watermark opt-in, file-naming.
- [`captions.md`](captions.md) — how to write a caption sidecar
  *when asked* (plain language, three blocks).
- [`per-figure-style-sheets/`](per-figure-style-sheets/) — one
  short `.md` per figure family (heatmaps, swarms, networks, …) with
  the canonical style for that family.
- [`fc_templates/`](fc_templates/README.md) — **canonical FC
  adjacency-matrix templates** (single, row-per-phase, mosaics).
  Read this before plotting any `imcoh_abs`, MSC, or `corr` matrix.
  Each template ships a `.md` style sheet + `.py` script that
  produces a vector PDF obeying every rule in this folder.

## "Read this before plotting" entry point

Before you write any figure code, run:

```bash
cat .agents/guides/05_plotting/README.md
```

Or invoke the matching skill / hook (see
[`hook-and-skill.md`](hook-and-skill.md)).

A `pre_tool_use` hook in `.claude/settings.json` triggers a reminder
on Write/Edit when the touched file path matches
`scripts/.*\.py|src/lrg_eegfc/visuals/.*\.py|figures.*` — agents are
prompted to consult this guide before producing or editing a figure.

## Per-figure-class templates

For recurring figure classes we keep ready-to-run, opinionated
templates in subfolders. **Always check for a matching template
before writing new code from scratch.**

| If the request is about… | Go to |
|---|---|
| FC / adjacency / connectivity matrix / `imcoh_abs` / MSC / correlation matrix | [`fc_templates/`](fc_templates/README.md) |
| Network drawing / graph layout / spring/KK/SFDP / connectome / LRG-seeded | [`network_templates/`](network_templates/README.md) |
| (more classes coming: LRG dendrograms, distance heatmaps, swarms…) | per-figure-style-sheets/ for now |

## Library files agents should know about

| file | role |
|---|---|
| `lrgsglib/src/lrgsglib/plotlib/colorbars.py` | **canonical** `imshow_colorbar_caxdivider` |
| `lrgsglib/src/lrgsglib/plotlib/const_plotlib.py` | central matplotlib imports + colormap constants |
| `lrgsglib/src/lrgsglib/plotlib/colormaps.py` | reusable colormaps (incl. project-specific) |
| `lrgsglib/src/lrgsglib/plotlib/color.py` | colour utilities |
| `lrgsglib/src/lrgsglib/plotlib/lrg.py` | LRG-specific plot helpers |
| `lrgsglib/src/lrgsglib/plotlib/mathplot.py` | math-style plot helpers (log distributions, etc.) |
| `src/lrg_eegfc/visuals/correlation.py` | FC heatmap + network panel |
| `src/lrg_eegfc/visuals/lrg.py` | LRG dendrogram + entropy |
| `src/lrg_eegfc/visuals/spatial.py` | nilearn glass-brain / 3D plots |
| `src/lrg_eegfc/visuals/cross_condition.py` | cohort-level cross-condition figures |
| `src/lrg_eegfc/visuals/plotting.py` | low-level CLI plot utilities |
| `src/lrg_eegfc/cli/plot.py` | CLI front-end (`lrg-eegfc plot …`) |

## Spirit

- Don't reinvent. If a helper exists, import it.
- Don't paste matplotlib boilerplate. Wrap repeated patterns into a
  function and put it in `visuals/` or `plotlib/`.
- Don't decorate a single axis with information that belongs on the
  figure. Symmetry is more important than convenience.
- Captions in companion `.md` files, never `suptitle`.
