---
name: 05_plotting / output-and-rasterization
type: style_sheet
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Output, rasterisation, watermarks, file naming

## Format — PDF only, full vector

The default output is **vector PDF**, fully vector for every artist
(including `imshow` pixels). PNG siblings are not written. This is
already a project rule (`feedback_no_png_duplicates.md`,
`04_rules/never-always-list.md`).

```python
out_path = output_dir / f"{patient}_{band}_{phase}_pattern.pdf"
fig.savefig(out_path, bbox_inches="tight")
plt.close(fig)
```

**Vector is sharper at any zoom, the file size for typical FC
matrices (`N ≤ 130`) is small (tens of KB), and rasterising loses
crispness for no real win.**

## Rasterisation — do not rasterise

This project's rule is **never rasterise**. Earlier guides told
agents to set `im.set_rasterized(True)` for "heavy" artists; that
rule is **withdrawn**. Concretely:

- **FC adjacency matrices, distance matrices, audit heatmaps** —
  vector. Do **not** call `im.set_rasterized(True)`.
- **Scatter plots, swarms, line plots** — vector.
- **`pcolormesh`** — vector.
- **Dendrograms, network drawings** — vector.

If a specific figure ever produces a genuinely unmanageable PDF
(say a `4096 × 4096` correlation cube), bring it up explicitly and
we will decide together. Until then, the rule is: **never call
`set_rasterized(True)`**.

## Watermarks / provenance footer — opt-in only

A small grey provenance string in the bottom-right corner is
sometimes useful (slide deck, screenshot detached from filename) but
it is **off by default**. Every figure helper in this project that
emits a footer must:

1. Default to `watermark=False`.
2. Only write the footer when explicitly requested via
   `watermark=True` (or its equivalent kwarg name).

Rationale: the file name already carries
`<patient> · <band> · <phase> · <fc_method>`; the footer duplicates
that information and clutters the figure.

If you really want one (e.g. for a slide), opt in:

```python
plot_single_adjacency(M, ..., watermark=True)
```

Or call the helper directly:

```python
from lrg_eegfc.visuals.layout import add_provenance_footer
add_provenance_footer(fig, "single_adjacency · Pat_02 · β · rest_pre · imcoh_abs")
```

## File naming

- Per-patient outputs: `<Patient>/<band>_<phase>_<kind>.pdf`
- Cohort outputs: `cohort_<kind>.pdf`
- Audit outputs: `<audit_NN>_<kind>.pdf` or
  `<datafolder>/<kind>.pdf` where datafolder is the
  audit's namespaced folder.

The file name is the canonical place for provenance metadata
(patient, band, phase, fc_method, template). Make it descriptive
enough that no on-figure watermark is needed.

## Caption .md alongside a figure — only when requested

Caption sidecars (`<figure-name>.md` next to the PDF) are NOT a
default — write one only when the user asks for it, or when a folder
explicitly accumulates them (e.g. `data/audit/raw_fc_phase_distance/`
keeps captions as part of its handoff format).

When asked to write a caption, see [`captions.md`](captions.md) for
the short style guide (plain language, no jargon, sentence-first).

## Always close figures

```python
plt.close(fig)
```

Memory leaks build up fast in batch loops; always close. Project rule
in `CLAUDE.md` §"Technical invariants" #6.

## Always create the parent dir

```python
output_path.parent.mkdir(parents=True, exist_ok=True)
```

Project rule in `CLAUDE.md` §"Technical invariants" #7.
