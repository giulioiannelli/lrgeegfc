---
name: network_templates / row_per_phase
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - row_per_phase.py
  - row_per_band.md
  - single_network.md
  - README.md
---

# `row_per_phase` — one patient × one band × 4 phases

A 1 × 4 row of FC network drawings; columns iterate over phases
(rest$_{\mathrm{pre}}$, task$_{\mathrm{learn}}$, task$_{\mathrm{test}}$,
rest$_{\mathrm{post}}$); patient, band, and `fc_method` fixed.

## What this figure shows

Phase-induced reorganisation of the FC network within a single band.
Comparing panels left → right exposes how connectivity restructures
through learning and back to baseline.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `layout` | layout algorithm (LRG-seeded excluded) | **`spring`** |
| `node_color` | `shaft / uniform` | `shaft` |
| `coloring` | `probe / signed / weight` (auto-default) | `auto` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/row_per_phase/` |

## Style choices

- **Same layout across panels** — within a single (patient, band) the
  layout for each phase is computed independently.  The natural reading
  is shape-change: nodes that move close together task-learn → task-test
  show emergent coupling; nodes that drift apart in rest_post are
  trace-relaxed.
- **Edge recipe / probe coloring** identical to `single_network`.
- **Panel size** 3.5 × 3.5; **figure size** 14 × 3.5 (4 panels).
- **Per-panel titles** are phase TeX glyphs; no `fig.suptitle`.
- **Tight layout** packs panels.  No shared cb (every panel is
  independently drawn).
- **File name:**
  `<Pat_NN>_<band>_<fc_method>_<layout>_<node_color>_<coloring>.pdf`.

## How to call

```bash
python .agents/guides/05_plotting/network_templates/row_per_phase.py \
    --patient Pat_05 --band beta --layout kk
# add --fc-method corr to switch to signed-FC diverging coloring
```

## Anti-patterns

- ❌ Using a different layout per panel — destroys the comparison.
- ❌ Forcing same node positions across phases by computing the
  layout on phase-averaged FC — masks the reorganisation we're trying
  to show.  (For "anchor" figures where stable positions help, use
  `lrg_seeded` with task-test community labels and cross-render
  manually.)
- ❌ LRG-seeded layouts in this template — community labels change
  across phases by construction; pinning them to one phase introduces
  bias.  Use `lrg_seeded` for that purpose; `row_per_phase` is for
  layout-driven (not LRG-driven) phase comparison.
