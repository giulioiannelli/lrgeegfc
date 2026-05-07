---
name: network_templates / row_per_band
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - row_per_band.py
  - row_per_phase.md
  - single_network.md
  - README.md
---

# `row_per_band` — one patient × one phase × 6 bands

A 1 × 6 row of FC network drawings; columns iterate over bands
(δ, θ, α, β, γ$_l$, γ$_h$); patient, phase, and `fc_method` fixed.

## What this figure shows

Spectral structure of the FC network at a single phase.  Comparing
panels left → right shows how connectivity reorganises across
frequency bins.  γ$_l$/γ$_h$ are typically much weaker than α/β under
`imcoh_abs`, so γ panels appear sparser even with the same γ-power
edge scaling — that's information, not a bug.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `layout` | layout algorithm (LRG-seeded excluded) | **`spring`** |
| `node_color` | `shaft / uniform` | `shaft` |
| `coloring` | `probe / signed / weight` (auto-default) | `auto` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/row_per_band/` |

## Style choices

- **Same layout algorithm across panels**, computed per band on that
  band's FC matrix.  Nodes will move; that movement is the figure's
  signal.
- **Edge γ-power normalised per panel** to its own `|w|_max` — γ$_h$
  panels won't be uniformly black just because they have lower
  absolute amplitudes.  This means edge widths/alphas are not
  comparable across panels in absolute terms; only the *relative
  pattern* within each panel is.
- **Panel size** 3.0 × 3.0; **figure size** 18 × 3.0 (6 panels).
- **Per-panel titles** are band TeX glyphs.  No `fig.suptitle`.
- **File name:**
  `<Pat_NN>_<phase>_<fc_method>_<layout>_<node_color>_<coloring>.pdf`.

## How to call

```bash
python .agents/guides/05_plotting/network_templates/row_per_band.py \
    --patient Pat_05 --phase rest_pre --layout kk
```

## Anti-patterns

- ❌ Sharing a global edge-norm across panels — γ$_h$ goes black,
  α/β saturate.  Per-panel norm is the right default.
- ❌ Using `--layout circular_by_shaft` here — that layout is
  band-independent (only depends on probe ids), so all six panels
  would look identical.
