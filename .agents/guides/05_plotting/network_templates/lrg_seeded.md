---
name: network_templates / lrg_seeded
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - lrg_seeded.py
  - single_network.md
  - mosaic_layout_compare.md
  - README.md
  - .agents/guides/02_methods/lrg-framework-guide.md
---

# `lrg_seeded` — LRG-driven network drawing (the cluster figure)

A single-axis network drawing where the layout is *forced* to separate
LRG modules.  Two backends share the same visual contract:

- **`gt` (default)** — graph-tool SFDP with `groups=` set from the
  LRG dendrogram cut at `n_communities`.  Same-group vertices get an
  extra attractive force; `sfdp_gamma` controls module-separation
  strength.
- **`nx`** — NetworkX Kamada-Kawai with the distance matrix inflated
  by `separation` on inter-community edges.  Slower but no graph-tool
  dependency required at runtime if the user is on a NetworkX-only env.

## What this figure shows

LRG modules pulled apart in the plane.  Edges drawn at full density
with the canonical γ-power recipe.  Default node coloring is
**community** (one tab20 color per LRG cluster); switch to
`--node-color shaft` to confirm modules aren't trivially probe-driven.

If a putative module collapses back into one shaft when re-coloured
by shaft, that module is sEEG geometry, not connectivity.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_05` |
| `band` | one of `BRAIN_BANDS` | `beta` |
| `phase` | one of `PHASE_LABELS` | `rest_pre` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `n_communities` | LRG dendrogram cut | `10` |
| `backend` | `gt / nx` | **`gt`** |
| `separation` | NetworkX KK inter-comm distance multiplier | `5.0` |
| `sfdp_gamma` | graph-tool SFDP module-separation γ | `0.1` |
| `node_color` | `community / shaft / uniform` | **`community`** |
| `coloring` | `probe / signed / weight` (auto-default) | `auto` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/lrg_seeded/` |

## Style choices

- **LRG cut as a fixed prior.**  The same `n_communities` is used
  whether `backend=gt` (passed as `groups=`) or `backend=nx` (used
  to construct the inter-community distance mask).  Changing
  `n_communities` *will* change the layout — that's the point.
- **Default `n_communities=10`** matches the cohort-LRG convention
  used elsewhere (e.g. dendrogram color thresholds).  Sweep
  `--n-communities {5, 10, 15, 25}` to inspect multi-scale module
  structure.
- **Edge recipe identical to `single_network`.**  γ-power on
  width/alpha; same-probe red, cross-probe gray (or RdBu_r for signed
  FC).  Module separation comes from the layout, not the edge style.
- **Node size 24** (slightly larger than `single_network`'s 22) so the
  cluster colours dominate visually.
- **Figure size** 6 × 6.  Single panel.
- **No `fig.suptitle`.**  Provenance lives in the file name.
- **File name:**
  `<Pat_NN>_<band>_<phase>_<fc_method>_n<N>_<backend>_<node_color>_<coloring>.pdf`.

## How to call

```bash
# Default — graph-tool SFDP, community-coloured, 10 LRG modules
python .agents/guides/05_plotting/network_templates/lrg_seeded.py \
    --patient Pat_05 --band beta --phase rest_post

# NetworkX KK fallback with stronger module separation
python .agents/guides/05_plotting/network_templates/lrg_seeded.py \
    --patient Pat_05 --band beta --phase rest_post \
    --backend nx --separation 10.0

# Sanity check: re-colour the same layout by shaft to expose probe bias
python .agents/guides/05_plotting/network_templates/lrg_seeded.py \
    --patient Pat_05 --band beta --phase rest_post \
    --node-color shaft

# Multi-scale sweep (manual; one figure per cut):
for n in 5 10 15 25; do
    python .agents/guides/05_plotting/network_templates/lrg_seeded.py \
        --patient Pat_05 --band beta --phase rest_post --n-communities $n
done
```

## Anti-patterns

- ❌ Comparing `lrg_seeded` panels across phases without re-checking
  the LRG cut — community labels change between phases by construction
  (the dendrograms are different), so panels are not directly
  comparable.  For phase comparison use `row_per_phase` with a
  layout-driven (not LRG-driven) backend.
- ❌ Drawing this with `node_color="uniform"` — defeats the entire
  purpose; use `community` or `shaft`.
- ❌ Setting `n_communities` higher than the natural number of stable
  modules in the LRG dendrogram (visible as merge plateaus).  At
  `n=30+` the SFDP layout pulls every contact into its own
  micro-module and the figure becomes incoherent.
- ❌ Using SFDP `gamma > 1.0` — the layout collapses into disjoint
  clusters with no visible inter-module edges; readers can't tell
  cross-module strength.
