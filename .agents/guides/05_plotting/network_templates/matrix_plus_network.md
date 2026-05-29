---
name: network_templates / matrix_plus_network
type: figure_template
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-26
updated: 2026-05-26
pointers:
  - matrix_plus_network.py
  - circular_dendrogram_network.md
  - chord_highlighted_edges.md
  - single_network.md
  - README.md
  - src/lrg_eegfc/visuals/network_templates.py
  - src/lrg_eegfc/visuals/fc_templates.py
  - data/outputs/figures/presentation_figures/Pat_02/fig1_msc_alpha.pdf
---

# `matrix_plus_network` — FC matrix + cmap-linked network drawing

A row of **(matrix, network)** where the network's edges are coloured
by the SAME ``cmap + norm`` as the matrix imshow.  Reader's eye reads
both panels in the same colour language: dark cells map to dark edges,
bright cells to bright edges.

Crystallised from the MSC-era figure
``data/outputs/figures/presentation_figures/Pat_02/fig1_msc_alpha.pdf``
(two-row rsPre / rsPost with magma matrices + spring networks), but
the concept is generic — works for any magnitude FC (``imcoh_abs``,
``imcoh_sq``, ``msc``) and adapts to signed FC by switching to the
diverging family.

## Three independent layers (modular)

1.  **Matrix** — :func:`plot_fc_adjacency`.  Tick mode
    (``generic`` / ``index`` / ``chnames``), optional ``LogNorm``,
    canonical colorbar label per ``fc_method`` / ``band``.
2.  **Layout** — any name in :data:`LAYOUT_REGISTRY`.  Default
    ``spring`` with the legacy ``msc_era`` preset (``k=0.1``,
    ``iterations=100``).  Computed ONCE on the cross-phase average FC
    when ``shared_layout=True`` (default) so node positions are
    identical across phase rows — direct side-by-side comparison.
3.  **Edge style** — :func:`draw_gamma_edges` with γ-power on
    width/alpha and a ``coloring`` mode.  Default ``cmap`` paints
    each edge as ``cmap(norm(|w|))``, locked to the matrix.

## The three knobs that make structure visible

The legacy MSC-presentation figure had visible chains and clusters;
naïve defaults (single_network's ``γ=6``, fully-converged spring with
``k=k_base/√N``, per-row layouts) produce a ball even on the same
matrix.  Three knobs flip the picture:

- **``γ=2`` (not 6)** — ``γ=6`` is the right choice for ``single_network``
  where the goal is "show me the strongest 1% on a faint backdrop"; here
  we want a mid-tail visible.  ``γ=2`` lets the top quartile of edges
  read crisply while the bottom half fades to invisibility.
- **``width_range=(0, 4)``, ``alpha_range=(0, 1)``** — anchored at
  zero, not at the small ``width_min=0.15`` / ``alpha_min=0.03`` of
  ``single_network``.  Without a floor, the 6786 weakest edges
  genuinely disappear instead of forming a uniform grey wash that
  reads as "ball".
- **``shared_layout`` from the cross-phase mean** — spring layout
  computed once on ``mean(A_per_phase)``; both rows use those
  positions.  Without this, each row's spring relaxation lands in a
  different local minimum and side-by-side comparison is meaningless.

The legacy MSC-era spring preset (``k=0.1``, ``iter=100``) is the
fourth knob; small ``k`` lets the heavy-tailed strong edges form
tight clusters, and ``iter=100`` stops before relaxation equilibrates
everything to a ball.  Pair only with heavy-tailed FC (MSC).  For
the more uniform ImCoh distribution the same recipe is less dramatic
— strong-edge clustering is genuinely weaker in that data.

## Phases — single or multi-row

``--phases`` is ``nargs='+'``.  One phase → one row.  Multiple phases
→ stacked rows (matches the reference figure with
``rest_pre rest_post``).  ``shared_scale=True`` (default) uses a
cohort-wide ``vmax`` across the selected phases so the panels are
directly comparable; ``--no-shared-scale`` gives each phase its own
norm.

## Inputs

| arg | meaning | default |
|---|---|---|
| `patient` | one of `PATIENTS_*` | `Pat_02` |
| `band` | one of `BRAIN_BANDS` | `alpha` |
| `phases` | one or more phases (`nargs='+'`) | `rest_pre rest_post` |
| `fc_method` | `corr / msc / imcoh / imcoh_abs / imcoh_sq` | `imcoh_abs` |
| `cmap` | colormap shared by matrix + edges | `magma` (mag.) / `RdBu_r` (signed) |
| `vmin` / `vmax` / `log_scale` | imshow controls | `0 / None / False` |
| `tick_labels` | `generic / index / chnames` | `index` |
| `no_shared_scale` | per-phase vmax | `False` |
| `layout` | `spring / kk / spectral / laplacian_pca / lrg_kk / lrg_sfdp / arf / …` | `spring` |
| `spring_k_scale` | `kbase_sqrtN / inv_sqrtN / fixed` (auto-k rule) | `kbase_sqrtN` |
| `spring_k_base` | k_base in ``k = k_base/√N`` | `5.0` |
| `spring_k_fixed` | k when ``--spring-k-scale=fixed`` | `0.5` |
| `k` | raw spring k (overrides ``--spring-k-scale``) | `None` |
| `iterations` | spring solver iterations | `600` |
| `n_communities` | LRG cut for LRG-seeded layouts / community node-color | `10` |
| `coloring` | `cmap / probe / shaft / signed / weight` | `cmap` |
| `node_color` | `shaft / community / uniform` | `shaft` |
| `gamma` | γ in `t = (\|w\|/\|w\|_max)**γ` | `6.0` |
| `width_min` / `width_max` | edge width range | `0.15 / 4.0` |
| `alpha_min` / `alpha_max` | edge alpha range | `0.03 / 0.9` |
| `min_alpha` | drop edges below this alpha for file weight | `0.0` |
| `node_size` | scatter node size | `18` |
| `out_dir` | output directory | `data/outputs/figures/network_templates/matrix_plus_network/` |

## Style choices

- **Default `cmap="magma"` for magnitude FC** matches the reference
  figure and reads on white backgrounds.  ``RdBu_r`` is auto-selected
  for signed FC (``corr``, ``imcoh``).
- **`coloring="cmap"`** is the figure's defining choice.  Other modes
  (``probe``, ``shaft``, ``signed``, ``weight``) are available for
  comparison but break the matrix↔network colour link.
- **Spring `k = k_base/√N`** (canonical FR scaling) reads gracefully
  across N=80–130.  NetworkX's documented default ``1/√N`` (``--spring-k-scale
  inv_sqrtN``) tends to hairball for dense FC.  Sweep with
  ``--spring-k-scale fixed --spring-k-fixed 0.05`` for diagnostics.
- **`iterations=600`** is high enough for convergence in our dense
  N≈100 regime.  Drop to 200 for faster sweeps; raise to 2000+ only
  if you suspect non-convergence.
- **`tick_labels="index"`** matches the reference figure; switch to
  ``chnames`` for probe-family midpoints when the audience needs to
  read individual shafts.
- **Shared scale** is on by default so the two phases share a colorbar
  range — anything else makes the reader miscompare brightnesses.
  Turn off with ``--no-shared-scale`` only for diagnostic per-phase
  saturation.
- **No `fig.suptitle`.**  Per-row ax-level titles carry patient / band /
  phase / fc_method on the matrix and layout + coloring on the
  network; the file name carries everything else.

## How to call

```bash
# Reproduce the reference figure (Pat_02 alpha rsPre/rsPost MSC)
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_02 --band alpha --phases rest_pre rest_post \
    --fc-method msc

# Single panel, current era (imcoh_abs)
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre

# Switch to Kamada-Kawai layout, keep cmap-linked edges
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre --layout kk

# Diagnostic: probe-coloured edges (same-probe red / cross-probe gray)
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre \
    --coloring probe --node-color shaft

# LogNorm so weak edges read above the floor
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre --log-scale

# LRG-seeded SFDP layout (graph-tool) + community node colouring
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre \
    --layout lrg_sfdp --node-color community

# Spring k sweep — fixed k for visual comparison
python .agents/guides/05_plotting/network_templates/matrix_plus_network.py \
    --patient Pat_05 --band beta --phases rest_pre \
    --spring-k-scale fixed --spring-k-fixed 0.05
```

## Layout catalog

| Name | Backend | Determinism | When to use |
|---|---|---|---|
| `spring` | nx | stochastic (seed=42) | default; canonical `k=k_base/√N` |
| `kk` | nx | deterministic | always-reads; slow at N≈115 |
| `spectral` | nx | deterministic | reveals bipartition; flat for nearly fully-connected |
| `laplacian_pca` | nx | deterministic | richer than spectral when no clear bipartition |
| `community_grouped` | nx | stochastic (sub-spring) | Louvain modules; sub-spring unstable |
| `backbone_guided` | nx | stochastic | layout on disparity-filter backbone, draw full graph |
| `circular_by_shaft` | nx | deterministic | sanity-check baseline |
| `lrg_kk` | nx | deterministic given dendrogram | LRG modules forced apart via distance inflation |
| `lrg_sfdp` | gt | deterministic | LRG modules forced apart via `groups=`; scales best |
| `sfdp` | gt | deterministic | dense without module hints |
| `arf` | gt | deterministic post-init | alternative to SFDP for very dense graphs |

Pick by your figure's purpose — see `README.md` §3 for a longer
discussion.

## Edge-coloring catalog

| Mode | Source | Use |
|---|---|---|
| `cmap` | `cmap(norm(\|w\|))` shared with matrix | **default**; matrix↔network colour link |
| `probe` | same-probe red, cross-probe gray | probe-bias diagnostic |
| `shaft` | same-probe painted in shaft's tab20, cross-probe gray | shaft-aware backdrop |
| `signed` | `RdBu_r` on signed weight (TwoSlopeNorm) | signed FC (`corr`, `imcoh`) |
| `weight` | sequential gray on `\|w\|` | cohort comparisons where probe geometry differs |

When ``coloring`` is not ``cmap``, the matrix's colorbar still uses
the same cmap, but the network ink is independent — the figure
becomes "matrix + diagnostic network" rather than
"matrix-linked-to-network".

## Anti-patterns

- ❌ Setting ``--cmap`` to a diverging palette (e.g. ``RdBu_r``) for
  magnitude FC — the ``vmin=0`` lower bound means half the colorbar
  is unused.  Stick to sequential cmaps (``magma``, ``viridis``,
  ``inferno``, ``plasma``) for magnitude.
- ❌ ``--no-shared-scale`` on a 2-phase comparison figure — the two
  rows then live in different colour scales and side-by-side
  brightness comparisons silently lie.
- ❌ ``--coloring cmap`` with ``--log-scale`` and tiny ``vmin``
  (default 0) — the LogNorm fails for non-positive entries.  Either
  pass ``--vmin 1e-3`` or let the helper auto-pick from the minimum
  positive entry.
- ❌ Plain ``--layout spring`` with very high or very low ``--k``
  values on a dense ImCoh graph — produces a hairball or a starburst,
  neither informative.  Trust the ``kbase_sqrtN`` auto-tuner unless
  you have a specific diagnostic reason.
- ❌ ``--layout lrg_sfdp`` (or any LRG-seeded layout) without a
  cached LRG result — the helper raises ``FileNotFoundError``.
  Compute LRG for the (patient, anchor_phase, band, fc_method) first.
- ❌ Forgetting ``--fc-method`` for the MSC-era reproduction — the
  reference figure used ``msc``; today's pipeline defaults to
  ``imcoh_abs`` which has a different colour distribution.
- ❌ Locking ``--spring-k-fixed`` for a figure that ships in the
  paper — magic constants don't survive a patient swap.  Use the
  auto-tuner for any non-diagnostic figure.

## Library entry points

```python
from lrg_eegfc.visuals.network_templates import (
    plot_fc_matrix_and_network,         # one row (matrix, network)
    plot_fc_matrix_and_network_rows,    # row-per-phase mosaic
    spring_auto_k,                      # canonical k chooser
    draw_gamma_edges,                   # underlying edge renderer
    LAYOUT_REGISTRY,                    # layout names
)
```

See also: ``single_network.md`` for the network-only template and
``circular_dendrogram_network.md`` / ``chord_highlighted_edges.md`` for
the curvy-chord family.
