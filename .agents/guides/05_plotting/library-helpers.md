---
name: 05_plotting / library-helpers
type: reference
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
---

# Library helpers — what to import before writing your own

Before writing any plotting code, scan this index. If a helper
already exists, **import it; do not re-implement**.

## Canonical (lrgsglib.plotlib)

```python
from lrgsglib.plotlib import (
    imshow_colorbar_caxdivider,   # default colorbar (see colorbars.md)
    # colormap registry helpers
    # mathplot helpers (log distributions etc.)
)
```

| helper | file | role |
|---|---|---|
| `imshow_colorbar_caxdivider(im, ax, ...)` | `colorbars.py` | properly spaced colorbar; default for any imshow |
| `plot_log_distribution(data, ...)` | `mathplot.py` | log-binned histogram with PDF overlay |
| named colormaps | `colormaps.py` | project-specific colormaps |
| colour utilities (hex / rgb / lighten) | `color.py` | colour manipulation |
| tilings / lattices / chladni / 3d / ax_patches | `tilings.py`, `lattices.py`, `chladni.py`, `plot3d.py`, `ax_patches.py` | shape primitives |

## Project-specific (lrg_eegfc.visuals)

```python
from lrg_eegfc.visuals.correlation import plot_correlation_and_network
from lrg_eegfc.visuals.msc import plot_msc_and_network, plot_msc_summary
from lrg_eegfc.visuals.lrg import (
    plot_lrg_full_panel, plot_lrg_dendrogram,
    plot_lrg_entropy_curves, plot_ultrametric_heatmap,
)
from lrg_eegfc.visuals.spatial import (
    plot_spatial_network_3d,        # plotly
    view_brain_connectome,          # nilearn glass-brain
)
from lrg_eegfc.visuals.reorganization import (
    plot_phase_reorganization,
    plot_reorganization_distance_matrix,
)
from lrg_eegfc.visuals.metastable import (
    compute_clustering_across_tau, create_sankey_diagram,
)
from lrg_eegfc.visuals.compare import plot_fc_comparison
from lrg_eegfc.visuals.cross_condition import generate_cross_condition_figure
from lrg_eegfc.visuals.plotting import (
    plot_correlation_matrix,        # low-level CLI helper
    plot_entropy, plot_dendrogram, plot_graph,
)
```

## Open coding tasks (helpers worth promoting)

These appear ≥ 2 times across `scripts/` or `visuals/` and should be
promoted into a single library helper. Each entry: where it lives now
→ proposed library home.

1. **`figure_legend(fig, handles, where, ...)`**
   - Currently scattered: `cross_condition.py:176/237/398/725`,
     `scripts/01_compute/audit/queries/q_missing_figures.py:104`,
     etc.
   - Proposed: `lrg_eegfc.visuals.layout.figure_legend(...)`.
   - See [`legends.md`](legends.md) for the canonical signature.

2. **Provenance footer `fig.text(0.99, 0.01, ...)`** — opt-in only
   - Currently inlined across `audit_25/26/28` and the queries.
   - Helper exists: `lrg_eegfc.visuals.layout.add_provenance_footer(fig, label)`.
   - **Default is no footer.** Call this helper only when the figure
     will be detached from its file name (slide deck, screenshot).
     Plotting functions that expose a footer must do so through a
     `watermark=False` kwarg.

3. **Common-V* indices `_common_giant_indices`**
   - Defined in `audit_25.py`, `audit_26.py`, `audit_28.py`,
     `q_dP_dS_examples.py` — four copies.
   - Proposed: `lrg_eegfc.utils.fc.common_giant_indices(adj_phases)`.

4. **`_triu_offdiag` (upper-triangle vector excluding diag)**
   - Multiple copies across audits.
   - Proposed: `lrg_eegfc.utils.fc.triu_offdiag(A)` or use
     `A[np.triu_indices_from(A, k=1)]` directly with no helper.

5. **Distance functions `d_P / d_S / d_F`**
   - Three copies across `audit_25/26/28`.
   - Proposed: `lrg_eegfc.utils.metrics.fc_distances` module
     exporting `d_P, d_S, d_F, all_distances`.

6. **`imshow_colorbar_caxdivider` duplicate**
   - Lives in BOTH `lrgsglib.plotlib.colorbars` (canonical, more
     featured) and `lrg_eegfc.visuals.correlation` (duplicate copy).
   - Action: remove the `correlation.py` copy, update its 3 callers
     to import from `lrgsglib.plotlib`. Library-first rule violated
     by the duplicate.

These tasks are listed for the next person who edits the relevant
file. Don't promote helpers preemptively — only when the second use
appears (per `04_rules/coding-rules.md`).

## Import-first checklist (before writing any new helper)

1. `rg "def <name>" src/lrg_eegfc/ lrgsglib/src/lrgsglib/`
2. Skim `lrgsglib.plotlib.__init__` — most plot primitives surface there.
3. Skim `lrg_eegfc.visuals.__init__` for project-specific helpers.
4. If your helper would be the second copy of an existing function,
   promote and import — do not paste.
