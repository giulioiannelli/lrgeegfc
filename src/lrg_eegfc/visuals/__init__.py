"""Visualization functions for LRG EEG FC analysis."""

# Canonical colorbar helper (re-export from lrgsglib so callers have one
# source). Use for any single-imshow axis. Cannot serve a colorbar shared
# across multiple columns in the same row — for that case keep the
# explicit `make_axes_locatable` / `fig.add_axes([...])` pattern (rule
# locked 2026-05-28; see CLAUDE.md plotting rule 2 and the
# `feedback_imshow_colorbar_caxdivider_scope` memory).
from lrgsglib.plotlib.colorbars import imshow_colorbar_caxdivider

# Canonical band colour palette (spectrum: slow=red -> fast=blue). Figures must
# call `band_color(band)` instead of hardcoding per-band colours, so the whole
# codebase restyles from one constant. See `visuals.styles`.
from .styles import (
    use_lrg_style,
    band_color,
    band_colors_list,
    band_cmap,
    BAND_COLORS,
)

from .correlation import (
    plot_correlation_heatmap,
    plot_correlation_and_network,
    plot_marchenko_pastur_comparison,
    plot_percolation_curves,
)
from .msc import (
    plot_msc_heatmap,
    plot_msc_and_network,
    plot_msc_comparison_dense_vs_validated,
)
from .lrg import (
    plot_lrg_entropy_curves,
    plot_lrg_dendrogram,
    plot_ultrametric_heatmap,
)
from .lrg_panels import plot_lrg_full_panel
from .spatial import (
    load_spatial_metadata,
    prepare_spatial_coordinates,
    build_edge_traces,
    plot_spatial_network_3d,
    plot_spatial_network_3d_mpl,
    plot_spatial_clusters_comparison,
    view_brain_connectome,
    plot_brain_connectome,
)
from .fc_templates import (
    plot_fc_adjacency,
    plot_fc_adjacency_row,
    plot_fc_adjacency_grid,
    fc_method_colorbar_label,
    probe_groups,
    draw_probe_outlines,
)
from .network_layouts import (
    compute_network_layout,
    compute_percolation_threshold,
)
from .network_drawing import (
    EDGE_GAMMA,
    scale_edge_weights,
    draw_network_edges,
    render_sbm_panel,
    render_lrg_panel,
)
from .network_templates import (
    plot_fc_network,
    plot_fc_network_row,
    plot_fc_network_grid,
    plot_fc_network_lrg,
    plot_layout_gallery,
    compute_layout,
    LAYOUT_REGISTRY,
    DEFAULT_GALLERY,
    draw_gamma_edges,
    draw_nodes,
    shaft_colors,
    community_colors,
    load_probe_labels,
    nx_to_gt,
    matrix_to_gt,
    SAME_PROBE_RGB,
    CROSS_PROBE_RGB,
)

__all__ = [
    # Correlation visualizations
    "plot_correlation_heatmap",
    "plot_correlation_and_network",
    "plot_marchenko_pastur_comparison",
    "plot_percolation_curves",
    # MSC visualizations
    "plot_msc_heatmap",
    "plot_msc_and_network",
    "plot_msc_comparison_dense_vs_validated",
    # LRG visualizations
    "plot_lrg_entropy_curves",
    "plot_lrg_dendrogram",
    "plot_ultrametric_heatmap",
    "plot_lrg_full_panel",
    # Spatial 3D visualizations
    "load_spatial_metadata",
    "prepare_spatial_coordinates",
    "build_edge_traces",
    "plot_spatial_network_3d",
    "plot_spatial_network_3d_mpl",
    "plot_spatial_clusters_comparison",
    # Brain connectome visualizations (nilearn)
    "view_brain_connectome",
    "plot_brain_connectome",
    # FC adjacency-matrix templates
    "plot_fc_adjacency",
    "plot_fc_adjacency_row",
    "plot_fc_adjacency_grid",
    "fc_method_colorbar_label",
    "probe_groups",
    "draw_probe_outlines",
    # Network-drawing templates
    "plot_fc_network",
    "plot_fc_network_row",
    "plot_fc_network_grid",
    "plot_fc_network_lrg",
    "plot_layout_gallery",
    "compute_layout",
    "LAYOUT_REGISTRY",
    "DEFAULT_GALLERY",
    "draw_gamma_edges",
    "draw_nodes",
    "shaft_colors",
    "community_colors",
    "load_probe_labels",
    "nx_to_gt",
    "matrix_to_gt",
    "SAME_PROBE_RGB",
    "CROSS_PROBE_RGB",
    # Network layouts (promoted from figures_for_notes/_shared.py 2026-05-29)
    "compute_network_layout",
    "compute_percolation_threshold",
    # Network drawing (promoted from figures_for_notes/_shared.py 2026-05-29)
    "EDGE_GAMMA",
    "scale_edge_weights",
    "draw_network_edges",
    "render_sbm_panel",
    "render_lrg_panel",
]
