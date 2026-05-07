"""Visualization functions for LRG EEG FC analysis."""

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
    plot_lrg_full_panel,
)
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
]
