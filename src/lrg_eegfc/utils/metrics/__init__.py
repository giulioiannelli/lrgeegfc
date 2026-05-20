"""Metrics and comparison utilities."""

from .vi import compute_vi, conditional_entropy
from .comparison import (
    compute_phase_distance_matrix,
    compute_cross_patient_consistency,
    rank_distance_measures,
)
from .compare import (
    UltrametricComparison,
    compare_ultrametric_matrices,
    compare_fc_methods,
    compare_phases,
    aggregate_comparisons,
    batch_compare_fc_methods,
    batch_compare_phases,
)
from .reorganization import (
    build_metric_specs,
    compute_metric_matrix,
    compute_cluster_labels,
    cluster_swap_coefficient,
    compute_cluster_swap_matrix,
    compute_ari_matrix,
)
from .tree import (
    dmax_from_Z,
    tree_internal_nodes,
    jaccard_leafsets,
    h_log_grid,
    fcluster_at_h_rel,
)
from .tree_distance import (
    kc_vectors,
    kc_distance,
    matching_cluster_distance,
    weighted_rf_distance,
)
from .functional_tree_distance import (
    FunctionalTreeDistanceResult,
    tau_star_from_C,
    tau_interval,
    delta_curves,
    compute_functional_tree_distance,
)
from .spatial import (
    SUMMARY_NAMES as SPATIAL_SUMMARY_NAMES,
    cluster_spatial_scale,
)
from .hypothesis import (
    wilcoxon_z,
    rank_biserial,
    boot_ci_mean,
    bh_fdr,
    cluster_stats,
)
from .spectral import (
    principal_angles,
    chordal_distance,
    chordal_from_angles,
    grassmann_to_coord_subspace,
    chordal_full_vs_resect,
)

__all__ = [
    "compute_vi",
    "conditional_entropy",
    "compute_phase_distance_matrix",
    "compute_cross_patient_consistency",
    "rank_distance_measures",
    "UltrametricComparison",
    "compare_ultrametric_matrices",
    "compare_fc_methods",
    "compare_phases",
    "aggregate_comparisons",
    "batch_compare_fc_methods",
    "batch_compare_phases",
    "build_metric_specs",
    "compute_metric_matrix",
    "compute_cluster_labels",
    "cluster_swap_coefficient",
    "compute_cluster_swap_matrix",
    "compute_ari_matrix",
    "dmax_from_Z",
    "tree_internal_nodes",
    "jaccard_leafsets",
    "h_log_grid",
    "fcluster_at_h_rel",
    "kc_vectors",
    "kc_distance",
    "matching_cluster_distance",
    "weighted_rf_distance",
    "FunctionalTreeDistanceResult",
    "tau_star_from_C",
    "tau_interval",
    "delta_curves",
    "compute_functional_tree_distance",
    "SPATIAL_SUMMARY_NAMES",
    "cluster_spatial_scale",
    "wilcoxon_z",
    "rank_biserial",
    "boot_ci_mean",
    "bh_fdr",
    "cluster_stats",
    "principal_angles",
    "chordal_distance",
    "chordal_from_angles",
    "grassmann_to_coord_subspace",
    "chordal_full_vs_resect",
]
