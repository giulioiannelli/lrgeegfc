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
]
