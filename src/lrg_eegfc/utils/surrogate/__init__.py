"""Graph-level surrogate utilities (matched-strength, etc.).

Distinct from `lrg_eegfc.utils.fc.msc.surrogates` which operates on raw
timeseries (Welch / phase randomization). This subpackage operates on
already-computed FC adjacency matrices.
"""
from .matched_strength import (
    SURROGATE_LRG_CACHE,
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_adjacency,
    cophenetic_condensed_from_eigs,
    coupled_surrogate_cophenet,
    load_or_compute_eigs_at_path,
    load_or_compute_surrogate_eigs,
    strength_preserving_shuffle,
    surrogate_cache_path,
    verify_strengths,
)
from .timeseries_nulls import (
    band_bin_frequencies,
    block_partition_indices,
    blocks_to_segment_sets,
    coherency_from_csd,
    csd_from_segment_subset,
    global_rescale_invariance_residual,
    imcoh_abs_from_coherency,
    lag_randomized_coherency,
    phase_randomized_coherency,
    segment_fft,
    segment_shifted_csd,
)
from .coherency_surrogate import (
    complex_coherency_band,
    complex_coherency_bands,
    coupling_randomized_coherency,
    deviation_rotated_coherency,
    haar_orthogonal,
    haar_unitary,
    shared_backbone_deviation,
)

__all__ = [
    "SURROGATE_LRG_CACHE",
    "adjacency_from_laplacian_eigs",
    "cophenetic_condensed_from_adjacency",
    "cophenetic_condensed_from_eigs",
    "coupled_surrogate_cophenet",
    "load_or_compute_eigs_at_path",
    "load_or_compute_surrogate_eigs",
    "strength_preserving_shuffle",
    "surrogate_cache_path",
    "verify_strengths",
    # on-manifold coherency surrogates (CRC)
    "complex_coherency_band",
    "complex_coherency_bands",
    "coupling_randomized_coherency",
    "haar_orthogonal",
    "haar_unitary",
    # coordinated cross-phase null (SB-CRC)
    "shared_backbone_deviation",
    "deviation_rotated_coherency",
    # timeseries / coherency-level nulls (the W0-B ladder)
    "band_bin_frequencies",
    "block_partition_indices",
    "blocks_to_segment_sets",
    "coherency_from_csd",
    "csd_from_segment_subset",
    "global_rescale_invariance_residual",
    "imcoh_abs_from_coherency",
    "lag_randomized_coherency",
    "phase_randomized_coherency",
    "segment_fft",
    "segment_shifted_csd",
]
