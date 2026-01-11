"""Workflow orchestration subpackage."""

from .core import BandComputationResult, compute_band_connectivity
from .corr import CorrResult, compute_corr_matrix, load_corr_matrix, get_corr_cache_path, compute_corr_for_patient
from .msc import MSCResult, compute_msc_matrix, load_msc_matrix, get_msc_cache_path, compute_msc_for_patient
from .lrg import LRGResult, compute_lrg_analysis, load_lrg_result, get_lrg_cache_path, compute_lrg_for_patient
from .cleaning import (
    CleanedCorrResult,
    clean_correlation_matrix_full,
    load_cleaned_corr_matrix,
    get_cleaned_corr_cache_path,
    compute_cleaned_corr_for_patient,
)

__all__ = [
    "BandComputationResult",
    "compute_band_connectivity",
    "CorrResult",
    "compute_corr_matrix",
    "load_corr_matrix",
    "get_corr_cache_path",
    "compute_corr_for_patient",
    "MSCResult",
    "compute_msc_matrix",
    "load_msc_matrix",
    "get_msc_cache_path",
    "compute_msc_for_patient",
    "LRGResult",
    "compute_lrg_analysis",
    "load_lrg_result",
    "get_lrg_cache_path",
    "compute_lrg_for_patient",
    "CleanedCorrResult",
    "clean_correlation_matrix_full",
    "load_cleaned_corr_matrix",
    "get_cleaned_corr_cache_path",
    "compute_cleaned_corr_for_patient",
]
