"""MSC-based functional connectivity workflow.

This module provides a simple workflow for computing and caching MSC-based
functional connectivity matrices. It mirrors the correlation workflow but
uses magnitude-squared coherence instead.
"""

from __future__ import annotations

import gc
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import networkx as nx
import numpy as np

from lrg_eegfc.config.const import BRAIN_BANDS, DEFAULT_NPERSEG, PHASE_LABELS
from lrg_eegfc.config.paths import MSC_CACHE, MSC_DEV_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc import (
    coherence_fc_pipeline,
    surrogate_msc_null,
    soft_sparsify_surrogate,
    fdr_sparsify_surrogate,
    disparity_filter,
    hybrid_sparsify,
    ecm_sparsify,
    ecm_sparsify_adaptive,
)
from lrg_eegfc.utils.io import load_timeseries, load_patient_dataset_robust

DEFAULT_MSC_CACHE_ROOT = MSC_CACHE
DEFAULT_MSC_DEV_CACHE_ROOT = MSC_DEV_CACHE

__all__ = ["MSCResult", "compute_msc_matrix", "load_msc_matrix", "get_msc_cache_path"]


def _resolve_cache_root(cache_root: Path, filter_time: Optional[int]) -> Path:
    if filter_time is not None and filter_time > 0 and cache_root == DEFAULT_MSC_CACHE_ROOT:
        return DEFAULT_MSC_DEV_CACHE_ROOT
    return cache_root


@dataclass
class MSCResult:
    """Result from MSC computation.

    Attributes
    ----------
    adjacency_matrix : np.ndarray
        MSC adjacency matrix (N x N)
    graph : nx.Graph
        NetworkX graph from adjacency matrix
    patient : str
        Patient identifier
    phase : str
        Recording phase
    band : str
        Frequency band
    n_channels : int
        Number of channels
    mean_msc : float
        Mean MSC value (excluding diagonal)
    sparsify : str
        Sparsification method ("none", "soft", "fdr", "disparity", "hybrid",
        "ecm", "ecm_adaptive")
    n_surrogates : int
        Number of surrogates used (0 if no surrogates needed)
    nperseg : int
        Window length for Welch's method
    fdr_q : float or None
        FDR q-value (only for sparsify="fdr")
    disparity_alpha : float or None
        Disparity alpha level (only for sparsify="disparity" or "hybrid")
    ecm_alpha : float or None
        ECM significance level (only for sparsify="ecm"; for "ecm_adaptive"
        this stores the *discovered* alpha)
    ecm_n_ensemble : int or None
        ECM ensemble size (only for sparsify="ecm" or "ecm_adaptive")
    ecm_weight_scale : int or None
        ECM weight scaling factor (only for sparsify="ecm" or "ecm_adaptive")
    ecm_alpha_min : float or None
        Lower bound of adaptive alpha search (only for sparsify="ecm_adaptive")
    ecm_alpha_max : float or None
        Upper bound of adaptive alpha search (only for sparsify="ecm_adaptive")
    """

    adjacency_matrix: np.ndarray
    graph: nx.Graph
    patient: str
    phase: str
    band: str
    n_channels: int
    mean_msc: float
    sparsify: str
    n_surrogates: int
    nperseg: int
    fdr_q: Optional[float] = None
    disparity_alpha: Optional[float] = None
    ecm_alpha: Optional[float] = None
    ecm_n_ensemble: Optional[int] = None
    ecm_weight_scale: Optional[int] = None
    ecm_alpha_min: Optional[float] = None
    ecm_alpha_max: Optional[float] = None


def get_msc_cache_path(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = DEFAULT_MSC_CACHE_ROOT,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = DEFAULT_NPERSEG,
    filter_time: Optional[int] = None,
    fdr_q: Optional[float] = None,
    disparity_alpha: Optional[float] = None,
    ecm_alpha: Optional[float] = None,
    ecm_n_ensemble: Optional[int] = None,
    ecm_weight_scale: Optional[int] = None,
    ecm_alpha_min: Optional[float] = None,
    ecm_alpha_max: Optional[float] = None,
) -> Path:
    """Get cache file path for MSC matrix.

    Parameters
    ----------
    patient : str
        Patient identifier
    phase : str
        Recording phase
    band : str
        Frequency band
    cache_root : Path, optional
        Root directory for cache files
    sparsify : str, optional
        Sparsification method
    n_surrogates : int, optional
        Number of surrogates (0 if not needed)
    nperseg : int, optional
        Window length for Welch's method
    fdr_q : float, optional
        FDR q-value (only for sparsify="fdr")
    disparity_alpha : float, optional
        Disparity significance level (only for sparsify="disparity" or "hybrid")

    Returns
    -------
    Path
        Path to cache file
    """
    cache_root = _resolve_cache_root(cache_root, filter_time)
    cache_dir = cache_root / patient
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Include sparsification parameters in filename
    if sparsify == "none":
        suffix = f"sparsify-none_nperseg-{nperseg}"
    elif sparsify == "soft":
        suffix = f"sparsify-soft_nsurr-{n_surrogates}_nperseg-{nperseg}"
    elif sparsify == "fdr":
        suffix = f"sparsify-fdr_nsurr-{n_surrogates}_q-{fdr_q}_nperseg-{nperseg}"
    elif sparsify == "disparity":
        suffix = f"sparsify-disparity_alpha-{disparity_alpha}_nperseg-{nperseg}"
    elif sparsify == "hybrid":
        suffix = f"sparsify-hybrid_nsurr-{n_surrogates}_alpha-{disparity_alpha}_nperseg-{nperseg}"
    elif sparsify == "ecm":
        suffix = f"sparsify-ecm_alpha-{ecm_alpha}_nens-{ecm_n_ensemble}_wscale-{ecm_weight_scale}_nperseg-{nperseg}"
    elif sparsify == "ecm_adaptive":
        suffix = f"sparsify-ecm_adaptive_alpharange-{ecm_alpha_min}-{ecm_alpha_max}_nens-{ecm_n_ensemble}_wscale-{ecm_weight_scale}_nperseg-{nperseg}"
    else:
        suffix = f"sparsify-{sparsify}_nsurr-{n_surrogates}_nperseg-{nperseg}"

    if filter_time is not None and filter_time > 0:
        suffix = f"{suffix}_ftime-{filter_time}"

    return cache_dir / f"{band}_{phase}_msc_{suffix}.npy"


def load_msc_matrix(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = DEFAULT_MSC_CACHE_ROOT,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = DEFAULT_NPERSEG,
    filter_time: Optional[int] = None,
    fdr_q: Optional[float] = None,
    disparity_alpha: Optional[float] = None,
    ecm_alpha: Optional[float] = None,
    ecm_n_ensemble: Optional[int] = None,
    ecm_weight_scale: Optional[int] = None,
    ecm_alpha_min: Optional[float] = None,
    ecm_alpha_max: Optional[float] = None,
) -> Optional[np.ndarray]:
    """Load cached MSC matrix if it exists.

    Parameters
    ----------
    patient : str
        Patient identifier
    phase : str
        Recording phase
    band : str
        Frequency band
    cache_root : Path, optional
        Root directory for cache files
    sparsify : str, optional
        Sparsification method
    n_surrogates : int, optional
        Number of surrogates
    nperseg : int, optional
        Window length for Welch's method
    fdr_q : float, optional
        FDR q-value (only for sparsify="fdr")
    disparity_alpha : float, optional
        Disparity significance level

    Returns
    -------
    np.ndarray or None
        Cached MSC matrix, or None if not cached
    """
    cache_path = get_msc_cache_path(
        patient, phase, band, cache_root,
        sparsify, n_surrogates, nperseg, filter_time,
        fdr_q=fdr_q, disparity_alpha=disparity_alpha,
        ecm_alpha=ecm_alpha, ecm_n_ensemble=ecm_n_ensemble,
        ecm_weight_scale=ecm_weight_scale,
        ecm_alpha_min=ecm_alpha_min, ecm_alpha_max=ecm_alpha_max,
    )

    if cache_path.exists():
        return np.load(cache_path)
    return None


def compute_msc_matrix(
    patient: str,
    phase: str,
    band: str,
    dataset_root: Path = SEEG_DATAPATH,
    cache_root: Path = DEFAULT_MSC_CACHE_ROOT,
    *,
    use_cache: bool = True,
    overwrite_cache: bool = False,
    sample_rate: float = 2048.0,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = DEFAULT_NPERSEG,
    noverlap: Optional[int] = None,
    batch_size: int = 64,
    n_workers: Optional[int] = None,
    filter_time: Optional[int] = None,
    verbose: bool = False,
    fdr_q: float = 0.05,
    disparity_alpha: float = 0.05,
    ecm_alpha: float = 0.05,
    ecm_n_ensemble: int = 100,
    ecm_weight_scale: int = 1000,
    ecm_alpha_min: float = 0.01,
    ecm_alpha_max: float = 0.50,
) -> MSCResult:
    """Compute MSC-based functional connectivity matrix.

    This function computes magnitude-squared coherence (MSC) for a given
    patient, phase, and frequency band. Results are automatically cached
    for faster subsequent access.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02")
    phase : str
        Recording phase (e.g., "rest_pre", "task_learn")
    band : str
        Frequency band (must exist in BRAIN_BANDS)
    dataset_root : Path, optional
        Root directory containing patient data
    cache_root : Path, optional
        Root directory for cache files
    use_cache : bool, optional
        If True, load from cache if available (default: True)
    overwrite_cache : bool, optional
        If True, recompute even if cached (default: False)
    sample_rate : float, optional
        Sampling rate in Hz (default: 2048.0)
    sparsify : str, optional
        Sparsification method: "none" (dense MSC) or "soft" (surrogate-based
        soft sparsification) (default: "none")
    n_surrogates : int, optional
        Number of circular shift surrogates for validation. Only used if
        sparsify="soft". Set to 0 for no validation (default: 0)
    nperseg : int, optional
        Window length for Welch's method (default: 1024)
    noverlap : int, optional
        Overlap between windows (default: nperseg // 2)
    batch_size : int, optional
        Number of Welch segments to process per batch (reduces Python/FFT overhead;
        default: 64)
    n_workers : int, optional
        Number of parallel workers for surrogate computation. If None, uses all CPU
        cores. Set to 1 to disable parallelism (default: None)
    filter_time : int, optional
        Limit to first N samples (for testing)
    verbose : bool, optional
        Print progress information (default: False)

    Returns
    -------
    MSCResult
        Result object containing adjacency matrix and metadata

    Raises
    ------
    KeyError
        If band not defined in BRAIN_BANDS
    FileNotFoundError
        If patient/phase data doesn't exist

    Examples
    --------
    >>> # Dense MSC (no validation)
    >>> result = compute_msc_matrix("Pat_02", "rest_pre", "beta")
    >>> print(result.adjacency_matrix.shape)
    (117, 117)
    >>> print(f"Mean MSC: {result.mean_msc:.4f}")
    Mean MSC: 0.0523
    >>> print(f"Sparsify: {result.sparsify}, Surrogates: {result.n_surrogates}")
    Sparsify: none, Surrogates: 0

    >>> # Validated MSC with surrogates
    >>> result = compute_msc_matrix("Pat_02", "rest_pre", "beta",
    ...                             sparsify="soft", n_surrogates=200)
    >>> print(f"Sparsify: {result.sparsify}, Surrogates: {result.n_surrogates}")
    Sparsify: soft, Surrogates: 200
    """

    # Validate band
    if band not in BRAIN_BANDS:
        available = ", ".join(sorted(BRAIN_BANDS))
        raise KeyError(f"Band '{band}' not defined. Available: {available}")

    cache_root = _resolve_cache_root(cache_root, filter_time)

    _NEEDS_SURROGATES = {"soft", "fdr", "hybrid"}

    # Check cache
    _ECM_METHODS = {"ecm", "ecm_adaptive"}
    cache_path = get_msc_cache_path(
        patient, phase, band, cache_root,
        sparsify, n_surrogates, nperseg, filter_time,
        fdr_q=fdr_q if sparsify == "fdr" else None,
        disparity_alpha=disparity_alpha if sparsify in ("disparity", "hybrid") else None,
        ecm_alpha=ecm_alpha if sparsify == "ecm" else None,
        ecm_n_ensemble=ecm_n_ensemble if sparsify in _ECM_METHODS else None,
        ecm_weight_scale=ecm_weight_scale if sparsify in _ECM_METHODS else None,
        ecm_alpha_min=ecm_alpha_min if sparsify == "ecm_adaptive" else None,
        ecm_alpha_max=ecm_alpha_max if sparsify == "ecm_adaptive" else None,
    )

    if use_cache and not overwrite_cache and cache_path.exists():
        if verbose:
            print(f"  loading cached msc: {cache_path.name}")
        adj_matrix = np.load(cache_path)

        # Reconstruct result
        graph = nx.from_numpy_array(adj_matrix)
        n_channels = adj_matrix.shape[0]

        # Compute mean MSC (upper triangle, excluding diagonal)
        triu_indices = np.triu_indices_from(adj_matrix, k=1)
        mean_msc = float(np.mean(adj_matrix[triu_indices]))

        return MSCResult(
            adjacency_matrix=adj_matrix,
            graph=graph,
            patient=patient,
            phase=phase,
            band=band,
            n_channels=n_channels,
            mean_msc=mean_msc,
            sparsify=sparsify,
            n_surrogates=n_surrogates if sparsify in _NEEDS_SURROGATES else 0,
            nperseg=nperseg,
            fdr_q=fdr_q if sparsify == "fdr" else None,
            disparity_alpha=disparity_alpha if sparsify in ("disparity", "hybrid") else None,
            ecm_alpha=ecm_alpha if sparsify == "ecm" else None,
            ecm_n_ensemble=ecm_n_ensemble if sparsify in _ECM_METHODS else None,
            ecm_weight_scale=ecm_weight_scale if sparsify in _ECM_METHODS else None,
            ecm_alpha_min=ecm_alpha_min if sparsify == "ecm_adaptive" else None,
            ecm_alpha_max=ecm_alpha_max if sparsify == "ecm_adaptive" else None,
        )

    # Compute MSC from timeseries
    if verbose:
        print(f"Computing MSC for {patient} {phase} {band}...")
        print(f"  Sparsify: {sparsify}, n_surrogates: {n_surrogates}")
        print(f"  Loading timeseries and extracting sampling rate...")

    # Load data with robust loader to extract sampling rate
    dataset = load_patient_dataset_robust(patient, dataset_root, phases=[phase])

    if phase not in dataset:
        raise FileNotFoundError(f"Phase '{phase}' not found for patient '{patient}'")

    recording = dataset[phase]
    data = recording.timeseries

    # Extract sampling rate from parameters, or use provided default
    if 'fs' in recording.parameters and recording.parameters['fs'] is not None:
        actual_fs = float(recording.parameters['fs'])
        if verbose:
            print(f"  Using sampling rate from data: {actual_fs} Hz")
        # Override sample_rate with actual value from data
        sample_rate = actual_fs
    else:
        if verbose:
            print(f"  WARNING: No sampling rate found in data, using provided value: {sample_rate} Hz")

    if filter_time is not None and filter_time > 0:
        data = data[:, :filter_time]

    if verbose:
        print(f"  Data shape: {data.shape}")

    # Compute MSC for single band
    bands_dict = {band: BRAIN_BANDS[band]}

    # For soft sparsification, we first need the dense MSC
    # Check if dense MSC is already cached (to avoid recomputing)
    dense_cache_path = get_msc_cache_path(
        patient, phase, band, cache_root,
        sparsify="none", n_surrogates=0, nperseg=nperseg, filter_time=filter_time,
    )

    if dense_cache_path.exists() and not overwrite_cache:
        if verbose:
            print(f"  Loading cached dense MSC: {dense_cache_path.name}")
        dense_msc = np.load(dense_cache_path)
    else:
        if verbose:
            print(f"  Computing dense MSC with Welch's method...")
        # Compute dense MSC (no sparsification)
        dense_matrices = coherence_fc_pipeline(
            data,
            fs=sample_rate,
            bands=bands_dict,
            sparsify="none",
            nperseg=nperseg,
            noverlap=noverlap,
            batch_size=batch_size,
            zero_diagonal=True,
            verbose=verbose,
        )
        dense_msc = dense_matrices[band]

        # Cache the dense MSC for future use
        if verbose:
            print(f"  Saving dense MSC to cache: {dense_cache_path.name}")
        np.save(dense_cache_path, dense_msc)

    # Apply sparsification if requested
    if sparsify in _NEEDS_SURROGATES and n_surrogates > 0:
        if verbose:
            print(f"  Computing {n_surrogates} surrogates for {sparsify} sparsification...")

        W_null = surrogate_msc_null(
            data, sample_rate, bands_dict, n_surrogates,
            nperseg=nperseg, noverlap=noverlap, n_workers=n_workers,
        )

        if verbose:
            print(f"  Applying {sparsify} sparsification...")

        if sparsify == "soft":
            adj_matrix = soft_sparsify_surrogate(dense_msc, W_null[band])
        elif sparsify == "fdr":
            adj_matrix = fdr_sparsify_surrogate(dense_msc, W_null[band], q=fdr_q)
        elif sparsify == "hybrid":
            adj_matrix = hybrid_sparsify(dense_msc, W_null[band], alpha=disparity_alpha)

        np.fill_diagonal(adj_matrix, 0.0)

        del W_null
        gc.collect()

    elif sparsify == "disparity":
        if verbose:
            print(f"  Applying disparity filter (alpha={disparity_alpha})...")
        adj_matrix = disparity_filter(dense_msc, alpha=disparity_alpha)

    elif sparsify == "ecm":
        if verbose:
            print(f"  Applying ECM filter (alpha={ecm_alpha}, n_ensemble={ecm_n_ensemble}, scale={ecm_weight_scale})...")
        adj_matrix = ecm_sparsify(
            dense_msc,
            alpha=ecm_alpha,
            n_ensemble=ecm_n_ensemble,
            weight_scale=ecm_weight_scale,
        )

    elif sparsify == "ecm_adaptive":
        if verbose:
            print(f"  Applying adaptive ECM filter (alpha_range=[{ecm_alpha_min}, {ecm_alpha_max}], "
                  f"n_ensemble={ecm_n_ensemble}, scale={ecm_weight_scale})...")
        adj_matrix, alpha_used = ecm_sparsify_adaptive(
            dense_msc,
            alpha_min=ecm_alpha_min,
            alpha_max=ecm_alpha_max,
            n_ensemble=ecm_n_ensemble,
            weight_scale=ecm_weight_scale,
        )
        # Store the discovered alpha in ecm_alpha for the result
        ecm_alpha = alpha_used
        if verbose:
            print(f"  Adaptive ECM selected alpha={alpha_used:.4f}")

    else:
        adj_matrix = dense_msc

    # Free intermediate data no longer needed
    del data, dataset, recording, dense_msc
    gc.collect()

    # Cache the final result (sparsified or dense)
    if verbose:
        print(f"  Saving to cache: {cache_path.name}")
    np.save(cache_path, adj_matrix)

    # Build graph
    graph = nx.from_numpy_array(adj_matrix)
    n_channels = adj_matrix.shape[0]

    # Compute statistics
    triu_indices = np.triu_indices_from(adj_matrix, k=1)
    mean_msc = float(np.mean(adj_matrix[triu_indices]))

    if verbose:
        print(f"  Done! Mean MSC: {mean_msc:.4f}")

    return MSCResult(
        adjacency_matrix=adj_matrix,
        graph=graph,
        patient=patient,
        phase=phase,
        band=band,
        n_channels=n_channels,
        mean_msc=mean_msc,
        sparsify=sparsify,
        n_surrogates=n_surrogates if sparsify in _NEEDS_SURROGATES else 0,
        nperseg=nperseg,
        fdr_q=fdr_q if sparsify == "fdr" else None,
        disparity_alpha=disparity_alpha if sparsify in ("disparity", "hybrid") else None,
        ecm_alpha=ecm_alpha if sparsify in _ECM_METHODS else None,
        ecm_n_ensemble=ecm_n_ensemble if sparsify in _ECM_METHODS else None,
        ecm_weight_scale=ecm_weight_scale if sparsify in _ECM_METHODS else None,
        ecm_alpha_min=ecm_alpha_min if sparsify == "ecm_adaptive" else None,
        ecm_alpha_max=ecm_alpha_max if sparsify == "ecm_adaptive" else None,
    )


def compute_msc_for_patient(
    patient: str,
    bands: Optional[list] = None,
    phases: Optional[list] = None,
    return_results: bool = False,
    **kwargs
) -> Dict[str, Dict[str, Optional[MSCResult]]]:
    """Compute MSC matrices for all band/phase combinations for a patient.

    Parameters
    ----------
    patient : str
        Patient identifier
    bands : list, optional
        List of band names (default: all BRAIN_BANDS)
    phases : list, optional
        List of phase names (default: all phases)
    return_results : bool, optional
        If True, return full MSCResult objects (high memory usage).
        If False (default), return lightweight status dict with just success/failure
        markers. Use False for batch processing to avoid memory accumulation.
    **kwargs
        Additional arguments passed to compute_msc_matrix()

    Returns
    -------
    Dict[str, Dict[str, MSCResult | bool | None]]
        Nested dict indexed by [band][phase]. If return_results=True, contains
        MSCResult objects. If return_results=False, contains True for success,
        None for failure.

    Examples
    --------
    >>> # Compute all MSC matrices for Pat_02
    >>> results = compute_msc_for_patient("Pat_02")
    >>>
    >>> # Access specific result
    >>> beta_rsPre = results["beta"]["rest_pre"]
    >>> print(beta_rsPre.mean_msc)
    """
    if bands is None:
        bands = list(BRAIN_BANDS.keys())
    if phases is None:
        phases = list(PHASE_LABELS)

    results = {band: {} for band in bands}

    for band in bands:
        for phase in phases:
            try:
                result = compute_msc_matrix(patient, phase, band, **kwargs)
                if return_results:
                    results[band][phase] = result
                else:
                    # Don't store the full result to save memory
                    # Just mark as successful (True-like for backwards compat)
                    results[band][phase] = True
                    del result
            except Exception as e:
                print(f"WARNING: Failed to compute {patient} {phase} {band}: {e}")
                results[band][phase] = None

            # Force garbage collection after each computation to free memory
            gc.collect()

    return results
