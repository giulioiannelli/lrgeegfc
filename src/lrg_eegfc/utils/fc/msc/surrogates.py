"""Surrogate data generation for coherence null model estimation."""

from __future__ import annotations

import gc
import os
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory
from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from .msc import compute_msc_welch, band_average_msc


__all__ = [
    'circular_shift_surrogates',
    'surrogate_msc_null',
]


def _compute_single_surrogate_msc(
    X: NDArray,
    fs: float,
    bands: Dict[str, Tuple[float, float]],
    nperseg: int,
    noverlap: int | None,
    seed: int,
) -> Dict[str, NDArray]:
    """
    Worker function to compute MSC for a single surrogate (non-shared memory version).

    This function is designed to be called in a separate process.
    Uses memory-efficient circular shifting.
    """
    rng = np.random.default_rng(seed)
    N, L = X.shape

    # Generate random shifts
    shifts = rng.integers(0, L, size=N)

    # Memory-efficient circular shift (no large indices array)
    X_surr = _circular_shift_memory_efficient(X, shifts)
    del shifts

    # Compute MSC for this surrogate
    freqs, Coh = compute_msc_welch(X_surr, fs, nperseg=nperseg, noverlap=noverlap)
    del X_surr

    # Band average
    result = band_average_msc(Coh, freqs, bands)
    del Coh, freqs

    return result


def _circular_shift_memory_efficient(X: NDArray, shifts: NDArray) -> NDArray:
    """
    Apply circular shifts to each channel without creating a large indices array.

    Uses np.roll channel-by-channel which is much more memory efficient than
    creating a full (N, L) indices array for advanced indexing.

    Parameters
    ----------
    X : NDArray
        Input data of shape (N, L)
    shifts : NDArray
        Array of shift amounts of shape (N,)

    Returns
    -------
    NDArray
        Shifted data of shape (N, L)
    """
    N, L = X.shape
    X_surr = np.empty_like(X)
    for i in range(N):
        X_surr[i] = np.roll(X[i], -int(shifts[i]))
    return X_surr


def _compute_surrogate_msc_shm(
    shm_name: str,
    shape: Tuple[int, int],
    dtype: np.dtype,
    fs: float,
    bands: Dict[str, Tuple[float, float]],
    nperseg: int,
    noverlap: int | None,
    seed: int,
    max_samples: int | None = None,
) -> Dict[str, NDArray]:
    """
    Worker function to compute MSC for a single surrogate using shared memory.

    Instead of receiving a copy of the data array, this function attaches to
    an existing shared memory block and creates a numpy view of the data.
    This avoids the memory overhead of pickling large arrays.

    Uses memory-efficient circular shifting that avoids creating large
    intermediate index arrays.

    Parameters
    ----------
    shm_name : str
        Name of the shared memory block containing the data
    shape : Tuple[int, int]
        Shape of the data array (N, L)
    dtype : np.dtype
        Data type of the array
    fs : float
        Sampling frequency in Hz
    bands : Dict[str, Tuple[float, float]]
        Dictionary mapping band names to (fmin, fmax) tuples
    nperseg : int
        Length of each segment for coherence computation
    noverlap : int | None
        Number of overlapping samples
    seed : int
        Random seed for reproducibility
    max_samples : int | None
        If set, truncate data to this many samples for faster computation.
        Surrogates only need enough data for statistical threshold estimation.

    Returns
    -------
    Dict[str, NDArray]
        Dictionary mapping band names to MSC matrices
    """
    # Attach to shared memory (read-only access)
    shm = shared_memory.SharedMemory(name=shm_name)
    try:
        # Create numpy array view of the shared memory
        X_full = np.ndarray(shape, dtype=dtype, buffer=shm.buf)

        # Truncate for faster surrogate computation if requested
        N, L = shape
        if max_samples is not None and L > max_samples:
            X = X_full[:, :max_samples]
            L = max_samples
        else:
            X = X_full

        rng = np.random.default_rng(seed)

        # Generate random shifts
        shifts = rng.integers(0, L, size=N)

        # Memory-efficient circular shift (no large indices array)
        X_surr = _circular_shift_memory_efficient(X, shifts)
        del shifts

        # Compute MSC for this surrogate
        freqs, Coh = compute_msc_welch(X_surr, fs, nperseg=nperseg, noverlap=noverlap)
        del X_surr
        gc.collect()

        # Band average
        result = band_average_msc(Coh, freqs, bands)

        # Explicitly free intermediate arrays to reduce worker memory
        del Coh, freqs
        gc.collect()

        return result
    finally:
        # Close (but don't unlink) the shared memory in the worker
        shm.close()


def circular_shift_surrogates(
    X: NDArray,
    n_surrogates: int,
    rng: np.random.Generator | None = None,
) -> NDArray:
    """
    Generate surrogate data via circular shifts of each channel.

    This method preserves each channel's spectrum and amplitude distribution but
    removes cross-channel coupling by shifting each channel by a random amount.

    Parameters
    ----------
    X : NDArray
        Input data of shape (N, L) where N is the number of channels and L is the
        number of samples
    n_surrogates : int
        Number of surrogate datasets to generate
    rng : np.random.Generator, optional
        Random number generator for reproducibility. If None, uses default RNG.

    Returns
    -------
    X_surr : NDArray
        Surrogate data of shape (n_surrogates, N, L) where X_surr[r, i, :] is the
        r-th surrogate of channel i

    Notes
    -----
    For each surrogate r and channel i:

    x_i^(r)(t) = x_i((t + Δ_i^(r)) mod L)

    where Δ_i^(r) is a random integer in [0, L-1], independent per channel.
    This preserves each channel's spectrum and amplitude distribution but removes
    cross-channel coupling.
    """
    if rng is None:
        rng = np.random.default_rng()

    N, L = X.shape
    X_surr = np.zeros((n_surrogates, N, L), dtype=X.dtype)

    # Vectorized circular shifts using advanced indexing
    base_indices = np.arange(L)
    for r in range(n_surrogates):
        # Generate all shifts for this surrogate at once
        shifts = rng.integers(0, L, size=N)
        # Create index array: (indices - shifts) mod L for each channel
        indices = (base_indices[None, :] - shifts[:, None]) % L
        X_surr[r] = np.take_along_axis(X, indices, axis=1)

    return X_surr


def surrogate_msc_null(
    X: NDArray,
    fs: float,
    bands: Dict[str, Tuple[float, float]],
    n_surrogates: int,
    nperseg: int = 256,
    noverlap: int | None = None,
    rng: np.random.Generator | None = None,
    n_workers: int | None = None,
    max_samples: int | None = 500_000,
) -> Dict[str, NDArray]:
    """
    Compute null distribution of band-averaged MSC using circular shift surrogates.

    Parameters
    ----------
    X : NDArray
        Input data of shape (N, L)
    fs : float
        Sampling frequency in Hz
    bands : Dict[str, Tuple[float, float]]
        Dictionary mapping band names to (fmin, fmax) tuples
    n_surrogates : int
        Number of surrogates to generate
    nperseg : int, optional
        Length of each segment for coherence computation (default: 256)
    noverlap : int, optional
        Number of overlapping samples (default: nperseg // 2)
    rng : np.random.Generator, optional
        Random number generator for reproducibility
    n_workers : int, optional
        Number of parallel workers. If None, uses number of CPU cores.
        Set to 1 to disable parallelism.
    max_samples : int, optional
        Maximum samples to use for surrogate computation (default: 500,000).
        Surrogates only need enough data for statistical threshold estimation,
        not the full timeseries. Set to None to use all samples.

    Returns
    -------
    W_null : Dict[str, NDArray]
        Dictionary mapping band names to null MSC matrices of shape (n_surrogates, N, N)
        where W_null[band][r, i, j] is the MSC between channels i and j in band for
        surrogate r

    Notes
    -----
    This function uses shared memory to avoid copying the data array to each worker
    process. This dramatically reduces memory usage when processing large datasets
    with multiple workers.

    Using max_samples (default 500k) speeds up computation ~4x for 2M+ sample datasets
    while maintaining statistical validity for null distribution estimation.
    """
    if rng is None:
        rng = np.random.default_rng()

    if n_workers is None:
        n_workers = os.cpu_count() or 1

    N, L = X.shape

    # Initialize null distribution dict
    W_null = {band_name: np.zeros((n_surrogates, N, N)) for band_name in bands}

    # Generate independent seeds for each surrogate (for reproducibility)
    seeds = rng.integers(0, 2**31, size=n_surrogates)

    if n_workers == 1:
        # Sequential execution (useful for debugging)
        for r in tqdm(range(n_surrogates), desc="Surrogates", unit="surr"):
            W_bands_r = _compute_single_surrogate_msc(
                X, fs, bands, nperseg, noverlap, int(seeds[r])
            )
            for band_name in bands:
                W_null[band_name][r] = W_bands_r[band_name]
            del W_bands_r
        gc.collect()
    else:
        # Parallel execution with shared memory
        # Ensure array is contiguous for shared memory
        X_contiguous = np.ascontiguousarray(X)

        # Create shared memory block with unique name
        shm_name = f"msc_surr_{uuid.uuid4().hex[:8]}"
        shm = shared_memory.SharedMemory(name=shm_name, create=True, size=X_contiguous.nbytes)

        try:
            # Copy data into shared memory
            shm_array = np.ndarray(X_contiguous.shape, dtype=X_contiguous.dtype, buffer=shm.buf)
            np.copyto(shm_array, X_contiguous)

            # Free the local copy now that it's in shared memory
            del X_contiguous
            gc.collect()

            # Parallel execution using shared memory
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                # Submit all surrogate computations
                futures = {
                    executor.submit(
                        _compute_surrogate_msc_shm,
                        shm_name,
                        shm_array.shape,
                        shm_array.dtype,
                        fs,
                        bands,
                        nperseg,
                        noverlap,
                        int(seeds[r]),
                        max_samples,
                    ): r
                    for r in range(n_surrogates)
                }

                # Collect results as they complete with progress bar
                for future in tqdm(
                    as_completed(futures),
                    total=n_surrogates,
                    desc=f"Surrogates ({n_workers} workers)",
                    unit="surr",
                ):
                    r = futures[future]
                    W_bands_r = future.result()
                    for band_name in bands:
                        W_null[band_name][r] = W_bands_r[band_name]
                    # Free the result dict to avoid accumulation
                    del W_bands_r

                # Clear the futures dict to release references
                futures.clear()

            # Executor is now closed, workers should be terminated
        finally:
            # Clean up shared memory
            shm.close()
            shm.unlink()

        # Force garbage collection after parallel execution
        gc.collect()

    return W_null
