"""Surrogate data generation for coherence null model estimation."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    Worker function to compute MSC for a single surrogate.

    This function is designed to be called in a separate process.
    """
    rng = np.random.default_rng(seed)
    N, L = X.shape

    # Generate surrogate via vectorized circular shifts
    shifts = rng.integers(0, L, size=N)
    indices = (np.arange(L)[None, :] - shifts[:, None]) % L
    X_surr = np.take_along_axis(X, indices, axis=1)

    # Compute MSC for this surrogate
    freqs, Coh = compute_msc_welch(X_surr, fs, nperseg=nperseg, noverlap=noverlap)

    # Band average and return
    return band_average_msc(Coh, freqs, bands)


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

    Returns
    -------
    W_null : Dict[str, NDArray]
        Dictionary mapping band names to null MSC matrices of shape (n_surrogates, N, N)
        where W_null[band][r, i, j] is the MSC between channels i and j in band for
        surrogate r

    Notes
    -----
    This function generates circular shift surrogates in parallel using ProcessPoolExecutor.
    Each surrogate MSC computation is independent (embarrassingly parallel).

    Memory usage per worker: O(N × L)
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
    else:
        # Parallel execution
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Submit all surrogate computations
            futures = {
                executor.submit(
                    _compute_single_surrogate_msc,
                    X, fs, bands, nperseg, noverlap, int(seeds[r])
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

    return W_null
