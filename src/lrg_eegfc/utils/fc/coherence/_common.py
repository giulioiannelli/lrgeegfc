"""Shared spectral estimation primitives for coherence-based FC.

This module factors out the Welch/cross-spectral-density (CSD) computation
and the band-averaging step from the legacy ``msc.py`` implementation so
that MSC, ImCoh, and future coherence-family metrics share a single,
well-tested backend.

References
----------
- Welch, P. D. (1967) "The use of fast Fourier transform for the
  estimation of power spectra: A method based on time averaging over
  short, modified periodograms." IEEE Trans. Audio Electroacoust.
  15(2):70-73.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.signal import get_window
from numpy.lib.stride_tricks import sliding_window_view


__all__ = [
    "welch_csd",
    "band_average",
]


def welch_csd(
    X: NDArray,
    fs: float,
    nperseg: int = 256,
    noverlap: int | None = None,
    batch_size: int = 64,
) -> Tuple[NDArray, NDArray]:
    """Vectorized Welch estimator of the (complex) cross-spectral density.

    Computes the full ``(N, N, F)`` CSD matrix in a single pass using
    Hanning-windowed segments and an FFT-batch accumulation strategy.

    Parameters
    ----------
    X : NDArray
        Time series of shape ``(N, L)`` (N channels, L samples).
    fs : float
        Sampling rate in Hz.
    nperseg : int
        Segment length in samples.
    noverlap : int, optional
        Samples of overlap between successive segments.  Defaults to
        ``nperseg // 2``.
    batch_size : int
        Number of segments per FFT batch.

    Returns
    -------
    freqs : NDArray
        One-sided frequency grid of length ``F``.
    CSD : NDArray
        Complex CSD of shape ``(N, N, F)``, scaled by ``1 / (fs * W *
        n_segments)`` where ``W = sum(window^2)``.

    Notes
    -----
    The returned CSD is Hermitian in the channel axes: ``CSD[j, i, :] =
    conj(CSD[i, j, :])``.  Diagonal entries (``i == j``) are real and
    equal to the one-sided auto-power-spectral densities.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    N, L = X.shape
    if noverlap is None:
        noverlap = nperseg // 2

    step = nperseg - noverlap
    segments = sliding_window_view(X, window_shape=nperseg, axis=1)[:, ::step, :]
    n_segments = segments.shape[1]
    if n_segments == 0:
        raise ValueError("nperseg larger than signal length")

    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    F = len(freqs)

    window = get_window("hann", nperseg).astype(X.dtype, copy=False)
    scale = 1.0 / (fs * (window * window).sum() * n_segments)

    CSD = np.zeros((N, N, F), dtype=complex)

    for start in range(0, n_segments, batch_size):
        stop = start + batch_size
        block = segments[:, start:stop, :] * window[None, None, :]
        fft_block = np.fft.rfft(block, n=nperseg, axis=2)
        CSD += np.einsum(
            "nbf,mbf->nmf", fft_block, np.conj(fft_block), optimize="greedy"
        )

    CSD *= scale
    return freqs, CSD


def band_average(
    Coh: NDArray,
    freqs: NDArray,
    bands: Dict[str, Tuple[float, float]],
) -> Dict[str, NDArray]:
    """Average a ``(N, N, F)`` coherence array over named frequency bands.

    Parameters
    ----------
    Coh : NDArray
        Coherence-like tensor of shape ``(N, N, F)``.
    freqs : NDArray
        Frequency grid of length ``F`` in Hz.
    bands : dict of {str: (fmin, fmax)}
        Band definitions in Hz.

    Returns
    -------
    dict of {str: NDArray}
        For each band, the ``(N, N)`` mean of ``Coh[:, :, mask]`` where
        ``mask = (freqs >= fmin) & (freqs <= fmax)``.  Bands with no
        frequencies in range return an all-zero matrix.
    """
    N = Coh.shape[0]
    W_bands: Dict[str, NDArray] = {}
    for band_name, (fmin, fmax) in bands.items():
        band_mask = (freqs >= fmin) & (freqs <= fmax)
        if not np.any(band_mask):
            W_bands[band_name] = np.zeros((N, N))
            continue
        W_bands[band_name] = np.mean(Coh[:, :, band_mask], axis=2)
    return W_bands
