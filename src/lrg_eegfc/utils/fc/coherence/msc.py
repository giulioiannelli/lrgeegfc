"""Magnitude-squared coherence (MSC) from the Welch cross-spectral density.

MSC between channels i and j at frequency f is::

    MSC_ij(f) = |S_ij(f)|^2 / (S_ii(f) * S_jj(f))

which lies in ``[0, 1]`` and is symmetric (``MSC_ij = MSC_ji``). MSC is
sensitive to volume conduction because instantaneous mixing (zero-phase
coupling) increases ``|S_ij|``; for volume-conduction-immune alternatives
see :mod:`lrg_eegfc.utils.fc.coherence.imcoh`.

References
----------
- Carter, G. C. (1987) "Coherence and time delay estimation." Proc. IEEE
  75(2):236-255.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from ._common import welch_csd


__all__ = ["compute_msc"]


def compute_msc(
    X: NDArray,
    fs: float,
    nperseg: int = 256,
    noverlap: int | None = None,
    batch_size: int = 64,
) -> Tuple[NDArray, NDArray]:
    """Compute magnitude-squared coherence via vectorized Welch.

    Parameters
    ----------
    X : NDArray
        Time series ``(N, L)``.
    fs : float
        Sampling rate in Hz.
    nperseg, noverlap, batch_size
        Passed through to :func:`welch_csd`.

    Returns
    -------
    freqs : NDArray
        Frequency grid ``(F,)``.
    MSC : NDArray
        ``(N, N, F)`` MSC in ``[0, 1]``; diagonal is set to 1.
    """
    N = X.shape[0]
    freqs, CSD = welch_csd(
        X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
    )
    F = freqs.shape[0]
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)  # (N, F)
    denom = PSD[:, None, :] * PSD[None, :, :]
    CSD_mag_sq = np.abs(CSD) ** 2
    MSC = np.divide(
        CSD_mag_sq, denom, out=np.zeros((N, N, F)), where=denom > 0
    )
    np.einsum("iif->if", MSC)[...] = 1.0
    return freqs, MSC
