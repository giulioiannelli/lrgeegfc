"""Coherence-family FC metrics built on a shared Welch CSD backend.

Public API
----------
- :func:`welch_csd` : Welch cross-spectral density (complex ``(N, N, F)``).
- :func:`band_average` : Average a coherence tensor over named frequency bands.
- :func:`compute_msc` : Magnitude-squared coherence ``in [0, 1]``.
- :func:`compute_imcoh` : Signed imaginary coherency ``in [-1, 1]``
  (Nolte 2004).
- :func:`compute_coherence` : Unified dispatcher keyed on ``coh_method``.
"""
from __future__ import annotations

from typing import Tuple

from numpy.typing import NDArray

from ._common import welch_csd, band_average
from .msc import compute_msc
from .imcoh import compute_imcoh


__all__ = [
    "welch_csd",
    "band_average",
    "compute_msc",
    "compute_imcoh",
    "compute_coherence",
]


def compute_coherence(
    X: NDArray,
    fs: float,
    coh_method: str = "msc",
    nperseg: int = 256,
    noverlap: int | None = None,
    batch_size: int = 64,
) -> Tuple[NDArray, NDArray]:
    """Dispatch to :func:`compute_msc` or :func:`compute_imcoh`.

    Parameters
    ----------
    X : NDArray
        Time series ``(N, L)``.
    fs : float
        Sampling rate in Hz.
    coh_method : {"msc", "imcoh"}
        Which coherence variant to compute.  Magnitude and squared-magnitude
        ImCoh transforms are applied at load time by
        :func:`lrg_eegfc.workflow.fc.load_fc_matrix`, not here.
    nperseg, noverlap, batch_size
        Welch CSD parameters.

    Returns
    -------
    (freqs, Coh)
        Frequency grid and ``(N, N, F)`` coherence tensor.
    """
    if coh_method == "msc":
        return compute_msc(
            X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
        )
    if coh_method == "imcoh":
        return compute_imcoh(
            X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
        )
    raise ValueError(
        f"Unknown coh_method {coh_method!r}; choose 'msc' or 'imcoh'"
    )
