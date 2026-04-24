"""Backwards-compatibility shim for the legacy ``compute_msc_welch`` API.

The actual math now lives in :mod:`lrg_eegfc.utils.fc.coherence`.  This
module re-exports the shared Welch CSD primitives and dispatches the
``metric`` keyword (``"msc" | "imcoh" | "wpli"``) to the appropriate
implementation.  ``msc`` and ``imcoh`` are handled by the new coherence
subpackage; ``wpli`` retains its original per-segment Vinck-2011
accumulator because it cannot be derived from the batch-averaged CSD.

All ~20 downstream callers continue to work unchanged.  New code should
import from :mod:`lrg_eegfc.utils.fc.coherence` directly.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.signal import get_window
from numpy.lib.stride_tricks import sliding_window_view

from ..coherence import compute_imcoh, compute_msc
from ..coherence._common import band_average as _band_average


__all__ = [
    "compute_msc_welch",
    "band_average_msc",
]


#: Back-compat alias: ``band_average_msc`` is a thin wrapper over the
#: generic band-averaging helper now exposed from ``coherence._common``.
band_average_msc = _band_average


def _compute_wpli_welch(
    X: NDArray,
    fs: float,
    nperseg: int,
    noverlap: int | None,
    batch_size: int,
) -> Tuple[NDArray, NDArray]:
    """Weighted Phase-Lag Index (Vinck et al. 2011).

    wPLI requires per-segment cross-spectra so it cannot share the
    batch-averaged CSD path; this function preserves the original
    segment-wise accumulator.
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

    wpli_num = np.zeros((N, N, F))
    wpli_den = np.zeros((N, N, F))

    for start in range(0, n_segments, batch_size):
        stop = start + batch_size
        block = segments[:, start:stop, :] * window[None, None, :]
        fft_block = np.fft.rfft(block, n=nperseg, axis=2)
        B_actual = fft_block.shape[1]
        for b in range(B_actual):
            seg = fft_block[:, b, :]
            cross_seg = seg[:, None, :] * np.conj(seg[None, :, :])
            im_cross = np.imag(cross_seg)
            wpli_num += np.abs(im_cross) * np.sign(im_cross)
            wpli_den += np.abs(im_cross)

    Coh = np.divide(
        np.abs(wpli_num), wpli_den, out=np.zeros((N, N, F)), where=wpli_den > 0
    )
    np.einsum("iif->if", Coh)[...] = 0.0
    return freqs, Coh


def compute_msc_welch(
    X: NDArray,
    fs: float,
    nperseg: int = 256,
    noverlap: int | None = None,
    batch_size: int = 64,
    metric: str = "msc",
) -> Tuple[NDArray, NDArray]:
    """Back-compat entry point that dispatches to the new coherence subpackage.

    Parameters
    ----------
    metric : {"msc", "imcoh", "wpli"}
        - ``"msc"``   → :func:`coherence.compute_msc` (range ``[0, 1]``).
        - ``"imcoh"`` → :func:`coherence.compute_imcoh` (signed, ``[-1, 1]``,
          Nolte 2004). Magnitude and squared-magnitude transforms are
          derived at load time by the FC loader.
        - ``"wpli"``  → local Vinck 2011 implementation (range ``[0, 1]``).

    Returns
    -------
    (freqs, Coh)
        See the individual metric modules for exact semantics.
    """
    if metric == "msc":
        return compute_msc(
            X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
        )
    if metric == "imcoh":
        return compute_imcoh(
            X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
        )
    if metric == "wpli":
        return _compute_wpli_welch(
            X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
        )
    raise ValueError(
        f"Unknown metric {metric!r}; choose 'msc', 'imcoh', or 'wpli'"
    )
