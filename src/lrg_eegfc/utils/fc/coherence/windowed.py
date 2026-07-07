"""Window-level |ImCoh| from segment FFTs (time-resolved functional connectivity).

This module lets a caller compute ``|ImCoh|`` on arbitrary *subsets* of a single
window's Welch segments (e.g. even/odd interleaved halves for split-half
reliability, or jackknife) **without recomputing the per-segment FFTs**. It is the
time-resolved backend for sliding-window connectivity / replay-state analyses.

The recipe — per-frequency-bin ``|Im(coherency)|`` first, then average over the
band — matches the canonical :func:`lrg_eegfc.workflow.fc.load_fc_matrix`
``imcoh_abs`` transform exactly. Verified continuity anchor: a full-window pass
through these helpers reproduces the cached freq-resolved ``imcoh_abs`` matrix to
Spearman/Pearson ``1.0`` (``audit_124``). So a windowed measure built on this
backend reduces to the canonical full-phase FC in the full-window limit.

References
----------
- Nolte et al. (2004) Clin. Neurophysiol. 115(10):2292-2307 (imaginary coherency).
- Ewald et al. (2012) NeuroImage 60(1):476-488; Bastos & Schoffelen (2016)
  Front. Syst. Neurosci. 9:175 (``|ImCoh|`` magnitude convention).
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from numpy.typing import NDArray
from scipy.signal import get_window

__all__ = [
    "segment_ffts",
    "imcoh_abs_cube",
    "imcoh_signed_cube",
    "band_abs_average",
    "band_signed_average",
    "windowed_imcoh_abs",
]


def segment_ffts(
    x: NDArray,
    fs: float,
    nperseg: int,
    noverlap: Optional[int] = None,
    fmax_keep: Optional[float] = None,
) -> Tuple[NDArray, NDArray]:
    """Per-segment Hann-windowed rFFT of a time window ``x = (N, L)``.

    No Welch averaging is performed — the caller forms the cross-spectral density
    from whichever segment subset it wants. Frequencies above ``fmax_keep`` are
    discarded to save memory (default: keep all up to Nyquist).

    Returns
    -------
    freqs : (F,) ndarray
    ffts : (N, n_seg, F) complex ndarray
    """
    if noverlap is None:
        noverlap = nperseg // 2
    step = nperseg - noverlap
    seg = sliding_window_view(x, window_shape=nperseg, axis=1)[:, ::step, :]
    win = get_window("hann", nperseg).astype(float)
    seg = seg * win[None, None, :]
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    if fmax_keep is not None:
        keep = freqs <= fmax_keep
        freqs = freqs[keep]
        ff = np.fft.rfft(seg, n=nperseg, axis=2)[:, :, keep]
    else:
        ff = np.fft.rfft(seg, n=nperseg, axis=2)
    return freqs, ff


def imcoh_abs_cube(ff: NDArray, seg_idx: Optional[NDArray] = None) -> Optional[NDArray]:
    """``|ImCoh|(f)`` cube ``(N, N, F)`` from a subset of segment FFTs.

    ``ff`` is the ``(N, n_seg, F)`` output of :func:`segment_ffts`. ``seg_idx``
    selects which segments to average over (default: all). Returns ``None`` if
    fewer than 2 segments are selected (coherency is degenerate with <2 segments).
    """
    if seg_idx is not None:
        ff = ff[:, seg_idx, :]
    nseg = ff.shape[1]
    if nseg < 2:
        return None
    CSD = np.einsum("nbf,mbf->nmf", ff, np.conj(ff), optimize="greedy") / nseg
    PSD = np.real(np.einsum("nnf->nf", CSD))                       # (N, F)
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    imcoh = np.divide(np.imag(CSD), denom, out=np.zeros_like(denom), where=denom > 0)
    return np.abs(imcoh)                                           # (N, N, F)


def imcoh_signed_cube(ff: NDArray, seg_idx: Optional[NDArray] = None) -> Optional[NDArray]:
    """SIGNED ``ImCoh(f)`` cube ``(N, N, F)`` — the lead/lag direction kept.

    Identical to :func:`imcoh_abs_cube` but WITHOUT the final ``np.abs`` — the
    returned cube is real and antisymmetric in the node axes per frequency bin
    (``Im(CSD)_{ji} = -Im(CSD)_{ij}``), so it carries the phase-lead/lag direction
    that ``|ImCoh|`` discards. For directional (magnetic / signed-Laplacian)
    analyses. ``np.abs`` of this equals :func:`imcoh_abs_cube` bin-for-bin.
    """
    if seg_idx is not None:
        ff = ff[:, seg_idx, :]
    nseg = ff.shape[1]
    if nseg < 2:
        return None
    CSD = np.einsum("nbf,mbf->nmf", ff, np.conj(ff), optimize="greedy") / nseg
    PSD = np.real(np.einsum("nnf->nf", CSD))
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    return np.divide(np.imag(CSD), denom, out=np.zeros_like(denom), where=denom > 0)


def band_signed_average(
    cube: NDArray,
    freqs: NDArray,
    lo: float,
    hi: float,
    cleanup: bool = True,
) -> Optional[NDArray]:
    """Band-average a SIGNED ``ImCoh(f)`` cube to one antisymmetric ``(N, N)``.

    The net directional flow in the band. With ``cleanup`` (default) the result is
    antisymmetrized (``0.5*(A - Aᵀ)``) and its diagonal zeroed. Returns ``None`` if
    no bins fall in ``[lo, hi]``. NB: band-averaging signed values can cancel
    sub-band sign flips — this is the *net* in-band direction (lossy by design;
    same convention as the signed-flow / Hodge audits).
    """
    mask = (freqs >= lo) & (freqs <= hi)
    if not mask.any():
        return None
    A = cube[:, :, mask].mean(axis=-1)
    if cleanup:
        A = 0.5 * (A - A.T)
        np.fill_diagonal(A, 0.0)
    return A


def band_abs_average(
    cube: NDArray,
    freqs: NDArray,
    lo: float,
    hi: float,
    cleanup: bool = True,
) -> Optional[NDArray]:
    """Band-average a ``|ImCoh|(f)`` cube to a single ``(N, N)`` matrix.

    With ``cleanup`` (default) the result is symmetrized and its diagonal zeroed —
    the connectivity convention used everywhere downstream. Returns ``None`` if no
    frequency bins fall in ``[lo, hi]``.
    """
    mask = (freqs >= lo) & (freqs <= hi)
    if not mask.any():
        return None
    W = cube[:, :, mask].mean(axis=-1)
    if cleanup:
        W = 0.5 * (W + W.T)
        np.fill_diagonal(W, 0.0)
    return W


def windowed_imcoh_abs(
    x: NDArray,
    fs: float,
    lo: float,
    hi: float,
    nperseg: int,
    noverlap: Optional[int] = None,
) -> Optional[NDArray]:
    """Convenience: full ``|ImCoh|_band`` ``(N, N)`` for one time window ``x``.

    Equivalent to ``segment_ffts`` -> ``imcoh_abs_cube`` (all segments) ->
    ``band_abs_average``.
    """
    freqs, ff = segment_ffts(x, fs, nperseg, noverlap=noverlap, fmax_keep=hi + 5.0)
    cube = imcoh_abs_cube(ff)
    if cube is None:
        return None
    return band_abs_average(cube, freqs, lo, hi)
