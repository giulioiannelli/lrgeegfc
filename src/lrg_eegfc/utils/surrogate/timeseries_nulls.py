"""Timeseries- and coherency-level surrogate nulls for coherence-based FC.

Where :mod:`lrg_eegfc.utils.surrogate.matched_strength` shuffles the **finished**
``N x N`` adjacency -- and therefore cannot test the connectivity estimator, the
band transform, the frequency-band split, session nonstationarity, drift, or the
phase segmentation -- everything here acts **upstream** of the adjacency: on the
per-segment Fourier coefficients, or on the complex band coherency ``C(f)``.

This module deliberately does **not** rebuild the archived matrix-level
coherency-rotation nulls (CRC / SB-CRC, see the note at the top of
:mod:`lrg_eegfc.utils.surrogate.coherency_surrogate`). Those rotated a finished
matrix by an opaque Haar transform. The constructions here are elementary
signal-processing operations with an explicit invariance statement each.

The four families
-----------------
``lag_randomized_coherency``            (N1)
    One circular time shift ``Delta_c`` per channel, applied as its exact
    frequency-domain equivalent: ``C_ij(f) -> C_ij(f) * exp(-2*pi*i*f*(Delta_i -
    Delta_j)/fs)``. Preserves ``|C_ij(f)|`` and every channel's power spectrum
    **exactly**; randomizes only the cross-spectral phase, i.e. the lag.

``segment_shifted_csd``                 (N1b)
    A genuine time-domain circular shift, quantized to the Welch segment
    lattice: channel ``c``'s Welch segment index is rolled by ``m_c``, so the
    surrogate is built from the *same* segments in a rotated order. Each
    channel's spectrum is preserved exactly; cross-channel *temporal
    correspondence* is destroyed. This is the classical **shift predictor** of
    Perkel, Gerstein & Moore (1967) carried into the spectral domain.

``phase_randomized_coherency``          (N2)
    Independent uniform phase per (channel, frequency):
    ``C_ij(f) -> C_ij(f) * exp(i(theta_i(f) - theta_j(f)))``. The frequency-domain
    form of the Theiler et al. (1992) Fourier-transform surrogate restricted to
    the cross-spectrum. Preserves ``|C_ij(f)|`` and per-channel spectra exactly;
    destroys lag **and** cross-frequency phase consistency.

``block_partition_indices``             (N3 / N4)
    Reassignment of contiguous session blocks to pseudo-phases, in a free and an
    order-preserving variant. Combined with :func:`csd_from_segment_subset` this
    yields the only null in the module that reaches session nonstationarity,
    slow drift, artifact epochs and the phase segmentation itself.

The invariance that makes N1/N2 well-posed downstream
------------------------------------------------------
Randomizing the cross-spectral phase drives ``<|Im C|>_f`` toward
``(2/pi) <|C|>_f`` (for a phase uniform on the circle, ``E|sin| = 2/pi``). On
sEEG, where volume conduction makes ``|Im C| << |C|``, that is a several-fold
**inflation** of the edge weights -- the same arithmetic that invalidated the
``complex_eigbasis`` CRC variant. It does **not** invalidate N1/N2 here, because
the readout is invariant to a global positive rescale of ``W``:

* the ``mst_union_top_fraction`` backbone is rank-based, so ``W -> cW`` keeps the
  identical edge set;
* the dimensionless scale ``s = tau * lambda_max`` means ``tau L`` is invariant
  under ``W -> cW`` (``L -> cL``, ``lambda_max -> c lambda_max``, ``tau -> tau/c``);
* ``rho_sym`` is a Spearman statistic on cophenetic differences.

Use :func:`global_rescale_invariance_residual` to verify this on real data rather
than trusting the argument. What N1/N2 therefore test is precise and worth
stating in those terms: *does the coherence-magnitude structure* ``|C|`` *alone,
with the lag structure destroyed, reproduce the cross-phase trace?*

What none of these can reject
-----------------------------
N1/N2 hold ``|C_ij(f)|`` fixed at its observed, phase-specific value, so they
cannot reject an explanation in which the trace lives in cross-phase changes of
coherence **magnitude**. N1b and N3/N4 do reach that, at the price of being
coarser. Read the ladder, not any single rung.

References
----------
Perkel, D. H., Gerstein, G. L. & Moore, G. P. (1967) "Neuronal spike trains and
stochastic point processes. II. Simultaneous spike trains." *Biophys. J.*
7(4):419-440. -- the shift predictor.

Theiler, J., Eubank, S., Longtin, A., Galdrikian, B. & Farmer, J. D. (1992)
"Testing for nonlinearity in time series: the method of surrogate data."
*Physica D* 58(1-4):77-94. -- Fourier-transform (phase-randomized) surrogates.

Nolte, G. et al. (2004) *Clin. Neurophysiol.* 115(10):2292-2307. -- imaginary
coherency. Ewald, A. et al. (2012) *NeuroImage* 60(1):476-488. -- ``|ImCoh|``.

Welch, P. D. (1967) *IEEE Trans. Audio Electroacoust.* 15(2):70-73.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.signal import get_window
from numpy.lib.stride_tricks import sliding_window_view

__all__ = [
    "segment_fft",
    "csd_from_segment_subset",
    "coherency_from_csd",
    "imcoh_abs_from_coherency",
    "band_bin_frequencies",
    "lag_randomized_coherency",
    "phase_randomized_coherency",
    "segment_shifted_csd",
    "block_partition_indices",
    "blocks_to_segment_sets",
    "global_rescale_invariance_residual",
]


# --------------------------------------------------------------------------- #
# Segment-level spectral primitives
#
# The Welch CSD is a *sum over segments*. Materializing the per-segment Fourier
# coefficients once turns every downstream surrogate that reorders, subsets, or
# rephases segments into cheap linear algebra instead of a fresh FFT pass. This
# is the enabling step for N1b and N3.
# --------------------------------------------------------------------------- #
def segment_fft(
    X: NDArray,
    fs: float,
    nperseg: int,
    band: tuple[float, float] | None = None,
    noverlap: int | None = None,
    dtype: type = np.complex64,
) -> tuple[NDArray, NDArray, float]:
    """Per-segment Hann-windowed rFFT coefficients, optionally band-restricted.

    Segmentation is **identical** to :func:`lrg_eegfc.utils.fc.coherence._common.welch_csd`
    (Hann window, ``noverlap = nperseg // 2`` by default, segments taken by a
    strided view), so summing the returned outer products over all segments
    reproduces that function's CSD exactly up to the scale factor returned here.

    Parameters
    ----------
    X : (N, L) ndarray
        Channels-first timeseries.
    fs : float
        Sampling rate (Hz).
    nperseg : int
        Welch segment length in samples.
    band : (fmin, fmax) or None
        If given, keep only bins with ``fmin <= f <= fmax`` (inclusive, matching
        the band mask used throughout the FC code).
    noverlap : int or None
        Defaults to ``nperseg // 2``.
    dtype : numpy complex dtype
        ``complex64`` by default -- the coefficients feed a normalized ratio, and
        halving the footprint is what keeps a whole session resident.

    Returns
    -------
    (freqs, F, scale)
        ``freqs`` (F_b,) the kept frequencies; ``F`` (N, n_seg, F_b) coefficients;
        ``scale`` the Welch normalization ``1 / (fs * sum(w^2))`` **without** the
        ``1 / n_segments`` factor (which depends on the subset and is applied by
        :func:`csd_from_segment_subset`).
    """
    N, L = X.shape
    if noverlap is None:
        noverlap = nperseg // 2
    step = nperseg - noverlap
    segments = sliding_window_view(X, window_shape=nperseg, axis=1)[:, ::step, :]
    n_seg = segments.shape[1]
    if n_seg == 0:
        raise ValueError("nperseg larger than signal length")

    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    window = get_window("hann", nperseg)
    scale = 1.0 / (fs * float((window * window).sum()))
    # Match the window to X's dtype: a float64 window silently upcasts a float32
    # recording mid-multiply and doubles the transient footprint. The scale
    # factor is computed in float64 above, so precision of the normalization is
    # unaffected.
    window = window.astype(X.dtype, copy=False)

    mask = slice(None)
    if band is not None:
        m = (freqs >= band[0]) & (freqs <= band[1])
        if not m.any():
            raise ValueError(f"no rFFT bins inside band {band} at nperseg={nperseg}")
        mask = m
        freqs = freqs[m]

    out = np.empty((N, n_seg, freqs.size), dtype=dtype)
    for start in range(0, n_seg, 64):                     # batch to bound peak RAM
        stop = min(start + 64, n_seg)
        blk = segments[:, start:stop, :] * window[None, None, :]
        out[:, start:stop, :] = np.fft.rfft(blk, n=nperseg, axis=2)[:, :, mask]
    return freqs, out, scale


def csd_from_segment_subset(
    F: NDArray, idx: NDArray | slice | None = None, scale: float = 1.0
) -> NDArray:
    """Welch CSD ``(F_b, N, N)`` from a subset of the per-segment coefficients.

    ``S(f) = scale / n_sub * sum_{k in idx} F[:, k, f] F[:, k, f]^H``. Evaluated
    as one complex GEMM per frequency, which is where the speed comes from.

    ``idx`` may also be an ``(N, n_sub)`` integer array giving a **per-channel**
    segment selection -- that is how :func:`segment_shifted_csd` realizes a
    channel-wise circular shift on the segment lattice.
    """
    if idx is None:
        idx = slice(None)
    if isinstance(idx, np.ndarray) and idx.ndim == 2:
        N, n_sub = idx.shape
        G = np.take_along_axis(F, idx[:, :, None], axis=1)          # (N, n_sub, F_b)
    else:
        G = F[:, idx, :]
        N, n_sub = G.shape[0], G.shape[1]
    if n_sub == 0:
        raise ValueError("empty segment subset")
    Gt = np.ascontiguousarray(np.moveaxis(G, 2, 0))                 # (F_b, N, n_sub)
    S = Gt @ np.conj(np.swapaxes(Gt, 1, 2))                         # (F_b, N, N)
    return S * (scale / n_sub)


def coherency_from_csd(S: NDArray) -> NDArray:
    """``C(f) = D^{-1/2} S(f) D^{-1/2}``, ``D = diag(S)``. Input/output ``(F_b, N, N)``.

    Same normalization as
    :func:`lrg_eegfc.utils.surrogate.coherency_surrogate.complex_coherency_band`
    and :func:`lrg_eegfc.utils.fc.coherence.imcoh.compute_imcoh`, so the derived
    ``<|Im C|>_f`` reproduces the cached ``imcoh_abs`` adjacency.
    """
    P = np.real(np.diagonal(S, axis1=1, axis2=2))                   # (F_b, N)
    den = np.sqrt(P[:, :, None] * P[:, None, :])
    return np.divide(S, den, out=np.zeros_like(S), where=den > 0)


def imcoh_abs_from_coherency(C: NDArray) -> NDArray:
    """``W = <|Im C(f)|>_f`` -- the ``imcoh_abs`` adjacency, cleaned for LRG.

    Symmetric, zero-diagonal, clipped to ``[0, 1]``: identical post-processing to
    the ``load_phase`` used by the gate scripts, so observed and surrogate
    adjacencies enter the Laplacian through the same door.
    """
    W = np.abs(np.imag(C)).mean(axis=0).astype(np.float64)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def band_bin_frequencies(fs: float, nperseg: int, band: tuple[float, float]) -> NDArray:
    """In-band rFFT bin frequencies -- the ``f`` axis the phase ramps act on."""
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    return freqs[(freqs >= band[0]) & (freqs <= band[1])]


# --------------------------------------------------------------------------- #
# N1 -- lag randomization by per-channel circular shift (frequency domain)
# --------------------------------------------------------------------------- #
def lag_randomized_coherency(
    C: NDArray,
    freqs: NDArray,
    fs: float,
    rng: np.random.Generator,
    nperseg: int | None = None,
) -> NDArray:
    """One N1 surrogate coherency: per-channel circular shift as a phase ramp.

    Draws ``Delta_c`` uniformly on ``[0, nperseg)`` samples and applies

        ``C_ij(f) -> C_ij(f) * exp(-2*pi*i*f*Delta_i/fs) * exp(+2*pi*i*f*Delta_j/fs)``

    which is the exact frequency-domain image of a circular time shift (the DFT
    shift theorem). ``|C_ij(f)|`` is preserved to machine precision and the
    diagonal stays real and unit, so the surrogate is a genuine coherency.

    Drawing on ``[0, nperseg)`` rather than ``[0, L)`` is not an approximation:
    the ramp ``exp(-2*pi*i*f*Delta/fs)`` evaluated on the rFFT grid
    ``f = k*fs/nperseg`` is **periodic in Delta with period nperseg**, so shifts
    differing by a multiple of ``nperseg`` are the identical surrogate.

    Parameters
    ----------
    C : (F_b, N, N) complex
        Observed band coherency.
    freqs : (F_b,) ndarray
        The matching bin frequencies (see :func:`band_bin_frequencies`).
    fs : float
        Sampling rate (Hz).
    rng : numpy Generator
    nperseg : int or None
        Shift support. Defaults to ``round(fs / (freqs[1] - freqs[0]))`` inferred
        from the bin spacing.

    Returns
    -------
    (F_b, N, N) complex surrogate coherency.
    """
    N = C.shape[1]
    if nperseg is None:
        if freqs.size < 2:
            raise ValueError("cannot infer nperseg from a single frequency bin")
        nperseg = int(round(fs / (freqs[1] - freqs[0])))
    delta = rng.integers(0, nperseg, size=N)
    ramp = np.exp(-2j * np.pi * np.outer(freqs, delta) / fs)        # (F_b, N)
    return C * ramp[:, :, None] * np.conj(ramp)[:, None, :]


# --------------------------------------------------------------------------- #
# N2 -- independent per-(channel, frequency) phase randomization
# --------------------------------------------------------------------------- #
def phase_randomized_coherency(C: NDArray, rng: np.random.Generator) -> NDArray:
    """One N2 surrogate coherency: iid uniform phase per (channel, frequency).

    ``C_ij(f) -> C_ij(f) * exp(i(theta_i(f) - theta_j(f)))`` with
    ``theta ~ U[0, 2*pi)`` independent over channels and frequencies -- the
    cross-spectral restriction of the Theiler et al. (1992) FT surrogate.
    ``|C_ij(f)|`` and every channel's power spectrum are preserved exactly.

    Relative to :func:`lag_randomized_coherency`, N2 has the **same first moment**
    for ``<|Im C|>_f`` (both send the phase to uniform) but far smaller
    realization variance, because N1's phase is a single coherent ramp shared
    across the whole band while N2's is independent per bin and averages down.
    N2 is therefore the low-variance limit of N1, not an independent rung -- say
    so when reporting them together.
    """
    Fb, N, _ = C.shape
    theta = rng.uniform(0.0, 2.0 * np.pi, size=(Fb, N))
    ph = np.exp(1j * theta)
    return C * ph[:, :, None] * np.conj(ph)[:, None, :]


# --------------------------------------------------------------------------- #
# N1b -- genuine time-domain circular shift on the Welch segment lattice
# --------------------------------------------------------------------------- #
def segment_shifted_csd(
    F: NDArray, rng: np.random.Generator, scale: float = 1.0,
    min_shift: int = 1,
) -> NDArray:
    """One N1b surrogate CSD: per-channel circular shift of the segment index.

    Channel ``c``'s segment sequence is rolled by ``m_c`` (uniform on
    ``[min_shift, n_seg)``), i.e. the surrogate is assembled from exactly the same
    Welch segments, in a channel-wise rotated order. Consequences, all exact:

    * each channel's Welch power spectrum is **unchanged** (same multiset of
      segments);
    * cross-channel *temporal correspondence* is destroyed, so genuine coupling
      -- magnitude included -- averages down toward the ``~1/sqrt(n_seg)`` noise
      floor.

    This is the spectral form of the shift predictor (Perkel, Gerstein & Moore
    1967). Unlike N1/N2 it is **not** magnitude-preserving: it is an independence
    null, and the report must label it as such rather than as "a circular shift"
    tout court, because for a Welch-averaged estimator those are different
    animals.

    Parameters
    ----------
    F : (N, n_seg, F_b) complex
        Per-segment coefficients from :func:`segment_fft`.
    rng : numpy Generator
    scale : float
        Welch scale from :func:`segment_fft`.
    min_shift : int
        Smallest allowed roll. ``1`` is the natural floor (``0`` is the observed
        data); with 50%-overlapping segments a roll of 1 still shares half the
        samples, so callers wanting a clean separation should pass a larger
        value.

    Returns
    -------
    (F_b, N, N) complex surrogate CSD.
    """
    N, n_seg, _ = F.shape
    if n_seg <= min_shift:
        raise ValueError(f"n_seg={n_seg} too small for min_shift={min_shift}")
    m = rng.integers(min_shift, n_seg, size=N)
    idx = (np.arange(n_seg)[None, :] + m[:, None]) % n_seg          # (N, n_seg)
    return csd_from_segment_subset(F, idx.astype(np.intp), scale)


# --------------------------------------------------------------------------- #
# N3 / N4 -- block phase-label permutation
# --------------------------------------------------------------------------- #
def block_partition_indices(
    n_blocks: int,
    sizes: dict[str, int],
    rng: np.random.Generator,
    mode: str = "free",
    rotation: int | None = None,
) -> dict[str, NDArray]:
    """Assign contiguous session blocks to pseudo-phases of prescribed sizes.

    Parameters
    ----------
    n_blocks : int
        Number of contiguous blocks the concatenated session was cut into.
    sizes : dict label -> int
        Blocks per pseudo-phase, **in the true temporal order of the labels**.
        ``sum(sizes.values())`` must equal ``n_blocks``, so every surrogate
        reproduces the true per-phase durations (up to block quantization) and
        the comparison to the observed statistic is never confounded by differing
        segment counts / spectral degrees of freedom.
    rng : numpy Generator
    mode : {"free", "order_preserving", "identity"}
        ``free``
            Uniformly random assignment of blocks to labels. Destroys temporal
            order entirely, hence any session structure including drift and the
            clustering of artifact epochs.
        ``order_preserving``
            **Circular rotation of the phase labels along the session.** The
            session is treated as a ring, rotated by ``rotation`` blocks, and then
            cut at the true boundaries. Every pseudo-phase stays a contiguous
            interval of the true length, the local temporal order is preserved
            everywhere except at the single seam, and any smooth session-level
            trend is *retained* rather than removed. Only the placement of the
            cut points moves.

            The realization set is finite and small (``n_blocks`` rotations), so
            the honest way to use it is to **enumerate all rotations** and form an
            exact test, rather than to sample ``R`` of them with replacement. The
            runner does this.

            This is **not** the retired drift null. That one asked whether a
            monotone session trend could manufacture a *distance asymmetry*, and
            was degenerate because a directional task makes drift and trace the
            same object. This rung asks a different, non-degenerate question:
            holding the session, its drift, and the phase durations fixed, does
            the *specific placement* of the real phase boundaries matter?
        ``identity``
            The true partition. Used to compute the observed statistic through the
            identical code path, so the null is attached to the right observed
            value.
    rotation : int or None
        Rotation offset for ``order_preserving``. If ``None`` one is drawn from
        ``rng``; pass an explicit value to enumerate the group exhaustively.

    Returns
    -------
    dict label -> int array of block indices.
    """
    labels = list(sizes)
    total = int(sum(sizes[k] for k in labels))
    if total != n_blocks:
        raise ValueError(f"sizes sum to {total}, expected n_blocks={n_blocks}")

    if mode == "identity":
        order = np.arange(n_blocks)
    elif mode == "free":
        order = rng.permutation(n_blocks)
    elif mode == "order_preserving":
        r = int(rng.integers(0, n_blocks)) if rotation is None else int(rotation) % n_blocks
        order = np.roll(np.arange(n_blocks), -r)
    else:
        raise ValueError(f"unknown mode {mode!r}")

    out, k = {}, 0
    for lab in labels:
        n = int(sizes[lab])
        out[str(lab)] = order[k:k + n]
        k += n
    return out


def blocks_to_segment_sets(
    block_of_segment: NDArray, block_idx: dict[str, NDArray]
) -> dict[str, NDArray]:
    """Map a block-level assignment to the Welch segment indices it selects.

    ``block_of_segment`` is the ``(n_seg,)`` block label of each segment, with
    ``-1`` marking segments to discard (those straddling a block or phase
    boundary -- the runner marks them so no segment ever mixes two blocks).
    """
    n_blk = int(block_of_segment.max()) + 1
    by_block = [np.flatnonzero(block_of_segment == b) for b in range(n_blk)]
    out = {}
    for lab, idx in block_idx.items():
        parts = [by_block[int(b)] for b in idx]
        out[lab] = (np.concatenate(parts) if parts else
                    np.array([], dtype=np.intp)).astype(np.intp)
    return out


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def global_rescale_invariance_residual(readout, W_by_phase: dict, c: float = 7.0) -> float:
    """Max absolute change in ``readout`` when every phase's ``W`` is scaled by ``c``.

    The N1/N2 surrogates inflate ``<|Im C|>`` by a roughly constant factor
    (``2/pi * <|C|> / <|Im C|>``), which is only harmless if the readout is
    invariant to a global positive rescale. This checks that claim numerically
    instead of asserting it. ``readout`` maps ``{phase: W}`` to an array.
    """
    a = np.asarray(readout(W_by_phase), dtype=float)
    b = np.asarray(readout({k: c * v for k, v in W_by_phase.items()}), dtype=float)
    return float(np.nanmax(np.abs(a - b)))
