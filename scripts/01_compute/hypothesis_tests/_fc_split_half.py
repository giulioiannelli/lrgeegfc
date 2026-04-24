"""Compute |ImCoh| per band on the two halves of a time-series recording.

Pure helper used by H2e (session-drift noise floor). No I/O. Given a
channel × time matrix, returns the six-band |ImCoh| adjacency matrices
for the first half and the second half, computed with a halved Welch
nperseg so that the spectral sample count stays comparable to the
full-duration run.

See plan §2 in
`/home/giulio/.claude/plans/i-just-added-home-giulio-documents-resea-stateless-swing.md`
for design rationale; see `.agents/guides/02_methods/H2_METRICS.md` for
the FC recipe (`imcoh_abs` = band-averaged |signed ImCoh|).
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch


__all__ = ["compute_imcoh_abs_halves"]


def _band_abs(Coh_signed: np.ndarray, freqs: np.ndarray,
              flo: float, fhi: float) -> np.ndarray | None:
    """Band-abs-average a signed ImCoh cube → (N, N) non-negative matrix.

    ``imcoh_abs`` convention: take |·| at each freq bin, then average
    over bins inside the band. Matches the loader transform applied to
    canonical freq-resolved caches.
    """
    mask = (freqs >= flo) & (freqs <= fhi)
    if not mask.any():
        return None
    return np.abs(Coh_signed[:, :, mask]).mean(axis=-1)


def compute_imcoh_abs_halves(
    X: np.ndarray,
    fs: float,
    nperseg: int,
    bands: dict[str, tuple[float, float]],
) -> dict[tuple[str, str], np.ndarray]:
    """Return {(band, 'A'|'B'): |ImCoh|_band} for both halves of X.

    ``X`` must be (N_channels, T_samples); it is split along axis 1 at
    ``T // 2``. ``nperseg`` should typically be halved vs full-duration
    (we call ``nperseg_for_fs(fs) // 2`` from the driver) so each half
    retains 2-s windows worth of segments.
    """
    T = X.shape[1]
    halves = {"A": X[:, : T // 2], "B": X[:, T // 2:]}
    out: dict[tuple[str, str], np.ndarray] = {}
    for tag, Xh in halves.items():
        freqs, Coh = compute_msc_welch(Xh, fs, nperseg=nperseg, metric="imcoh")
        for band, (flo, fhi) in bands.items():
            A_band = _band_abs(Coh, freqs, flo, fhi)
            if A_band is None:
                continue
            # Hermitian cleanup; zero diagonal (connectivity convention).
            A_band = 0.5 * (A_band + A_band.T)
            np.fill_diagonal(A_band, 0.0)
            out[(band, tag)] = A_band.astype(np.float64)
        del Coh
    return out
