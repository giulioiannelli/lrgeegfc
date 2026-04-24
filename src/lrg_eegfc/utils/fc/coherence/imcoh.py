"""Signed imaginary part of coherency (Nolte et al. 2004).

Imaginary coherency ImC is defined from the complex cross-spectral
density ``S_ij(f)`` as::

    ImC_ij(f) = Im(S_ij(f)) / sqrt(S_ii(f) * S_jj(f))

and lies in ``[-1, 1]``.  Because instantaneous (zero-phase) mixing
produces a purely real CSD, ImC vanishes exactly under volume
conduction; non-zero values indicate a non-zero-phase-lag interaction.
The sign of ``ImC_ij`` encodes lead/lag: ``ImC_ji = -ImC_ij``
(skew-symmetric).

This module stores **only the signed form** on disk.  The downstream
loader (:func:`lrg_eegfc.workflow.fc.load_fc_matrix`) derives the
magnitude ``|ImC|`` (Ewald 2012 / Bastos-Schoffelen 2016 convention) and
the squared magnitude ``|ImC|^2`` on the fly via ``np.abs(raw)`` and
``raw ** 2``; they are **never** recomputed from the time series or
cached separately.

References
----------
- Nolte, G., Bai, O., Wheaton, L., Mari, Z., Vorbach, S., Hallett, M.
  (2004) "Identifying true brain interaction from EEG data using the
  imaginary part of coherency." Clin. Neurophysiol. 115(10):2292-2307.
- Ewald, A., Marzetti, L., Zappasodi, F., Meinecke, F. C., Nolte, G.
  (2012) "Estimating true brain connectivity from EEG/MEG data invariant
  to linear and static transformations in sensor space." NeuroImage
  60(1):476-488.
- Bastos, A. M., Schoffelen, J.-M. (2016) "A tutorial review of
  functional connectivity analysis methods and their interpretational
  pitfalls." Front. Syst. Neurosci. 9:175.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from ._common import welch_csd


__all__ = ["compute_imcoh"]


def compute_imcoh(
    X: NDArray,
    fs: float,
    nperseg: int = 256,
    noverlap: int | None = None,
    batch_size: int = 64,
) -> Tuple[NDArray, NDArray]:
    """Compute the signed imaginary part of coherency.

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
    ImCoh : NDArray
        Signed ``(N, N, F)`` imaginary coherency in ``[-1, 1]``.
        Skew-symmetric: ``ImCoh[j, i, :] = -ImCoh[i, j, :]``; diagonal is
        zero.  Derived magnitude/square (``|ImCoh|`` or ``|ImCoh|^2``)
        are applied at load time, not here.
    """
    N = X.shape[0]
    freqs, CSD = welch_csd(
        X, fs, nperseg=nperseg, noverlap=noverlap, batch_size=batch_size
    )
    F = freqs.shape[0]
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)  # (N, F)
    denom_sqrt = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    # Canonical signed ImCoh (Nolte 2004): Im(S_ij) / sqrt(S_ii * S_jj)
    ImCoh = np.divide(
        np.imag(CSD),
        denom_sqrt,
        out=np.zeros((N, N, F)),
        where=denom_sqrt > 0,
    )
    np.einsum("iif->if", ImCoh)[...] = 0.0  # self-imcoh = 0
    return freqs, ImCoh
