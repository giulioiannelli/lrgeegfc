"""Split-half connectivity matrices from a single continuous recording.

A recording is cut in two contiguous halves and each half is turned into a
band-resolved adjacency with the **same** estimator that produced the
full-duration matrix. Two things need this:

- a **within-condition baseline** for any cross-phase statistic (the two halves
  of ``rest_pre`` supply the ``A`` / ``B`` arms of the symmetric split-half
  estimator, so a cross-phase correlation is referenced to the estimator's own
  test-retest rather than to zero);
- the **test-retest reliability of the estimator itself** -- a property of the
  connectivity measure that is independent of any downstream hypothesis.

The Welch segment length should normally be halved relative to the
full-duration run so each half keeps a comparable number of segments (the
callers pass ``nperseg_for_fs(fs) // 2``).

Transform convention (order of operations matters, by Jensen's inequality):
the requested magnitude transform is applied **per frequency bin first**, then
averaged over the band -- identical to
:func:`lrg_eegfc.workflow.fc.load_fc_matrix`, so a half matrix is directly
comparable to a cached full-phase one.
"""
from __future__ import annotations

from typing import Dict, Mapping, Tuple

import numpy as np
from numpy.typing import NDArray

from .coherence.imcoh import compute_imcoh

__all__ = ["band_transform_signed", "imcoh_split_half_adjacencies"]

#: Magnitude transforms of a signed imaginary coherency, keyed as in
#: ``lrg_eegfc.workflow.fc.FC_METHOD_SHORTCUTS``.
_TRANSFORMS = {
    "identity": lambda x: x,
    "abs": np.abs,
    "sq": lambda x: x ** 2,
}


def band_transform_signed(
    Coh_signed: NDArray,
    freqs: NDArray,
    band: Tuple[float, float],
    transform: str = "abs",
) -> NDArray | None:
    """Band-average a signed ``(N, N, F)`` coherency under a magnitude transform.

    Returns a symmetric, zero-diagonal ``(N, N)`` matrix, or ``None`` if the
    band contains no frequency bin. ``transform`` is one of ``"identity"``,
    ``"abs"`` (``⟨|Im C|⟩_f``, Ewald 2012 / Bastos-Schoffelen 2016) or ``"sq"``
    (``⟨(Im C)²⟩_f``, Ewald 2012).
    """
    try:
        fn = _TRANSFORMS[transform]
    except KeyError:                                   # pragma: no cover
        raise ValueError(
            f"Unknown transform {transform!r}; expected one of {list(_TRANSFORMS)}"
        ) from None
    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not mask.any():
        return None
    A = fn(Coh_signed[:, :, mask]).mean(axis=-1)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return np.asarray(A, float)


def imcoh_split_half_adjacencies(
    X: NDArray,
    fs: float,
    nperseg: int,
    bands: Mapping[str, Tuple[float, float]],
    transform: str = "abs",
) -> Dict[Tuple[str, str], NDArray]:
    """``{(band, 'A'|'B'): adjacency}`` for the two contiguous halves of ``X``.

    ``X`` is ``(N_channels, T_samples)`` and is split along axis 1 at ``T // 2``;
    ``'A'`` is the first half, ``'B'`` the second. One Welch pass per half is
    shared across all bands.

    Generalises the previous script-local ``compute_imcoh_abs_halves`` over the
    magnitude transform, so the ``abs`` / ``sq`` substrate comparison uses one
    code path rather than two.
    """
    T = X.shape[1]
    out: Dict[Tuple[str, str], NDArray] = {}
    for tag, Xh in (("A", X[:, : T // 2]), ("B", X[:, T // 2:])):
        freqs, Coh = compute_imcoh(np.ascontiguousarray(Xh), fs, nperseg=nperseg)
        for band, lohi in bands.items():
            A = band_transform_signed(Coh, freqs, lohi, transform)
            if A is not None:
                out[(band, tag)] = A
        del Coh
    return out
