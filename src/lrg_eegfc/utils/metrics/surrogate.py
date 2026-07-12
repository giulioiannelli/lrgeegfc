"""Strength-preserving surrogate graphs (matched-strength null).

The canonical FC cohort null: a **4-cycle +/-delta edge shuffle** that leaves
every node's weighted strength exactly invariant while randomising topology.
Picking two disjoint edges ``(a,b),(c,d)`` and shifting weight
``w_ab += dl, w_cd += dl, w_ad -= dl, w_cb -= dl`` conserves the strengths of
``a,b,c,d`` for any ``dl``; ``dl`` is drawn uniformly in the feasibility window
that keeps all four weights in ``[0, w_max]``.

The jitted loop draws its ``(samples, fracs)`` arrays *outside* numba so the
result is reproducible from a numpy ``Generator`` and **bit-identical** to the
original pure-Python implementation in ``audit_63``/``audit_150`` (byte-copied
here so promotion to the library does not perturb any established anchor).

See ``feedback_matched_strength_mandatory`` and ``feedback_numba_for_surrogates``.
"""
from __future__ import annotations

import numba
import numpy as np
from numpy.typing import NDArray

__all__ = ["matched_strength_shuffle"]


@numba.njit
def _swap_loop(W, s, fr, w_max):
    for i in range(s.shape[0]):
        a = s[i, 0]; b = s[i, 1]; c = s[i, 2]; d = s[i, 3]
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]; w2 = W[c, d]; w3 = W[a, d]; w4 = W[c, b]
        lo = -w1
        if -w2 > lo: lo = -w2
        if w3 - w_max > lo: lo = w3 - w_max
        if w4 - w_max > lo: lo = w4 - w_max
        hi = w_max - w1
        if w_max - w2 < hi: hi = w_max - w2
        if w3 < hi: hi = w3
        if w4 < hi: hi = w4
        if lo >= hi:
            continue
        dl = lo + fr[i] * (hi - lo)
        W[a, b] = w1 + dl; W[b, a] = w1 + dl
        W[c, d] = w2 + dl; W[d, c] = w2 + dl
        W[a, d] = w3 - dl; W[d, a] = w3 - dl
        W[c, b] = w4 - dl; W[b, c] = w4 - dl
    return W


def matched_strength_shuffle(W: NDArray, n_swaps: int, rng: np.random.Generator,
                             w_max: float = 1.0) -> NDArray:
    """One strength-preserving surrogate of ``W``.

    Parameters
    ----------
    W : (N, N) symmetric weighted adjacency in ``[0, w_max]``, zero diagonal.
    n_swaps : number of 4-cycle attempts (``SWAP_FACTOR * N(N-1)/2`` canonically).
    rng : numpy ``Generator``; both the vertex quadruples and the fractional
        shifts are drawn from it *before* the jitted loop (reproducible + fast).
    w_max : upper weight bound (``1.0`` for ``imcoh_abs``).

    Returns a **new** matrix (input not modified). Bit-identical to the
    pure-Python audit_63/audit_150 shuffle for the same ``rng`` state.
    """
    W = np.asarray(W, float).copy()
    N = W.shape[0]
    s = rng.integers(0, N, size=(n_swaps, 4)).astype(np.int64)
    fr = rng.uniform(0.0, 1.0, size=n_swaps)
    return _swap_loop(W, s, fr, float(w_max))
