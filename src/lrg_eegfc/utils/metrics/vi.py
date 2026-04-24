"""Variation of Information (Meilă 2007) between two partitions.

Canonical implementation used across every cross-phase analysis script in
this project. See `.agents/guides/02_methods/H2_METRICS.md` §1 for the
mathematical definition and for how VI fits into the H2 metric suite.

Prior to 2026-04-24 this function was re-implemented in ~17 scripts with
minor variations; all scripts now import from here.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


__all__ = ["compute_vi", "conditional_entropy"]


def _joint_histogram(labels1: np.ndarray, labels2: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Shared machinery for VI and conditional entropy.

    Returns (p1, p2, joint_counts, n) where ``joint_counts[i, j]`` counts
    items in cluster ``uniq1[i]`` of ``labels1`` and cluster ``uniq2[j]``
    of ``labels2``.
    """
    n = len(labels1)
    uniq1, inv1 = np.unique(labels1, return_inverse=True)
    uniq2, inv2 = np.unique(labels2, return_inverse=True)
    joint = np.zeros((len(uniq1), len(uniq2)), dtype=np.int64)
    np.add.at(joint, (inv1, inv2), 1)
    p1 = joint.sum(axis=1) / n
    p2 = joint.sum(axis=0) / n
    return p1, p2, joint, n


def conditional_entropy(labels_y: ArrayLike, labels_x: ArrayLike) -> float:
    """Conditional entropy H(Y | X) of partition Y given partition X.

    Uses natural log (nats), matching ``compute_vi``'s convention. H(Y|X)
    = 0 iff X is a refinement of Y (knowing X determines Y). Upper bound
    is ``log(n)``. This is the asymmetric half of VI:

        VI(Y, X) = H(Y | X) + H(X | Y)

    and in general ``H(Y|X) != H(X|Y)``. Used by the H2a′ directed
    partition-level memory test — see §1 of
    ``.agents/guides/02_methods/H2_METRICS.md``.

    Parameters
    ----------
    labels_y, labels_x
        Cluster labels for the same N items (same length required).

    Returns
    -------
    h_y_given_x : float
        Non-negative; clipped at 0 to absorb tiny floating-point drift.
    """
    labels_y = np.asarray(labels_y)
    labels_x = np.asarray(labels_x)
    n = len(labels_y)
    if n == 0:
        return 0.0
    if len(labels_x) != n:
        raise ValueError(
            f"conditional_entropy: partitions must have equal length "
            f"({n} vs {len(labels_x)})."
        )
    # H(Y|X) = H(Y) - I(Y;X); compute directly from the joint histogram.
    p_y, p_x, joint, _ = _joint_histogram(labels_y, labels_x)
    h_y = -np.sum(p_y * np.log(p_y))
    with np.errstate(divide="ignore", invalid="ignore"):
        pj = joint / n
        mask = pj > 0
        mi = np.sum(
            pj[mask] * (np.log(pj[mask]) - np.log(np.outer(p_y, p_x)[mask]))
        )
    return float(max(h_y - mi, 0.0))


def compute_vi(labels1: ArrayLike, labels2: ArrayLike) -> float:
    """Variation of Information: ``VI = H(P) + H(Q) - 2·I(P;Q)``.

    A proper metric on the space of partitions (Meilă 2007). Returns 0 iff
    the two partitions are identical; upper bound is ``log(n)`` where
    ``n`` is the number of items. Symmetric.

    Parameters
    ----------
    labels1, labels2
        Cluster labels for the same N items, any hashable (or int) type.
        Must have the same length.

    Returns
    -------
    vi : float
        Non-negative; small values indicate similar partitions. Clipped at
        ``0.0`` to guard against tiny negative floating-point artefacts
        from the ``H + H - 2·MI`` formulation.
    """
    labels1 = np.asarray(labels1)
    labels2 = np.asarray(labels2)
    n = len(labels1)
    if n == 0:
        return 0.0
    if len(labels2) != n:
        raise ValueError(
            f"compute_vi: partitions must have equal length "
            f"({n} vs {len(labels2)})."
        )

    uniq1, counts1 = np.unique(labels1, return_counts=True)
    uniq2, counts2 = np.unique(labels2, return_counts=True)
    p1 = counts1 / n
    p2 = counts2 / n
    h1 = -np.sum(p1 * np.log(p1))
    h2 = -np.sum(p2 * np.log(p2))

    mi = 0.0
    for c1, n1 in zip(uniq1, counts1):
        mask1 = labels1 == c1
        for c2, n2 in zip(uniq2, counts2):
            n12 = int((mask1 & (labels2 == c2)).sum())
            if n12 > 0:
                mi += (n12 / n) * np.log(n12 * n / (n1 * n2))

    return float(max(h1 + h2 - 2.0 * mi, 0.0))
