"""Smoke tests for lrg_eegfc.utils.metrics.tree.

Covers helpers that have been promoted from scripts into the library
(2026-04 onward). Each test fixes a small synthetic input so the
expected value is deterministic.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.cluster.hierarchy import cophenet as scipy_cophenet
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

from lrg_eegfc.utils.metrics.tree import cophenet_matrix, dmax_from_Z, induced_linkage


def _toy_linkage(seed: int = 0, n: int = 8) -> np.ndarray:
    """Random small UPGMA linkage on Euclidean distances of n 2D points."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 2))
    return linkage(pdist(X), method="average")


def test_cophenet_matrix_default_is_square_with_zero_diag():
    Z = _toy_linkage()
    n = Z.shape[0] + 1
    M = cophenet_matrix(Z)
    assert M.shape == (n, n)
    # symmetric
    assert np.allclose(M, M.T)
    # zero diagonal
    assert np.allclose(np.diag(M), 0.0)


def test_cophenet_matrix_condensed_matches_scipy():
    Z = _toy_linkage(seed=7)
    M_condensed = cophenet_matrix(Z, condensed=True)
    expected = scipy_cophenet(Z)
    if isinstance(expected, tuple):
        expected = expected[1]
    assert np.allclose(M_condensed, expected)


def test_cophenet_matrix_full_squareform_matches_scipy():
    Z = _toy_linkage(seed=1)
    M_full = cophenet_matrix(Z, condensed=False)
    expected_condensed = scipy_cophenet(Z)
    if isinstance(expected_condensed, tuple):
        expected_condensed = expected_condensed[1]
    assert np.allclose(M_full, squareform(expected_condensed))


def test_cophenet_matrix_top_height_matches_dmax():
    Z = _toy_linkage(seed=3)
    M = cophenet_matrix(Z)
    # The max cophenetic distance equals the top merge height of the linkage.
    assert M.max() == pytest.approx(dmax_from_Z(Z), rel=1e-9)


def test_cophenet_matrix_consistent_with_induced_linkage():
    Z = _toy_linkage(seed=11, n=10)
    M = cophenet_matrix(Z)
    # induced_linkage is the canonical user of squareform(cophenet(Z));
    # the helper should produce exactly the same intermediate.
    leaves = [0, 2, 4, 6]
    Z_sub_via_helper_path = M[np.ix_(leaves, leaves)]
    Z_sub = induced_linkage(Z, leaves)
    # Sub-linkage should encode the same pairwise heights as the helper's
    # restricted matrix (UPGMA preserves heights on sub-cophenetic).
    sub_full = cophenet_matrix(Z_sub)
    assert np.allclose(sub_full, Z_sub_via_helper_path)
