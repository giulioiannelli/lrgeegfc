"""Sparsification methods for coherence matrices.

Provides multiple approaches to extract statistically significant edges
from MSC (magnitude-squared coherence) matrices:

- **soft**: surrogate-based soft weighting (heuristic, ``(1-p)*W``)
- **fdr**: surrogate-based FDR-corrected thresholding (Benjamini-Hochberg)
- **disparity**: Serrano-Boguna-Vespignani disparity filter (network backbone)
- **hybrid**: surrogate excess coherence + disparity filter
- **ecm**: Enhanced Configuration Model via NEMtropy (maximum-entropy null)
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
from numpy.typing import NDArray


__all__ = [
    'soft_sparsify_surrogate',
    'fdr_sparsify_surrogate',
    'disparity_filter',
    'hybrid_sparsify',
    'ecm_sparsify',
]


# ---------------------------------------------------------------------------
# Helper: empirical p-values from surrogate null
# ---------------------------------------------------------------------------

def _empirical_pvalues(W: NDArray, W_null: NDArray) -> NDArray:
    """Compute per-edge empirical p-values from a surrogate null distribution.

    Parameters
    ----------
    W : NDArray, shape (N, N)
        Observed matrix.
    W_null : NDArray, shape (n_surrogates, N, N)
        Null distribution.

    Returns
    -------
    p_values : NDArray, shape (N, N)
        p_ij = (# surrogates with W_null[r,i,j] >= W[i,j]) / n_surrogates.
    """
    n_surrogates = W_null.shape[0]
    return np.sum(W_null >= W[None, :, :], axis=0) / n_surrogates


# ---------------------------------------------------------------------------
# Method 1: soft (original)
# ---------------------------------------------------------------------------

def soft_sparsify_surrogate(
    W: NDArray,
    W_null: NDArray,
) -> NDArray:
    """
    Apply soft sparsification to a coherence matrix using surrogate null distribution.

    This function computes empirical p-values for each edge based on a surrogate null
    distribution and uses them to downweight non-significant edges while preserving
    the weight geometry.

    Parameters
    ----------
    W : NDArray
        Observed band-averaged MSC matrix of shape (N, N)
    W_null : NDArray
        Null distribution from surrogates of shape (n_surrogates, N, N)

    Returns
    -------
    A : NDArray
        Soft-sparsified adjacency matrix of shape (N, N)

    Notes
    -----
    The sparsification is computed as:

    1. For each edge (i, j), compute empirical p-value:
       p_ij = (#surrogates r with W_null[r, i, j] >= W[i, j]) / n_surrogates

    2. Compute soft weight:
       g_ij = 1 - p_ij

    3. Apply soft sparsification:
       A_ij = g_ij * W_ij

    This approach preserves the full weighted matrix but downweights non-significant
    edges without binarization or strict backbone extraction.
    """
    p_values = _empirical_pvalues(W, W_null)

    # Soft sparsification: g_ij = 1 - p_ij, then apply to weights
    A = (1.0 - p_values) * W

    # Ensure diagonal is zero (no self-loops)
    np.fill_diagonal(A, 0.0)

    return A


# ---------------------------------------------------------------------------
# Method 2: FDR-corrected thresholding
# ---------------------------------------------------------------------------

def fdr_sparsify_surrogate(
    W: NDArray,
    W_null: NDArray,
    q: float = 0.05,
) -> NDArray:
    """
    FDR-corrected thresholding using surrogate null distribution.

    Computes empirical p-values per edge, applies Benjamini-Hochberg FDR
    correction across all edges, and retains edges whose corrected p-value
    falls below the threshold *q*.  Retained edges keep their **original**
    weight; all others are set to zero.

    Parameters
    ----------
    W : NDArray
        Observed band-averaged MSC matrix of shape (N, N).
    W_null : NDArray
        Null distribution from surrogates of shape (n_surrogates, N, N).
    q : float, optional
        False discovery rate level (default: 0.05).

    Returns
    -------
    A : NDArray
        FDR-thresholded adjacency matrix of shape (N, N).

    Notes
    -----
    1. Empirical p-value per edge from surrogate null.
    2. Upper-triangle p-values sorted ascending; Benjamini-Hochberg step-up:
       reject H_0 for edge with rank *k* if ``p_(k) <= k * q / m``.
    3. Surviving edges keep original weight ``W_ij``; others zeroed.
    """
    N = W.shape[0]
    p_values = _empirical_pvalues(W, W_null)

    # Upper triangle only (symmetric matrix, avoid double-counting)
    triu_i, triu_j = np.triu_indices(N, k=1)
    p_upper = p_values[triu_i, triu_j]
    m = len(p_upper)

    # Benjamini-Hochberg procedure
    sorted_idx = np.argsort(p_upper)
    sorted_p = p_upper[sorted_idx]
    thresholds = np.arange(1, m + 1) * q / m

    reject_mask_sorted = sorted_p <= thresholds
    if np.any(reject_mask_sorted):
        # Largest k such that p_(k) <= k*q/m; reject all ranks 1..k
        max_k = np.max(np.where(reject_mask_sorted)[0])
        reject_sorted = np.zeros(m, dtype=bool)
        reject_sorted[: max_k + 1] = True
    else:
        reject_sorted = np.zeros(m, dtype=bool)

    # Map back to original edge ordering
    reject = np.zeros(m, dtype=bool)
    reject[sorted_idx] = reject_sorted

    # Build symmetric output
    A = np.zeros_like(W)
    passed_i = triu_i[reject]
    passed_j = triu_j[reject]
    A[passed_i, passed_j] = W[passed_i, passed_j]
    A[passed_j, passed_i] = W[passed_j, passed_i]

    np.fill_diagonal(A, 0.0)
    return A


# ---------------------------------------------------------------------------
# Method 3: Disparity filter (Serrano, Boguna, Vespignani 2009)
# ---------------------------------------------------------------------------

def disparity_filter(
    W: NDArray,
    alpha: float = 0.05,
) -> NDArray:
    """
    Multiscale backbone extraction via the disparity filter.

    For each node *i* with strength ``s_i`` and degree ``k_i``, the normalised
    weight of edge (i, j) is ``p_ij = w_ij / s_i``.  The null hypothesis is
    that the ``k_i`` edge-weights of node *i* are drawn uniformly.  Under that
    null the probability of observing a normalised weight >= ``p_ij`` is:

        alpha_ij = (1 - p_ij)^{k_i - 1}

    An edge is retained if it is significant (``alpha_ij < alpha``) from
    **either** endpoint.

    Parameters
    ----------
    W : NDArray
        Weighted adjacency matrix of shape (N, N).  Must be non-negative
        and symmetric.
    alpha : float, optional
        Significance level (default: 0.05).

    Returns
    -------
    A : NDArray
        Backbone adjacency matrix of shape (N, N).

    References
    ----------
    Serrano, M. A., Boguna, M. & Vespignani, A.  "Extracting the multiscale
    backbone of complex weighted networks."  *PNAS* **106**, 6483-6488 (2009).
    """
    N = W.shape[0]
    W_pos = np.maximum(W, 0.0).copy()
    np.fill_diagonal(W_pos, 0.0)

    strength = W_pos.sum(axis=1)        # (N,)
    degree = (W_pos > 0).sum(axis=1)    # (N,)

    # Normalised weights: p_ij = w_ij / s_i
    with np.errstate(divide='ignore', invalid='ignore'):
        P = W_pos / strength[:, None]
    P = np.nan_to_num(P, nan=0.0, posinf=0.0, neginf=0.0)

    # Disparity significance: alpha_ij = (1 - p_ij)^(k_i - 1)
    exponent = np.maximum(degree[:, None] - 1, 0)
    alpha_values = np.power(1.0 - P, exponent)

    # Nodes with degree <= 1: their single edge is trivially significant
    alpha_values[degree <= 1, :] = 0.0

    # Keep edge if significant from EITHER endpoint
    keep = (alpha_values < alpha) | (alpha_values.T < alpha)

    A = np.where(keep, W_pos, 0.0)
    np.fill_diagonal(A, 0.0)
    return A


# ---------------------------------------------------------------------------
# Method 4: Hybrid (surrogate excess + disparity filter)
# ---------------------------------------------------------------------------

def hybrid_sparsify(
    W: NDArray,
    W_null: NDArray,
    alpha: float = 0.05,
) -> NDArray:
    """
    Hybrid surrogate + disparity filter sparsification.

    Combines surrogate-based signal-level validation with the network-
    topological disparity filter:

    1. Compute per-edge surrogate mean ``<W_null>_ij``.
    2. Excess coherence: ``W_excess = max(W - <W_null>, 0)``.
    3. Apply :func:`disparity_filter` to ``W_excess``.

    The surrogate subtraction removes the spectral bias (signal physics);
    the disparity filter removes topologically trivial edges (network physics).

    Parameters
    ----------
    W : NDArray
        Observed band-averaged MSC matrix of shape (N, N).
    W_null : NDArray
        Null distribution from surrogates of shape (n_surrogates, N, N).
    alpha : float, optional
        Significance level for the disparity filter (default: 0.05).

    Returns
    -------
    A : NDArray
        Hybrid-sparsified adjacency matrix of shape (N, N).
    """
    # Surrogate mean per edge
    W_null_mean = np.mean(W_null, axis=0)

    # Excess coherence (clip at zero)
    W_excess = np.maximum(W - W_null_mean, 0.0)
    np.fill_diagonal(W_excess, 0.0)

    # Disparity filter on the excess
    A = disparity_filter(W_excess, alpha=alpha)
    return A


# ---------------------------------------------------------------------------
# Method 5: Enhanced Configuration Model (NEMtropy CReMa)
# ---------------------------------------------------------------------------

def ecm_sparsify(
    W: NDArray,
    alpha: float = 0.05,
    n_ensemble: int = 100,
    weight_scale: int = 1000,
) -> NDArray:
    """
    Enhanced Configuration Model sparsification using NEMtropy.

    Uses the CReMa (Conditional Reconstruction and Estimation Model Assembly)
    as a maximum-entropy null model that preserves both the **degree sequence**
    and **strength sequence** of the observed network.  An ensemble of weighted
    graphs is sampled from the fitted model, and per-edge z-scores are used to
    identify edges significantly stronger than expected.

    Parameters
    ----------
    W : NDArray
        Weighted adjacency matrix of shape (N, N).  Must be non-negative
        and symmetric.  Typically a dense MSC matrix with values in [0, 1].
    alpha : float, optional
        Significance level for the one-sided z-test (default: 0.05).
    n_ensemble : int, optional
        Number of weighted graphs sampled from the CReMa ensemble
        (default: 100).  More samples give more stable z-scores.
    weight_scale : int, optional
        Scaling factor to convert float weights to integers, since the CReMa
        model requires integer-valued adjacency matrices (default: 1000).
        ``W_int = round(W * weight_scale)``.

    Returns
    -------
    A : NDArray
        Sparsified adjacency matrix of shape (N, N).  Significant edges
        keep their **original** (unscaled) weight; all others are zero.

    Notes
    -----
    1. Float weights are scaled to integers: ``W_int = round(W * scale)``.
    2. The CReMa model is fitted via ``NEMtropy.UndirectedGraph.solve_tool()``
       with ``model='crema'``.  CReMa preserves both degree and strength
       sequences, generating a maximum-entropy ensemble of weighted graphs.
    3. An ensemble of ``n_ensemble`` weighted graphs is sampled; each sample
       is a weighted edge list (source, target, weight).
    4. Per-edge z-scores: ``z_ij = (w_obs - <w_null>) / std(w_null)``.
    5. Edges with ``z_ij > z_{1-alpha}`` (one-sided) are deemed significant
       and retain the original weight ``W_ij``.

    References
    ----------
    Squartini, T. & Garlaschelli, D.  "Analytical maximum-likelihood method
    to detect patterns in real networks."  *New J. Phys.* **13**, 083001 (2011).

    Parisi, F. et al.  "A faster horse on a safer trail: generalized inference
    for the efficient reconstruction of weighted networks."  *New J. Phys.*
    **22**, 053053 (2020).
    """
    from NEMtropy import UndirectedGraph
    from scipy.stats import norm

    N = W.shape[0]
    W_pos = np.maximum(W, 0.0).copy()
    np.fill_diagonal(W_pos, 0.0)

    # Bail out early if the graph is empty
    if W_pos.max() == 0.0:
        return np.zeros_like(W)

    # Scale to integers (CReMa requires integer weights)
    W_int = np.round(W_pos * weight_scale).astype(int)

    # Fit the CReMa model (preserves degree + strength sequences)
    G = UndirectedGraph(adjacency=W_int)
    G.solve_tool(model="crema", method="newton")

    # Sample weighted ensemble to a temporary directory
    triu_i, triu_j = np.triu_indices(N, k=1)
    n_edges = len(triu_i)

    with tempfile.TemporaryDirectory() as tmpdir:
        # NEMtropy requires trailing slash on output_dir
        G.ensemble_sampler(n=n_ensemble, output_dir=tmpdir + "/", cpu_n=1)

        # Load weighted ensemble samples
        # Each file has columns: source, target, weight
        files = sorted(
            f for f in os.listdir(tmpdir)
            if f.endswith(".txt")
        )
        ensemble_weights = np.zeros((len(files), n_edges))
        for k, fname in enumerate(files):
            data = np.loadtxt(os.path.join(tmpdir, fname))
            # Reconstruct symmetric matrix from edge list
            sources = data[:, 0].astype(int)
            targets = data[:, 1].astype(int)
            weights = data[:, 2]
            A_sample = np.zeros((N, N))
            A_sample[sources, targets] = weights
            A_sample[targets, sources] = weights
            ensemble_weights[k] = A_sample[triu_i, triu_j]

    # Observed upper triangle (scaled integer space)
    obs = W_int[triu_i, triu_j].astype(float)

    # Per-edge null statistics
    null_mean = ensemble_weights.mean(axis=0)
    null_std = ensemble_weights.std(axis=0, ddof=1)

    # z-scores (handle zero variance: obs == mean -> z=0; obs > mean -> +inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (obs - null_mean) / null_std
    z = np.nan_to_num(z, nan=0.0, posinf=np.inf, neginf=-np.inf)

    # One-sided significance threshold
    z_threshold = norm.ppf(1.0 - alpha)
    significant = z > z_threshold

    # Build symmetric output with ORIGINAL (unscaled) weights
    A = np.zeros_like(W)
    sig_i = triu_i[significant]
    sig_j = triu_j[significant]
    A[sig_i, sig_j] = W[sig_i, sig_j]
    A[sig_j, sig_i] = W[sig_j, sig_i]

    np.fill_diagonal(A, 0.0)
    return A
