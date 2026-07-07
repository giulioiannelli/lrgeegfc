"""Diffusion geometry from Laplacian eigenpairs (Coifman-Lafon 2006).

The tau-resolved **diffusion distance** on a graph, built from the eigenpairs
(lambda_a, u_a) of the combinatorial Laplacian ``L_hat = D_hat - W``:

    D_tau(i, j)^2 = sum_{a>=2} exp(-2 * tau * lambda_a) * (u_a(i) - u_a(j))^2

This is the eigenvalue-WEIGHTED node geometry of the heat kernel ``e^{-tau L}``:
it weights every eigenvector by a smooth, tau-tunable function of its eigenvalue,
continuously interpolating between the small-tau strength/adjacency regime and the
large-tau Fiedler collapse. It is *distinct* from the LRG communication distance
``1 / rho_hat(tau)`` (which is inverted and then made ultrametric by UPGMA linkage);
here no linkage and no inversion are applied. The trivial constant mode (a = 1,
lambda_1 = 0) contributes zero to every distance and is dropped.

General graph/network primitive - carries no manuscript-local scope.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "diffusion_gram",
    "diffusion_distance_matrix",
    "diffusion_distance_condensed",
    "effective_modes",
]


def _weights(eigvals: np.ndarray, tau: float):
    """Non-trivial eigenvalues and their diffusion weights exp(-2 tau lambda)."""
    lam = np.asarray(eigvals, dtype=float)[1:]     # drop constant mode a=1
    return lam, np.exp(-2.0 * float(tau) * lam)


def diffusion_gram(eigvals: np.ndarray, eigvecs: np.ndarray, tau: float) -> np.ndarray:
    """Gram matrix ``G_ij = <Psi_tau(i), Psi_tau(j)>`` of diffusion coordinates."""
    _, w = _weights(eigvals, tau)
    U = np.asarray(eigvecs, dtype=float)[:, 1:]
    psi = U * np.sqrt(w)[None, :]
    return psi @ psi.T


def diffusion_distance_matrix(eigvals: np.ndarray, eigvecs: np.ndarray,
                              tau: float) -> np.ndarray:
    """Square ``(N, N)`` matrix of diffusion distances ``D_tau(i, j)``."""
    gram = diffusion_gram(eigvals, eigvecs, tau)
    g = np.diag(gram)
    d2 = g[:, None] + g[None, :] - 2.0 * gram
    np.fill_diagonal(d2, 0.0)
    return np.sqrt(np.clip(d2, 0.0, None))


def diffusion_distance_condensed(eigvals: np.ndarray, eigvecs: np.ndarray,
                                 tau: float) -> np.ndarray:
    """Condensed upper-triangle vector (scipy ``squareform`` order) of ``D_tau``."""
    dmat = diffusion_distance_matrix(eigvals, eigvecs, tau)
    iu = np.triu_indices(dmat.shape[0], k=1)
    return dmat[iu]


def effective_modes(eigvals: np.ndarray, tau: float) -> float:
    """Participation ratio of the diffusion weights = effective number of active
    modes. Runs from ``N - 1`` (tau -> 0) down to ``~1`` (Fiedler collapse)."""
    _, w = _weights(eigvals, tau)
    return float(w.sum() ** 2 / np.sum(w ** 2))
