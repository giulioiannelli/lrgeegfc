"""Functional tree distance: τ-resolved phase-pair similarity on LRG flows.

Scope report:
    .agents/guides/task-persistence-investigation/2026-04-28_functional-tree-distance.md

For each (patient, band, phase-pair (A, B)) the LRG flow is treated as a
τ-indexed family of matrices.  Two parallel views are computed in lock-step:

    δ_S(τ; A, B) = 1 − Spearman( triu(D^A(τ)),  triu(D^B(τ)) )      # rank, on D = 1/ρ
    δ_P(τ; A, B) = 1 − Pearson ( triu(ρ^A(τ)),  triu(ρ^B(τ)) )      # value, on ρ

`δ_S` reads dendrogram-topology agreement (rank-invariant; robust to the
1/ρ blow-up at small τ).  `δ_P` reads height/scale agreement on the
bounded similarity ρ ∈ [0, 1/N].

Snapshot scalar: δ̂_X = δ_X(τ_min) where
    τ_min = max(1/λ_max^A, 1/λ_max^B)        (finest τ at which both phases
                                              are at-or-coarser than their
                                              fastest mode).
Functional integral:
    Δ_X = (log τ_max − log τ_min)^{−1} · ∫_{τ_min}^{τ_max} δ_X(τ) d(log τ)
where τ_max = min(τ*_A, τ*_B), τ*_φ = argmax_τ C_φ(τ) (the LRG single-peak
characteristic scale; smoothed via Savitzky-Golay before argmax).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple

import networkx as nx
import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import expm
from scipy.signal import savgol_filter
from scipy.stats import pearsonr, spearmanr

from lrgsglib.utils.lrg.infocomm import entropy as lrg_entropy
from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

__all__ = [
    "FunctionalTreeDistanceResult",
    "tau_star_from_C",
    "tau_interval",
    "delta_curves",
    "compute_functional_tree_distance",
    "delta_integrals_from_adjacencies",
]


@dataclass
class FunctionalTreeDistanceResult:
    """Per (patient, band, phase-pair) functional-tree-distance output."""

    patient: str
    band: str
    phase_A: str
    phase_B: str
    fc_method: str

    tau_grid: np.ndarray  # (n_tau,) log-spaced over [tau_min, tau_max]
    delta_S: np.ndarray   # (n_tau,) Spearman distance on D = 1/ρ
    delta_P: np.ndarray   # (n_tau,) Pearson distance on ρ

    delta_hat_S: float    # snapshot at tau_min
    delta_hat_P: float
    Delta_S: float        # log-tau-normalised integral over [tau_min, tau_max]
    Delta_P: float

    tau_min: float
    tau_max: float
    tau_star_A: float
    tau_star_B: float
    lambda_max_A: float
    lambda_max_B: float
    n_nodes: int

    interval_compatible: bool = True
    n_tau: int = 100
    savgol_window: int = 5

    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def tau_star_from_C(
    entropy_tau: np.ndarray,
    entropy_C: np.ndarray,
    savgol_window: int = 5,
    polyorder: int = 2,
) -> Tuple[float, int, bool]:
    """Find τ* = argmax_τ C(τ) on a smoothed C array.

    `entropy_C` has length `len(entropy_tau) - 1` (np.diff convention from
    `lrgsglib.utils.lrg.infocomm.entropy`).  τ* is reported as the geometric
    midpoint of the τ-interval enclosing the smoothed argmax.

    Returns
    -------
    tau_star : float
        Geometric midpoint of `entropy_tau[i_star], entropy_tau[i_star + 1]`.
    i_star : int
        argmax index into the smoothed C array.
    interior : bool
        True if the peak is strictly interior (`0 < i_star < len(C) − 1`).
    """
    C = np.asarray(entropy_C)
    if savgol_window is not None and savgol_window > polyorder + 1:
        if savgol_window > len(C):
            savgol_window = len(C) if len(C) % 2 == 1 else len(C) - 1
        if savgol_window > polyorder + 1 and savgol_window % 2 == 1:
            C_smooth = savgol_filter(C, window_length=savgol_window, polyorder=polyorder)
        else:
            C_smooth = C
    else:
        C_smooth = C

    i_star = int(np.argmax(C_smooth))
    interior = 0 < i_star < (len(C_smooth) - 1)
    tau = np.asarray(entropy_tau)
    tau_star = float(np.sqrt(tau[i_star] * tau[i_star + 1]))
    return tau_star, i_star, interior


def tau_interval(
    L_A: np.ndarray,
    L_B: np.ndarray,
    entropy_tau_A: np.ndarray,
    entropy_C_A: np.ndarray,
    entropy_tau_B: np.ndarray,
    entropy_C_B: np.ndarray,
    savgol_window: int = 5,
) -> dict:
    """Compute the meaningful τ-interval `[τ_min, τ_max]` for phase-pair (A, B).

    `τ_min = max(1/λ_max^A, 1/λ_max^B)`, `τ_max = min(τ*_A, τ*_B)`.
    """
    spec_A = np.linalg.eigvalsh(L_A)
    spec_B = np.linalg.eigvalsh(L_B)
    lambda_max_A = float(spec_A[-1])
    lambda_max_B = float(spec_B[-1])

    tau_star_A, i_A, interior_A = tau_star_from_C(
        entropy_tau_A, entropy_C_A, savgol_window=savgol_window
    )
    tau_star_B, i_B, interior_B = tau_star_from_C(
        entropy_tau_B, entropy_C_B, savgol_window=savgol_window
    )

    tau_min = max(1.0 / lambda_max_A, 1.0 / lambda_max_B)
    tau_max = min(tau_star_A, tau_star_B)

    return {
        "tau_min": tau_min,
        "tau_max": tau_max,
        "lambda_max_A": lambda_max_A,
        "lambda_max_B": lambda_max_B,
        "tau_star_A": tau_star_A,
        "tau_star_B": tau_star_B,
        "interior_A": interior_A,
        "interior_B": interior_B,
        "interval_compatible": tau_min < tau_max,
    }


def _rho_from_L(L: np.ndarray, tau: float) -> np.ndarray:
    """Heat-kernel density operator at diffusion time τ."""
    K = expm(-tau * L)
    Z = np.trace(K)
    return K / Z


def _D_from_rho(rho: np.ndarray) -> np.ndarray:
    """Communication distance D = 1/ρ, max-symmetrised, zero diagonal."""
    with np.errstate(divide="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    return D


def _upper_triangle(M: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(M.shape[0], k=1)
    return M[iu]


def delta_curves(
    L_A: np.ndarray,
    L_B: np.ndarray,
    tau_grid: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Per-τ Spearman (on `1/ρ`) and Pearson (on `ρ`) correlation distances.

    Both phases must have the same node set / Laplacian dimension.
    """
    if L_A.shape != L_B.shape:
        raise ValueError(
            f"L_A and L_B must have matching shape; got {L_A.shape} vs {L_B.shape}"
        )

    n_tau = len(tau_grid)
    delta_S = np.zeros(n_tau)
    delta_P = np.zeros(n_tau)

    for i, tau in enumerate(tau_grid):
        rho_A = _rho_from_L(L_A, tau)
        rho_B = _rho_from_L(L_B, tau)
        D_A = _D_from_rho(rho_A)
        D_B = _D_from_rho(rho_B)

        triu_DA = _upper_triangle(D_A)
        triu_DB = _upper_triangle(D_B)
        triu_rA = _upper_triangle(rho_A)
        triu_rB = _upper_triangle(rho_B)

        rho_s, _ = spearmanr(triu_DA, triu_DB)
        r_p, _ = pearsonr(triu_rA, triu_rB)

        delta_S[i] = 1.0 - rho_s
        delta_P[i] = 1.0 - r_p

    return delta_S, delta_P


# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------

def delta_integrals_from_adjacencies(
    adj_A: np.ndarray,
    adj_B: np.ndarray,
    tau_grid: np.ndarray,
) -> dict:
    """Compute (δ_S, δ_P) curves and log-τ-normalised integrals on a fixed τ-grid.

    Adjacency matrices must already be on a common node set.  Used for
    within-baseline split-half nulls where the τ-grid is anchored to the
    full-data computation rather than re-derived from the halves.
    """
    L_A = _laplacian(np.asarray(adj_A))
    L_B = _laplacian(np.asarray(adj_B))
    delta_S, delta_P = delta_curves(L_A, L_B, tau_grid)
    log_tau = np.log(tau_grid)
    L = float(log_tau[-1] - log_tau[0]) if len(tau_grid) >= 2 else 0.0
    if L > 0:
        Delta_S = float(trapezoid(delta_S, log_tau) / L)
        Delta_P = float(trapezoid(delta_P, log_tau) / L)
    else:
        Delta_S = float(delta_S[0])
        Delta_P = float(delta_P[0])
    return {
        "tau_grid": tau_grid,
        "delta_S": delta_S,
        "delta_P": delta_P,
        "Delta_S": Delta_S,
        "Delta_P": Delta_P,
    }


def _common_giant_component_indices(adj_phases: dict) -> np.ndarray:
    """Indices that lie in every phase's giant component.

    `adj_phases` maps phase name → adjacency matrix.  All adjacencies must
    share the same node ordering (i.e. same N rows/cols index the same
    contacts).  Returns sorted integer indices.
    """
    sets = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(adj)
        components = list(nx.connected_components(G))
        if not components:
            sets.append(set())
            continue
        giant = max(components, key=len)
        sets.append(set(giant))
    common = set.intersection(*sets) if sets else set()
    return np.array(sorted(common), dtype=int)


def _laplacian(adj: np.ndarray) -> np.ndarray:
    """Combinatorial graph Laplacian L = D − A on a dense adjacency."""
    deg = adj.sum(axis=1)
    return np.diag(deg) - adj


def compute_functional_tree_distance(
    patient: str,
    band: str,
    phase_A: str,
    phase_B: str,
    fc_method: str = "imcoh_abs",
    n_tau: int = 100,
    savgol_window: int = 5,
    *,
    entropy_steps: int = 400,
    entropy_t1: float = -3.0,
    entropy_t2: float = 5.0,
    common_giant_phases: Optional[Tuple[str, ...]] = None,
) -> FunctionalTreeDistanceResult:
    """Compute the functional-tree-distance pair (δ_S, δ_P) curves and scalars.

    Parameters
    ----------
    patient, band, phase_A, phase_B :
        Standard cohort identifiers.
    fc_method :
        Defaults to ``"imcoh_abs"`` per the current era.
    n_tau :
        Number of log-spaced τ points in the integration interval.
    savgol_window :
        Smoothing window for the C(τ) argmax; odd integer ≥ 5.
    common_giant_phases :
        Optional phase-set used to define the *common* giant component
        (intersection across phases).  Defaults to the union of
        ``(phase_A, phase_B)`` only — keeps the analysis primitive.  Pass
        ``("rest_pre", "task_learn", "task_test", "rest_post")`` to
        intersect across all four phases.
    entropy_steps, entropy_t1, entropy_t2 :
        τ-grid parameters for the entropy curves used to find τ*.
    """
    from lrg_eegfc.workflow.fc import load_fc_matrix

    if common_giant_phases is None:
        common_giant_phases = (phase_A, phase_B)

    adj_phases: dict[str, np.ndarray] = {}
    for phase in common_giant_phases:
        adj = load_fc_matrix(patient, phase, band, fc_method)
        if adj is None:
            raise FileNotFoundError(
                f"FC matrix missing for {patient}/{phase}/{band}/{fc_method}"
            )
        adj_phases[phase] = np.asarray(adj)

    common_idx = _common_giant_component_indices(adj_phases)
    if common_idx.size < 4:
        raise ValueError(
            f"Common giant component too small ({common_idx.size}) for "
            f"{patient}/{band}; bailing."
        )

    A_full = adj_phases[phase_A]
    B_full = adj_phases[phase_B]
    A = A_full[np.ix_(common_idx, common_idx)]
    B = B_full[np.ix_(common_idx, common_idx)]

    L_A = _laplacian(A)
    L_B = _laplacian(B)

    G_A = nx.from_numpy_array(A)
    G_B = nx.from_numpy_array(B)

    sm1_A, dS_A, _, t_A = lrg_entropy(
        G_A, steps=entropy_steps, t1=entropy_t1, t2=entropy_t2
    )
    sm1_B, dS_B, _, t_B = lrg_entropy(
        G_B, steps=entropy_steps, t1=entropy_t1, t2=entropy_t2
    )

    interval = tau_interval(
        L_A=L_A,
        L_B=L_B,
        entropy_tau_A=t_A,
        entropy_C_A=dS_A,
        entropy_tau_B=t_B,
        entropy_C_B=dS_B,
        savgol_window=savgol_window,
    )

    tau_min = interval["tau_min"]
    tau_max = interval["tau_max"]
    interval_ok = bool(interval["interval_compatible"])

    if interval_ok:
        tau_grid = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)
    else:
        # Snapshot-only mode: a single τ at the (would-be) lower bound.
        tau_grid = np.array([tau_min])

    delta_S, delta_P = delta_curves(L_A, L_B, tau_grid)

    delta_hat_S = float(delta_S[0])
    delta_hat_P = float(delta_P[0])

    if interval_ok and len(tau_grid) >= 2:
        log_tau = np.log(tau_grid)
        L = float(log_tau[-1] - log_tau[0])
        Delta_S = float(trapezoid(delta_S, log_tau) / L)
        Delta_P = float(trapezoid(delta_P, log_tau) / L)
    else:
        Delta_S = float("nan")
        Delta_P = float("nan")

    return FunctionalTreeDistanceResult(
        patient=patient,
        band=band,
        phase_A=phase_A,
        phase_B=phase_B,
        fc_method=fc_method,
        tau_grid=tau_grid,
        delta_S=delta_S,
        delta_P=delta_P,
        delta_hat_S=delta_hat_S,
        delta_hat_P=delta_hat_P,
        Delta_S=Delta_S,
        Delta_P=Delta_P,
        tau_min=tau_min,
        tau_max=tau_max,
        tau_star_A=float(interval["tau_star_A"]),
        tau_star_B=float(interval["tau_star_B"]),
        lambda_max_A=float(interval["lambda_max_A"]),
        lambda_max_B=float(interval["lambda_max_B"]),
        n_nodes=int(common_idx.size),
        interval_compatible=interval_ok,
        n_tau=int(len(tau_grid)),
        savgol_window=int(savgol_window),
        extra={
            "common_giant_phases": list(common_giant_phases),
            "interior_A": bool(interval["interior_A"]),
            "interior_B": bool(interval["interior_B"]),
            "common_node_indices": common_idx.tolist(),
        },
    )
