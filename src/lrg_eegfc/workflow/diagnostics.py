"""MSLCD cross-condition diagnostics.

Compute scalar diagnostics from the LRG multiscale community detection
pipeline for every (patient, phase, band) triplet.  Covers spectral
properties, entropic susceptibility, partition richness, coarsening
trajectories, metastability statistics, and community-size balance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from scipy.signal import find_peaks
from scipy.spatial.distance import squareform

from lrg_eegfc.config.paths import MSC_CACHE
from lrg_eegfc.visuals.lrg import compute_partition_stability_index
from lrg_eegfc.workflow.msc import load_msc_matrix

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_TAU_META = 40          # log-spaced tau for data-driven metastability
N_TAU_ENTROPY = 300      # points for S(tau) / C(tau) curve
PAD_FACTOR = 0.1         # extend tau grid beyond resolution window (decades)
LOG_DROP_MAX = 0.5       # sensible PSI peak log-drop threshold
MIN_N_SPACING = 5        # min spacing between sensible n* values


# ---------------------------------------------------------------------------
# Dataclass for the full diagnostic record
# ---------------------------------------------------------------------------
@dataclass
class TripletDiagnostics:
    """All scalar diagnostics for one (patient, phase, band) triplet."""

    patient: str
    phase: str
    band: str

    # A1 — spectral
    N: int = 0
    lambda_max: float = np.nan
    lambda_2: float = np.nan
    W: float = np.nan            # resolution window width (decades)
    gap_ratio: float = np.nan    # lambda_2 / lambda_max
    tau_min: float = np.nan
    tau_max: float = np.nan

    # A2 — entropic susceptibility
    N_peaks: int = 0
    tau_star1: float = np.nan
    tau_star1_norm: float = np.nan
    C_max: float = np.nan
    W_peak: float = np.nan       # FWHM in log10 tau

    # A3 — partition richness at tau_min
    N_sens: int = 0
    sensible_ns: str = ""        # comma-separated list
    n_global: int = 0            # argmax PSI
    psi_max: float = np.nan

    # A4 — coarsening trajectory endpoints
    n_start: int = 0
    n_end: int = 0
    coarsening_ratio: float = np.nan

    # A5 — metastability
    mu_mean: float = np.nan
    mu_max: float = np.nan
    f_meta_01: float = np.nan
    f_meta_02: float = np.nan
    sigma_mu: float = np.nan

    # A6 — community size balance
    H_size: float = np.nan
    C_max_frac: float = np.nan
    n_singleton: int = 0

    def as_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _entropy_at_tau(eigenvalues: np.ndarray, tau: float) -> float:
    """Von Neumann entropy S(tau) via Boltzmann distribution."""
    log_boltz = -tau * eigenvalues
    log_Z = np.logaddexp.reduce(log_boltz)
    log_p = log_boltz - log_Z
    p = np.exp(log_p)
    return float(-np.sum(p * np.where(p > 1e-300, log_p, 0)))


def _compute_propagator(eigenvalues: np.ndarray,
                        eigenvectors: np.ndarray,
                        tau: float) -> np.ndarray:
    """Heat-kernel propagator K(tau) = V exp(-tau Lambda) V^T."""
    exp_vals = np.exp(-tau * eigenvalues)
    return (eigenvectors * exp_vals[None, :]) @ eigenvectors.T


def _ultrametric_from_propagator(K: np.ndarray) -> np.ndarray:
    """D_ij = 1 / K_ij  (off-diagonal)."""
    N = K.shape[0]
    D = np.zeros_like(K)
    mask = ~np.eye(N, dtype=bool)
    D[mask] = 1.0 / np.where(K[mask] > 1e-30, K[mask], 1e-30)
    return D


def sensible_psi_peaks(psi_vals: np.ndarray,
                       n_comms: np.ndarray,
                       log_drop_max: float = LOG_DROP_MAX,
                       min_n_spacing: int = MIN_N_SPACING):
    """Find PSI peaks via mean-floor + log-drop chain rule.

    Returns (kept_indices, n_sensible).
    """
    if len(psi_vals) == 0:
        return [], 0
    peak_idx, _ = find_peaks(psi_vals)
    if len(peak_idx) == 0:
        return [int(np.argmax(psi_vals))], 1
    psi_mean = float(psi_vals.mean())
    peak_idx = peak_idx[psi_vals[peak_idx] > psi_mean]
    if len(peak_idx) == 0:
        return [int(np.argmax(psi_vals))], 1
    order = np.argsort(psi_vals[peak_idx])[::-1]
    peak_idx = peak_idx[order]
    kept = [int(peak_idx[0])]
    kept_n = {int(n_comms[peak_idx[0]])}
    for ci in peak_idx[1:]:
        n_val = int(n_comms[ci])
        if any(abs(n_val - kn) < min_n_spacing for kn in kept_n):
            continue
        log_drop = (np.log10(psi_vals[kept[-1]])
                    - np.log10(max(psi_vals[ci], 1e-30)))
        if log_drop <= log_drop_max:
            kept.append(int(ci))
            kept_n.add(n_val)
        else:
            break
    return kept, len(kept)


def _align_cluster_labels(labels_list: list[np.ndarray]) -> list[np.ndarray]:
    """Hungarian-match cluster labels across consecutive tau scales."""
    aligned = [labels_list[0].copy()]
    for s in range(1, len(labels_list)):
        prev, curr = aligned[-1], labels_list[s]
        u_prev, u_curr = np.unique(prev), np.unique(curr)
        overlap = np.zeros((len(u_prev), len(u_curr)))
        for i, lp in enumerate(u_prev):
            for j, lc in enumerate(u_curr):
                overlap[i, j] = np.sum((prev == lp) & (curr == lc))
        dim = max(len(u_prev), len(u_curr))
        cost = np.zeros((dim, dim))
        cost[:len(u_prev), :len(u_curr)] = -overlap
        row_ind, col_ind = linear_sum_assignment(cost)
        label_map = {}
        for r, c in zip(row_ind, col_ind):
            if r < len(u_prev) and c < len(u_curr):
                label_map[u_curr[c]] = u_prev[r]
        used = set(label_map.values())
        next_lbl = max(max(used) + 1, max(u_prev) + 1) if used else 1
        for lc in u_curr:
            if lc not in label_map:
                label_map[lc] = next_lbl
                next_lbl += 1
        aligned.append(np.array([label_map[l] for l in curr]))
    return aligned


# ---------------------------------------------------------------------------
# A1 — Spectral properties
# ---------------------------------------------------------------------------

def compute_spectral_properties(eigenvalues: np.ndarray) -> dict:
    """Compute spectral diagnostics from Laplacian eigenvalues."""
    N = len(eigenvalues)
    lam2 = float(eigenvalues[1])
    lam_max = float(eigenvalues[-1])
    tau_min = 1.0 / lam_max if lam_max > 0 else np.nan
    tau_max = 1.0 / lam2 if lam2 > 1e-10 else np.nan
    W = np.log10(lam_max / lam2) if lam2 > 1e-10 else np.nan
    gap_ratio = lam2 / lam_max if lam_max > 0 else np.nan

    return dict(N=N, lambda_max=lam_max, lambda_2=lam2,
                W=W, gap_ratio=gap_ratio, tau_min=tau_min, tau_max=tau_max)


# ---------------------------------------------------------------------------
# A2 — Entropic susceptibility
# ---------------------------------------------------------------------------

def compute_entropy_curve(eigenvalues: np.ndarray,
                          n_points: int = N_TAU_ENTROPY,
                          pad_decades: float = PAD_FACTOR):
    """Compute S(tau), C(tau) on extended grid.

    Returns dict with keys: tau, log10_tau, S, C, tau_min, tau_max,
    lambda_2, lambda_max.
    """
    eigenvalues = np.maximum(eigenvalues, 0.0)
    lam2 = eigenvalues[1]
    lam_max = eigenvalues[-1]
    tau_min = 1.0 / lam_max
    tau_max = 1.0 / lam2

    tau_grid = np.logspace(
        np.log10(tau_min) - pad_decades,
        np.log10(tau_max) + pad_decades,
        n_points,
    )
    S = np.array([_entropy_at_tau(eigenvalues, t) for t in tau_grid])
    log10_tau = np.log10(tau_grid)
    C = -np.gradient(S, log10_tau)

    return dict(tau=tau_grid, log10_tau=log10_tau, S=S, C=C,
                tau_min=tau_min, tau_max=tau_max,
                lambda_2=lam2, lambda_max=lam_max)


def compute_susceptibility_diagnostics(eigenvalues: np.ndarray) -> dict:
    """Compute A2 diagnostics from entropic susceptibility."""
    ec = compute_entropy_curve(eigenvalues)
    tau_min, tau_max = ec["tau_min"], ec["tau_max"]

    # Restrict to resolution window for peak detection
    in_window = (ec["tau"] >= tau_min) & (ec["tau"] <= tau_max)
    C_in = ec["C"].copy()
    C_in[~in_window] = 0.0

    C_max_val = float(C_in.max()) if C_in.max() > 0 else np.nan
    prominence = 0.05 * C_max_val if not np.isnan(C_max_val) else 0.01
    peaks_idx, props = find_peaks(C_in, prominence=max(prominence, 1e-6))

    N_peaks = len(peaks_idx)
    if N_peaks == 0:
        # Fallback: use global max in window
        idx_max = np.argmax(C_in)
        tau_star1 = float(ec["tau"][idx_max])
        tau_star1_norm = tau_star1 / tau_min if tau_min > 0 else np.nan
        W_peak = np.nan
    else:
        # Dominant peak = highest C
        dom_idx = peaks_idx[np.argmax(ec["C"][peaks_idx])]
        tau_star1 = float(ec["tau"][dom_idx])
        tau_star1_norm = tau_star1 / tau_min if tau_min > 0 else np.nan

        # FWHM of dominant peak
        half_max = C_max_val / 2.0
        log_tau = ec["log10_tau"]
        dom_log = log_tau[dom_idx]
        above_half = ec["C"] >= half_max
        # Find left and right boundaries
        left = dom_idx
        while left > 0 and above_half[left]:
            left -= 1
        right = dom_idx
        while right < len(above_half) - 1 and above_half[right]:
            right += 1
        W_peak = float(log_tau[right] - log_tau[left])

    return dict(N_peaks=N_peaks, tau_star1=tau_star1,
                tau_star1_norm=tau_star1_norm, C_max=C_max_val,
                W_peak=W_peak)


# ---------------------------------------------------------------------------
# A3 — Partition richness at tau_min
# ---------------------------------------------------------------------------

def compute_partition_richness(eigenvalues: np.ndarray,
                               eigenvectors: np.ndarray) -> dict:
    """Compute PSI profile at tau_min and return sensible-peak diagnostics."""
    lam_max = eigenvalues[-1]
    tau_min = 1.0 / lam_max
    N = len(eigenvalues)

    K = _compute_propagator(eigenvalues, eigenvectors, tau_min)
    D = _ultrametric_from_propagator(K)
    D_cond = squareform(D, checks=False)
    Z = linkage(D_cond, method="average")

    psi_vals, n_comms = compute_partition_stability_index(Z)
    if len(psi_vals) == 0:
        return dict(N_sens=0, sensible_ns="", n_global=1, psi_max=np.nan)

    # Global PSI maximum
    idx_global = int(np.argmax(psi_vals))
    n_global = int(n_comms[idx_global])
    psi_max = float(psi_vals[idx_global])

    # Sensible peaks
    kept_idx, n_sens = sensible_psi_peaks(psi_vals, n_comms)
    sensible_n_vals = sorted([int(n_comms[k]) for k in kept_idx])

    return dict(
        N_sens=n_sens,
        sensible_ns=",".join(map(str, sensible_n_vals)),
        n_global=n_global,
        psi_max=psi_max,
        _linkage_at_tau_min=Z,         # pass through for A6
        _n_comms=n_comms,
        _psi_vals=psi_vals,
    )


# ---------------------------------------------------------------------------
# A4 + A5 — Coarsening trajectory + metastability (data-driven)
# ---------------------------------------------------------------------------

def compute_coarsening_and_metastability(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    n_tau: int = N_TAU_META,
) -> dict:
    """Data-driven coarsening + metastability in a single pass.

    At each of *n_tau* log-spaced tau in [1/lambda_max, 1/lambda_2]:
      1. propagator -> ultrametric -> linkage -> PSI
      2. sensible peaks -> n_raw(tau) = max(sensible), capped at N/2
      3. enforce monotonic non-increasing: n_max(tau_k) = min(n_raw_k, n_max_{k-1})
      4. cut dendrogram at n_max(tau)
    Then: align labels via Hungarian, compute mu_i.

    The monotonic constraint ensures that coarsening only removes communities
    as the diffusion time grows; without it, marginal PSI peaks that flicker
    in/out cause wild jumps in n_max.
    """
    N = len(eigenvalues)
    lam2 = eigenvalues[1]
    lam_max = eigenvalues[-1]
    tau_min = 1.0 / lam_max
    tau_max = 1.0 / lam2
    nmax_cap = N // 2

    tau_grid = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    # --- First pass: raw n_max at each tau (independent PSI) ---
    nmax_raw = np.zeros(n_tau, dtype=int)
    linkages = []  # store for second pass

    for i, tau in enumerate(tau_grid):
        K = _compute_propagator(eigenvalues, eigenvectors, tau)
        D = _ultrametric_from_propagator(K)
        D_cond = squareform(D, checks=False)
        Z = linkage(D_cond, method="average")
        linkages.append(Z)
        psi_vals, n_comms = compute_partition_stability_index(Z)

        peaks_idx, n_sensible = sensible_psi_peaks(psi_vals, n_comms)
        if n_sensible > 0:
            sensible_n = [int(n_comms[pi]) for pi in peaks_idx
                          if int(n_comms[pi]) <= nmax_cap]
            nmax_raw[i] = max(sensible_n) if sensible_n else 1
        else:
            nmax_raw[i] = 1

    # --- Enforce monotonic non-increasing ---
    nmax_traj = nmax_raw.copy()
    for i in range(1, n_tau):
        nmax_traj[i] = min(nmax_traj[i], nmax_traj[i - 1])

    # --- Second pass: cut dendrograms with monotonic n_max ---
    communities = []
    for i in range(n_tau):
        labels = fcluster(linkages[i], max(nmax_traj[i], 1), criterion="maxclust")
        communities.append(labels)

    # Align labels and compute metastability
    communities_aligned = _align_cluster_labels(communities)

    # Count only steps where n_max actually changes (coarsening events).
    # Plateau steps (identical n_max) are excluded from the denominator so
    # that increasing n_tau doesn't dilute the metastability signal.
    change_mask = np.diff(nmax_traj) != 0  # True where coarsening happens
    n_coarsening = int(change_mask.sum())

    mu = np.zeros(N)
    if n_coarsening > 0:
        for s in range(n_tau - 1):
            if change_mask[s]:
                mu += (communities_aligned[s] != communities_aligned[s + 1]).astype(float)
        mu /= n_coarsening

    # A4 scalars
    n_start = int(nmax_traj[0])
    n_end = int(nmax_traj[-1])
    coarsening_ratio = n_start / n_end if n_end > 0 else np.nan

    # A5 scalars
    mu_mean = float(mu.mean())
    mu_max_val = float(mu.max())
    f_meta_01 = float(np.mean(mu > 0.1))
    f_meta_02 = float(np.mean(mu > 0.2))
    sigma_mu = float(mu.std())

    return dict(
        # A4
        n_start=n_start, n_end=n_end, coarsening_ratio=coarsening_ratio,
        tau_grid=tau_grid, nmax_trajectory=nmax_traj,
        nmax_raw=nmax_raw,  # unconstrained for comparison
        # A5
        mu_mean=mu_mean, mu_max=mu_max_val,
        f_meta_01=f_meta_01, f_meta_02=f_meta_02, sigma_mu=sigma_mu,
        mu=mu,
    )


# ---------------------------------------------------------------------------
# A6 — Community size balance
# ---------------------------------------------------------------------------

def compute_community_balance(linkage_matrix: np.ndarray,
                              n_global: int,
                              N: int) -> dict:
    """Size-entropy, largest-community fraction, singletons."""
    if n_global < 2 or N < 2:
        return dict(H_size=np.nan, C_max_frac=np.nan, n_singleton=0)

    labels = fcluster(linkage_matrix, n_global, criterion="maxclust")
    _, counts = np.unique(labels, return_counts=True)
    fracs = counts / N
    # Normalized size entropy: H / ln(n*)
    H_raw = float(-np.sum(fracs * np.log(fracs + 1e-30)))
    H_size = H_raw / np.log(n_global) if n_global > 1 else 0.0
    C_max_frac = float(counts.max() / N)
    n_singleton = int(np.sum(counts == 1))

    return dict(H_size=H_size, C_max_frac=C_max_frac, n_singleton=n_singleton)


# ---------------------------------------------------------------------------
# Master: compute everything for one (patient, phase, band) triplet
# ---------------------------------------------------------------------------

def compute_triplet_diagnostics(
    patient: str,
    phase: str,
    band: str,
    msc_cache: Path = MSC_CACHE,
    nperseg: int = 4096,
    n_tau: int = N_TAU_META,
    verbose: bool = False,
) -> Optional[TripletDiagnostics]:
    """Compute all MSLCD diagnostics for one triplet.

    Returns None if the MSC matrix is missing or the graph is disconnected.
    """
    A = load_msc_matrix(patient, phase, band,
                        cache_root=msc_cache, sparsify="none",
                        n_surrogates=0, nperseg=nperseg)
    if A is None:
        if verbose:
            print(f"  SKIP {patient}/{phase}/{band}: no MSC cache")
        return None

    np.fill_diagonal(A, 0)
    N = A.shape[0]

    # Laplacian eigendecomposition
    D_diag = A.sum(axis=1)
    L = np.diag(D_diag) - A
    eigenvalues, eigenvectors = np.linalg.eigh(L)
    eigenvalues = np.maximum(eigenvalues, 0.0)

    if eigenvalues[1] < 1e-10:
        if verbose:
            print(f"  SKIP {patient}/{phase}/{band}: disconnected graph")
        return None

    diag = TripletDiagnostics(patient=patient, phase=phase, band=band)

    # A1 — spectral
    sp = compute_spectral_properties(eigenvalues)
    diag.N = sp["N"]
    diag.lambda_max = sp["lambda_max"]
    diag.lambda_2 = sp["lambda_2"]
    diag.W = sp["W"]
    diag.gap_ratio = sp["gap_ratio"]
    diag.tau_min = sp["tau_min"]
    diag.tau_max = sp["tau_max"]

    # A2 — susceptibility
    sus = compute_susceptibility_diagnostics(eigenvalues)
    diag.N_peaks = sus["N_peaks"]
    diag.tau_star1 = sus["tau_star1"]
    diag.tau_star1_norm = sus["tau_star1_norm"]
    diag.C_max = sus["C_max"]
    diag.W_peak = sus["W_peak"]

    # A3 — partition richness
    pr = compute_partition_richness(eigenvalues, eigenvectors)
    diag.N_sens = pr["N_sens"]
    diag.sensible_ns = pr["sensible_ns"]
    diag.n_global = pr["n_global"]
    diag.psi_max = pr["psi_max"]

    # A4 + A5 — coarsening + metastability
    cm = compute_coarsening_and_metastability(eigenvalues, eigenvectors,
                                              n_tau=n_tau)
    diag.n_start = cm["n_start"]
    diag.n_end = cm["n_end"]
    diag.coarsening_ratio = cm["coarsening_ratio"]
    diag.mu_mean = cm["mu_mean"]
    diag.mu_max = cm["mu_max"]
    diag.f_meta_01 = cm["f_meta_01"]
    diag.f_meta_02 = cm["f_meta_02"]
    diag.sigma_mu = cm["sigma_mu"]

    # A6 — community balance (at n_global from A3)
    Z_tau_min = pr.get("_linkage_at_tau_min")
    if Z_tau_min is not None and diag.n_global >= 2:
        cb = compute_community_balance(Z_tau_min, diag.n_global, N)
        diag.H_size = cb["H_size"]
        diag.C_max_frac = cb["C_max_frac"]
        diag.n_singleton = cb["n_singleton"]

    return diag, cm["tau_grid"], cm["nmax_trajectory"], eigenvalues


# ---------------------------------------------------------------------------
# Part C — Threshold analysis
# ---------------------------------------------------------------------------

THRESHOLD_PERCENTILES = [0, 50, 70, 80, 85, 90, 95, 97, 99]


def threshold_matrix(A: np.ndarray, percentile: float) -> np.ndarray:
    """Zero out edges below the given weight percentile.

    Uses the upper-triangle non-zero values as the reference distribution.
    """
    if percentile <= 0:
        return A.copy()
    triu_vals = A[np.triu_indices_from(A, k=1)]
    nonzero = triu_vals[triu_vals > 0]
    if len(nonzero) == 0:
        return A.copy()
    thresh_val = np.percentile(nonzero, percentile)
    A_t = A.copy()
    A_t[A_t < thresh_val] = 0.0
    return A_t


def compute_threshold_analysis(
    patient: str,
    phase: str,
    bands: list[str],
    msc_cache: Path = MSC_CACHE,
    nperseg: int = 4096,
    thresholds: list[int] | None = None,
    verbose: bool = False,
) -> dict:
    """Compute threshold analysis for a patient/phase across all bands.

    For each band, loads the full MSC matrix and progressively removes the
    bottom X% of edge weights.  At each threshold, computes:
      - Laplacian eigenvalues
      - Entropy curve S(tau), susceptibility C(tau)
      - Connectivity metrics (n_components, LCC fraction, edge fraction,
        mean degree)

    Returns dict suitable for ``np.savez_compressed``.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components as cc_sparse

    if thresholds is None:
        thresholds = THRESHOLD_PERCENTILES

    result: dict = {}
    result["percentiles"] = np.array(thresholds)

    for band in bands:
        A = load_msc_matrix(
            patient, phase, band,
            cache_root=msc_cache, sparsify="none",
            n_surrogates=0, nperseg=nperseg,
        )
        if A is None:
            if verbose:
                print(f"  SKIP {band}: no MSC cache")
            continue

        np.fill_diagonal(A, 0)
        N = A.shape[0]

        # Reference tau grid from unthresholded Laplacian.
        # Ensure the grid extends to at least log10(tau) = 2 (tau = 100)
        # so that the C(tau) tail decays fully.
        D0 = A.sum(axis=1)
        L0 = np.diag(D0) - A
        ref_eig = np.maximum(np.linalg.eigvalsh(L0), 0.0)
        ref_tau_min = 1.0 / ref_eig[-1]
        log10_left = np.log10(ref_tau_min) - PAD_FACTOR
        log10_right = max(np.log10(1.0 / max(ref_eig[1], 1e-30)) + 0.5, 3.0)
        tau_grid = np.logspace(log10_left, log10_right, N_TAU_ENTROPY)
        log10_tau = np.log10(tau_grid)

        result[f"{band}__tau"] = tau_grid
        result[f"{band}__log10_tau"] = log10_tau

        total_edges = float((A > 0).sum()) / 2.0
        # Columns: n_components, lcc_frac, lcc_density, transitivity
        connectivity = np.zeros((len(thresholds), 4))

        for ti, pct in enumerate(thresholds):
            A_t = threshold_matrix(A, pct)

            # Connected components (on full thresholded graph)
            n_comp, labels = cc_sparse(csr_matrix(A_t > 0))
            _, counts = np.unique(labels, return_counts=True)
            lcc_label = np.argmax(counts)
            lcc_mask = labels == lcc_label
            lcc_size = int(counts.max())
            lcc_frac = float(lcc_size) / N

            # Extract LCC subgraph — used for spectral analysis AND
            # topology metrics.  Computing on the full graph after
            # disconnection injects spurious zero eigenvalues.
            A_lcc = A_t[np.ix_(lcc_mask, lcc_mask)]

            # LCC density: edges_in_lcc / max_possible_edges_in_lcc
            A_lcc_bin = (A_lcc > 0).astype(float)
            n_lcc_edges = float(A_lcc_bin.sum()) / 2.0
            max_lcc_edges = lcc_size * (lcc_size - 1) / 2.0
            lcc_density = (
                n_lcc_edges / max_lcc_edges if max_lcc_edges > 0 else 0.0
            )

            # Transitivity (global clustering coefficient of LCC):
            # 3 * n_triangles / n_connected_triples
            A2 = A_lcc_bin @ A_lcc_bin
            n_triangles = float(np.trace(A2 @ A_lcc_bin)) / 6.0
            k = A_lcc_bin.sum(axis=1)
            n_triples = float(np.sum(k * (k - 1))) / 2.0
            transitivity = (
                3.0 * n_triangles / n_triples if n_triples > 0 else 0.0
            )

            connectivity[ti] = [n_comp, lcc_frac, lcc_density, transitivity]

            # Spectral analysis on LCC
            D_lcc = A_lcc.sum(axis=1)
            L_lcc = np.diag(D_lcc) - A_lcc
            eig_t = np.maximum(np.linalg.eigvalsh(L_lcc), 0.0)

            result[f"{band}__eigenvalues__{pct}"] = eig_t
            result[f"{band}__lcc_size__{pct}"] = lcc_size

            # Entropy / susceptibility on reference tau grid (LCC spectrum)
            S = np.array([_entropy_at_tau(eig_t, t) for t in tau_grid])
            C = -np.gradient(S, log10_tau)
            result[f"{band}__C__{pct}"] = C
            result[f"{band}__S__{pct}"] = S

            if verbose:
                print(
                    f"  {band} pct={pct:2d}: "
                    f"lcc={lcc_size}/{N}, lam2={eig_t[1]:.4e}, "
                    f"n_comp={n_comp}, dens={lcc_density:.3f}, "
                    f"trans={transitivity:.3f}"
                )

        result[f"{band}__connectivity"] = connectivity

    return result
