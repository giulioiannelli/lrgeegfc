#!/usr/bin/env python3
"""Sweep of FUNDAMENTALLY DIFFERENT similarity metrics.

Previous sweep (sweep_tree_metrics.py) tested 17 cophenetic-distance-based
metrics — all failed (max 2/6 unanimous bands).

This sweep tries 4 new families:
  A) Entropy / complexity curve similarity (global hierarchical shape)
  B) Merge-height distribution similarity (dendrogram shape, not node-specific)
  C) Direct MSC matrix comparison (skip the dendrogram entirely)
  D) Laplacian spectral comparison (eigenvalue fingerprint)

For each metric, computes:
  gap = mean_similarity(within_pairs) - mean_similarity(cross_pairs)
  gap > 0 → within-condition more similar → REORGANIZATION detected

Winning metric: gap sign unanimous across ALL patients per band.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr, wasserstein_distance, ks_2samp
from scipy.interpolate import interp1d

from lrg_eegfc.config import BRAIN_BANDS_NAMES
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import LRG_CACHE, MSC_CACHE
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ===================================================================
# Data loading
# ===================================================================
def load_all():
    """Load LRG results + MSC matrices for all (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                # Load LRG result
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    print(f"  SKIP {pat}/{ph}/{band}: no LRG cache")
                    continue

                Z = lrg.linkage_matrix
                c_condensed = cophenet(Z)
                n_nodes = lrg.n_nodes

                # Entropy curves (stored in LRG result)
                ent_tau = lrg.entropy_tau
                ent_S = lrg.entropy_1_minus_S  # 1-S values
                ent_C = lrg.entropy_C           # complexity

                # Merge heights from linkage
                merge_heights = Z[:, 2]

                # Load MSC matrix
                msc = load_msc_matrix(pat, ph, band, cache_root=MSC_CACHE)
                msc_mat = None
                if msc is not None:
                    msc_mat = msc if isinstance(msc, np.ndarray) else msc.adjacency_matrix
                    # Extract upper triangle (excluding diagonal)
                    idx = np.triu_indices(msc_mat.shape[0], k=1)
                    msc_upper = msc_mat[idx]
                else:
                    msc_upper = None

                # Laplacian eigenvalues from MSC
                eig_vals = None
                if msc_mat is not None:
                    try:
                        D = np.diag(msc_mat.sum(axis=1))
                        L = D - msc_mat
                        eig_vals = np.sort(np.linalg.eigvalsh(L))
                    except Exception:
                        pass

                data[(pat, ph, band)] = {
                    "Z": Z,
                    "c_cond": c_condensed,
                    "n": n_nodes,
                    "ent_tau": ent_tau,
                    "ent_S": ent_S,
                    "ent_C": ent_C,
                    "merge_heights": merge_heights,
                    "msc_upper": msc_upper,
                    "msc_mat": msc_mat,
                    "eig_vals": eig_vals,
                }
    return data


# ===================================================================
# FAMILY A: Entropy / complexity curve metrics
# ===================================================================

def _interp_to_common_grid(tau1, s1, tau2, s2, n_points=200):
    """Interpolate two curves to a common tau grid."""
    lo = max(tau1.min(), tau2.min())
    hi = min(tau1.max(), tau2.max())
    if lo >= hi:
        return None, None, None
    common_tau = np.linspace(lo, hi, n_points)
    f1 = interp1d(tau1, s1, kind='linear', bounds_error=False, fill_value='extrapolate')
    f2 = interp1d(tau2, s2, kind='linear', bounds_error=False, fill_value='extrapolate')
    return common_tau, f1(common_tau), f2(common_tau)


def m_entropy_pearson(d1, d2):
    """Pearson correlation of entropy curves S(tau)."""
    _, s1, s2 = _interp_to_common_grid(d1["ent_tau"], d1["ent_S"],
                                        d2["ent_tau"], d2["ent_S"])
    if s1 is None:
        return None
    r, _ = pearsonr(s1, s2)
    return r


def m_entropy_l2_sim(d1, d2):
    """1 / (1 + L2 distance of entropy curves)."""
    _, s1, s2 = _interp_to_common_grid(d1["ent_tau"], d1["ent_S"],
                                        d2["ent_tau"], d2["ent_S"])
    if s1 is None:
        return None
    dist = np.sqrt(np.mean((s1 - s2) ** 2))
    return 1.0 / (1.0 + dist)


def m_entropy_cosine(d1, d2):
    """Cosine similarity of entropy curves."""
    _, s1, s2 = _interp_to_common_grid(d1["ent_tau"], d1["ent_S"],
                                        d2["ent_tau"], d2["ent_S"])
    if s1 is None:
        return None
    dot = np.dot(s1, s2)
    norm = np.linalg.norm(s1) * np.linalg.norm(s2)
    if norm < 1e-12:
        return 0.0
    return dot / norm


def m_entropy_sup_sim(d1, d2):
    """1 / (1 + sup|S1 - S2|). Max deviation of entropy curves."""
    _, s1, s2 = _interp_to_common_grid(d1["ent_tau"], d1["ent_S"],
                                        d2["ent_tau"], d2["ent_S"])
    if s1 is None:
        return None
    sup = np.max(np.abs(s1 - s2))
    return 1.0 / (1.0 + sup)


def m_complexity_pearson(d1, d2):
    """Pearson correlation of complexity curves C(tau)."""
    tau1 = d1["ent_tau"][:len(d1["ent_C"])]
    tau2 = d2["ent_tau"][:len(d2["ent_C"])]
    _, c1, c2 = _interp_to_common_grid(tau1, d1["ent_C"], tau2, d2["ent_C"])
    if c1 is None:
        return None
    r, _ = pearsonr(c1, c2)
    return r


def m_complexity_l2_sim(d1, d2):
    """1 / (1 + L2 distance of complexity curves)."""
    tau1 = d1["ent_tau"][:len(d1["ent_C"])]
    tau2 = d2["ent_tau"][:len(d2["ent_C"])]
    _, c1, c2 = _interp_to_common_grid(tau1, d1["ent_C"], tau2, d2["ent_C"])
    if c1 is None:
        return None
    dist = np.sqrt(np.mean((c1 - c2) ** 2))
    return 1.0 / (1.0 + dist)


def m_complexity_peak_sim(d1, d2):
    """Similarity of complexity peak locations.

    Peak of C(tau) = main phase transition in the hierarchy.
    1 / (1 + |log(tau*_1) - log(tau*_2)|)
    """
    tau1 = d1["ent_tau"][:len(d1["ent_C"])]
    tau2 = d2["ent_tau"][:len(d2["ent_C"])]
    if len(d1["ent_C"]) == 0 or len(d2["ent_C"]) == 0:
        return None
    peak1 = tau1[np.argmax(d1["ent_C"])]
    peak2 = tau2[np.argmax(d2["ent_C"])]
    if peak1 <= 0 or peak2 <= 0:
        return None
    dist = abs(np.log10(peak1) - np.log10(peak2))
    return 1.0 / (1.0 + dist)


def m_entropy_area_sim(d1, d2):
    """Similarity of area under entropy curve.

    Area under S(tau) captures overall hierarchy complexity.
    """
    _, s1, s2 = _interp_to_common_grid(d1["ent_tau"], d1["ent_S"],
                                        d2["ent_tau"], d2["ent_S"])
    if s1 is None:
        return None
    a1 = np.trapz(s1)
    a2 = np.trapz(s2)
    denom = max(abs(a1), abs(a2), 1e-12)
    return 1.0 - abs(a1 - a2) / denom


# ===================================================================
# FAMILY B: Merge height distribution metrics
# ===================================================================

def m_merge_wasserstein_sim(d1, d2):
    """1 / (1 + Wasserstein of merge height distributions)."""
    h1 = np.sort(d1["merge_heights"])
    h2 = np.sort(d2["merge_heights"])
    w = wasserstein_distance(h1, h2)
    return 1.0 / (1.0 + w)


def m_merge_log_wasserstein_sim(d1, d2):
    """1 / (1 + Wasserstein of LOG merge height distributions)."""
    h1 = np.log(np.clip(np.sort(d1["merge_heights"]), 1e-12, None))
    h2 = np.log(np.clip(np.sort(d2["merge_heights"]), 1e-12, None))
    w = wasserstein_distance(h1, h2)
    return 1.0 / (1.0 + w)


def m_merge_ks_sim(d1, d2):
    """1 - KS statistic of merge height distributions."""
    stat, _ = ks_2samp(d1["merge_heights"], d2["merge_heights"])
    return 1.0 - stat


def m_merge_quantile_sim(d1, d2):
    """Quantile profile similarity of merge heights."""
    qs = np.linspace(0.05, 0.95, 19)
    q1 = np.quantile(d1["merge_heights"], qs)
    q2 = np.quantile(d2["merge_heights"], qs)
    r = max(np.ptp(q1), np.ptp(q2), 1e-12)
    rmse = np.sqrt(np.mean((q1 - q2) ** 2)) / r
    return 1.0 - min(rmse, 1.0)


def m_merge_pearson(d1, d2):
    """Pearson of sorted merge heights (as 1D signatures)."""
    h1 = np.sort(d1["merge_heights"])
    h2 = np.sort(d2["merge_heights"])
    # Pad shorter to same length
    n = max(len(h1), len(h2))
    if len(h1) < n:
        h1 = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(h1)), h1)
    if len(h2) < n:
        h2 = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(h2)), h2)
    return pearsonr(h1, h2)[0]


def m_merge_cosine(d1, d2):
    """Cosine similarity of sorted merge heights."""
    h1 = np.sort(d1["merge_heights"])
    h2 = np.sort(d2["merge_heights"])
    n = max(len(h1), len(h2))
    if len(h1) < n:
        h1 = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(h1)), h1)
    if len(h2) < n:
        h2 = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(h2)), h2)
    dot = np.dot(h1, h2)
    norm = np.linalg.norm(h1) * np.linalg.norm(h2)
    if norm < 1e-12:
        return 0.0
    return dot / norm


# ===================================================================
# FAMILY C: Direct MSC matrix comparison
# ===================================================================

def m_msc_pearson(d1, d2):
    """Pearson of vectorized MSC upper triangles."""
    if d1["msc_upper"] is None or d2["msc_upper"] is None:
        return None
    return pearsonr(d1["msc_upper"], d2["msc_upper"])[0]


def m_msc_spearman(d1, d2):
    """Spearman of vectorized MSC upper triangles."""
    if d1["msc_upper"] is None or d2["msc_upper"] is None:
        return None
    return spearmanr(d1["msc_upper"], d2["msc_upper"])[0]


def m_msc_cosine(d1, d2):
    """Cosine similarity of MSC upper triangle vectors."""
    if d1["msc_upper"] is None or d2["msc_upper"] is None:
        return None
    u1, u2 = d1["msc_upper"], d2["msc_upper"]
    dot = np.dot(u1, u2)
    norm = np.linalg.norm(u1) * np.linalg.norm(u2)
    if norm < 1e-12:
        return 0.0
    return dot / norm


def m_msc_frobenius_sim(d1, d2):
    """1 / (1 + normalized Frobenius distance of MSC matrices)."""
    if d1["msc_mat"] is None or d2["msc_mat"] is None:
        return None
    diff = np.linalg.norm(d1["msc_mat"] - d2["msc_mat"], 'fro')
    denom = np.linalg.norm(d1["msc_mat"], 'fro') + np.linalg.norm(d2["msc_mat"], 'fro')
    if denom < 1e-12:
        return 0.0
    return 1.0 - diff / denom


def m_msc_rank_sim(d1, d2):
    """Jaccard overlap of top-20% strongest MSC edges."""
    if d1["msc_upper"] is None or d2["msc_upper"] is None:
        return None
    n = len(d1["msc_upper"])
    k = max(int(n * 0.2), 10)
    idx1 = set(np.argsort(d1["msc_upper"])[-k:])
    idx2 = set(np.argsort(d2["msc_upper"])[-k:])
    return len(idx1 & idx2) / len(idx1 | idx2)


# ===================================================================
# FAMILY D: Spectral comparison
# ===================================================================

def m_spectral_l2_sim(d1, d2):
    """1 / (1 + L2 distance of Laplacian eigenvalue sequences)."""
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    e1, e2 = d1["eig_vals"], d2["eig_vals"]
    # Pad to same length
    n = max(len(e1), len(e2))
    e1p = np.zeros(n); e1p[:len(e1)] = e1
    e2p = np.zeros(n); e2p[:len(e2)] = e2
    dist = np.sqrt(np.mean((e1p - e2p) ** 2))
    return 1.0 / (1.0 + dist)


def m_spectral_pearson(d1, d2):
    """Pearson correlation of Laplacian eigenvalues."""
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    e1, e2 = d1["eig_vals"], d2["eig_vals"]
    n = max(len(e1), len(e2))
    e1p = np.zeros(n); e1p[:len(e1)] = e1
    e2p = np.zeros(n); e2p[:len(e2)] = e2
    return pearsonr(e1p, e2p)[0]


def m_spectral_cosine(d1, d2):
    """Cosine similarity of Laplacian eigenvalues."""
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    e1, e2 = d1["eig_vals"], d2["eig_vals"]
    n = max(len(e1), len(e2))
    e1p = np.zeros(n); e1p[:len(e1)] = e1
    e2p = np.zeros(n); e2p[:len(e2)] = e2
    dot = np.dot(e1p, e2p)
    norm = np.linalg.norm(e1p) * np.linalg.norm(e2p)
    if norm < 1e-12:
        return 0.0
    return dot / norm


def m_spectral_wasserstein_sim(d1, d2):
    """1 / (1 + Wasserstein of eigenvalue distributions)."""
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    w = wasserstein_distance(d1["eig_vals"], d2["eig_vals"])
    return 1.0 / (1.0 + w)


def m_spectral_energy_sim(d1, d2):
    """Similarity of graph energies (sum of absolute eigenvalues).

    Graph energy captures global network complexity.
    """
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    e1 = np.sum(np.abs(d1["eig_vals"]))
    e2 = np.sum(np.abs(d2["eig_vals"]))
    denom = max(e1, e2, 1e-12)
    return 1.0 - abs(e1 - e2) / denom


# ===================================================================
# FAMILY E: Hybrid / novel approaches
# ===================================================================

def m_coph_log_l2_sim(d1, d2):
    """1 / (1 + L2 distance of log cophenetic distances).

    Unlike correlations, this is sensitive to absolute scale shifts.
    Best cophenetic metric from prev sweep was LogHierEmbedCos — try L2.
    """
    lc1 = np.sort(np.log(np.clip(d1["c_cond"], 1e-12, None)))
    lc2 = np.sort(np.log(np.clip(d2["c_cond"], 1e-12, None)))
    dist = np.sqrt(np.mean((lc1 - lc2) ** 2))
    return 1.0 / (1.0 + dist)


def m_entropy_deriv_pearson(d1, d2):
    """Pearson correlation of entropy derivative (dS/d_tau).

    Derivative captures where the hierarchy transitions happen.
    More sensitive to structural changes than the cumulative S(tau).
    """
    tau1 = d1["ent_tau"]
    tau2 = d2["ent_tau"]
    s1 = d1["ent_S"]
    s2 = d2["ent_S"]
    # Compute derivative
    ds1 = np.diff(s1) / np.diff(tau1)
    ds2 = np.diff(s2) / np.diff(tau2)
    # Interpolate to common grid
    t1_mid = 0.5 * (tau1[:-1] + tau1[1:])
    t2_mid = 0.5 * (tau2[:-1] + tau2[1:])
    _, ds1_int, ds2_int = _interp_to_common_grid(t1_mid, ds1, t2_mid, ds2)
    if ds1_int is None:
        return None
    return pearsonr(ds1_int, ds2_int)[0]


def m_entropy_log_pearson(d1, d2):
    """Pearson of entropy curves on LOG tau scale.

    Use log(tau) as x-axis, which compresses the large-tau region
    and expands the small-tau region. More sensitive to fine scales.
    """
    log_tau1 = np.log10(np.clip(d1["ent_tau"], 1e-12, None))
    log_tau2 = np.log10(np.clip(d2["ent_tau"], 1e-12, None))
    _, s1, s2 = _interp_to_common_grid(log_tau1, d1["ent_S"],
                                        log_tau2, d2["ent_S"])
    if s1 is None:
        return None
    return pearsonr(s1, s2)[0]


def m_msc_eigen_profile_sim(d1, d2):
    """Cosine similarity of cumulative eigenvalue profiles.

    CDF of normalized eigenvalues captures the spectrum's shape
    independent of scale.
    """
    if d1["eig_vals"] is None or d2["eig_vals"] is None:
        return None
    e1 = np.sort(d1["eig_vals"])
    e2 = np.sort(d2["eig_vals"])
    # Normalize eigenvalues by their sum (spectral distribution)
    s1 = e1[1:].sum()  # skip zero eigenvalue
    s2 = e2[1:].sum()
    if s1 < 1e-12 or s2 < 1e-12:
        return 0.0
    # Cumulative distribution
    c1 = np.cumsum(e1[1:]) / s1
    c2 = np.cumsum(e2[1:]) / s2
    # Interpolate to same length
    n = max(len(c1), len(c2))
    c1i = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(c1)), c1)
    c2i = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(c2)), c2)
    dot = np.dot(c1i, c2i)
    norm = np.linalg.norm(c1i) * np.linalg.norm(c2i)
    if norm < 1e-12:
        return 0.0
    return dot / norm


# ===================================================================
# Registry (all 26 metrics)
# ===================================================================
METRICS = {
    # A: Entropy / complexity curves
    "EntPearson":       m_entropy_pearson,
    "EntL2Sim":         m_entropy_l2_sim,
    "EntCosine":        m_entropy_cosine,
    "EntSupSim":        m_entropy_sup_sim,
    "EntAreaSim":       m_entropy_area_sim,
    "EntLogPearson":    m_entropy_log_pearson,
    "EntDerivPearson":  m_entropy_deriv_pearson,
    "ComplPearson":     m_complexity_pearson,
    "ComplL2Sim":       m_complexity_l2_sim,
    "ComplPeakSim":     m_complexity_peak_sim,
    # B: Merge height distribution
    "MrgWasserstSim":   m_merge_wasserstein_sim,
    "MrgLogWassSim":    m_merge_log_wasserstein_sim,
    "MrgKSSim":         m_merge_ks_sim,
    "MrgQuantileSim":   m_merge_quantile_sim,
    "MrgPearson":       m_merge_pearson,
    "MrgCosine":        m_merge_cosine,
    # C: Direct MSC matrix
    "MSC_Pearson":      m_msc_pearson,
    "MSC_Spearman":     m_msc_spearman,
    "MSC_Cosine":       m_msc_cosine,
    "MSC_FrobSim":      m_msc_frobenius_sim,
    "MSC_RankSim":      m_msc_rank_sim,
    # D: Spectral
    "SpecL2Sim":        m_spectral_l2_sim,
    "SpecPearson":      m_spectral_pearson,
    "SpecCosine":       m_spectral_cosine,
    "SpecWassSim":      m_spectral_wasserstein_sim,
    "SpecEnergySim":    m_spectral_energy_sim,
    # E: Hybrid
    "CophLogL2Sim":     m_coph_log_l2_sim,
    "SpecEigProfCos":   m_msc_eigen_profile_sim,
}


# ===================================================================
# Computation
# ===================================================================
def compute_gaps(data):
    gaps = {}
    sims = {}

    for metric_name, metric_fn in METRICS.items():
        print(f"  {metric_name}...", end="", flush=True)
        n_ok = 0
        for pat in PATIENTS:
            phases = PATIENT_PHASES[pat]
            for band in BANDS:
                available = [ph for ph in phases if (pat, ph, band) in data]
                pairs = list(combinations(available, 2))

                within_sims = []
                cross_sims = []
                for p1, p2 in pairs:
                    d1 = data[(pat, p1, band)]
                    d2 = data[(pat, p2, band)]
                    try:
                        s = metric_fn(d1, d2)
                        if s is None or np.isnan(s):
                            continue
                    except Exception:
                        continue

                    n_ok += 1
                    cat = classify_pair(p1, p2)
                    sims[(metric_name, pat, band, (p1, p2))] = s
                    if cat == "within":
                        within_sims.append(s)
                    else:
                        cross_sims.append(s)

                if within_sims and cross_sims:
                    gap = np.mean(within_sims) - np.mean(cross_sims)
                    gaps[(metric_name, pat, band)] = gap

        print(f" ({n_ok} comparisons)")

    return gaps, sims


def score_metrics(gaps):
    results = []

    for metric_name in METRICS:
        band_scores = []
        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS:
                key = (metric_name, pat, band)
                if key in gaps:
                    pat_gaps.append((pat, gaps[key]))

            if len(pat_gaps) < 3:
                continue

            gap_vals = [g for _, g in pat_gaps]
            all_pos = all(g > 0 for g in gap_vals)
            all_neg = all(g < 0 for g in gap_vals)
            unanimous = all_pos or all_neg
            mean_gap = np.mean(gap_vals)
            std_gap = np.std(gap_vals)
            cv = std_gap / abs(mean_gap) if abs(mean_gap) > 1e-6 else 999

            band_scores.append({
                "band": band,
                "unanimous": unanimous,
                "direction": "reorg" if all_pos else ("persist" if all_neg else "mixed"),
                "mean_gap": mean_gap,
                "std_gap": std_gap,
                "cv": cv,
                "n_pos": sum(g > 0 for g in gap_vals),
                "n_neg": sum(g < 0 for g in gap_vals),
                "n_patients": len(gap_vals),
                "pat_gaps": pat_gaps,
            })

        n_unanimous = sum(b["unanimous"] for b in band_scores)
        n_strong = sum(b["unanimous"] and abs(b["mean_gap"]) > 0.01 for b in band_scores)
        mean_cv = np.mean([b["cv"] for b in band_scores]) if band_scores else 999

        results.append({
            "metric": metric_name,
            "n_unanimous": n_unanimous,
            "n_strong_unanimous": n_strong,
            "mean_cv": mean_cv,
            "band_scores": band_scores,
        })

    results.sort(key=lambda r: (-r["n_unanimous"], -r["n_strong_unanimous"], r["mean_cv"]))
    return results


def main():
    print("Loading all data (LRG + MSC + spectra)...")
    data = load_all()
    print(f"  Loaded {len(data)} entries\n")

    print("Computing 28 metrics across 4 families...")
    gaps, sims = compute_gaps(data)

    print("\n\nScoring metrics...")
    results = score_metrics(gaps)

    # ---- Print ranking ----
    print("\n" + "=" * 90)
    print("METRIC RANKING — cross-patient consistency (gap sign unanimous)")
    print("=" * 90)
    print(f"{'Rank':<5} {'Metric':<20} {'Family':>8} {'Unan':>5} {'Strong':>7} {'MeanCV':>8}")
    print("-" * 60)

    family_map = {}
    for k in METRICS:
        if k.startswith("Ent") or k.startswith("Compl"):
            family_map[k] = "entropy"
        elif k.startswith("Mrg"):
            family_map[k] = "merge"
        elif k.startswith("MSC"):
            family_map[k] = "MSC"
        elif k.startswith("Spec"):
            family_map[k] = "spectral"
        else:
            family_map[k] = "hybrid"

    for i, r in enumerate(results):
        fam = family_map.get(r["metric"], "?")
        print(f"{i+1:<5} {r['metric']:<20} {fam:>8} {r['n_unanimous']:>5}/6 "
              f"{r['n_strong_unanimous']:>7}/6 {r['mean_cv']:>8.2f}")

    # ---- Detailed top 8 ----
    print("\n\n" + "=" * 90)
    print("TOP 8 METRICS — per-band detail")
    print("=" * 90)
    for r in results[:8]:
        print(f"\n--- {r['metric']} ({r['n_unanimous']}/6 unanimous, "
              f"{r['n_strong_unanimous']}/6 strong) ---")
        print(f"  {'Band':<12} {'Dir':>7} {'Mean':>9} {'Std':>8} "
              f"{'Pos':>4}/{'Neg':>4} {'Unan':>5}")
        for bs in r["band_scores"]:
            mark = "***" if bs["unanimous"] else ""
            print(f"  {bs['band']:<12} {bs['direction']:>7} {bs['mean_gap']:>+9.5f} "
                  f"{bs['std_gap']:>8.5f} {bs['n_pos']:>4}/{bs['n_neg']:>4} "
                  f"{'YES' if bs['unanimous'] else 'no':>5} {mark}")

    # ---- Per-patient detail for top 3 ----
    print("\n\n" + "=" * 90)
    print("TOP 3 — per-patient gaps")
    print("=" * 90)
    for r in results[:3]:
        metric = r["metric"]
        print(f"\n--- {metric} ---")
        print(f"  {'Band':<12}", end="")
        for pat in PATIENTS:
            print(f"  {pat:>8}", end="")
        print("  Dir")
        for bs in r["band_scores"]:
            band = bs["band"]
            print(f"  {band:<12}", end="")
            pat_dict = dict(bs["pat_gaps"])
            for pat in PATIENTS:
                if pat in pat_dict:
                    g = pat_dict[pat]
                    print(f"  {g:>+8.5f}", end="")
                else:
                    print(f"  {'N/A':>8}", end="")
            print(f"  {bs['direction']:>7} {'***' if bs['unanimous'] else ''}")

    # ---- Summary by family ----
    print("\n\n" + "=" * 90)
    print("FAMILY SUMMARY — best metric per family")
    print("=" * 90)
    families = {}
    for r in results:
        fam = family_map.get(r["metric"], "?")
        if fam not in families:
            families[fam] = r
    for fam, r in families.items():
        print(f"  {fam:<10} → {r['metric']:<20} {r['n_unanimous']}/6 unanimous, "
              f"{r['n_strong_unanimous']}/6 strong")

    # ==================================================================
    # ADDITIONAL ANALYSIS: Exclude Pat_07 (3 phases, asymmetric)
    # ==================================================================
    PATIENTS_4PH = PATIENTS_4PHASE
    print("\n\n" + "=" * 90)
    print("ANALYSIS 2: EXCLUDING Pat_07 (only 4-phase patients)")
    print("=" * 90)

    results2 = []
    for metric_name in METRICS:
        band_scores = []
        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS_4PH:
                key = (metric_name, pat, band)
                if key in gaps:
                    pat_gaps.append(gaps[key])
            if len(pat_gaps) < 3:
                continue
            all_pos = all(g > 0 for g in pat_gaps)
            all_neg = all(g < 0 for g in pat_gaps)
            unanimous = all_pos or all_neg
            mean_gap = np.mean(pat_gaps)
            band_scores.append({
                "band": band, "unanimous": unanimous,
                "direction": "reorg" if all_pos else ("persist" if all_neg else "mixed"),
                "mean_gap": mean_gap,
            })
        n_unan = sum(b["unanimous"] for b in band_scores)
        n_strong = sum(b["unanimous"] and abs(b["mean_gap"]) > 0.01 for b in band_scores)
        results2.append({"metric": metric_name, "n_unanimous": n_unan,
                         "n_strong": n_strong, "band_scores": band_scores})
    results2.sort(key=lambda r: (-r["n_unanimous"], -r["n_strong"]))

    print(f"{'Rank':<5} {'Metric':<20} {'Unan':>5} {'Strong':>7}")
    print("-" * 40)
    for i, r in enumerate(results2[:15]):
        print(f"{i+1:<5} {r['metric']:<20} {r['n_unanimous']:>5}/6 {r['n_strong']:>7}/6")

    # Detail for best without Pat_07
    for r in results2[:3]:
        print(f"\n  --- {r['metric']} (excl Pat_07): {r['n_unanimous']}/6 ---")
        for bs in r["band_scores"]:
            print(f"    {bs['band']:<12} {bs['direction']:>7} gap={bs['mean_gap']:>+.5f} "
                  f"{'***' if bs['unanimous'] else ''}")

    # ==================================================================
    # ANALYSIS 3: Specific pair comparisons
    # For each metric and band, check if specific orderings hold
    # across ALL 4-phase patients:
    #   Q1: rest-rest < rest-task?  (rest_pre,rest_post more similar than rest_pre,task_learn)
    #   Q2: task-task < rest-task?  (task_learn,task_test more similar than task_learn,rest_post)
    # ==================================================================
    print("\n\n" + "=" * 90)
    print("ANALYSIS 3: SPECIFIC PAIR ORDERINGS (4-phase patients only)")
    print("Does sim(rest_pre,rest_post) > sim(rest_pre,task_learn) for ALL patients?")
    print("Does sim(task_learn,task_test) > sim(task_learn,rest_post) for ALL patients?")
    print("=" * 90)

    pair_results = []
    for metric_name in METRICS:
        q1_band_count = 0  # bands where Q1 unanimous
        q2_band_count = 0
        band_details = []
        for band in BANDS:
            q1_ok = []
            q2_ok = []
            for pat in PATIENTS_4PH:
                rr = sims.get((metric_name, pat, band, ("rest_pre", "rest_post")))
                rt = sims.get((metric_name, pat, band, ("rest_pre", "task_learn")))
                tt = sims.get((metric_name, pat, band, ("task_learn", "task_test")))
                tr = sims.get((metric_name, pat, band, ("task_learn", "rest_post")))
                # Also check reverse pair names
                if rr is None:
                    rr = sims.get((metric_name, pat, band, ("rest_post", "rest_pre")))
                if rt is None:
                    rt = sims.get((metric_name, pat, band, ("task_learn", "rest_pre")))
                if tt is None:
                    tt = sims.get((metric_name, pat, band, ("task_test", "task_learn")))
                if tr is None:
                    tr = sims.get((metric_name, pat, band, ("rest_post", "task_learn")))

                if rr is not None and rt is not None:
                    q1_ok.append(rr > rt)
                if tt is not None and tr is not None:
                    q2_ok.append(tt > tr)

            q1_unan = len(q1_ok) >= 3 and all(q1_ok)
            q2_unan = len(q2_ok) >= 3 and all(q2_ok)
            if q1_unan:
                q1_band_count += 1
            if q2_unan:
                q2_band_count += 1
            band_details.append({
                "band": band,
                "q1_unan": q1_unan, "q1_count": sum(q1_ok) if q1_ok else 0,
                "q1_total": len(q1_ok),
                "q2_unan": q2_unan, "q2_count": sum(q2_ok) if q2_ok else 0,
                "q2_total": len(q2_ok),
            })
        pair_results.append({
            "metric": metric_name,
            "q1_bands": q1_band_count, "q2_bands": q2_band_count,
            "total": q1_band_count + q2_band_count,
            "band_details": band_details,
        })

    pair_results.sort(key=lambda r: -r["total"])
    print(f"\n{'Rank':<5} {'Metric':<20} {'Q1(rr>rt)':>10} {'Q2(tt>tr)':>10} {'Total':>6}")
    print("-" * 55)
    for i, r in enumerate(pair_results[:15]):
        print(f"{i+1:<5} {r['metric']:<20} {r['q1_bands']:>10}/6 "
              f"{r['q2_bands']:>10}/6 {r['total']:>6}/12")

    # Detail for top 3
    for r in pair_results[:3]:
        print(f"\n  --- {r['metric']} ---")
        print(f"    {'Band':<12} Q1(rr>rt) Q2(tt>tr)")
        for bd in r["band_details"]:
            q1s = f"{bd['q1_count']}/{bd['q1_total']}" + (" ***" if bd['q1_unan'] else "")
            q2s = f"{bd['q2_count']}/{bd['q2_total']}" + (" ***" if bd['q2_unan'] else "")
            print(f"    {bd['band']:<12} {q1s:<12} {q2s:<12}")

    # ==================================================================
    # ANALYSIS 4: Check actual similarity values for entropy metrics
    # (diagnostics — are they all ~1.0?)
    # ==================================================================
    print("\n\n" + "=" * 90)
    print("DIAGNOSTIC: Sample similarity values for key metrics")
    print("=" * 90)
    for metric_name in ["EntPearson", "EntCosine", "ComplPearson", "MSC_Pearson",
                         "SpecPearson", "MrgPearson"]:
        print(f"\n  {metric_name}:")
        for band in ["delta", "theta", "alpha"]:
            vals = []
            for pat in PATIENTS_4PH:
                for p1, p2 in combinations(PATIENT_PHASES[pat], 2):
                    key = (metric_name, pat, band, (p1, p2))
                    if key in sims:
                        vals.append(sims[key])
            if vals:
                print(f"    {band}: min={min(vals):.6f} max={max(vals):.6f} "
                      f"range={max(vals)-min(vals):.6f} mean={np.mean(vals):.6f}")


if __name__ == "__main__":
    main()
