#!/usr/bin/env python3
"""Focused sweep of ultrametric D(tau) comparison metrics.

Compares the LRG ultrametric distance matrices between phases.
ALL 6 patients included. Both standard MSC and CReMa-validated MSC.

Key design choices:
  - Distances are in log-scale (from diffusion exp(-tau*L))
  - Cophenetic distances have only N-1 unique values (discrete)
  - Same nodes across phases for a given patient
  - MSC Spearman correlation as control metric

Metrics focus on:
  A) Log-space comparisons of D(tau)
  B) Scale-specific comparisons (fine vs coarse structure)
  C) Node-level stability features
  D) MSC matrix control
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
from scipy.cluster.hierarchy import cophenet, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr, wasserstein_distance
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import networkx as nx

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, LRG_CREMA_CACHE, MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result, LRGResult
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrgsglib.utils.lrg.spectral import compute_laplacian_properties, get_graph_lspectrum
from lrgsglib.utils.lrg.clustering import compute_normalized_linkage, compute_optimal_threshold
from lrgsglib.utils.lrg.infocomm import entropy, extract_ultrametric_matrix

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
LRG_CACHE = _LRG_CACHE

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


# ===================================================================
# Load CReMa-validated MSC and compute LRG
# ===================================================================
def load_crema_msc(pat, ph, band):
    """Load CReMa-validated MSC matrix from cache."""
    pattern = f"{band}_{ph}_msc_sparsify-ecm_adaptive_alpharange-0.01-0.5_nens-100_wscale-1000_nperseg-4096.npy"
    path = MSC_CACHE / pat / pattern
    if path.exists():
        return np.load(path)
    return None


CREMA_LRG_CACHE = {}  # in-memory cache for CReMa LRG results


def get_giant_component(adj):
    """Extract giant component from adjacency matrix."""
    G = nx.from_numpy_array(adj)
    if not nx.is_connected(G):
        cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(cc).copy()
        G = nx.convert_node_labels_to_integers(G)
    return G


def compute_lrg_from_adj(adj, pat, ph, band):
    """Compute LRG analysis directly from adjacency matrix (no fc_method check)."""
    G = get_giant_component(adj)
    n_nodes = G.number_of_nodes()
    if n_nodes < 5:
        return None

    # Entropy
    sm1, spec, _, tau = entropy(G, steps=400, t1=-3.0, t2=5.0)
    ent_C = np.abs(np.diff(sm1) / np.diff(np.log10(tau)))

    # Diffusion distances
    _, _, _, Trho, _ = compute_laplacian_properties(G)
    dists = squareform(Trho)

    # Linkage
    Z, labels, tmax = compute_normalized_linkage(dists, G)

    # Optimal threshold
    threshold, *_ = compute_optimal_threshold(Z)

    # Ultrametric
    U_sq = extract_ultrametric_matrix(Z, n_nodes)
    U_cond = squareform(U_sq)

    return LRGResult(
        ultrametric_matrix=U_cond,
        linkage_matrix=Z,
        entropy_tau=tau,
        entropy_1_minus_S=sm1,
        entropy_C=ent_C,
        optimal_threshold=threshold,
        patient=pat,
        phase=ph,
        band=band,
        fc_method="msc_crema",
        n_nodes=n_nodes,
    )


def ensure_crema_lrg(pat, ph, band):
    """Compute LRG from CReMa MSC matrix (cached in memory)."""
    cache_key = (pat, ph, band)
    if cache_key in CREMA_LRG_CACHE:
        return CREMA_LRG_CACHE[cache_key]

    # Also check disk cache
    cache_path = LRG_CACHE / pat / f"{band}_{ph}_lrg_msc_crema.npz"
    if cache_path.exists():
        lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
        # Won't match — need custom loading
        pass

    msc = load_crema_msc(pat, ph, band)
    if msc is None:
        CREMA_LRG_CACHE[cache_key] = None
        return None

    try:
        lrg = compute_lrg_from_adj(msc, pat, ph, band)
        CREMA_LRG_CACHE[cache_key] = lrg
        return lrg
    except Exception as e:
        print(f"    FAIL {pat}/{ph}/{band} CReMa LRG: {e}")
        CREMA_LRG_CACHE[cache_key] = None
        return None


# ===================================================================
# Data loading
# ===================================================================
def load_all(fc_methods=("msc", "msc_crema")):
    """Load ultrametric data for all patients/phases/bands."""
    data = {}

    for fc in fc_methods:
        print(f"\n  Loading fc_method={fc}...")
        count = 0
        for pat, phases in PATIENT_PHASES.items():
            for band in BANDS:
                for ph in phases:
                    if fc == "msc_crema":
                        lrg = ensure_crema_lrg(pat, ph, band)
                    else:
                        lrg = load_lrg_result(pat, ph, band, fc, cache_root=LRG_CACHE)

                    if lrg is None:
                        continue

                    Z = lrg.linkage_matrix
                    U_cond = lrg.ultrametric_matrix
                    # Make sure it's condensed
                    if U_cond.ndim == 2:
                        U_cond = squareform(U_cond)
                    U_sq = squareform(U_cond)
                    n = lrg.n_nodes
                    opt_thresh = lrg.optimal_threshold

                    # Log-transformed ultrametric (clip zeros)
                    U_log = np.log(np.clip(U_cond, 1e-15, None))

                    # MSC matrix for control (only for standard MSC)
                    msc_upper = None
                    if fc == "msc":
                        msc = load_msc_matrix(pat, ph, band, cache_root=MSC_CACHE)
                        if msc is not None:
                            mat = msc if isinstance(msc, np.ndarray) else msc.adjacency_matrix
                            idx = np.triu_indices(mat.shape[0], k=1)
                            msc_upper = mat[idx]
                    elif fc == "msc_crema":
                        msc = load_crema_msc(pat, ph, band)
                        if msc is not None:
                            idx = np.triu_indices(msc.shape[0], k=1)
                            msc_upper = msc[idx]

                    data[(fc, pat, ph, band)] = {
                        "Z": Z,
                        "U_cond": U_cond,
                        "U_sq": U_sq,
                        "U_log": U_log,
                        "n": n,
                        "opt_thresh": opt_thresh,
                        "msc_upper": msc_upper,
                    }
                    count += 1
        print(f"    → {count} entries loaded")

    return data


# ===================================================================
# METRICS: all return SIMILARITY (higher = more similar)
# ===================================================================

# --- A: Log-space comparisons of D(tau) ---

def m_log_spearman(d1, d2):
    """Spearman of log(D). Same as Spearman of D (rank-invariant)."""
    return spearmanr(d1["U_cond"], d2["U_cond"])[0]


def m_log_pearson(d1, d2):
    """Pearson of log(D). Weights fine structure more."""
    return pearsonr(d1["U_log"], d2["U_log"])[0]


def m_log_cosine(d1, d2):
    """Cosine similarity in log space."""
    a, b = d1["U_log"], d2["U_log"]
    dot = np.dot(a, b)
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return dot / max(norm, 1e-12)


def m_log_l2_sim(d1, d2):
    """1/(1 + RMSE of log D). Sensitive to absolute scale."""
    diff = d1["U_log"] - d2["U_log"]
    rmse = np.sqrt(np.mean(diff ** 2))
    return 1.0 / (1.0 + rmse)


def m_log_sup_sim(d1, d2):
    """1/(1 + max|logD1 - logD2|). Worst-case deviation."""
    diff = np.abs(d1["U_log"] - d2["U_log"])
    return 1.0 / (1.0 + np.max(diff))


def m_raw_pearson(d1, d2):
    """Pearson of raw D (for comparison)."""
    return pearsonr(d1["U_cond"], d2["U_cond"])[0]


# --- B: Scale-specific comparisons ---

def _split_by_scale(U, thresholds=(0.1, 0.3, 0.6)):
    """Split ultrametric entries by value into scale bands."""
    bands = {}
    prev = 0
    for i, t in enumerate(thresholds):
        mask = (U >= prev) & (U < t)
        bands[f"scale_{i}"] = mask
        prev = t
    bands[f"scale_{len(thresholds)}"] = U >= thresholds[-1]
    return bands


def m_fine_scale_spearman(d1, d2):
    """Spearman of D for fine-scale pairs only (D < 0.1)."""
    mask = (d1["U_cond"] < 0.1) & (d2["U_cond"] < 0.1)
    if mask.sum() < 10:
        return None
    return spearmanr(d1["U_cond"][mask], d2["U_cond"][mask])[0]


def m_coarse_scale_spearman(d1, d2):
    """Spearman of D for coarse-scale pairs only (D > 0.5)."""
    mask = (d1["U_cond"] > 0.5) & (d2["U_cond"] > 0.5)
    if mask.sum() < 10:
        return None
    return spearmanr(d1["U_cond"][mask], d2["U_cond"][mask])[0]


def m_mid_scale_spearman(d1, d2):
    """Spearman of D for mid-scale pairs (0.1 < D < 0.5)."""
    mask = ((d1["U_cond"] > 0.1) & (d1["U_cond"] < 0.5) &
            (d2["U_cond"] > 0.1) & (d2["U_cond"] < 0.5))
    if mask.sum() < 10:
        return None
    return spearmanr(d1["U_cond"][mask], d2["U_cond"][mask])[0]


def m_log_fine_pearson(d1, d2):
    """Pearson of log(D) for fine-scale pairs (D < median)."""
    med = np.median(np.concatenate([d1["U_cond"], d2["U_cond"]]))
    mask = (d1["U_cond"] < med) & (d2["U_cond"] < med)
    if mask.sum() < 10:
        return None
    return pearsonr(d1["U_log"][mask], d2["U_log"][mask])[0]


def m_log_coarse_pearson(d1, d2):
    """Pearson of log(D) for coarse-scale pairs (D > median)."""
    med = np.median(np.concatenate([d1["U_cond"], d2["U_cond"]]))
    mask = (d1["U_cond"] > med) & (d2["U_cond"] > med)
    if mask.sum() < 10:
        return None
    return pearsonr(d1["U_log"][mask], d2["U_log"][mask])[0]


# --- C: Node-level stability features ---

def m_node_profile_log_pearson(d1, d2):
    """Pearson of per-node mean log(D) profiles.

    For node i: mu_i = mean_j(log D_ij).
    Compare mu vectors across phases.
    """
    n = d1["n"]
    logD1 = np.log(np.clip(d1["U_sq"], 1e-15, None))
    logD2 = np.log(np.clip(d2["U_sq"], 1e-15, None))
    mu1 = logD1.mean(axis=1)
    mu2 = logD2.mean(axis=1)
    return pearsonr(mu1, mu2)[0]


def m_node_profile_log_spearman(d1, d2):
    """Spearman of per-node mean log(D) profiles."""
    logD1 = np.log(np.clip(d1["U_sq"], 1e-15, None))
    logD2 = np.log(np.clip(d2["U_sq"], 1e-15, None))
    mu1 = logD1.mean(axis=1)
    mu2 = logD2.mean(axis=1)
    return spearmanr(mu1, mu2)[0]


def m_node_rank_stability(d1, d2, k=10):
    """Mean Jaccard overlap of k-nearest neighbors per node.

    For each node i, find k closest nodes in D. Compare
    neighborhood sets across phases.
    """
    n = d1["n"]
    k_use = min(k, n - 1)
    total_overlap = 0
    for i in range(n):
        nn1 = set(np.argsort(d1["U_sq"][i])[:k_use + 1]) - {i}
        nn2 = set(np.argsort(d2["U_sq"][i])[:k_use + 1]) - {i}
        if nn1 or nn2:
            total_overlap += len(nn1 & nn2) / len(nn1 | nn2)
    return total_overlap / n


def m_node_rank_stability_k5(d1, d2):
    return m_node_rank_stability(d1, d2, k=5)


def m_node_rank_stability_k20(d1, d2):
    return m_node_rank_stability(d1, d2, k=20)


# --- D: Partition-based (at optimal threshold) ---

def m_ari_optimal(d1, d2):
    """ARI of partitions at optimal threshold (each uses own threshold)."""
    c1 = fcluster(d1["Z"], d1["opt_thresh"], criterion='distance')
    c2 = fcluster(d2["Z"], d2["opt_thresh"], criterion='distance')
    return adjusted_rand_score(c1, c2)


def m_nmi_optimal(d1, d2):
    """NMI of partitions at optimal threshold."""
    c1 = fcluster(d1["Z"], d1["opt_thresh"], criterion='distance')
    c2 = fcluster(d2["Z"], d2["opt_thresh"], criterion='distance')
    return normalized_mutual_info_score(c1, c2)


def m_ari_fixed_3(d1, d2):
    """ARI cutting at 3 clusters."""
    c1 = fcluster(d1["Z"], 3, criterion='maxclust')
    c2 = fcluster(d2["Z"], 3, criterion='maxclust')
    return adjusted_rand_score(c1, c2)


def m_ari_fixed_5(d1, d2):
    """ARI cutting at 5 clusters."""
    c1 = fcluster(d1["Z"], 5, criterion='maxclust')
    c2 = fcluster(d2["Z"], 5, criterion='maxclust')
    return adjusted_rand_score(c1, c2)


def m_ari_fixed_10(d1, d2):
    """ARI cutting at 10 clusters."""
    c1 = fcluster(d1["Z"], 10, criterion='maxclust')
    c2 = fcluster(d2["Z"], 10, criterion='maxclust')
    return adjusted_rand_score(c1, c2)


# --- E: Distribution-based ---

def m_merge_log_wasserstein(d1, d2):
    """1/(1 + Wasserstein of log merge heights)."""
    h1 = np.log(np.clip(d1["Z"][:, 2], 1e-15, None))
    h2 = np.log(np.clip(d2["Z"][:, 2], 1e-15, None))
    return 1.0 / (1.0 + wasserstein_distance(h1, h2))


def m_log_quantile_sim(d1, d2):
    """Quantile profile similarity of log(D)."""
    qs = np.linspace(0.05, 0.95, 19)
    q1 = np.quantile(d1["U_log"], qs)
    q2 = np.quantile(d2["U_log"], qs)
    r = max(np.ptp(q1), np.ptp(q2), 1e-12)
    rmse = np.sqrt(np.mean((q1 - q2) ** 2)) / r
    return 1.0 - min(rmse, 1.0)


# --- F: MSC control ---

def m_msc_spearman(d1, d2):
    """Spearman of vectorized MSC matrices (control)."""
    if d1["msc_upper"] is None or d2["msc_upper"] is None:
        return None
    return spearmanr(d1["msc_upper"], d2["msc_upper"])[0]


# ===================================================================
# Registry
# ===================================================================
METRICS = {
    # A: Log-space D(tau)
    "logPearson":       m_log_pearson,
    "logCosine":        m_log_cosine,
    "logL2Sim":         m_log_l2_sim,
    "logSupSim":        m_log_sup_sim,
    "rawPearson":       m_raw_pearson,
    "Spearman":         m_log_spearman,
    # B: Scale-specific
    "fineSpearman":     m_fine_scale_spearman,
    "midSpearman":      m_mid_scale_spearman,
    "coarseSpearman":   m_coarse_scale_spearman,
    "logFinePearson":   m_log_fine_pearson,
    "logCoarsePearson": m_log_coarse_pearson,
    # C: Node-level
    "nodeLogPearson":   m_node_profile_log_pearson,
    "nodeLogSpearman":  m_node_profile_log_spearman,
    "kNN5_Jaccard":     m_node_rank_stability_k5,
    "kNN10_Jaccard":    m_node_rank_stability,
    "kNN20_Jaccard":    m_node_rank_stability_k20,
    # D: Partition-based
    "ARI_optimal":      m_ari_optimal,
    "NMI_optimal":      m_nmi_optimal,
    "ARI_k3":           m_ari_fixed_3,
    "ARI_k5":           m_ari_fixed_5,
    "ARI_k10":          m_ari_fixed_10,
    # E: Distribution
    "mrgLogWass":       m_merge_log_wasserstein,
    "logQuantileSim":   m_log_quantile_sim,
    # F: Control
    "MSC_Spearman":     m_msc_spearman,
}


# ===================================================================
# Computation
# ===================================================================
def compute_all(data, fc_method):
    """Compute gaps and raw similarities for one FC method."""
    gaps = {}
    sims = {}

    for metric_name, metric_fn in METRICS.items():
        print(f"  {metric_name}...", end="", flush=True)
        n_ok = 0
        for pat in PATIENTS:
            phases = PATIENT_PHASES[pat]
            for band in BANDS:
                available = [ph for ph in phases
                             if (fc_method, pat, ph, band) in data]
                pairs = list(combinations(available, 2))

                within_sims = []
                cross_sims = []
                for p1, p2 in pairs:
                    d1 = data[(fc_method, pat, p1, band)]
                    d2 = data[(fc_method, pat, p2, band)]
                    try:
                        s = metric_fn(d1, d2)
                        if s is None or np.isnan(s):
                            continue
                    except Exception:
                        continue

                    n_ok += 1
                    cat = classify_pair(p1, p2)
                    sims[(metric_name, pat, band, p1, p2)] = s
                    if cat == "within":
                        within_sims.append(s)
                    else:
                        cross_sims.append(s)

                if within_sims and cross_sims:
                    gap = np.mean(within_sims) - np.mean(cross_sims)
                    gaps[(metric_name, pat, band)] = gap

        print(f" ({n_ok})")

    return gaps, sims


def score_and_report(gaps, sims, fc_label):
    """Score metrics and print report."""
    print(f"\n{'='*90}")
    print(f"RESULTS: {fc_label}")
    print(f"{'='*90}")

    results = []
    for metric_name in METRICS:
        band_scores = []
        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS:
                key = (metric_name, pat, band)
                if key in gaps:
                    pat_gaps.append((pat, gaps[key]))

            if len(pat_gaps) < 2:
                continue

            gap_vals = [g for _, g in pat_gaps]
            all_pos = all(g > 0 for g in gap_vals)
            all_neg = all(g < 0 for g in gap_vals)
            unanimous = all_pos or all_neg
            mean_gap = np.mean(gap_vals)

            band_scores.append({
                "band": band,
                "unanimous": unanimous,
                "direction": "reorg" if all_pos else ("persist" if all_neg else "mixed"),
                "mean_gap": mean_gap,
                "n_patients": len(gap_vals),
                "pat_gaps": pat_gaps,
            })

        n_unan = sum(b["unanimous"] for b in band_scores)
        n_strong = sum(b["unanimous"] and abs(b["mean_gap"]) > 0.01
                       for b in band_scores)
        results.append({
            "metric": metric_name,
            "n_unanimous": n_unan,
            "n_strong": n_strong,
            "band_scores": band_scores,
        })

    results.sort(key=lambda r: (-r["n_unanimous"], -r["n_strong"]))

    # Print ranking
    print(f"\n{'Rank':<5} {'Metric':<20} {'Unan':>5} {'Strong':>7}")
    print("-" * 40)
    for i, r in enumerate(results):
        print(f"{i+1:<5} {r['metric']:<20} {r['n_unanimous']:>5}/6 "
              f"{r['n_strong']:>7}/6")

    # Detailed top 5
    print(f"\n--- Top 5 detail ---")
    for r in results[:5]:
        print(f"\n  {r['metric']} ({r['n_unanimous']}/6 unanimous)")
        for bs in r["band_scores"]:
            mark = "***" if bs["unanimous"] else ""
            print(f"    {bs['band']:<12} {bs['direction']:>7} "
                  f"gap={bs['mean_gap']:>+.5f} (N={bs['n_patients']}) {mark}")

    # Per-patient detail for top 3
    print(f"\n--- Per-patient gaps (top 3) ---")
    for r in results[:3]:
        print(f"\n  {r['metric']}:")
        header = f"    {'Band':<12}"
        for pat in PATIENTS:
            header += f" {pat.replace('Pat_0', 'P'):>7}"
        header += "  Dir"
        print(header)
        for bs in r["band_scores"]:
            row = f"    {bs['band']:<12}"
            pd = dict(bs["pat_gaps"])
            for pat in PATIENTS:
                if pat in pd:
                    row += f" {pd[pat]:>+7.4f}"
                else:
                    row += f" {'N/A':>7}"
            row += f"  {bs['direction']:>7}"
            if bs["unanimous"]:
                row += " ***"
            print(row)

    # Q2 test: sim(task_learn,task_test) > sim(task_learn,rest_post)
    PATIENTS_4PH = PATIENTS_4PHASE
    print(f"\n--- Q2 test: task-task > task-rest (4-phase patients) ---")
    for metric_name in list(METRICS.keys())[:10]:  # check first 10
        q2_bands = 0
        for band in BANDS:
            q2_ok = []
            for pat in PATIENTS_4PH:
                tt = sims.get((metric_name, pat, band, "task_learn", "task_test"))
                tr = sims.get((metric_name, pat, band, "task_learn", "rest_post"))
                if tt is not None and tr is not None:
                    q2_ok.append(tt > tr)
            if len(q2_ok) >= 3 and all(q2_ok):
                q2_bands += 1
        if q2_bands > 0:
            print(f"  {metric_name:<20} Q2 unanimous: {q2_bands}/6 bands")

    return results


# ===================================================================
# Main
# ===================================================================
def main():
    print("=" * 90)
    print("Ultrametric D(tau) metric sweep — MSC + CReMa")
    print("=" * 90)

    print("\nLoading data...")
    data = load_all(fc_methods=("msc", "msc_crema"))

    for fc, label in [("msc", "Standard MSC"), ("msc_crema", "CReMa-validated MSC")]:
        print(f"\n\nComputing metrics for {label}...")
        gaps, sims = compute_all(data, fc)
        score_and_report(gaps, sims, label)


if __name__ == "__main__":
    main()
