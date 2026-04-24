#!/usr/bin/env python3
"""Graph Signal Processing analysis of epileptic node characterization.

For each patient, builds a graph Laplacian from the MSC FC matrix,
projects raw EEG time series onto the graph Fourier basis, and computes
per-node spectral metrics. Compares epileptic vs non-epileptic contacts
within the same electrode probes (mixed-probe design).

Metrics:
  - total_variation: weighted mean absolute difference from neighbours
  - hf_ratio: high-graph-freq energy / low-graph-freq energy
  - dirichlet_energy: quadratic form x^T L x decomposed per node
"""

import numpy as np
import pandas as pd
import re
import warnings
from pathlib import Path
from scipy.stats import mannwhitneyu

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
warnings.filterwarnings("ignore")

from lrg_eegfc.config.paths import MSC_CACHE, SEEG_DATAPATH
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.utils.io import load_epileptic_nodes, load_timeseries
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_ch(pat: str) -> list[str]:
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [str(l).strip('"').split(",")[0].strip().replace(" ", "") for l in df.iloc[:, 0]]


def probe(label: str) -> str:
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


# Patient -> sampling rate
PAT_FS = {
    "Pat_02": 2048,
    "Pat_03": 1024,
    "Pat_05": 2048,
    "Pat_07": 2048,
    "Pat_08": 2048,
}

PATIENTS = list(PAT_FS.keys())
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PHASES = ["rest_pre"]

# Collect all results for a summary table
all_results = []

for pat in PATIENTS:
    epi_set = set(load_epileptic_nodes(pat))
    ch = load_ch(pat)
    N = len(ch)
    epi_mask = np.array([l in epi_set for l in ch])

    # Build mixed-probe pools
    probes = {}
    for i, label in enumerate(ch):
        p = probe(label)
        if p not in probes:
            probes[p] = {"epi": [], "non": []}
        if epi_mask[i]:
            probes[p]["epi"].append(i)
        else:
            probes[p]["non"].append(i)

    mixed = {p: v for p, v in probes.items() if v["epi"] and v["non"]}
    if not mixed:
        print(f"\n{'='*70}")
        print(f"{pat}: no mixed probes (skip)")
        continue

    epi_pool = sum((v["epi"] for v in mixed.values()), [])
    non_pool = sum((v["non"] for v in mixed.values()), [])

    print(f"\n{'='*70}")
    print(f"{pat}  (N={N}, epi_nodes={sum(epi_mask)}, mixed_probes={len(mixed)})")
    print(f"  Mixed-probe epi contacts: {len(epi_pool)}, non-epi: {len(non_pool)}")

    fs = PAT_FS[pat]
    nperseg = nperseg_for_fs(fs)
    print(f"  fs={fs} Hz, nperseg={nperseg}")

    for phase in PHASES:
        # Load raw time series
        try:
            ts = load_timeseries(pat, phase, SEEG_DATAPATH)
        except Exception as e:
            print(f"  Cannot load timeseries for {phase}: {e}")
            continue

        print(f"\n  Phase: {phase}, timeseries shape: {ts.shape}")
        if ts.shape[0] != N:
            print(f"  Shape mismatch: {ts.shape[0]} vs {N} channels, skipping")
            continue

        for band in BANDS:
            A = load_msc_matrix(
                pat, phase, band,
                cache_root=MSC_CACHE,
                sparsify="none", n_surrogates=0, nperseg=nperseg,
            )
            if A is None:
                # Try with nperseg=4096 as fallback
                A = load_msc_matrix(
                    pat, phase, band,
                    cache_root=MSC_CACHE,
                    sparsify="none", n_surrogates=0, nperseg=4096,
                )
            if A is None or A.shape[0] != N:
                print(f"    {band}: no MSC matrix or shape mismatch, skipping")
                continue

            # Build graph Laplacian
            A = A.copy()
            np.fill_diagonal(A, 0)
            D = A.sum(axis=1)
            L = np.diag(D) - A

            # Eigen-decomposition
            ev, U = np.linalg.eigh(L)
            ev = np.maximum(ev, 0.0)
            if ev[1] < 1e-10:
                print(f"    {band}: graph disconnected (lambda_1 ~ 0), skipping")
                continue

            # Subsample time points for efficiency
            n_samples = min(1000, ts.shape[1])
            step = max(1, ts.shape[1] // n_samples)
            ts_sub = ts[:, ::step][:, :n_samples]

            # Graph Fourier Transform: (N_freq, n_samples)
            gft = U.T @ ts_sub

            n_freq = N - 1  # skip zero eigenvalue
            mid = n_freq // 2

            # --- Metric 1: Total Variation (vectorized) ---
            # TV_i = sum_j A_ij * mean_t |x_i - x_j|
            # Vectorize: for each pair (i,j) compute mean |x_i - x_j| * A_ij
            tv = np.zeros(N)
            for i in range(N):
                diff = np.abs(ts_sub[i:i+1, :] - ts_sub)  # (N, T)
                tv[i] = (A[i, :] * diff.mean(axis=1)).sum()
            tv_norm = tv / (D + 1e-30)

            # --- Metric 2: High-freq / Low-freq graph spectral energy ratio ---
            # Node i's energy at graph freq k = U[i,k]^2 * mean_t(gft[k,t]^2)
            gft_power = np.mean(gft[1:]**2, axis=1)  # (N-1,) skip k=0
            node_spec_energy = U[:, 1:]**2 * gft_power[np.newaxis, :]  # (N, N-1)
            low_energy = node_spec_energy[:, :mid].sum(axis=1)
            high_energy = node_spec_energy[:, mid:].sum(axis=1)
            hf_ratio = high_energy / (low_energy + 1e-30)

            # --- Metric 3: Dirichlet energy per node ---
            # E_i = sum_j L_ij * cov(x_i, x_j)
            cov_mat = ts_sub @ ts_sub.T / n_samples
            dirichlet_per_node = np.diag(L @ cov_mat)

            # --- Metric 4: Graph smoothness (Laplacian quadratic form per node) ---
            # s_i = sum_j A_ij (x_i - x_j)^2  averaged over time
            smooth = np.zeros(N)
            for i in range(N):
                diff2 = (ts_sub[i:i+1, :] - ts_sub)**2  # (N, T)
                smooth[i] = (A[i, :] * diff2.mean(axis=1)).sum()
            smooth_norm = smooth / (D + 1e-30)

            # --- Metric 5: Graph spectral centroid per node ---
            # Centroid_i = sum_k lambda_k * E_ik / sum_k E_ik
            # where E_ik = U[i,k]^2 * mean(gft[k]^2)
            total_energy_per_node = node_spec_energy.sum(axis=1)
            spectral_centroid = (node_spec_energy * ev[1:][np.newaxis, :]).sum(axis=1) / (total_energy_per_node + 1e-30)

            metrics = {
                "total_variation": tv_norm,
                "hf_ratio": hf_ratio,
                "dirichlet_energy": dirichlet_per_node,
                "laplacian_smooth": smooth_norm,
                "spectral_centroid": spectral_centroid,
            }

            print(f"\n    {band}:")
            for mname, mvals in metrics.items():
                ve = mvals[epi_pool]
                vn = mvals[non_pool]
                _, pg = mannwhitneyu(ve, vn, alternative="greater")
                _, pl = mannwhitneyu(ve, vn, alternative="less")
                p_best = min(pg, pl)
                direction = "epi>" if pg < pl else "epi<"
                sig = "***" if p_best < 0.001 else "**" if p_best < 0.01 else "*" if p_best < 0.05 else ""
                eff = (ve.mean() - vn.mean()) / (np.std(np.concatenate([ve, vn])) + 1e-30)
                print(
                    f"      {mname:20s}: epi={ve.mean():.4f}  non={vn.mean():.4f}  "
                    f"{direction}  p={p_best:.4f}  d={eff:+.2f}  {sig}"
                )
                all_results.append({
                    "patient": pat,
                    "phase": phase,
                    "band": band,
                    "metric": mname,
                    "epi_mean": ve.mean(),
                    "non_mean": vn.mean(),
                    "direction": direction,
                    "p_value": p_best,
                    "effect_size_d": eff,
                    "n_epi": len(epi_pool),
                    "n_non": len(non_pool),
                })

# ---------------------------------------------------------------------------
# Summary across patients
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY: Significant results (p < 0.05)")
print("=" * 70)

df = pd.DataFrame(all_results)
sig = df[df["p_value"] < 0.05].sort_values(["metric", "band", "patient"])
if len(sig) == 0:
    print("  No significant results found.")
else:
    for _, row in sig.iterrows():
        print(
            f"  {row['patient']:8s} {row['band']:12s} {row['metric']:20s} "
            f"{row['direction']}  p={row['p_value']:.4f}  d={row['effect_size_d']:+.2f}"
        )

# Cross-patient consistency
print("\n" + "=" * 70)
print("CONSISTENCY: Metrics with same direction across >= 3 patients")
print("=" * 70)

for band in BANDS:
    for metric in df["metric"].unique():
        sub = df[(df["band"] == band) & (df["metric"] == metric)]
        if len(sub) < 3:
            continue
        dirs = sub["direction"].value_counts()
        dominant = dirs.index[0]
        count = dirs.iloc[0]
        if count >= 3:
            mean_p = sub["p_value"].mean()
            mean_d = sub["effect_size_d"].mean()
            n_sig = (sub["p_value"] < 0.05).sum()
            print(
                f"  {band:12s} {metric:20s}: {dominant} in {count}/{len(sub)} patients, "
                f"mean_p={mean_p:.3f}, mean_d={mean_d:+.2f}, sig_in={n_sig}/{len(sub)}"
            )

# Cross-patient aggregated test (pool all patients' effect sizes)
print("\n" + "=" * 70)
print("AGGREGATED: Sign-consistency test per metric per band")
print("=" * 70)
from scipy.stats import binomtest, combine_pvalues

for metric in df["metric"].unique():
    print(f"\n  {metric}:")
    for band in BANDS:
        sub = df[(df["band"] == band) & (df["metric"] == metric)]
        if len(sub) < 3:
            continue
        n_greater = (sub["direction"] == "epi>").sum()
        n_total = len(sub)
        # Use Fisher's method to combine p-values for same-direction tests
        # Pick the dominant direction
        if n_greater >= n_total / 2:
            # Test epi > non across patients
            ps = []
            for _, row in sub.iterrows():
                # We need directional p-values; reconstruct from stored
                ps.append(row["p_value"] if row["direction"] == "epi>" else 1 - row["p_value"])
            direction_label = "epi>"
        else:
            ps = []
            for _, row in sub.iterrows():
                ps.append(row["p_value"] if row["direction"] == "epi<" else 1 - row["p_value"])
            direction_label = "epi<"

        # Clip to avoid 0/1 for combine_pvalues
        ps = np.clip(ps, 1e-10, 1 - 1e-10)
        try:
            _, combined_p = combine_pvalues(ps, method="fisher")
        except Exception:
            combined_p = np.nan

        sig_mark = "***" if combined_p < 0.001 else "**" if combined_p < 0.01 else "*" if combined_p < 0.05 else ""
        mean_d = sub["effect_size_d"].mean()
        print(
            f"    {band:12s}: {direction_label} {n_greater}/{n_total} patients, "
            f"Fisher_p={combined_p:.4f} {sig_mark}  mean_d={mean_d:+.2f}"
        )

print("\nDone.")
