#!/usr/bin/env python3
"""Cross-frequency analysis: do epileptic contacts show different
relationships across frequency bands?

Uses EXISTING cached MSC matrices (sparsify=none). For each node, computes:
1. Cross-band connectivity profile correlation (Spearman of MSC rows across bands)
2. Strength rank variability (CV of rank across bands)
3. Band entropy (uniformity of strength distribution across bands)
4. Delta/beta ratio (pathological coupling marker)
5. Per-band-pair correlations (to find which pair differs most)

Within-probe design: only compares epi vs non-epi contacts on the same probe.
"""

import numpy as np
import pandas as pd
import re, warnings
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, rankdata

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
warnings.filterwarnings("ignore")

from lrg_eegfc.config.paths import MSC_CACHE, SEEG_DATAPATH
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.utils.io import load_epileptic_nodes


def load_ch(pat):
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "") for l in df.iloc[:, 0]
    ]


def probe(label):
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


patients = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
bands = ["delta", "theta", "alpha", "beta"]
phases = ["rest_pre", "rest_post"]

# Pat_03 recorded at 1024 Hz -> nperseg=2048; others at 2048 Hz -> nperseg=4096
nperseg_map = {p: 4096 for p in patients}
nperseg_map["Pat_03"] = 2048

summary_rows = []

for pat in patients:
    epi_set = set(load_epileptic_nodes(pat))
    ch = load_ch(pat)
    N = len(ch)
    epi_mask = np.array([l in epi_set for l in ch])

    probes_dict = {}
    for i, label in enumerate(ch):
        p = probe(label)
        if p not in probes_dict:
            probes_dict[p] = {"epi": [], "non": []}
        if epi_mask[i]:
            probes_dict[p]["epi"].append(i)
        else:
            probes_dict[p]["non"].append(i)
    mixed = {p: v for p, v in probes_dict.items() if v["epi"] and v["non"]}
    if not mixed:
        print(f"{pat}: No mixed probes, skipping")
        continue
    epi_pool = sum((v["epi"] for v in mixed.values()), [])
    non_pool = sum((v["non"] for v in mixed.values()), [])

    print(f"\n{'='*70}")
    print(
        f"{pat}  (N={N}, epi_nodes={sum(epi_mask)}, mixed_probes={len(mixed)})"
    )
    print(
        f"  Within-probe: {len(epi_pool)} epi contacts vs {len(non_pool)} non-epi contacts"
    )

    nperseg = nperseg_map[pat]

    for phase in phases:
        # Load all band matrices
        As = {}
        for band in bands:
            A = load_msc_matrix(
                pat,
                phase,
                band,
                cache_root=MSC_CACHE,
                sparsify="none",
                n_surrogates=0,
                nperseg=nperseg,
            )
            if A is not None and A.shape[0] == N:
                A = A.copy()
                np.fill_diagonal(A, 0)
                As[band] = A

        if len(As) < 2:
            print(f"  {phase}: insufficient band data ({len(As)} bands)")
            continue

        # 1. Cross-band connectivity profile correlation per node
        band_pairs = [
            (b1, b2)
            for b1 in bands
            for b2 in bands
            if b1 < b2 and b1 in As and b2 in As
        ]

        cross_band_corr = np.zeros(N)
        for i in range(N):
            rhos = []
            for b1, b2 in band_pairs:
                r, _ = spearmanr(As[b1][i], As[b2][i])
                if not np.isnan(r):
                    rhos.append(r)
            cross_band_corr[i] = np.mean(rhos) if rhos else 0

        # 2. Strength rank variability across bands
        strength_ranks = {}
        for band in As:
            s = As[band].sum(axis=1)
            strength_ranks[band] = rankdata(s)

        rank_cv = np.zeros(N)
        for i in range(N):
            ranks_i = [strength_ranks[b][i] for b in strength_ranks]
            rank_cv[i] = np.std(ranks_i) / (np.mean(ranks_i) + 1e-30)

        # 3. Band entropy of strength distribution
        strengths_arr = np.array(
            [As[b].sum(axis=1) for b in sorted(As.keys())]
        )  # (n_bands, N)
        str_norm = strengths_arr / (
            strengths_arr.sum(axis=0, keepdims=True) + 1e-30
        )
        band_entropy = -np.sum(str_norm * np.log(str_norm + 1e-30), axis=0)

        # 4. Delta/beta ratio per node
        if "delta" in As and "beta" in As:
            db_ratio = As["delta"].sum(axis=1) / (
                As["beta"].sum(axis=1) + 1e-30
            )
        else:
            db_ratio = np.zeros(N)

        # 5. Per-band-pair cross-band correlation
        pair_corrs = {}
        for b1, b2 in band_pairs:
            pc = np.zeros(N)
            for i in range(N):
                r, _ = spearmanr(As[b1][i], As[b2][i])
                pc[i] = r if not np.isnan(r) else 0
            pair_corrs[f"{b1}-{b2}"] = pc

        metrics = {
            "cross_band_corr": cross_band_corr,
            "rank_CV": rank_cv,
            "band_entropy": band_entropy,
            "delta_beta_ratio": db_ratio,
        }
        metrics.update(pair_corrs)

        print(
            f"\n  {phase} (within-probe, {len(epi_pool)} epi vs {len(non_pool)} non-epi):"
        )
        for mname, mvals in metrics.items():
            ve = mvals[epi_pool]
            vn = mvals[non_pool]
            if len(ve) < 2 or len(vn) < 2:
                continue
            _, pg = mannwhitneyu(ve, vn, alternative="greater")
            _, pl = mannwhitneyu(ve, vn, alternative="less")
            p_best = min(pg, pl)
            d = "epi>" if pg < pl else "epi<"
            sig = (
                "***"
                if p_best < 0.001
                else "**" if p_best < 0.01 else "*" if p_best < 0.05 else ""
            )
            print(
                f"    {mname:20s}: epi={ve.mean():.4f} non={vn.mean():.4f} "
                f"{d} p={p_best:.4f} {sig}"
            )
            summary_rows.append(
                {
                    "patient": pat,
                    "phase": phase,
                    "metric": mname,
                    "epi_mean": ve.mean(),
                    "non_mean": vn.mean(),
                    "direction": d,
                    "p_value": p_best,
                    "sig": sig,
                }
            )

# Summary: which metrics are consistent across patients?
print(f"\n\n{'='*70}")
print("CROSS-PATIENT CONSISTENCY SUMMARY")
print("=" * 70)
df = pd.DataFrame(summary_rows)
for metric in df["metric"].unique():
    sub = df[df["metric"] == metric]
    n_sig = (sub["p_value"] < 0.05).sum()
    n_total = len(sub)
    # Check direction consistency
    dirs = sub[sub["p_value"] < 0.05]["direction"]
    dir_counts = dirs.value_counts() if len(dirs) > 0 else pd.Series(dtype=int)
    dominant_dir = dir_counts.index[0] if len(dir_counts) > 0 else "N/A"
    consistency = dir_counts.iloc[0] / len(dirs) if len(dirs) > 0 else 0

    # Patients with sig results
    sig_patients = sub[sub["p_value"] < 0.05]["patient"].unique()

    print(
        f"  {metric:20s}: {n_sig}/{n_total} sig  "
        f"dir={dominant_dir}  consistency={consistency:.0%}  "
        f"patients={list(sig_patients)}"
    )

# Extra: per-patient summary of which metrics are significant in BOTH phases
print(f"\n{'='*70}")
print("PER-PATIENT: METRICS SIGNIFICANT IN BOTH rest_pre AND rest_post")
print("=" * 70)
for pat in patients:
    sub = df[df["patient"] == pat]
    if len(sub) == 0:
        print(f"  {pat}: skipped (no mixed probes)")
        continue
    both_sig = []
    for metric in sub["metric"].unique():
        msub = sub[sub["metric"] == metric]
        if all(msub["p_value"] < 0.05) and len(msub) == 2:
            dirs = msub["direction"].unique()
            both_sig.append(f"{metric}({dirs[0]})")
    if both_sig:
        print(f"  {pat}: {'  '.join(both_sig)}")
    else:
        print(f"  {pat}: none")

print("\nDone.")
