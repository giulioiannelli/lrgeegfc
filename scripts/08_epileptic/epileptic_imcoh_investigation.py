#!/usr/bin/env python3
"""Epileptic nodes in the ImCoh-LRG multiscale framework — deep investigation.

ImCoh (imaginary coherence) removes volume conduction, eliminating the
same-probe bias that dominated the MSC-based analysis (March 2026).
This enables, for the first time, clean community-level and ultrametric
analyses of epileptic node connectivity.

Questions:
  Q1. Are epileptic contacts still weaker in ImCoh? (within-probe control)
  Q2. Do epileptic nodes cluster in ImCoh LRG communities? (impossible with MSC)
  Q3. Are epileptic pairs closer in the ImCoh ultrametric? (impossible with MSC)
  Q4. Is the weakness scale-amplified through the ImCoh heat kernel?
  Q5. Are epileptic nodes more/less metastable across LRG scales?
  Q6. Is there a cross-patient consensus across phases and bands?

Run:  python scripts/py/epileptic_imcoh_investigation.py [-v]
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from scipy.stats import binomtest, mannwhitneyu, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import (
    FIGURES_ROOT, IMCOH_CACHE, IMCOH_LRG_CACHE, SEEG_DATAPATH,
)
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.workflow.lrg import LRGResult

# ── Config ─────────────────────────────────────────────────────────────

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)
N_COMMS = [3, 5, 10, 15, 20, 30]  # community granularities to test

DATASET_ROOT = SEEG_DATAPATH

OUT = FIGURES_ROOT / "epileptic_imcoh_analysis"
OUT.mkdir(parents=True, exist_ok=True)

FS_MAP = {"Pat_03": 1024.0}
np.random.seed(42)

PAT_COLORS = {
    "Pat_02": "#e74c3c", "Pat_03": "#3498db", "Pat_05": "#9b59b6",
    "Pat_07": "#2ecc71", "Pat_08": "#f39c12",
}
BAND_LS = {"delta": "-", "theta": "--", "alpha": "-.", "beta": ":",
           "low_gamma": (0, (3, 1, 1, 1)), "high_gamma": (0, (1, 1))}

LINES: list[str] = []


def R(s: str = "") -> None:
    print(s)
    LINES.append(s)


# ── Helpers ────────────────────────────────────────────────────────────

def get_nperseg(pat: str) -> int:
    fs = FS_MAP.get(pat, 2048.0)
    return int(fs * 2.0)


def load_ch(pat: str) -> list[str]:
    """Load channel labels, handling Pat_05's header quirk."""
    p = DATASET_ROOT / pat / "channel_labels.csv"
    if not p.exists():
        # try .txt
        p = DATASET_ROOT / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "")
        for l in df.iloc[:, 0]
    ]


def probe(label: str) -> str:
    """Extract probe letter(s) from contact label (e.g., 'A12' → 'A')."""
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def mixed_probes(ch: list[str], epi_mask: np.ndarray) -> dict:
    """Find probes with both epileptic and non-epileptic contacts."""
    P: dict[str, dict[str, list[int]]] = {}
    for i, l in enumerate(ch):
        pr = probe(l)
        if pr not in P:
            P[pr] = {"epi": [], "non": []}
        (P[pr]["epi"] if epi_mask[i] else P[pr]["non"]).append(i)
    return {p: v for p, v in P.items() if v["epi"] and v["non"]}


def load_imcoh(pat: str, phase: str, band: str, N: int) -> np.ndarray | None:
    """Load |ImCoh| (= fc_method='imcoh_abs') from cache.

    Post-reset (2026-04-15): the canonical signed ImCoh is cached
    freq-resolved; `load_fc_matrix("imcoh_abs")` applies `np.abs` at
    load and band-averages.
    """
    from lrg_eegfc.workflow.fc import load_fc_matrix as _load_fc_matrix
    A = _load_fc_matrix(pat, phase, band, "imcoh_abs")
    if A is None or A.shape[0] != N:
        return None
    A = A.copy()
    np.fill_diagonal(A, 0)
    return A


def load_lrg(pat: str, phase: str, band: str) -> LRGResult | None:
    """Load ImCoh LRG result."""
    fname = f"{band}_{phase}_lrg_imcoh-abs.npz"
    fpath = IMCOH_LRG_CACHE / pat / fname
    if not fpath.exists():
        return None
    d = np.load(fpath)
    return LRGResult(
        ultrametric_matrix=d["ultrametric_matrix"],
        linkage_matrix=d["linkage_matrix"],
        entropy_tau=d["entropy_tau"],
        entropy_1_minus_S=d["entropy_1_minus_S"],
        entropy_C=d["entropy_C"],
        optimal_threshold=float(d["optimal_threshold"]),
        patient=str(d["patient"]),
        phase=str(d["phase"]),
        band=str(d["band"]),
        fc_method=str(d["fc_method"]),
        n_nodes=int(d["n_nodes"]),
    )


def laplacian(A: np.ndarray):
    """Compute graph Laplacian eigendecomposition."""
    np.fill_diagonal(A, 0)
    D = A.sum(axis=1)
    L = np.diag(D) - A
    ev, U = np.linalg.eigh(L)
    ev = np.maximum(ev, 0.0)
    return ev, U, D


def communities_at_k(lrg: LRGResult, k: int) -> np.ndarray:
    """Cut the LRG dendrogram at k communities."""
    return fcluster(lrg.linkage_matrix, t=k, criterion="maxclust")


def ultrametric_from_linkage(lrg: LRGResult, n: int) -> np.ndarray:
    """Reconstruct full ultrametric distance matrix from linkage."""
    # The ultrametric_matrix is stored as condensed form
    return squareform(lrg.ultrametric_matrix[:n * (n - 1) // 2])


# ── Parse args ─────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true")
args = parser.parse_args()
VERBOSE = args.verbose

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 0 — Patient summary
# ═══════════════════════════════════════════════════════════════════════

R("=" * 90)
R("EPILEPTIC NODES IN THE ImCoh-LRG FRAMEWORK — DEEP INVESTIGATION")
R("=" * 90)
R()
R("ImCoh removes volume conduction → eliminates same-probe bias → enables")
R("clean community-level analysis for the first time.")
R()

patient_info = {}
for pat in PATIENTS:
    epi_all = load_epileptic_nodes(pat)
    ch = load_ch(pat)
    N = len(ch)
    epi_set = set(epi_all)
    epi_in_analysis = [l for l in ch if l in epi_set]
    epi_mask = np.array([l in epi_set for l in ch])
    mx = mixed_probes(ch, epi_mask)
    patient_info[pat] = {
        "ch": ch, "N": N, "epi_all": epi_all,
        "epi_in_analysis": epi_in_analysis,
        "epi_mask": epi_mask, "mixed_probes": mx,
    }
    R(f"  {pat}: {N} channels, {len(epi_in_analysis)}/{len(epi_all)} epileptic "
      f"in analysis, {len(mx)} mixed probes")

R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 1 — Q1: WITHIN-PROBE STRENGTH (ImCoh)
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 1: Within-probe ImCoh strength — epileptic vs non-epileptic")
R("-" * 90)
R()

weakness_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    ch, N, epi_mask, mx = info["ch"], info["N"], info["epi_mask"], info["mixed_probes"]
    if not mx:
        continue
    for phase in PHASES:
        for band in BANDS:
            A = load_imcoh(pat, phase, band, N)
            if A is None:
                continue
            for pname, pi in mx.items():
                se = A[pi["epi"]].sum(axis=1).mean()
                sn = A[pi["non"]].sum(axis=1).mean()
                weakness_rows.append(dict(
                    patient=pat, phase=phase, band=band, probe=pname,
                    s_epi=se, s_non=sn, ratio=se / (sn + 1e-30),
                ))

wdf = pd.DataFrame(weakness_rows)
if len(wdf) == 0:
    R("ERROR: No within-probe data collected. Check ImCoh cache paths.")
    sys.exit(1)
n_weaker = int((wdf["ratio"] < 1).sum())
p_sign = binomtest(n_weaker, len(wdf), 0.5).pvalue
R(f"Total within-probe comparisons: {len(wdf)}")
R(f"Epileptic WEAKER in {n_weaker}/{len(wdf)} ({n_weaker / len(wdf):.1%}), "
  f"sign-test p = {p_sign:.2e}")
R()

for pat in PATIENTS:
    sub = wdf[wdf["patient"] == pat]
    if len(sub) == 0:
        continue
    nw = int((sub["ratio"] < 1).sum())
    R(f"  {pat}: weaker in {nw}/{len(sub)} ({nw / len(sub):.0%}), "
      f"median ratio = {sub['ratio'].median():.3f}")

R()

# ── Figure 1: Within-probe strength ratio ─────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: by patient
ax = axes[0]
for i, pat in enumerate(PATIENTS):
    sub = wdf[wdf["patient"] == pat]
    if len(sub) == 0:
        continue
    jitter = np.random.uniform(-0.2, 0.2, len(sub))
    ax.scatter(i + jitter, sub["ratio"], s=12, alpha=0.35,
               c=PAT_COLORS[pat], edgecolors="none")
    ax.boxplot([sub["ratio"].values], positions=[i], widths=0.45,
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=PAT_COLORS[pat], alpha=0.25),
               medianprops=dict(color="k", lw=1.5))
ax.axhline(1, color="k", ls="--", lw=0.8)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels([f"{p}\n(n={len(wdf[wdf['patient']==p])})" for p in PATIENTS])
ax.set_ylabel("ImCoh strength ratio  (epileptic / non-epileptic)")
ax.set_title(f"Within-probe ImCoh strength\n"
             f"({n_weaker}/{len(wdf)} weaker, p = {p_sign:.1e})", fontsize=11)

# Panel B: by band
ax = axes[1]
band_positions = {b: i for i, b in enumerate(BANDS)}
for band in BANDS:
    sub = wdf[wdf["band"] == band]
    if len(sub) == 0:
        continue
    bp = band_positions[band]
    jitter = np.random.uniform(-0.2, 0.2, len(sub))
    colors = [PAT_COLORS[r["patient"]] for _, r in sub.iterrows()]
    ax.scatter(bp + jitter, sub["ratio"], s=12, alpha=0.35,
               c=colors, edgecolors="none")
    ax.boxplot([sub["ratio"].values], positions=[bp], widths=0.45,
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#3498db", alpha=0.2),
               medianprops=dict(color="k", lw=1.5))
ax.axhline(1, color="k", ls="--", lw=0.8)
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(BANDS, rotation=30, ha="right")
ax.set_ylabel("ImCoh strength ratio")
ax.set_title("By frequency band", fontsize=11)

fig.tight_layout()
fig.savefig(OUT / "fig1_imcoh_strength_ratio.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig1_imcoh_strength_ratio.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig1_imcoh_strength_ratio.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 2 — Q2: COMMUNITY MEMBERSHIP (ImCoh LRG)
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 2: Community membership — do epileptic nodes cluster together?")
R("-" * 90)
R()
R("With ImCoh, same-probe bias is removed. If epileptic nodes share")
R("pathological connectivity, they should co-cluster at some scale k.")
R()

comm_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    ch, N = info["ch"], info["N"]
    epi_mask = info["epi_mask"]
    epi_idx = np.where(epi_mask)[0]
    non_idx = np.where(~epi_mask)[0]
    n_epi = len(epi_idx)
    if n_epi < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            lrg = load_lrg(pat, phase, band)
            if lrg is None:
                continue

            for k in N_COMMS:
                if k >= N:
                    continue
                labels = communities_at_k(lrg, k)

                # Count intra-epileptic same-community pairs
                n_epi_pairs = n_epi * (n_epi - 1) // 2
                if n_epi_pairs == 0:
                    continue
                same_epi = 0
                for ii in range(n_epi):
                    for jj in range(ii + 1, n_epi):
                        if labels[epi_idx[ii]] == labels[epi_idx[jj]]:
                            same_epi += 1
                obs_frac = same_epi / n_epi_pairs

                # Expected fraction under random: for each community c,
                # P(both in c) = (n_c/N)^2, sum over c
                comm_sizes = np.bincount(labels)[1:]  # fcluster is 1-indexed
                expected_frac = np.sum((comm_sizes / N) ** 2)

                enrichment = obs_frac / (expected_frac + 1e-30)

                # Permutation test: how often does a random group of n_epi
                # nodes have ≥ same_epi same-community pairs?
                n_perm = 2000
                null_counts = np.zeros(n_perm)
                for perm_i in range(n_perm):
                    fake_epi = np.random.choice(N, n_epi, replace=False)
                    cnt = 0
                    for ii in range(n_epi):
                        for jj in range(ii + 1, n_epi):
                            if labels[fake_epi[ii]] == labels[fake_epi[jj]]:
                                cnt += 1
                    null_counts[perm_i] = cnt
                p_val = (np.sum(null_counts >= same_epi) + 1) / (n_perm + 1)

                comm_rows.append(dict(
                    patient=pat, phase=phase, band=band, k=k,
                    same_epi_pairs=same_epi, n_epi_pairs=n_epi_pairs,
                    obs_frac=obs_frac, expected_frac=expected_frac,
                    enrichment=enrichment, p_value=p_val,
                ))

cdf = pd.DataFrame(comm_rows)

# Summary: how many conditions show significant clustering?
R(f"Total (patient × phase × band × k) tests: {len(cdf)}")
n_sig = int((cdf["p_value"] < 0.05).sum())
n_sig_bonf = int((cdf["p_value"] < 0.05 / len(cdf)).sum())
R(f"Significant at p < 0.05: {n_sig}/{len(cdf)} ({n_sig / len(cdf):.1%})")
R(f"Significant after Bonferroni: {n_sig_bonf}/{len(cdf)}")
R()

# Per patient × k summary
R(f"{'patient':8s} {'k':>4s} {'mean_enrich':>12s} {'p<0.05':>7s} {'n_tests':>8s}")
R("-" * 45)
for pat in PATIENTS:
    for k in N_COMMS:
        sub = cdf[(cdf["patient"] == pat) & (cdf["k"] == k)]
        if len(sub) == 0:
            continue
        me = sub["enrichment"].mean()
        ns = int((sub["p_value"] < 0.05).sum())
        R(f"{pat:8s} {k:4d} {me:12.2f} {ns:7d} {len(sub):8d}")
R()

# ── Figure 2: Enrichment heatmap per patient ─────────────────────────

fig, axes = plt.subplots(1, len(PATIENTS), figsize=(3.5 * len(PATIENTS), 5),
                         sharey=True)
for ax, pat in zip(axes, PATIENTS):
    sub = cdf[cdf["patient"] == pat]
    if len(sub) == 0:
        ax.set_title(f"{pat}\n(no data)")
        continue
    # Average enrichment over phases and bands → (k × band) or (k × 1)
    pivot = sub.groupby(["k", "band"])["enrichment"].mean().reset_index()
    mat = pivot.pivot(index="k", columns="band", values="enrichment")
    mat = mat.reindex(index=N_COMMS, columns=BANDS)
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r",
                   vmin=0.3, vmax=3.0)
    ax.set_yticks(range(len(N_COMMS)))
    ax.set_yticklabels(N_COMMS)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=8)
    ax.set_title(f"{pat}", fontsize=11)
    if ax == axes[0]:
        ax.set_ylabel("n communities (k)")

    # Mark significant cells
    for ki, k in enumerate(N_COMMS):
        for bi, band in enumerate(BANDS):
            sub_cell = cdf[(cdf["patient"] == pat) & (cdf["k"] == k) &
                           (cdf["band"] == band)]
            n_sig_cell = int((sub_cell["p_value"] < 0.05).sum())
            if n_sig_cell > 0:
                ax.text(bi, ki, f"{n_sig_cell}", ha="center", va="center",
                        fontsize=7, fontweight="bold")

fig.colorbar(im, ax=axes, label="Enrichment (obs/expected)", shrink=0.6)
fig.suptitle("Epileptic node co-clustering in ImCoh LRG communities\n"
             "(numbers = phases with p < 0.05; color = mean enrichment)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 0.92, 0.92])
fig.savefig(OUT / "fig2_community_enrichment.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig2_community_enrichment.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig2_community_enrichment.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 3 — Q3: ULTRAMETRIC DISTANCE
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 3: Ultrametric distance — are epileptic pairs closer?")
R("-" * 90)
R()

um_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N = info["N"]
    epi_idx = np.where(info["epi_mask"])[0]
    non_idx = np.where(~info["epi_mask"])[0]
    n_epi = len(epi_idx)
    if n_epi < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            lrg = load_lrg(pat, phase, band)
            if lrg is None:
                continue
            # Reconstruct ultrametric distance matrix
            um = squareform(lrg.ultrametric_matrix)

            # Intra-epileptic distances
            epi_dists = []
            for ii in range(n_epi):
                for jj in range(ii + 1, n_epi):
                    epi_dists.append(um[epi_idx[ii], epi_idx[jj]])
            epi_dists = np.array(epi_dists)

            # Intra-non-epileptic distances (sample same size for fairness)
            all_dists = um[np.triu_indices(N, k=1)]
            n_non = len(non_idx)
            non_dists = []
            for ii in range(n_non):
                for jj in range(ii + 1, n_non):
                    non_dists.append(um[non_idx[ii], non_idx[jj]])
            non_dists = np.array(non_dists)

            # Cross distances (epi-to-non)
            cross_dists = []
            for ii in epi_idx:
                for jj in non_idx:
                    cross_dists.append(um[ii, jj])
            cross_dists = np.array(cross_dists)

            # Mann-Whitney: intra-epi vs all pairwise
            if len(epi_dists) > 1 and len(all_dists) > 1:
                _, p_vs_all = mannwhitneyu(epi_dists, all_dists, alternative="less")
                _, p_vs_non = mannwhitneyu(epi_dists, non_dists, alternative="less")
                _, p_vs_cross = mannwhitneyu(epi_dists, cross_dists, alternative="less")
            else:
                p_vs_all = p_vs_non = p_vs_cross = 1.0

            um_rows.append(dict(
                patient=pat, phase=phase, band=band,
                mean_epi=epi_dists.mean(), mean_non=non_dists.mean(),
                mean_cross=cross_dists.mean(), mean_all=all_dists.mean(),
                ratio_vs_all=epi_dists.mean() / (all_dists.mean() + 1e-30),
                ratio_vs_non=epi_dists.mean() / (non_dists.mean() + 1e-30),
                p_vs_all=p_vs_all, p_vs_non=p_vs_non, p_vs_cross=p_vs_cross,
            ))

umdf = pd.DataFrame(um_rows)

n_sig_all = int((umdf["p_vs_all"] < 0.05).sum())
n_sig_non = int((umdf["p_vs_non"] < 0.05).sum())
R(f"Tests: {len(umdf)}")
R(f"Intra-epi CLOSER than all pairs (p < 0.05): {n_sig_all}/{len(umdf)}")
R(f"Intra-epi CLOSER than intra-non-epi (p < 0.05): {n_sig_non}/{len(umdf)}")
R()

for pat in PATIENTS:
    sub = umdf[umdf["patient"] == pat]
    if len(sub) == 0:
        continue
    mr = sub["ratio_vs_non"].mean()
    ns = int((sub["p_vs_non"] < 0.05).sum())
    R(f"  {pat}: mean epi/non distance ratio = {mr:.3f}, "
      f"sig in {ns}/{len(sub)} conditions")
R()

# ── Figure 3: Ultrametric distance comparison ─────────────────────────

fig, axes = plt.subplots(1, len(PATIENTS), figsize=(3.5 * len(PATIENTS), 5),
                         sharey=True)
for ax, pat in zip(axes, PATIENTS):
    sub = umdf[umdf["patient"] == pat]
    if len(sub) == 0:
        ax.set_title(f"{pat}\n(no data)")
        continue
    # Show distribution of ratio_vs_non across conditions
    for band in BANDS:
        sub_b = sub[sub["band"] == band]
        if len(sub_b) == 0:
            continue
        jitter = np.random.uniform(-0.15, 0.15, len(sub_b))
        bi = BANDS.index(band)
        ax.scatter(bi + jitter, sub_b["ratio_vs_non"], s=20, alpha=0.5,
                   color=PAT_COLORS[pat], edgecolors="none")
    ax.boxplot(
        [sub[sub["band"] == b]["ratio_vs_non"].values for b in BANDS
         if len(sub[sub["band"] == b]) > 0],
        positions=[i for i, b in enumerate(BANDS)
                   if len(sub[sub["band"] == b]) > 0],
        widths=0.5, showfliers=False, patch_artist=True,
        boxprops=dict(facecolor=PAT_COLORS[pat], alpha=0.25),
        medianprops=dict(color="k", lw=1.5))
    ax.axhline(1, color="k", ls="--", lw=0.8)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=8)
    ax.set_title(f"{pat}", fontsize=11)
    if ax == axes[0]:
        ax.set_ylabel("Ultrametric distance ratio\n(epi-epi / non-non)")

fig.suptitle("Epileptic nodes in ImCoh ultrametric space\n"
             "(< 1 = epileptic pairs closer than non-epileptic pairs)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT / "fig3_ultrametric_distance.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig3_ultrametric_distance.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig3_ultrametric_distance.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 4 — Q4: SCALE AMPLIFICATION (ImCoh heat kernel)
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 4: Scale amplification — ImCoh heat kernel")
R("-" * 90)
R()

sa_curves = {}
sa_stats = []

for pat in PATIENTS:
    info = patient_info[pat]
    ch, N, epi_mask, mx = info["ch"], info["N"], info["epi_mask"], info["mixed_probes"]
    if not mx:
        continue
    epi_pool = sum((v["epi"] for v in mx.values()), [])
    non_pool = sum((v["non"] for v in mx.values()), [])

    for phase in ["rsPre", "rsPost"]:
        for band in BANDS:
            A = load_imcoh(pat, phase, band, N)
            if A is None:
                continue
            ev, U, _ = laplacian(A)
            if ev[1] < 1e-10:
                continue

            n_tau = 100
            tau_g = np.logspace(np.log10(1 / ev[-1]), np.log10(1 / ev[1]), n_tau)
            sr_e = np.zeros(n_tau)
            sr_n = np.zeros(n_tau)
            for ti, tau in enumerate(tau_g):
                K = U @ np.diag(np.exp(-tau * ev)) @ U.T
                sr = np.diag(K)
                sr_e[ti] = sr[epi_pool].mean()
                sr_n[ti] = sr[non_pool].mean()

            mean_sr = (sr_e + sr_n) / 2
            rel = (sr_e - sr_n) / (mean_sr + 1e-30)
            lt = np.log10(tau_g)
            rho, p_rho = spearmanr(lt, rel)
            sa_curves[(pat, phase, band)] = (lt, rel)
            sa_stats.append(dict(
                patient=pat, phase=phase, band=band,
                rho=rho, p=p_rho,
                fine=rel[:10].mean(), coarse=rel[-10:].mean(),
            ))
            if VERBOSE:
                sig = "***" if p_rho < .001 else "** " if p_rho < .01 else "*  " if p_rho < .05 else "   "
                R(f"  {pat} {phase:10s} {band:12s}: ρ = {rho:+.2f} (p={p_rho:.1e}) "
                  f"fine={rel[:10].mean():+.1%} coarse={rel[-10:].mean():+.1%} {sig}")

sadf = pd.DataFrame(sa_stats)
if len(sadf) > 0:
    n_pos = int(((sadf["rho"] > 0.5) & (sadf["p"] < 0.01)).sum())
    n_neg = int(((sadf["rho"] < -0.5) & (sadf["p"] < 0.01)).sum())
    R(f"Scale amplification tests: {len(sadf)}")
    R(f"  ρ > 0.5 & p < 0.01 (weakness grows with scale): {n_pos}/{len(sadf)}")
    R(f"  ρ < -0.5 & p < 0.01 (weakness shrinks): {n_neg}/{len(sadf)}")
    R()
    for pat in PATIENTS:
        sub = sadf[sadf["patient"] == pat]
        if len(sub) == 0:
            continue
        R(f"  {pat}: mean ρ = {sub['rho'].mean():+.2f}, "
          f"fine = {sub['fine'].mean():+.1%}, coarse = {sub['coarse'].mean():+.1%}")
    R()

# ── Figure 4: Scale amplification ────────────────────────────────────

if len(sa_curves) > 0:
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Main panel: overlay all rsPre curves
    ax_main = fig.add_subplot(gs[0, :])
    for (pat, phase, band), (lt, rel) in sa_curves.items():
        if phase != "rsPre":
            continue
        ax_main.plot(lt, rel * 100, ls=BAND_LS.get(band, "-"),
                     color=PAT_COLORS[pat], alpha=0.6, lw=1.2,
                     label=f"{pat} {band}")
    ax_main.axhline(0, color="k", ls="-", lw=0.5)
    ax_main.set_xlabel("log₁₀(τ)  [diffusion time]")
    ax_main.set_ylabel("Relative self-return difference (%)\n(epi − non-epi) / mean")
    ax_main.set_title("ImCoh scale amplification (rsPre)\n"
                       "Positive = epi retains more heat", fontsize=11)
    ax_main.legend(fontsize=6, ncol=5, loc="upper left")

    # Per-patient panels
    for pi, pat in enumerate(PATIENTS[:3]):
        ax = fig.add_subplot(gs[1, pi])
        for band in BANDS:
            key = (pat, "rsPre", band)
            if key not in sa_curves:
                continue
            lt, rel = sa_curves[key]
            ax.plot(lt, rel * 100, ls=BAND_LS.get(band, "-"),
                    color=PAT_COLORS[pat], lw=1.5, label=band)
        ax.axhline(0, color="k", ls="-", lw=0.5)
        ax.set_xlabel("log₁₀(τ)")
        ax.set_ylabel("Δ (%)")
        sub = sadf[(sadf["patient"] == pat) & (sadf["phase"] == "rsPre")]
        if len(sub) > 0:
            ax.set_title(f"{pat}  (mean ρ = {sub['rho'].mean():+.2f})", fontsize=10)
        else:
            ax.set_title(pat, fontsize=10)
        ax.legend(fontsize=7)

    fig.savefig(OUT / "fig4_scale_amplification.pdf", bbox_inches="tight", dpi=200)
    fig.savefig(OUT / "fig4_scale_amplification.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    R(f"  → {OUT / 'fig4_scale_amplification.pdf'}")
    R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 5 — Q5: METASTABILITY
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 5: Metastability — do epileptic nodes change community more often?")
R("-" * 90)
R()

meta_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N = info["N"]
    epi_idx = np.where(info["epi_mask"])[0]
    non_idx = np.where(~info["epi_mask"])[0]
    n_epi = len(epi_idx)
    if n_epi < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            lrg = load_lrg(pat, phase, band)
            if lrg is None:
                continue

            # Compute metastability: at each k, which nodes change community
            # relative to k-1?
            ks = sorted([k for k in N_COMMS if k < N], reverse=True)
            if len(ks) < 2:
                continue

            prev_labels = communities_at_k(lrg, ks[0])
            changes = np.zeros(N)
            n_transitions = 0
            for ki in range(1, len(ks)):
                curr_labels = communities_at_k(lrg, ks[ki])
                # A node "changes" if its community partners change
                # (relabeling doesn't matter — check co-membership)
                for i in range(N):
                    # Check if node i kept same neighbours
                    prev_mates = set(np.where(prev_labels == prev_labels[i])[0])
                    curr_mates = set(np.where(curr_labels == curr_labels[i])[0])
                    if prev_mates != curr_mates:
                        changes[i] += 1
                n_transitions += 1
                prev_labels = curr_labels

            mu = changes / (n_transitions + 1e-30)  # metastability index

            mu_epi = mu[epi_idx]
            mu_non = mu[non_idx]

            if len(mu_epi) > 0 and len(mu_non) > 0:
                _, p_gt = mannwhitneyu(mu_epi, mu_non, alternative="greater")
                _, p_lt = mannwhitneyu(mu_epi, mu_non, alternative="less")
            else:
                p_gt = p_lt = 1.0

            meta_rows.append(dict(
                patient=pat, phase=phase, band=band,
                mu_epi=mu_epi.mean(), mu_non=mu_non.mean(),
                ratio=mu_epi.mean() / (mu_non.mean() + 1e-30),
                p_more_meta=p_gt, p_less_meta=p_lt,
                direction="more" if mu_epi.mean() > mu_non.mean() else "less",
            ))

metadf = pd.DataFrame(meta_rows)
if len(metadf) > 0:
    n_more = int((metadf["direction"] == "more").sum())
    n_sig_more = int((metadf["p_more_meta"] < 0.05).sum())
    n_sig_less = int((metadf["p_less_meta"] < 0.05).sum())
    R(f"Metastability tests: {len(metadf)}")
    R(f"Epileptic MORE metastable: {n_more}/{len(metadf)} "
      f"({n_more / len(metadf):.0%})")
    R(f"Significantly more metastable (p < 0.05): {n_sig_more}")
    R(f"Significantly less metastable (p < 0.05): {n_sig_less}")
    R()

    for pat in PATIENTS:
        sub = metadf[metadf["patient"] == pat]
        if len(sub) == 0:
            continue
        nm = int((sub["direction"] == "more").sum())
        R(f"  {pat}: more metastable in {nm}/{len(sub)}, "
          f"mean μ_epi/μ_non = {sub['ratio'].mean():.3f}")
    R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 6 — CROSS-PATIENT CONSENSUS HEATMAPS
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 6: Cross-patient consensus — phase × band heatmaps")
R("-" * 90)
R()

# ── 6A: Strength weakness consensus ──────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# Panel 1: Per patient, fraction of conditions where epi is weaker
for ai, pat in enumerate(PATIENTS):
    ax = axes.flat[ai]
    sub = wdf[wdf["patient"] == pat]
    if len(sub) == 0:
        ax.set_title(f"{pat} (no data)")
        continue
    # phase × band heatmap of median ratio
    mat = np.full((len(PHASES), len(BANDS)), np.nan)
    for pi, phase in enumerate(PHASES):
        for bi, band in enumerate(BANDS):
            cell = sub[(sub["phase"] == phase) & (sub["band"] == band)]
            if len(cell) > 0:
                mat[pi, bi] = cell["ratio"].median()
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.5, vmax=1.5)
    ax.set_yticks(range(len(PHASES)))
    ax.set_yticklabels(PHASES, fontsize=8)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=7)
    ax.set_title(f"{pat}", fontsize=10)
    # Annotate with values
    for pi in range(len(PHASES)):
        for bi in range(len(BANDS)):
            if not np.isnan(mat[pi, bi]):
                ax.text(bi, pi, f"{mat[pi, bi]:.2f}", ha="center", va="center",
                        fontsize=6, color="k" if 0.7 < mat[pi, bi] < 1.3 else "w")

# Panel 6: consensus (n patients with ratio < 1)
ax = axes.flat[5]
consensus = np.zeros((len(PHASES), len(BANDS)))
for pi, phase in enumerate(PHASES):
    for bi, band in enumerate(BANDS):
        n_weaker_cell = 0
        n_tested = 0
        for pat in PATIENTS:
            sub = wdf[(wdf["patient"] == pat) & (wdf["phase"] == phase) &
                      (wdf["band"] == band)]
            if len(sub) > 0:
                n_tested += 1
                if sub["ratio"].median() < 1:
                    n_weaker_cell += 1
        consensus[pi, bi] = n_weaker_cell / max(n_tested, 1)

im = ax.imshow(consensus, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=1)
ax.set_yticks(range(len(PHASES)))
ax.set_yticklabels(PHASES, fontsize=8)
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=7)
ax.set_title("Consensus: fraction\npatients epi weaker", fontsize=10)
for pi in range(len(PHASES)):
    for bi in range(len(BANDS)):
        ax.text(bi, pi, f"{consensus[pi, bi]:.0%}", ha="center", va="center",
                fontsize=8, fontweight="bold")

fig.suptitle("ImCoh within-probe strength ratio (epileptic / non-epileptic)\n"
             "Blue < 1 = epileptic weaker; Red > 1 = epileptic stronger",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT / "fig5_consensus_strength.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig5_consensus_strength.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig5_consensus_strength.pdf'}")

# ── 6B: Community enrichment consensus ────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ki, k in enumerate([5, 10, 20]):  # three representative scales
    for row, phase_set in enumerate([["rsPre", "rsPost"],
                                      ["taskLearn", "taskTest"]]):
        ax = axes[row, ki]
        # For each band, count patients with enrichment > 1.5
        mat = np.zeros((len(PATIENTS), len(BANDS)))
        for pi, pat in enumerate(PATIENTS):
            for bi, band in enumerate(BANDS):
                sub = cdf[(cdf["patient"] == pat) & (cdf["k"] == k) &
                          (cdf["band"] == band) & (cdf["phase"].isin(phase_set))]
                if len(sub) > 0:
                    mat[pi, bi] = sub["enrichment"].mean()

        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.3, vmax=3.0)
        ax.set_yticks(range(len(PATIENTS)))
        ax.set_yticklabels(PATIENTS, fontsize=8)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=7)
        phase_label = "rest" if row == 0 else "task"
        ax.set_title(f"k={k}, {phase_label}", fontsize=10)
        for pi in range(len(PATIENTS)):
            for bi in range(len(BANDS)):
                if mat[pi, bi] > 0:
                    ax.text(bi, pi, f"{mat[pi, bi]:.1f}", ha="center",
                            va="center", fontsize=6)

fig.suptitle("Community co-clustering enrichment (obs / expected)\n"
             "> 1 = epileptic nodes cluster more than random",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT / "fig6_consensus_community.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig6_consensus_community.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig6_consensus_community.pdf'}")

# ── 6C: Ultrametric distance consensus ────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: ratio_vs_non by patient × band (averaged over phases)
ax = axes[0]
mat = np.full((len(PATIENTS), len(BANDS)), np.nan)
for pi, pat in enumerate(PATIENTS):
    for bi, band in enumerate(BANDS):
        sub = umdf[(umdf["patient"] == pat) & (umdf["band"] == band)]
        if len(sub) > 0:
            mat[pi, bi] = sub["ratio_vs_non"].mean()
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.6, vmax=1.4)
ax.set_yticks(range(len(PATIENTS)))
ax.set_yticklabels(PATIENTS)
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=8)
ax.set_title("Ultrametric distance ratio (epi/non)\n< 1 = epi pairs closer",
             fontsize=10)
for pi in range(len(PATIENTS)):
    for bi in range(len(BANDS)):
        if not np.isnan(mat[pi, bi]):
            ax.text(bi, pi, f"{mat[pi, bi]:.2f}", ha="center", va="center",
                    fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.8)

# Panel B: significance count by patient
ax = axes[1]
for pi, pat in enumerate(PATIENTS):
    sub = umdf[umdf["patient"] == pat]
    if len(sub) == 0:
        continue
    n_sig_p = int((sub["p_vs_non"] < 0.05).sum())
    ax.bar(pi, n_sig_p / max(len(sub), 1) * 100,
           color=PAT_COLORS[pat], alpha=0.7)
    ax.text(pi, n_sig_p / max(len(sub), 1) * 100 + 1,
            f"{n_sig_p}/{len(sub)}", ha="center", fontsize=9)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels(PATIENTS)
ax.set_ylabel("% conditions with p < 0.05")
ax.set_title("Fraction of conditions where epileptic\npairs are significantly closer",
             fontsize=10)
ax.set_ylim(0, 100)

fig.tight_layout()
fig.savefig(OUT / "fig7_consensus_ultrametric.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig7_consensus_ultrametric.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig7_consensus_ultrametric.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 7 — PROBE-BIAS CHECK
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 7: Probe-bias sanity check — ImCoh same-probe enrichment")
R("-" * 90)
R()
R("Verifying that ImCoh indeed reduces same-probe community enrichment.")
R()

probe_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    ch, N = info["ch"], info["N"]
    probes = [probe(l) for l in ch]
    probe_labels = np.array(probes)

    for phase in ["rsPre"]:
        for band in ["alpha", "beta"]:
            lrg = load_lrg(pat, phase, band)
            if lrg is None:
                continue
            for k in N_COMMS:
                if k >= N:
                    continue
                comm = communities_at_k(lrg, k)
                # Fraction of same-community pairs that are also same-probe
                same_comm_same_probe = 0
                same_comm_total = 0
                for i in range(N):
                    for j in range(i + 1, N):
                        if comm[i] == comm[j]:
                            same_comm_total += 1
                            if probe_labels[i] == probe_labels[j]:
                                same_comm_same_probe += 1
                obs_frac = same_comm_same_probe / max(same_comm_total, 1)

                # Expected: fraction of all pairs that are same-probe
                total_pairs = N * (N - 1) // 2
                same_probe_pairs = 0
                for p_name in set(probes):
                    n_p = probes.count(p_name)
                    same_probe_pairs += n_p * (n_p - 1) // 2
                exp_frac = same_probe_pairs / total_pairs

                probe_rows.append(dict(
                    patient=pat, band=band, k=k,
                    same_probe_in_comm=obs_frac,
                    expected=exp_frac,
                    enrichment=obs_frac / (exp_frac + 1e-30),
                ))

pdf = pd.DataFrame(probe_rows)
R(f"{'patient':8s} {'band':6s} {'k':>4s} {'same_probe%':>12s} {'expected%':>10s} {'enrichment':>11s}")
R("-" * 55)
for _, row in pdf.iterrows():
    flag = " ⚠" if row["enrichment"] > 2.0 else ""
    R(f"{row['patient']:8s} {row['band']:6s} {row['k']:4.0f} "
      f"{row['same_probe_in_comm']:12.1%} {row['expected']:10.1%} "
      f"{row['enrichment']:11.2f}{flag}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 8 — LAPLACIAN MARKERS (ImCoh)
# ═══════════════════════════════════════════════════════════════════════

R("SECTION 8: Laplacian spectral markers on ImCoh adjacency")
R("-" * 90)
R()

marker_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    ch, N, epi_mask, mx = info["ch"], info["N"], info["epi_mask"], info["mixed_probes"]
    if not mx:
        continue
    epi_pool = sum((v["epi"] for v in mx.values()), [])
    non_pool = sum((v["non"] for v in mx.values()), [])
    if len(epi_pool) < 2 or len(non_pool) < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            A = load_imcoh(pat, phase, band, N)
            if A is None:
                continue
            ev, U, D = laplacian(A)
            if ev[1] < 1e-10:
                continue

            # Marker 1: Node strength
            strength = D.copy()

            # Marker 2: L† diagonal (effective resistance to ground)
            lam_inv = np.zeros_like(ev)
            lam_inv[ev > 1e-10] = 1.0 / ev[ev > 1e-10]
            L_pinv_diag = np.sum(U ** 2 * lam_inv, axis=1)

            # Marker 3: Spectral entropy
            n_ev = min(30, N - 1)
            V2 = U[:, 1:n_ev + 1] ** 2
            V2n = V2 / (V2.sum(axis=1, keepdims=True) + 1e-30)
            spec_entropy = -np.sum(V2n * np.log(V2n + 1e-30), axis=1)

            # Marker 4: Return time
            return_time = np.sum(U[:, 1:] ** 2 / ev[1:], axis=1)

            # Marker 5: Fiedler norm
            D_inv_sqrt = np.diag(1.0 / np.sqrt(D + 1e-30))
            L_norm = D_inv_sqrt @ (np.diag(D) - A) @ D_inv_sqrt
            ev_n, U_n = np.linalg.eigh(L_norm)
            fiedler_norm = np.abs(U_n[:, 1])

            # Marker 6: Local spectral gap
            v2_sq = U[:, 1] ** 2
            v3_sq = U[:, 2] ** 2
            local_gap = v2_sq / (v3_sq + 1e-30)

            metrics = {
                "strength": strength,
                "L†_diag": L_pinv_diag,
                "spectral_entropy": spec_entropy,
                "return_time": return_time,
                "fiedler_norm": fiedler_norm,
                "local_gap": local_gap,
            }

            for mname, mvals in metrics.items():
                ve = mvals[epi_pool]
                vn = mvals[non_pool]
                _, p_gt = mannwhitneyu(ve, vn, alternative="greater")
                _, p_lt = mannwhitneyu(ve, vn, alternative="less")
                p_best = min(p_gt, p_lt)
                direction = "epi>" if p_gt < p_lt else "epi<"
                try:
                    from sklearn.metrics import roc_auc_score
                    y = np.zeros(N)
                    y[epi_pool] = 1
                    pool_idx = np.array(epi_pool + non_pool)
                    auc = roc_auc_score(y[pool_idx], mvals[pool_idx])
                    auc = max(auc, 1 - auc)
                except Exception:
                    auc = 0.5

                marker_rows.append(dict(
                    patient=pat, phase=phase, band=band,
                    metric=mname, direction=direction,
                    p=p_best, auc=auc,
                ))

mdf = pd.DataFrame(marker_rows)

R(f"{'metric':20s} {'epi_lower':>9s} {'epi_higher':>10s} {'p<0.05':>7s} "
  f"{'mean_AUC':>9s} {'dominant':>10s} {'consistency':>12s}")
R("-" * 85)
for mname in ["strength", "L†_diag", "spectral_entropy",
              "return_time", "fiedler_norm", "local_gap"]:
    sub = mdf[mdf["metric"] == mname]
    if len(sub) == 0:
        continue
    n_lower = int((sub["direction"] == "epi<").sum())
    n_higher = int((sub["direction"] == "epi>").sum())
    n_sig = int((sub["p"] < 0.05).sum())
    mean_auc = sub["auc"].mean()
    dominant = "epi<" if n_lower > n_higher else "epi>"
    consistency = max(n_lower, n_higher) / len(sub)
    R(f"{mname:20s} {n_lower:9d} {n_higher:10d} {n_sig:7d} {mean_auc:9.3f} "
      f"{dominant:>10s} {consistency:>12.0%}")
R()

# Per-patient consistency
R("Per-patient direction consistency (strength metric):")
for pat in PATIENTS:
    sub = mdf[(mdf["metric"] == "strength") & (mdf["patient"] == pat)]
    if len(sub) == 0:
        continue
    n_lower = int((sub["direction"] == "epi<").sum())
    R(f"  {pat}: epi weaker in {n_lower}/{len(sub)} ({n_lower / len(sub):.0%})")
R()

# ── Figure 8: Laplacian marker comparison ─────────────────────────────

metrics_order = ["strength", "L†_diag", "return_time",
                 "spectral_entropy", "fiedler_norm", "local_gap"]
metric_labels = ["Node\nstrength", "L† diagonal\n(eff. resist.)",
                 "Return\ntime", "Spectral\nentropy",
                 "Fiedler\n|v₂|", "Local\ngap v₂/v₃"]

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel A: AUC boxplot per metric
ax = axes[0]
auc_data = [mdf[mdf["metric"] == m]["auc"].values for m in metrics_order]
bp = ax.boxplot(auc_data, labels=metric_labels, widths=0.6,
                patch_artist=True, showfliers=False)
for patch in bp["boxes"]:
    patch.set_facecolor("#3498db")
    patch.set_alpha(0.4)
ax.axhline(0.5, color="gray", ls="--", lw=0.8, label="chance")
ax.set_ylabel("ROC AUC (within-probe, epi vs non-epi)")
ax.set_title("Discriminative power of ImCoh Laplacian metrics", fontsize=11)
ax.legend()
ax.tick_params(axis="x", labelsize=8)

# Panel B: Consistency per patient
ax = axes[1]
x = np.arange(len(metrics_order))
width = 0.15
for pi, pat in enumerate(PATIENTS):
    vals = []
    for mname in metrics_order:
        sub = mdf[(mdf["metric"] == mname) & (mdf["patient"] == pat)]
        if len(sub) > 0:
            dominant = sub["direction"].mode().iloc[0]
            c = (sub["direction"] == dominant).sum() / len(sub)
        else:
            c = 0
        vals.append(c * 100)
    ax.bar(x + pi * width, vals, width, label=pat,
           color=PAT_COLORS[pat], alpha=0.7)
ax.axhline(50, color="gray", ls="--", lw=0.8)
ax.set_xticks(x + 2 * width)
ax.set_xticklabels(metric_labels, fontsize=8)
ax.set_ylabel("Direction consistency (%)")
ax.set_title("Per-patient consistency", fontsize=11)
ax.legend(fontsize=8)
ax.set_ylim(0, 105)

fig.tight_layout()
fig.savefig(OUT / "fig8_laplacian_metrics.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "fig8_laplacian_metrics.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'fig8_laplacian_metrics.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  SECTION 9 — GRAND SUMMARY
# ═══════════════════════════════════════════════════════════════════════

R("=" * 90)
R("GRAND SUMMARY")
R("=" * 90)
R()

R("Q1. Within-probe ImCoh strength:")
if len(wdf) > 0:
    R(f"     Epileptic weaker in {n_weaker}/{len(wdf)} "
      f"({n_weaker / len(wdf):.0%}), p = {p_sign:.2e}")
    for pat in PATIENTS:
        sub = wdf[wdf["patient"] == pat]
        if len(sub) > 0:
            nw = int((sub["ratio"] < 1).sum())
            R(f"     {pat}: {nw}/{len(sub)} weaker ({nw / len(sub):.0%})")
R()

R("Q2. Community co-clustering:")
if len(cdf) > 0:
    n_sig_q2 = int((cdf["p_value"] < 0.05).sum())
    R(f"     {n_sig_q2}/{len(cdf)} conditions show significant clustering "
      f"(p < 0.05)")
    for pat in PATIENTS:
        sub = cdf[cdf["patient"] == pat]
        if len(sub) > 0:
            ns = int((sub["p_value"] < 0.05).sum())
            me = sub["enrichment"].mean()
            R(f"     {pat}: {ns}/{len(sub)} sig, mean enrichment = {me:.2f}")
R()

R("Q3. Ultrametric distance:")
if len(umdf) > 0:
    for pat in PATIENTS:
        sub = umdf[umdf["patient"] == pat]
        if len(sub) > 0:
            ns = int((sub["p_vs_non"] < 0.05).sum())
            mr = sub["ratio_vs_non"].mean()
            R(f"     {pat}: {ns}/{len(sub)} sig closer, mean ratio = {mr:.3f}")
R()

R("Q4. Scale amplification:")
if len(sadf) > 0:
    for pat in PATIENTS:
        sub = sadf[sadf["patient"] == pat]
        if len(sub) > 0:
            R(f"     {pat}: mean ρ = {sub['rho'].mean():+.2f}")
R()

R("Q5. Metastability:")
if len(metadf) > 0:
    for pat in PATIENTS:
        sub = metadf[metadf["patient"] == pat]
        if len(sub) > 0:
            nm = int((sub["direction"] == "more").sum())
            R(f"     {pat}: more metastable in {nm}/{len(sub)}")
R()

R("Q6. Cross-patient consensus:")
R("     See heatmap figures 5-7.")
R()

# ── Save report ───────────────────────────────────────────────────────

report_path = OUT / "report.txt"
report_path.write_text("\n".join(LINES))
R(f"Report saved: {report_path}")
R(f"Figures saved: {OUT}")
