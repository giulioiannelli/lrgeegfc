#!/usr/bin/env python3
"""Deep LRG & network-topology analysis of epileptic nodes under ImCoh.

Extends the initial investigation (epileptic_imcoh_investigation.py) with:

  LRG-specific:
    A1. Probe-matched permutation test on ultrametric closeness (spatial control)
    A2. Merge-height profiles (epi–epi vs epi–non vs non–non)
    A3. Ultrametric nearest-neighbor enrichment
    A4. Dendrogram visualisation with epileptic nodes coloured
    A5. Band-specific breakdown of Q3 (strongest effect finder)

  Standard network topology:
    B1. Weighted node strength
    B2. Weighted clustering coefficient
    B3. Betweenness centrality
    B4. Eigenvector centrality
    B5. Participation coefficient (community diversity)
    B6. Within-module degree z-score
    B7. Intra-epileptic vs cross-group edge weight analysis

  Publication figures (all panels, per-patient + consensus).

Run:  python scripts/08_epileptic/epileptic_imcoh_deep_lrg.py [-v]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import binomtest, mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score

# ── Project setup ─────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import (
    FIGURES_ROOT, IMCOH_CACHE, IMCOH_LRG_CACHE, SEEG_DATAPATH,
)
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

# ── Config ─────────────────────────────────────────────────────────────

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)
N_PERMS = 3000  # probe-matched permutation count

OUT = FIGURES_ROOT / "epileptic_imcoh_deep"
OUT.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical
np.random.seed(42)

PAT_COLORS = {
    "Pat_02": "#e74c3c", "Pat_03": "#3498db", "Pat_05": "#9b59b6",
    "Pat_07": "#2ecc71", "Pat_08": "#f39c12",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}

LINES: list[str] = []


def R(s: str = "") -> None:
    print(s)
    LINES.append(s)


# ── Helpers ────────────────────────────────────────────────────────────

def load_ch(pat: str) -> list[str]:
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not p.exists():
        p = SEEG_DATAPATH / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "")
        for l in df.iloc[:, 0]
    ]


def probe(label: str) -> str:
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def probe_indices(ch: list[str]) -> dict[str, list[int]]:
    """Map probe name → list of channel indices."""
    P: dict[str, list[int]] = {}
    for i, l in enumerate(ch):
        pr = probe(l)
        P.setdefault(pr, []).append(i)
    return P


def load_A(pat, phase, band, N):
    """Load |ImCoh| adjacency (zero diagonal).

    Uses ``fc_method='imcoh_abs'`` so the loader returns ``np.abs`` of the
    signed Nolte-2004 cache — no manual transform needed here.
    """
    A = load_fc_matrix(pat, phase, band, "imcoh_abs")
    if A is None or A.shape[0] != N:
        return None
    A = A.copy()
    np.fill_diagonal(A, 0)
    return A


def laplacian(A):
    np.fill_diagonal(A, 0)
    D = A.sum(axis=1)
    L = np.diag(D) - A
    ev, U = np.linalg.eigh(L)
    ev = np.maximum(ev, 0.0)
    return ev, U, D


# ── Network metrics ───────────────────────────────────────────────────

def weighted_clustering_coefficient(A):
    """Onnela et al. (2005) weighted clustering coefficient."""
    N = A.shape[0]
    A_third = np.cbrt(A)
    num = A_third @ A_third @ A_third
    diag_num = np.diag(num)
    k = (A > 0).sum(axis=1)
    denom = k * (k - 1)
    cc = np.zeros(N)
    mask = denom > 0
    cc[mask] = diag_num[mask] / denom[mask]
    return cc


def betweenness_centrality(A):
    """Approximate betweenness via networkx."""
    import networkx as nx
    G = nx.from_numpy_array(A)
    bc = nx.betweenness_centrality(G, weight="weight")
    return np.array([bc[i] for i in range(A.shape[0])])


def eigenvector_centrality(A):
    """Leading eigenvector of adjacency."""
    evals, evecs = np.linalg.eigh(A)
    # Largest eigenvalue is last (sorted ascending)
    ec = np.abs(evecs[:, -1])
    ec /= ec.max() + 1e-30
    return ec


def participation_coefficient(A, communities):
    """Guimerà & Amaral (2005) participation coefficient."""
    N = A.shape[0]
    ki = A.sum(axis=1)
    P = np.zeros(N)
    for c in np.unique(communities):
        mask = communities == c
        kis = A[:, mask].sum(axis=1)
        P += (kis / (ki + 1e-30)) ** 2
    return 1.0 - P


def within_module_degree_zscore(A, communities):
    """Guimerà & Amaral (2005) within-module degree z-score."""
    N = A.shape[0]
    z = np.zeros(N)
    for c in np.unique(communities):
        mask = communities == c
        idx = np.where(mask)[0]
        if len(idx) < 2:
            continue
        ki_c = A[np.ix_(idx, idx)].sum(axis=1)
        mu = ki_c.mean()
        sigma = ki_c.std()
        if sigma > 1e-10:
            z[idx] = (ki_c - mu) / sigma
    return z


# ── Parse args ─────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true")
args = parser.parse_args()
VERBOSE = args.verbose

# ── Load patient info ─────────────────────────────────────────────────

patient_info = {}
for pat in PATIENTS:
    epi_all = load_epileptic_nodes(pat)
    ch = load_ch(pat)
    N = len(ch)
    epi_set = set(epi_all)
    epi_in = [l for l in ch if l in epi_set]
    epi_mask = np.array([l in epi_set for l in ch])
    epi_idx = np.where(epi_mask)[0]
    non_idx = np.where(~epi_mask)[0]
    probes = probe_indices(ch)
    patient_info[pat] = dict(
        ch=ch, N=N, epi_mask=epi_mask, epi_idx=epi_idx, non_idx=non_idx,
        epi_in=epi_in, probes=probes,
    )

R("=" * 90)
R("EPILEPTIC NODES — DEEP LRG & NETWORK-TOPOLOGY ANALYSIS (ImCoh)")
R("=" * 90)
R()
for pat in PATIENTS:
    info = patient_info[pat]
    R(f"  {pat}: {info['N']} ch, {len(info['epi_in'])} epileptic in analysis, "
      f"{len(info['probes'])} probes")
R()

# ═══════════════════════════════════════════════════════════════════════
#  A1 — PROBE-MATCHED PERMUTATION TEST ON ULTRAMETRIC CLOSENESS
# ═══════════════════════════════════════════════════════════════════════

R("A1: Probe-matched permutation test — does ultrametric closeness survive?")
R("-" * 90)
R()
R("The control that killed the MSC analysis. For each patient we build a")
R("null by replacing each epileptic probe with a random probe of the same")
R(f"size ({N_PERMS} permutations) and compute intra-group ultrametric distance.")
R()

perm_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_idx, ch = info["N"], info["epi_idx"], info["ch"]
    n_epi = len(epi_idx)
    if n_epi < 2:
        continue

    # Build epileptic probe structure: {probe_name: count_of_epi_contacts}
    epi_probe_sizes = {}
    for i in epi_idx:
        pr = probe(ch[i])
        epi_probe_sizes[pr] = epi_probe_sizes.get(pr, 0) + 1

    all_probes = info["probes"]  # probe_name → [indices]

    for phase in PHASES:
        for band in BANDS:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None:
                continue
            um = squareform(lrg.ultrametric_matrix)

            # Observed intra-epileptic mean ultrametric distance
            epi_dists = []
            for ii in range(n_epi):
                for jj in range(ii + 1, n_epi):
                    epi_dists.append(um[epi_idx[ii], epi_idx[jj]])
            obs_mean = np.mean(epi_dists)

            # Probe-matched null
            null_means = np.zeros(N_PERMS)
            probe_names = list(all_probes.keys())
            for pi in range(N_PERMS):
                selected = []
                for ep_probe, ep_count in epi_probe_sizes.items():
                    # Pick a random probe with enough contacts
                    avail = [p for p in probe_names if len(all_probes[p]) >= ep_count]
                    if not avail:
                        avail = probe_names
                    rp = avail[np.random.randint(len(avail))]
                    contacts = all_probes[rp]
                    max_start = max(1, len(contacts) - ep_count + 1)
                    start = np.random.randint(max_start)
                    selected.extend(contacts[start:start + ep_count])

                # Compute intra-group mean ultrametric distance
                n_sel = len(selected)
                if n_sel < 2:
                    null_means[pi] = np.inf
                    continue
                dists = []
                for ii in range(n_sel):
                    for jj in range(ii + 1, n_sel):
                        dists.append(um[selected[ii], selected[jj]])
                null_means[pi] = np.mean(dists)

            p_val = (np.sum(null_means <= obs_mean) + 1) / (N_PERMS + 1)
            perm_rows.append(dict(
                patient=pat, phase=phase, band=band,
                obs_mean=obs_mean, null_mean=np.mean(null_means),
                null_std=np.std(null_means), p_value=p_val,
                z_score=(obs_mean - np.mean(null_means)) / (np.std(null_means) + 1e-30),
            ))
            if VERBOSE:
                sig = "***" if p_val < .001 else "** " if p_val < .01 else "*  " if p_val < .05 else "   "
                R(f"  {pat} {phase:10s} {band:12s}: obs={obs_mean:.4f} "
                  f"null={np.mean(null_means):.4f}±{np.std(null_means):.4f} "
                  f"p={p_val:.4f} {sig}")

permdf = pd.DataFrame(perm_rows)
n_sig = int((permdf["p_value"] < 0.05).sum())
n_sig_bonf = int((permdf["p_value"] < 0.05 / len(permdf)).sum())
R(f"Total tests: {len(permdf)}")
R(f"Significant at p < 0.05: {n_sig}/{len(permdf)} ({n_sig / max(len(permdf), 1):.0%})")
R(f"Significant after Bonferroni (α = {0.05 / max(len(permdf), 1):.2e}): {n_sig_bonf}")
R()
for pat in PATIENTS:
    sub = permdf[permdf["patient"] == pat]
    if len(sub) == 0:
        continue
    ns = int((sub["p_value"] < 0.05).sum())
    mz = sub["z_score"].mean()
    R(f"  {pat}: {ns}/{len(sub)} sig (p<0.05), mean z-score = {mz:.2f}")
R()

# ── Figure A1: Probe-matched null ─────────────────────────────────────

fig, axes = plt.subplots(1, len(PATIENTS), figsize=(3.5 * len(PATIENTS), 5),
                         sharey=True)
for ax, pat in zip(axes, PATIENTS):
    sub = permdf[permdf["patient"] == pat]
    if len(sub) == 0:
        ax.set_title(f"{pat}\n(no data)")
        continue
    ax.scatter(range(len(sub)), sub["z_score"], s=15, alpha=0.5,
               c=[BAND_COLORS[b] for b in sub["band"]], edgecolors="none")
    ax.axhline(0, color="k", ls="-", lw=0.5)
    ax.axhline(-1.96, color="gray", ls="--", lw=0.8, label="p=0.05")
    ns = int((sub["p_value"] < 0.05).sum())
    ax.set_title(f"{pat}\n{ns}/{len(sub)} sig", fontsize=10)
    if ax == axes[0]:
        ax.set_ylabel("z-score\n(negative = epi closer than null)")
    ax.set_xlabel("condition index")
fig.suptitle("A1: Probe-matched permutation test — ultrametric closeness\n"
             "(z < -1.96 = epileptic pairs closer than probe-matched null)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.88])
fig.savefig(OUT / "figA1_probe_matched_null.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figA1_probe_matched_null.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figA1_probe_matched_null.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  A2 — MERGE-HEIGHT PROFILES
# ═══════════════════════════════════════════════════════════════════════

R("A2: Merge-height profiles — when do epileptic nodes merge?")
R("-" * 90)
R()

merge_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx, non_idx = info["epi_idx"], info["non_idx"]
    if len(epi_idx) < 2:
        continue

    for phase in ["rest_pre", "rest_post"]:
        for band in ["alpha", "beta"]:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None:
                continue
            um = squareform(lrg.ultrametric_matrix)

            # Three groups of pairwise distances
            epi_epi, epi_non, non_non = [], [], []
            for i in range(N):
                for j in range(i + 1, N):
                    d = um[i, j]
                    if epi_mask[i] and epi_mask[j]:
                        epi_epi.append(d)
                    elif epi_mask[i] or epi_mask[j]:
                        epi_non.append(d)
                    else:
                        non_non.append(d)

            merge_rows.append(dict(
                patient=pat, phase=phase, band=band,
                epi_epi_mean=np.mean(epi_epi), epi_epi_med=np.median(epi_epi),
                epi_non_mean=np.mean(epi_non), epi_non_med=np.median(epi_non),
                non_non_mean=np.mean(non_non), non_non_med=np.median(non_non),
                n_epi_epi=len(epi_epi), n_epi_non=len(epi_non), n_non_non=len(non_non),
            ))

mergedf = pd.DataFrame(merge_rows)
R(f"{'patient':8s} {'phase':8s} {'band':6s}  {'epi-epi':>10s} {'epi-non':>10s} {'non-non':>10s}  ratio ee/nn")
R("-" * 75)
for _, row in mergedf.iterrows():
    ratio = row["epi_epi_med"] / (row["non_non_med"] + 1e-30)
    R(f"{row['patient']:8s} {row['phase']:8s} {row['band']:6s}  "
      f"{row['epi_epi_med']:10.4f} {row['epi_non_med']:10.4f} {row['non_non_med']:10.4f}  "
      f"{ratio:.3f}")
R()

# ── Figure A2: Merge-height distributions ─────────────────────────────

fig, axes = plt.subplots(2, len(PATIENTS), figsize=(3.5 * len(PATIENTS), 8),
                         sharey="row")
for col, pat in enumerate(PATIENTS):
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx = info["epi_idx"]
    if len(epi_idx) < 2:
        continue

    for row_idx, (phase, band) in enumerate([("rest_pre", "alpha"), ("rest_pre", "beta")]):
        ax = axes[row_idx, col]
        lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
        if lrg is None:
            continue
        um = squareform(lrg.ultrametric_matrix)

        epi_epi, epi_non, non_non = [], [], []
        for i in range(N):
            for j in range(i + 1, N):
                d = um[i, j]
                if epi_mask[i] and epi_mask[j]:
                    epi_epi.append(d)
                elif epi_mask[i] or epi_mask[j]:
                    epi_non.append(d)
                else:
                    non_non.append(d)

        bins = np.linspace(0, max(np.max(non_non), np.max(epi_epi)), 40)
        ax.hist(non_non, bins=bins, density=True, alpha=0.3, color="gray",
                label=f"non–non (n={len(non_non)})")
        ax.hist(epi_non, bins=bins, density=True, alpha=0.4, color="#3498db",
                label=f"epi–non (n={len(epi_non)})")
        ax.hist(epi_epi, bins=bins, density=True, alpha=0.6, color="#e74c3c",
                label=f"epi–epi (n={len(epi_epi)})")
        ax.axvline(np.median(epi_epi), color="#e74c3c", ls="--", lw=1.5)
        ax.axvline(np.median(non_non), color="gray", ls="--", lw=1.5)
        if col == 0:
            ax.set_ylabel(f"{band}\nDensity")
        if row_idx == 0:
            ax.set_title(pat, fontsize=11)
        if row_idx == 1:
            ax.set_xlabel("Ultrametric distance")
        ax.legend(fontsize=6, loc="upper right")

fig.suptitle("A2: Merge-height distributions in ImCoh dendrogram (rest_pre)\n"
             "Dashed = median; red left-shifted = epileptic pairs merge earlier",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT / "figA2_merge_heights.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figA2_merge_heights.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figA2_merge_heights.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  A3 — ULTRAMETRIC NEAREST-NEIGHBOR ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════

R("A3: Ultrametric nearest-neighbor — is the closest node also epileptic?")
R("-" * 90)
R()

nn_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx = info["epi_idx"]
    n_epi = len(epi_idx)
    if n_epi < 2:
        continue
    frac_epi = n_epi / N  # baseline probability

    for phase in PHASES:
        for band in BANDS:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None:
                continue
            um = squareform(lrg.ultrametric_matrix)

            # For each epileptic node, find its nearest neighbor
            nn_is_epi = 0
            for i in epi_idx:
                dists = um[i].copy()
                dists[i] = np.inf  # exclude self
                nn = np.argmin(dists)
                if epi_mask[nn]:
                    nn_is_epi += 1

            obs_frac = nn_is_epi / n_epi
            enrichment = obs_frac / (frac_epi + 1e-30)
            # Binomial test
            p_val = binomtest(nn_is_epi, n_epi, frac_epi, alternative="greater").pvalue

            nn_rows.append(dict(
                patient=pat, phase=phase, band=band,
                nn_is_epi=nn_is_epi, n_epi=n_epi,
                obs_frac=obs_frac, expected_frac=frac_epi,
                enrichment=enrichment, p_value=p_val,
            ))

nndf = pd.DataFrame(nn_rows)
n_sig_nn = int((nndf["p_value"] < 0.05).sum())
R(f"Tests: {len(nndf)}")
R(f"NN is epileptic more than expected (p < 0.05): {n_sig_nn}/{len(nndf)}")
R()
for pat in PATIENTS:
    sub = nndf[nndf["patient"] == pat]
    if len(sub) == 0:
        continue
    me = sub["enrichment"].mean()
    ns = int((sub["p_value"] < 0.05).sum())
    R(f"  {pat}: mean enrichment = {me:.2f}x, sig in {ns}/{len(sub)}")
R()

# ── Figure A3: NN enrichment ─────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: enrichment by patient × band
ax = axes[0]
mat = np.full((len(PATIENTS), len(BANDS)), np.nan)
for pi, pat in enumerate(PATIENTS):
    for bi, band in enumerate(BANDS):
        sub = nndf[(nndf["patient"] == pat) & (nndf["band"] == band)]
        if len(sub) > 0:
            mat[pi, bi] = sub["enrichment"].mean()
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0, vmax=6)
ax.set_yticks(range(len(PATIENTS)))
ax.set_yticklabels(PATIENTS)
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=8)
ax.set_title("NN enrichment (obs/expected)", fontsize=10)
for pi in range(len(PATIENTS)):
    for bi in range(len(BANDS)):
        if not np.isnan(mat[pi, bi]):
            ax.text(bi, pi, f"{mat[pi, bi]:.1f}", ha="center", va="center", fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.8)

# Panel B: fraction of sig conditions per patient
ax = axes[1]
for pi, pat in enumerate(PATIENTS):
    sub = nndf[nndf["patient"] == pat]
    if len(sub) == 0:
        continue
    ns = int((sub["p_value"] < 0.05).sum())
    ax.bar(pi, ns / len(sub) * 100, color=PAT_COLORS[pat], alpha=0.7)
    ax.text(pi, ns / len(sub) * 100 + 1, f"{ns}/{len(sub)}", ha="center", fontsize=9)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels(PATIENTS)
ax.set_ylabel("% conditions with p < 0.05")
ax.set_title("Nearest-neighbor enrichment significance", fontsize=10)
ax.set_ylim(0, 100)

fig.suptitle("A3: Ultrametric nearest-neighbor enrichment\n"
             "(Is the dendrogram-closest node to an epileptic node also epileptic?)",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.88])
fig.savefig(OUT / "figA3_nn_enrichment.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figA3_nn_enrichment.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figA3_nn_enrichment.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  A4 — DENDROGRAM VISUALISATION WITH EPILEPTIC NODES COLOURED
# ═══════════════════════════════════════════════════════════════════════

R("A4: Dendrograms with epileptic nodes highlighted (rest_pre, alpha & beta)")
R("-" * 90)
R()

for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask, ch = info["N"], info["epi_mask"], info["ch"]

    for band in ["alpha", "beta"]:
        lrg = load_lrg_result(pat, "rest_pre", band, "imcoh_abs")
        if lrg is None:
            continue

        fig, ax = plt.subplots(figsize=(max(14, N * 0.12), 6))
        Z = lrg.linkage_matrix

        # Color each leaf
        leaf_colors = {}
        for i in range(N):
            leaf_colors[i] = "#e74c3c" if epi_mask[i] else "#333333"

        def color_func(k):
            if k < N:
                return leaf_colors.get(k, "#333333")
            return "#aaaaaa"

        dn = dendrogram(
            Z, ax=ax, orientation="top",
            labels=[ch[i] for i in range(N)],
            leaf_rotation=90, leaf_font_size=5,
            link_color_func=color_func,
            above_threshold_color="#aaaaaa",
        )

        # Color x-axis labels
        xlabels = ax.get_xticklabels()
        for lbl in xlabels:
            txt = lbl.get_text()
            if txt in set(info["epi_in"]):
                lbl.set_color("#e74c3c")
                lbl.set_fontweight("bold")

        ax.set_ylabel("Ultrametric distance")
        ax.set_title(f"{pat} — rest_pre {band} — ImCoh LRG dendrogram\n"
                     f"Red labels = epileptic contacts ({len(info['epi_in'])}/{N})",
                     fontsize=11)
        fig.tight_layout()
        fname = f"figA4_dendro_{pat}_{band}.pdf"
        fig.savefig(OUT / fname, bbox_inches="tight", dpi=200)
        fig.savefig(OUT / fname.replace(".pdf", ".png"), bbox_inches="tight", dpi=150)
        plt.close(fig)
        R(f"  → {OUT / fname}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  A5 — BAND-SPECIFIC BREAKDOWN OF Q3
# ═══════════════════════════════════════════════════════════════════════

R("A5: Band-specific breakdown — which frequencies drive the closeness?")
R("-" * 90)
R()

q3_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx, non_idx = info["epi_idx"], info["non_idx"]
    if len(epi_idx) < 2:
        continue

    for band in BANDS:
        ratios = []
        for phase in PHASES:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None:
                continue
            um = squareform(lrg.ultrametric_matrix)
            epi_d = [um[epi_idx[i], epi_idx[j]]
                     for i in range(len(epi_idx)) for j in range(i + 1, len(epi_idx))]
            non_d = [um[non_idx[i], non_idx[j]]
                     for i in range(len(non_idx)) for j in range(i + 1, len(non_idx))]
            ratios.append(np.mean(epi_d) / (np.mean(non_d) + 1e-30))

        if ratios:
            q3_rows.append(dict(
                patient=pat, band=band,
                mean_ratio=np.mean(ratios), std_ratio=np.std(ratios),
                n_phases=len(ratios),
            ))

q3df = pd.DataFrame(q3_rows)

# ── Figure A5: Band breakdown ─────────────────────────────────────────

fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(len(BANDS))
width = 0.15
for pi, pat in enumerate(PATIENTS):
    sub = q3df[q3df["patient"] == pat]
    vals = [sub[sub["band"] == b]["mean_ratio"].values[0]
            if len(sub[sub["band"] == b]) > 0 else np.nan for b in BANDS]
    ax.bar(x + pi * width, vals, width, label=pat,
           color=PAT_COLORS[pat], alpha=0.7)
ax.axhline(1, color="k", ls="--", lw=0.8)
ax.set_xticks(x + 2 * width)
ax.set_xticklabels(BANDS)
ax.set_ylabel("Ultrametric distance ratio\n(epi-epi / non-non)")
ax.set_title("A5: Band-specific ultrametric closeness\n"
             "(< 1 = epileptic pairs closer; averaged over 4 phases)",
             fontsize=11)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figA5_band_breakdown.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figA5_band_breakdown.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figA5_band_breakdown.pdf'}")
R()

# Consensus: which bands have ALL patients < 1?
R("Band consensus (all patients ratio < 1):")
for band in BANDS:
    sub = q3df[q3df["band"] == band]
    n_below = int((sub["mean_ratio"] < 1).sum())
    mean_r = sub["mean_ratio"].mean()
    R(f"  {band:12s}: {n_below}/{len(sub)} patients below 1, mean ratio = {mean_r:.3f}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  B — STANDARD NETWORK TOPOLOGY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

R("=" * 90)
R("PART B: STANDARD NETWORK TOPOLOGY — EPILEPTIC NODE PROPERTIES")
R("=" * 90)
R()

metric_names = ["strength", "clustering", "betweenness", "eigvec_cent",
                "participation", "within_mod_z"]
metric_labels_short = ["Strength", "Clust. Coeff.", "Betweenness",
                       "Eigvec. Cent.", "Participation", "Within-mod z"]

topo_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx, non_idx = info["epi_idx"], info["non_idx"]
    if len(epi_idx) < 2 or len(non_idx) < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            A = load_A(pat, phase, band, N)
            if A is None:
                continue

            # Compute communities at k=10 for participation/within-mod-z
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None:
                continue
            comms = fcluster(lrg.linkage_matrix, t=10, criterion="maxclust")

            # Compute metrics
            strength = A.sum(axis=1)
            cc = weighted_clustering_coefficient(A)
            bc = betweenness_centrality(A)
            ec = eigenvector_centrality(A)
            pc = participation_coefficient(A, comms)
            wz = within_module_degree_zscore(A, comms)

            metrics = dict(
                strength=strength, clustering=cc, betweenness=bc,
                eigvec_cent=ec, participation=pc, within_mod_z=wz,
            )

            for mname, mvals in metrics.items():
                ve = mvals[epi_idx]
                vn = mvals[non_idx]
                _, p_gt = mannwhitneyu(ve, vn, alternative="greater")
                _, p_lt = mannwhitneyu(ve, vn, alternative="less")
                p_best = min(p_gt, p_lt)
                direction = "epi>" if p_gt < p_lt else "epi<"
                try:
                    y = np.zeros(N)
                    y[epi_idx] = 1
                    pool = np.concatenate([epi_idx, non_idx])
                    auc = roc_auc_score(y[pool], mvals[pool])
                    auc = max(auc, 1 - auc)
                except Exception:
                    auc = 0.5

                topo_rows.append(dict(
                    patient=pat, phase=phase, band=band,
                    metric=mname, direction=direction,
                    p=p_best, auc=auc,
                    mean_epi=ve.mean(), mean_non=vn.mean(),
                    med_epi=np.median(ve), med_non=np.median(vn),
                ))

topodf = pd.DataFrame(topo_rows)

R(f"{'metric':16s} {'epi<':>5s} {'epi>':>5s} {'p<.05':>6s} {'AUC':>6s} "
  f"{'dominant':>9s} {'consist':>8s}")
R("-" * 65)
for mname, mlabel in zip(metric_names, metric_labels_short):
    sub = topodf[topodf["metric"] == mname]
    if len(sub) == 0:
        continue
    nl = int((sub["direction"] == "epi<").sum())
    ng = int((sub["direction"] == "epi>").sum())
    ns = int((sub["p"] < 0.05).sum())
    ma = sub["auc"].mean()
    dom = "epi<" if nl > ng else "epi>"
    con = max(nl, ng) / len(sub)
    R(f"{mlabel:16s} {nl:5d} {ng:5d} {ns:6d} {ma:6.3f} {dom:>9s} {con:8.0%}")
R()

# Per-patient breakdown
R("Per-patient dominant direction:")
for mname, mlabel in zip(metric_names, metric_labels_short):
    parts = []
    for pat in PATIENTS:
        sub = topodf[(topodf["metric"] == mname) & (topodf["patient"] == pat)]
        if len(sub) == 0:
            parts.append(f"{pat}: ?")
            continue
        nl = int((sub["direction"] == "epi<").sum())
        ng = int((sub["direction"] == "epi>").sum())
        dom = ">" if ng > nl else "<" if nl > ng else "="
        parts.append(f"{pat}:{dom}")
    R(f"  {mlabel:16s}  {'  '.join(parts)}")
R()

# ── Figure B1: Network topology summary ───────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for mi, (mname, mlabel) in enumerate(zip(metric_names, metric_labels_short)):
    ax = axes.flat[mi]
    # Per-patient boxplot of AUC values
    bp_data = []
    bp_labels = []
    for pat in PATIENTS:
        sub = topodf[(topodf["metric"] == mname) & (topodf["patient"] == pat)]
        if len(sub) > 0:
            bp_data.append(sub["auc"].values)
            bp_labels.append(pat)

    if bp_data:
        bp = ax.boxplot(bp_data, labels=bp_labels, widths=0.5,
                        showfliers=False, patch_artist=True)
        for patch, pat in zip(bp["boxes"], bp_labels):
            patch.set_facecolor(PAT_COLORS[pat])
            patch.set_alpha(0.4)
    ax.axhline(0.5, color="gray", ls="--", lw=0.8)
    ax.set_title(mlabel, fontsize=10)
    ax.set_ylim(0.4, 1.0)
    if mi % 3 == 0:
        ax.set_ylabel("AUC")
    ax.tick_params(axis="x", rotation=30, labelsize=8)

fig.suptitle("B: Network topology — epileptic vs non-epileptic (ImCoh)\n"
             "AUC > 0.5 = metric distinguishes groups",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT / "figB1_network_topology.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figB1_network_topology.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figB1_network_topology.pdf'}")

# ═══════════════════════════════════════════════════════════════════════
#  B7 — INTRA-EPILEPTIC EDGE WEIGHT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

R()
R("B7: Intra-epileptic edge weight — are epi-epi ImCoh edges stronger?")
R("-" * 90)
R()

edge_rows = []
for pat in PATIENTS:
    info = patient_info[pat]
    N, epi_mask = info["N"], info["epi_mask"]
    epi_idx, non_idx = info["epi_idx"], info["non_idx"]
    if len(epi_idx) < 2:
        continue

    for phase in PHASES:
        for band in BANDS:
            A = load_A(pat, phase, band, N)
            if A is None:
                continue

            # Three groups of edge weights
            epi_epi_w = A[np.ix_(epi_idx, epi_idx)][np.triu_indices(len(epi_idx), k=1)]
            non_non_w = A[np.ix_(non_idx, non_idx)][np.triu_indices(len(non_idx), k=1)]
            cross_w = A[np.ix_(epi_idx, non_idx)].ravel()

            if len(epi_epi_w) > 0 and len(non_non_w) > 0:
                _, p_vs_non = mannwhitneyu(epi_epi_w, non_non_w, alternative="greater")
                _, p_vs_cross = mannwhitneyu(epi_epi_w, cross_w, alternative="greater")
            else:
                p_vs_non = p_vs_cross = 1.0

            edge_rows.append(dict(
                patient=pat, phase=phase, band=band,
                mean_epi_epi=epi_epi_w.mean(), mean_non_non=non_non_w.mean(),
                mean_cross=cross_w.mean(),
                ratio_vs_non=epi_epi_w.mean() / (non_non_w.mean() + 1e-30),
                ratio_vs_cross=epi_epi_w.mean() / (cross_w.mean() + 1e-30),
                p_vs_non=p_vs_non, p_vs_cross=p_vs_cross,
            ))

edgedf = pd.DataFrame(edge_rows)
n_stronger = int((edgedf["ratio_vs_non"] > 1).sum())
n_sig_edge = int((edgedf["p_vs_non"] < 0.05).sum())
R(f"Tests: {len(edgedf)}")
R(f"Epi-epi edges STRONGER than non-non: {n_stronger}/{len(edgedf)} ({n_stronger / max(len(edgedf),1):.0%})")
R(f"Significantly stronger (p < 0.05): {n_sig_edge}/{len(edgedf)}")
R()
for pat in PATIENTS:
    sub = edgedf[edgedf["patient"] == pat]
    if len(sub) == 0:
        continue
    mr = sub["ratio_vs_non"].mean()
    ns = int((sub["p_vs_non"] < 0.05).sum())
    R(f"  {pat}: mean epi-epi/non-non = {mr:.3f}, sig in {ns}/{len(sub)}")
R()

# ── Figure B7: Edge weight comparison ─────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: ratio by patient × band
ax = axes[0]
mat = np.full((len(PATIENTS), len(BANDS)), np.nan)
for pi, pat in enumerate(PATIENTS):
    for bi, band in enumerate(BANDS):
        sub = edgedf[(edgedf["patient"] == pat) & (edgedf["band"] == band)]
        if len(sub) > 0:
            mat[pi, bi] = sub["ratio_vs_non"].mean()
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0.7, vmax=1.5)
ax.set_yticks(range(len(PATIENTS)))
ax.set_yticklabels(PATIENTS)
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(BANDS, rotation=45, ha="right", fontsize=8)
ax.set_title("Epi-epi / non-non edge ratio\n(> 1 = epi-epi stronger)", fontsize=10)
for pi in range(len(PATIENTS)):
    for bi in range(len(BANDS)):
        if not np.isnan(mat[pi, bi]):
            ax.text(bi, pi, f"{mat[pi, bi]:.2f}", ha="center", va="center", fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.8)

# Panel B: three-way comparison for rest_pre/alpha
ax = axes[1]
for pi, pat in enumerate(PATIENTS):
    sub = edgedf[(edgedf["patient"] == pat) & (edgedf["band"] == "alpha") &
                 (edgedf["phase"] == "rest_pre")]
    if len(sub) == 0:
        continue
    vals = [sub["mean_epi_epi"].values[0], sub["mean_cross"].values[0],
            sub["mean_non_non"].values[0]]
    x_pos = np.array([0, 1, 2]) + pi * 0.12
    ax.plot(x_pos, vals, "o-", color=PAT_COLORS[pat], label=pat, markersize=6)
ax.set_xticks([0.24, 1.24, 2.24])
ax.set_xticklabels(["epi–epi", "epi–non", "non–non"])
ax.set_ylabel("Mean ImCoh edge weight")
ax.set_title("Edge weight by group pair\n(rest_pre, alpha)", fontsize=10)
ax.legend(fontsize=8)

fig.suptitle("B7: Intra-epileptic edge strength (ImCoh)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT / "figB7_edge_weights.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figB7_edge_weights.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'figB7_edge_weights.pdf'}")
R()

# ═══════════════════════════════════════════════════════════════════════
#  GRAND SUMMARY
# ═══════════════════════════════════════════════════════════════════════

R("=" * 90)
R("GRAND SUMMARY")
R("=" * 90)
R()

R("A1. Probe-matched null (ultrametric closeness):")
if len(permdf) > 0:
    R(f"     {n_sig}/{len(permdf)} survive at p < 0.05, "
      f"{n_sig_bonf} after Bonferroni")
    for pat in PATIENTS:
        sub = permdf[permdf["patient"] == pat]
        if len(sub) > 0:
            ns = int((sub["p_value"] < 0.05).sum())
            R(f"     {pat}: {ns}/{len(sub)} sig")
R()

R("A2. Merge heights (median epi-epi / median non-non):")
if len(mergedf) > 0:
    for _, row in mergedf.iterrows():
        ratio = row["epi_epi_med"] / (row["non_non_med"] + 1e-30)
        R(f"     {row['patient']} {row['band']:6s}: {ratio:.3f}")
R()

R("A3. Nearest-neighbor enrichment:")
if len(nndf) > 0:
    for pat in PATIENTS:
        sub = nndf[nndf["patient"] == pat]
        if len(sub) > 0:
            R(f"     {pat}: mean {sub['enrichment'].mean():.1f}x, "
              f"sig {int((sub['p_value']<0.05).sum())}/{len(sub)}")
R()

R("A5. Band breakdown (ratio < 1 = epi closer):")
if len(q3df) > 0:
    for band in BANDS:
        sub = q3df[q3df["band"] == band]
        if len(sub) > 0:
            R(f"     {band:12s}: mean ratio = {sub['mean_ratio'].mean():.3f}, "
              f"{int((sub['mean_ratio']<1).sum())}/5 patients below 1")
R()

R("B. Network topology (best discriminative metric):")
if len(topodf) > 0:
    best_auc = 0
    best_metric = ""
    for mname in metric_names:
        sub = topodf[topodf["metric"] == mname]
        if len(sub) > 0 and sub["auc"].mean() > best_auc:
            best_auc = sub["auc"].mean()
            best_metric = mname
    R(f"     Best: {best_metric} (mean AUC = {best_auc:.3f})")
    for mname, mlabel in zip(metric_names, metric_labels_short):
        sub = topodf[topodf["metric"] == mname]
        if len(sub) > 0:
            R(f"     {mlabel:16s}: AUC = {sub['auc'].mean():.3f}")
R()

R("B7. Intra-epileptic edge strength:")
if len(edgedf) > 0:
    R(f"     Epi-epi stronger than non-non: {n_stronger}/{len(edgedf)}")
    for pat in PATIENTS:
        sub = edgedf[edgedf["patient"] == pat]
        if len(sub) > 0:
            R(f"     {pat}: mean ratio = {sub['ratio_vs_non'].mean():.3f}")
R()

# ── Save ──────────────────────────────────────────────────────────────

report_path = OUT / "report.txt"
report_path.write_text("\n".join(LINES))
R(f"\nReport: {report_path}")
R(f"Figures: {OUT}")

# Save key DataFrames for downstream use
permdf.to_csv(OUT / "probe_matched_null.csv", index=False)
nndf.to_csv(OUT / "nn_enrichment.csv", index=False)
topodf.to_csv(OUT / "network_topology.csv", index=False)
edgedf.to_csv(OUT / "edge_weights.csv", index=False)
q3df.to_csv(OUT / "band_breakdown.csv", index=False)
