#!/usr/bin/env python3
"""H2d — pairwise coactivation persistence.

At each dendrogram cut k, define the co-clustering indicator:
    C_phase(i, j, k) = 1 if nodes i and j share a cluster in that phase's LRG.

Partition pairs by pre-task × task_test status:
  - 'task_induced' pairs: C_rpre = 0 AND C_ttest = 1  (newly co-active in task)
  - 'inert'        pairs: C_rpre = 0 AND C_ttest = 0  (baseline-inactive pairs)

Persistence rate:
    rho_task  = P(C_rpost = 1 | task-induced)
    rho_inert = P(C_rpost = 1 | inert)

Hypothesis (one-sided):
    rho_task > rho_inert  per patient per band at each k.

Cross-patient test: one-sample Wilcoxon signed-rank on Δρ = rho_task - rho_inert
(k-averaged, or at each k with cluster-based permutation over k).

Also reports the ABSOLUTE persistence rate rho_task per (patient, band) —
readers may want to see "X% of pairs newly co-clustered in task stay
co-clustered in post-rest".

Outputs:
    data/reports/imcoh_vi/h2d_coactivation_persistence.md
    data/reports/imcoh_vi/h2d_coactivation_persistence.csv
"""
from __future__ import annotations

from itertools import groupby
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

K_RANGE = list(range(2, 50))         # multiscale sweep
N_PERM = 5_000
CLUSTER_Z_THRESH = 1.96


# ─────────────────────────── helpers ───────────────────────────

from _shared import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr  # canonical


def co_cluster_mask(labels: np.ndarray) -> np.ndarray:
    """N×N symmetric boolean: True if i, j share a cluster."""
    return labels[:, None] == labels[None, :]


def load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.linkage_matrix)


def compute_persistence(Z_rpre, Z_ttest, Z_rpost, k: int) -> tuple[float, float] | None:
    """Return (rho_task, rho_inert) at cut k; None if shapes mismatch or no task-induced pairs."""
    if Z_rpre is None or Z_ttest is None or Z_rpost is None:
        return None
    n_rpre = Z_rpre.shape[0] + 1
    n_ttest = Z_ttest.shape[0] + 1
    n_rpost = Z_rpost.shape[0] + 1
    if not (n_rpre == n_ttest == n_rpost):
        return None
    lab_pre  = fcluster(Z_rpre,  k, criterion="maxclust")
    lab_tt   = fcluster(Z_ttest, k, criterion="maxclust")
    lab_post = fcluster(Z_rpost, k, criterion="maxclust")
    C_pre  = co_cluster_mask(lab_pre)
    C_tt   = co_cluster_mask(lab_tt)
    C_post = co_cluster_mask(lab_post)
    # Upper triangle
    iu = np.triu_indices(n_rpre, k=1)
    c_pre  = C_pre[iu]
    c_tt   = C_tt[iu]
    c_post = C_post[iu]
    task_induced = (~c_pre) & c_tt
    inert        = (~c_pre) & (~c_tt)
    if task_induced.sum() == 0 or inert.sum() == 0:
        return None
    rho_task  = float(c_post[task_induced].mean())
    rho_inert = float(c_post[inert].mean())
    return rho_task, rho_inert


# ── cluster-based permutation over k for Δρ ──
def cluster_stats(z: np.ndarray, thresh: float) -> list[tuple[int, int, float]]:
    sup = np.isfinite(z) & (z > thresh)
    out = []
    in_run = False
    start = 0
    mass = 0.0
    for i, s in enumerate(sup):
        if s and not in_run:
            in_run = True; start = i; mass = float(z[i])
        elif s:
            mass += float(z[i])
        elif in_run:
            out.append((start, i - 1, mass)); in_run = False
    if in_run:
        out.append((start, len(sup) - 1, mass))
    return out


def per_k_z(mat: np.ndarray) -> np.ndarray:
    abs_ranks = np.apply_along_axis(stats.rankdata, 0, np.abs(mat))
    signs = np.sign(mat)
    w_plus = (abs_ranks * (signs > 0)).sum(axis=0)
    n = (signs != 0).sum(axis=0)
    mu = n * (n + 1) / 4.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sigma > 0, (w_plus - mu) / sigma, np.nan)


def cluster_perm(mat: np.ndarray, k_values: np.ndarray,
                 n_perm: int = N_PERM, thresh: float = CLUSTER_Z_THRESH) -> list[dict]:
    rng = np.random.default_rng(0)
    abs_ranks = np.apply_along_axis(stats.rankdata, 0, np.abs(mat))
    obs_signs = np.sign(mat)
    n = (obs_signs != 0).sum(axis=0)
    mu = n * (n + 1) / 4.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)

    def _z(signs):
        w_plus = (abs_ranks * (signs > 0)).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(sigma > 0, (w_plus - mu) / sigma, np.nan)

    z_obs = _z(obs_signs)
    obs_clusters = cluster_stats(z_obs, thresh)
    if not obs_clusters:
        return []

    n_pat = mat.shape[0]
    flips = rng.choice([-1.0, 1.0], size=(n_perm, n_pat))
    max_null = np.zeros(n_perm)
    for i in range(n_perm):
        z_null = _z(obs_signs * flips[i, :, None])
        cs = cluster_stats(z_null, thresh)
        max_null[i] = max((c[2] for c in cs), default=0.0)

    out = []
    for s, e, m in obs_clusters:
        p = float((max_null >= m).mean())
        out.append({"k_start": int(k_values[s]), "k_end": int(k_values[e]),
                    "len": e - s + 1, "mass": m, "p": p})
    return out


# ─────────────────────────── main ───────────────────────────

def main() -> None:
    patients = list(PATIENTS_LIST)
    raw = []
    for pat in patients:
        for band in BRAIN_BANDS_NAMES:
            Z_rpre  = load_Z(pat, "rest_pre",  band)
            Z_ttest = load_Z(pat, "task_test", band)
            Z_rpost = load_Z(pat, "rest_post", band)
            for k in K_RANGE:
                out = compute_persistence(Z_rpre, Z_ttest, Z_rpost, k)
                if out is None:
                    continue
                rho_task, rho_inert = out
                raw.append({
                    "patient": pat, "band": band, "k": k,
                    "rho_task": rho_task,
                    "rho_inert": rho_inert,
                    "delta_rho": rho_task - rho_inert,
                })
    rawdf = pd.DataFrame(raw)
    rawdf.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2d_persistence_raw.csv", index=False)

    lines: list[str] = []
    ap = lines.append
    ap("# H2d — pairwise coactivation persistence")
    ap("")
    ap("**Metric per (patient, band, k):**")
    ap("- `ρ_task = P(pair co-clusters in rpost | pair was co-clustered in ttest but NOT in rpre)`  ")
    ap("- `ρ_inert = P(pair co-clusters in rpost | pair NOT co-clustered in either rpre or ttest)`  ")
    ap("- `Δρ = ρ_task − ρ_inert` (task-trace excess)")
    ap("")
    ap("**Hypothesis:** Δρ > 0 — pairs newly co-active in task persist into "
       "post-rest at a higher rate than pairs that were never co-active.")
    ap("")
    ap(f"k range: {K_RANGE[0]}..{K_RANGE[-1]} (multiscale sweep).")
    ap("")
    ap("Patients with rest/task node-count mismatch (Pat_10), missing phases "
       "(Pat_13 rpre, Pat_14 ttest) are excluded automatically.")
    ap("")

    # Per-patient k-averaged Δρ per band
    pat_avg = (rawdf
               .groupby(["patient", "band"])[["delta_rho", "rho_task", "rho_inert"]]
               .mean()
               .reset_index())

    # Cross-patient Wilcoxon per band on k-averaged Δρ
    records = []
    ap("## Patient-level k-averaged test (k = 2..29)")
    ap("")
    ap("| band | n | mean Δρ | 95% CI | mean ρ_task | mean ρ_inert | r_rb | z | p | q (BH) | LOO worst p | p no-Pat_03 |")
    ap("|------|--:|--------:|:------|------------:|-------------:|-----:|--:|--:|-------:|------------:|------------:|")
    band_p = {}
    band_rows = []
    for band in BRAIN_BANDS_NAMES:
        bdf = pat_avg[pat_avg["band"] == band]
        pats = bdf["patient"].tolist()
        vals = bdf["delta_rho"].to_numpy()
        if len(vals) < 3:
            band_rows.append({"band": band, "n": len(vals)})
            continue
        z, p = wilcoxon_z(vals)
        r_rb = rank_biserial(vals)
        mean, lo, hi = boot_ci_mean(vals)
        mean_task = bdf["rho_task"].mean()
        mean_inert = bdf["rho_inert"].mean()
        loo_ps = []
        for i in range(len(vals)):
            _, pi = wilcoxon_z(np.delete(vals, i))
            if np.isfinite(pi):
                loo_ps.append(pi)
        loo_worst = max(loo_ps) if loo_ps else np.nan
        vals_no_p03 = np.array([v for pa, v in zip(pats, vals) if pa != "Pat_03"])
        if len(vals_no_p03) >= 3:
            _, p_no_p03 = wilcoxon_z(vals_no_p03)
        else:
            p_no_p03 = np.nan
        band_rows.append({
            "band": band, "n": len(vals),
            "mean_delta": mean, "ci_lo": lo, "ci_hi": hi,
            "mean_task": mean_task, "mean_inert": mean_inert,
            "r_rb": r_rb, "z": z, "p": p,
            "loo_worst_p": loo_worst, "p_no_p03": p_no_p03,
            "patients": pats, "values": vals.tolist(),
        })
    ps = [r.get("p", np.nan) for r in band_rows]
    finite_mask = [np.isfinite(p) for p in ps]
    finite_p = [p for p, ok in zip(ps, finite_mask) if ok]
    qs_finite = bh_fdr(finite_p)
    j = 0
    for r, ok in zip(band_rows, finite_mask):
        r["q"] = qs_finite[j] if ok else np.nan
        if ok:
            j += 1
    for r in band_rows:
        if "mean_delta" not in r:
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | — | — | — | — | — | — | — | — | — | — |")
            continue
        sig = "★" if r["q"] < 0.05 else ""
        ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | {r['mean_delta']:+.4f} | "
           f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | {r['mean_task']:.3f} | "
           f"{r['mean_inert']:.3f} | {r['r_rb']:+.3f} | {r['z']:+.2f} | "
           f"{r['p']:.4f} | {r['q']:.4f}{sig} | {r['loo_worst_p']:.4f} | "
           f"{r['p_no_p03']:.4f} |")
        records.append({"hypothesis": "H2d", "band": r["band"], **{k: r[k] for k in
                        ("n", "mean_delta", "ci_lo", "ci_hi", "mean_task",
                         "mean_inert", "r_rb", "z", "p", "q", "loo_worst_p", "p_no_p03")}})
    ap("")
    ap("★ = q < 0.05 (FDR-BH across 6 bands).")
    ap("")

    # Cluster-based permutation across k per band
    ap("## Cluster-based permutation over k (sign-flip, 5,000 permutations)")
    ap("")
    ap("| band | cluster k-range | len | mass | p (cluster-corr) |")
    ap("|------|:----------------|----:|-----:|----------------:|")
    any_cluster = False
    for band in BRAIN_BANDS_NAMES:
        bdf = rawdf[rawdf["band"] == band]
        pats = sorted(bdf["patient"].unique())
        ks = sorted(bdf["k"].unique())
        mat = np.full((len(pats), len(ks)), np.nan)
        for ip, pat in enumerate(pats):
            sub = bdf[bdf["patient"] == pat].set_index("k")
            for ik, k in enumerate(ks):
                if k in sub.index:
                    mat[ip, ik] = sub.at[k, "delta_rho"]
        if np.isnan(mat).all(axis=1).any():
            keep = ~np.isnan(mat).all(axis=1)
            mat = mat[keep]
        ok_col = ~np.isnan(mat).any(axis=0)
        if ok_col.any():
            clusters = cluster_perm(mat[:, ok_col], np.array(ks)[ok_col])
        else:
            clusters = []
        for c in clusters:
            sig = "★" if c["p"] < 0.05 else ""
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | k={c['k_start']}–{c['k_end']} | "
               f"{c['len']} | {c['mass']:.1f} | {c['p']:.4f}{sig} |")
            any_cluster = True
    if not any_cluster:
        ap("| — | no supra-threshold clusters | — | — | — |")
    ap("")

    # Per-patient k-averaged table
    ap("## Per-patient k-averaged Δρ")
    ap("")
    pats_all = sorted(rawdf["patient"].unique())
    header = "| band | " + " | ".join(pats_all) + " |"
    ap(header)
    ap("|" + "|".join(["------"] + [":-----:"] * len(pats_all)) + "|")
    for band in BRAIN_BANDS_NAMES:
        bdf = pat_avg[pat_avg["band"] == band].set_index("patient")
        cells = [f"{bdf.at[p, 'delta_rho']:+.3f}" if p in bdf.index else "n/a"
                 for p in pats_all]
        ap(f"| {BRAIN_BAND_TEX_DICT[band]} | " + " | ".join(cells) + " |")
    ap("")

    # Absolute persistence rate
    ap("## Absolute persistence rate ρ_task (k-averaged)")
    ap("")
    ap("Higher ρ_task = a larger fraction of task-induced pairs survive into post-rest.")
    ap("")
    ap(header)
    ap("|" + "|".join(["------"] + [":-----:"] * len(pats_all)) + "|")
    for band in BRAIN_BANDS_NAMES:
        bdf = pat_avg[pat_avg["band"] == band].set_index("patient")
        cells = [f"{bdf.at[p, 'rho_task']:.2f}" if p in bdf.index else "n/a"
                 for p in pats_all]
        ap(f"| {BRAIN_BAND_TEX_DICT[band]} | " + " | ".join(cells) + " |")
    ap("")

    out = REPORTS_ROOT / "imcoh_vi" / "h2d_coactivation_persistence.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    pd.DataFrame(records).to_csv(
        REPORTS_ROOT / "imcoh_vi" / "h2d_coactivation_persistence.csv", index=False
    )


if __name__ == "__main__":
    main()
