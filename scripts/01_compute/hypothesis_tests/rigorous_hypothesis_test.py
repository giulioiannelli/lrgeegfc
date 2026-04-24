#!/usr/bin/env python3
"""Rigorous cross-patient hypothesis tests on VI(k) contrasts.

Inputs: data/reports/imcoh_vi/hypothesis_contrasts.csv

For each (hypothesis, band):

(A) Patient-level pooled test — k-averaged contrast, one per patient.
    - Primary: one-sided Wilcoxon signed-rank (H0: median ≤ 0, H1: median > 0).
    - Effect size: rank-biserial correlation r_rb.
    - 95% bootstrap CI on the patient-mean contrast (B=10,000 resamples).
    - 24 hypothesis×band cells → FDR-BH correction within each hypothesis.

(B) Cluster-based permutation over k (Maris & Oostenveld 2007).
    - At each k: one-sample Wilcoxon W_z = sign(median) · |z(W)| where
      z(W) is the standardised Wilcoxon statistic.
    - Primary threshold: |W_z| > 1.96 (alpha = 0.05, two-sided threshold
      applied one-sided for positive clusters).
    - Cluster mass = sum of W_z over contiguous supra-threshold k's.
    - Null: flip each patient's sign independently (sign-flip permutation,
      respects within-patient dependence, 5,000 permutations), recompute
      max cluster mass.  Cluster p = fraction of null permutations whose
      max cluster mass ≥ observed cluster mass.

(C) Sensitivity analyses:
    - Leave-one-out (LOO) on (A) — is any single patient driving the result?
    - With vs without Pat_03 (1024 Hz outlier / negative control).

(D) H4 frequency-gradient concordance:
    - Kendall's W on per-patient band-rankings (mean cross-phase VI).
    - Two-sided Friedman test on band effect.

Outputs:
    data/reports/imcoh_vi/rigorous_tests.md
    data/reports/imcoh_vi/rigorous_tests.csv (machine-readable table)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


HYPS = ["H1", "H2a", "H2b", "H3"]
HYP_LABEL = {
    "H1":  "H1 — task stability (TL↔TT < others)",
    "H2a": "H2a — task trace (rpre↔rpost > ttest↔rpost)",
    "H2b": "H2b — approach > exit (rpre↔ttest > ttest↔rpost)",
    "H3":  "H3 — within < cross",
}

N_BOOT = 10_000
N_PERM = 5_000
ALPHA = 0.05
CLUSTER_Z_THRESH = 1.96  # one-sided 0.025 tail


# ─────────────────────────── helpers ───────────────────────────

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr, cluster_stats  # canonical


def per_k_wilcoxon_z(mat: np.ndarray) -> np.ndarray:
    """Vectorised one-sample Wilcoxon z per column (no ties correction — OK
    with continuous VI contrasts where exact zeros essentially never occur).

    For (n_patients × n_k) contrast matrix, returns z per k.
    """
    abs_ranks = np.apply_along_axis(stats.rankdata, 0, np.abs(mat))
    signs = np.sign(mat)
    w_plus = (abs_ranks * (signs > 0)).sum(axis=0)
    n = (signs != 0).sum(axis=0)
    mu = n * (n + 1) / 4.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sigma > 0, (w_plus - mu) / sigma, np.nan)
    return z



def cluster_perm_test(mat: np.ndarray, k_values: np.ndarray,
                      n_perm: int = N_PERM,
                      thresh: float = CLUSTER_Z_THRESH,
                      seed: int = 0) -> list[dict]:
    """Sign-flip cluster-based permutation test across k.

    mat: (n_patients × n_k) per-patient contrast at each k.
    Returns list of observed clusters with p-value.

    Vectorised: |ranks| are computed once, then each permutation is a
    vector operation that flips signs and recomputes z across k.
    """
    rng = np.random.default_rng(seed)

    # Pre-compute ranks once — they only depend on |mat|, invariant under sign flips.
    abs_ranks = np.apply_along_axis(stats.rankdata, 0, np.abs(mat))
    obs_signs = np.sign(mat)

    n = (obs_signs != 0).sum(axis=0)
    mu = n * (n + 1) / 4.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)

    def _z_from_signs(signs: np.ndarray) -> np.ndarray:
        w_plus = (abs_ranks * (signs > 0)).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(sigma > 0, (w_plus - mu) / sigma, np.nan)

    z_obs = _z_from_signs(obs_signs)
    obs_clusters = cluster_stats(z_obs, thresh)
    if not obs_clusters:
        return []

    n_pat = mat.shape[0]
    flips = rng.choice([-1.0, 1.0], size=(n_perm, n_pat))
    max_null = np.zeros(n_perm)
    for i in range(n_perm):
        signs = obs_signs * flips[i, :, None]
        z_null = _z_from_signs(signs)
        clusters = cluster_stats(z_null, thresh)
        max_null[i] = max((c[2] for c in clusters), default=0.0)

    out = []
    for start, end, mass in obs_clusters:
        p = float((max_null >= mass).mean())
        out.append({
            "k_start": int(k_values[start]),
            "k_end":   int(k_values[end]),
            "len":     int(end - start + 1),
            "mass":    float(mass),
            "p":       p,
        })
    return out


# ─────────────────────────── main ───────────────────────────

def main() -> None:
    src = REPORTS_ROOT / "imcoh_vi" / "hypothesis_contrasts.csv"
    df = pd.read_csv(src)

    # Build (hyp, band, patient) → contrast-per-k vector.
    per_patient_k: dict[tuple[str, str, str], dict[int, float]] = {}
    for (hyp, band, pat), grp in df.groupby(["hypothesis", "band", "patient"]):
        per_patient_k[(hyp, band, pat)] = dict(zip(grp["k"], grp["contrast"]))

    records = []
    lines: list[str] = []
    ap = lines.append

    ap("# Rigorous cross-patient hypothesis tests on VI(k) contrasts")
    ap("")
    ap(f"Source: `data/reports/imcoh_vi/hypothesis_contrasts.csv`  ")
    ap(f"Patients per hypothesis: **n varies by hypothesis** (contributing = has all required phases).  ")
    ap(f"FDR correction: Benjamini-Hochberg within each hypothesis (m=6 bands).  ")
    ap(f"Cluster-based permutation: {N_PERM} sign-flip permutations, primary |z| > {CLUSTER_Z_THRESH} (α=0.05 one-sided).  ")
    ap(f"Bootstrap CI: {N_BOOT} resamples on the patient-mean k-averaged contrast.  ")
    ap("")
    ap("Conventions: **positive contrast = hypothesis supported**.  "
       "Effect size = rank-biserial r_rb ∈ [−1, +1].  r_rb = 1 means every patient's "
       "k-averaged contrast is positive and larger (in absolute value) than any negative "
       "contrast that exists; r_rb = 0 is the null.")
    ap("")

    for hyp in HYPS:
        ap(f"## {HYP_LABEL[hyp]}")
        ap("")
        # Pass 1: gather per-patient k-averaged contrast + per-k matrix.
        band_rows = []
        band_clusters: dict[str, list] = {}
        all_patients = sorted({p for (h, b, p) in per_patient_k if h == hyp})
        for band in BRAIN_BANDS_NAMES:
            # per-patient k-averaged contrast
            xs = []
            for pat in all_patients:
                d = per_patient_k.get((hyp, band, pat))
                if d:
                    xs.append((pat, np.mean(list(d.values()))))
            if len(xs) < 3:
                band_rows.append({"band": band, "n": len(xs)})
                continue
            pats, vals = zip(*xs)
            vals = np.asarray(vals)
            # Primary test
            z, p = wilcoxon_z(vals)
            r_rb = rank_biserial(vals)
            mean, lo, hi = boot_ci_mean(vals)
            # LOO sensitivity: min p across leaving one patient out
            loo_ps = []
            loo_max_shift = 0.0
            for i in range(len(vals)):
                sub = np.delete(vals, i)
                _, p_i = wilcoxon_z(sub)
                if np.isfinite(p_i):
                    loo_ps.append(p_i)
                mean_i = np.mean(sub)
                loo_max_shift = max(loo_max_shift, abs(mean_i - mean))
            loo_worst_p = max(loo_ps) if loo_ps else np.nan
            # Pat_03 sensitivity
            no_p03_vals = [v for p_, v in xs if p_ != "Pat_03"]
            if len(no_p03_vals) >= 3:
                _, p_no_p03 = wilcoxon_z(np.array(no_p03_vals))
            else:
                p_no_p03 = np.nan
            band_rows.append({
                "band": band,
                "n": len(vals),
                "mean": mean, "ci_lo": lo, "ci_hi": hi,
                "r_rb": r_rb,
                "z": z, "p": p,
                "loo_worst_p": loo_worst_p,
                "loo_max_shift": loo_max_shift,
                "p_no_p03": p_no_p03,
                "patients": list(pats),
            })

            # (B) Cluster-based permutation across k
            k_values = sorted({k for d in (per_patient_k.get((hyp, band, pp), {})
                                            for pp in pats) for k in d})
            k_values = np.array(sorted(k_values))
            mat = np.full((len(pats), len(k_values)), np.nan)
            for ip, pat in enumerate(pats):
                d = per_patient_k.get((hyp, band, pat), {})
                for ik, k in enumerate(k_values):
                    if k in d:
                        mat[ip, ik] = d[k]
            # Drop rows with all-NaN
            if np.isnan(mat).all(axis=1).any():
                keep = ~np.isnan(mat).all(axis=1)
                mat = mat[keep]
            # Use only columns where all remaining patients have data (required for sign test)
            ok_col = ~np.isnan(mat).any(axis=0)
            if ok_col.any():
                clusters = cluster_perm_test(mat[:, ok_col], k_values[ok_col])
            else:
                clusters = []
            band_clusters[band] = clusters

        # FDR within hypothesis for primary test
        pvals = [r.get("p", np.nan) for r in band_rows]
        finite_mask = [np.isfinite(p) for p in pvals]
        finite_p = [p for p, ok in zip(pvals, finite_mask) if ok]
        finite_q = bh_fdr(finite_p)
        qs = []
        j = 0
        for ok, p in zip(finite_mask, pvals):
            if ok:
                qs.append(finite_q[j]); j += 1
            else:
                qs.append(np.nan)
        for r, q in zip(band_rows, qs):
            r["q"] = q

        # ── Primary table
        ap("### Patient-level k-averaged test (Wilcoxon signed-rank, one-sided)")
        ap("")
        ap("| band | n | mean contrast | 95% CI | r_rb | z | p (raw) | q (BH) | LOO worst p | p without Pat_03 |")
        ap("|------|--:|--------------:|:------|-----:|--:|--------:|-------:|------------:|-----------------:|")
        for r in band_rows:
            if "mean" not in r:
                ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | — | — | — | — | — | — | — | — |")
                continue
            sig = "★" if r["q"] < ALPHA else ""
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | {r['mean']:+.4f} | [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | "
               f"{r['r_rb']:+.3f} | {r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{sig} | "
               f"{r['loo_worst_p']:.4f} | {r['p_no_p03']:.4f} |")
            records.append({
                "hypothesis": hyp, "band": r["band"], "n": r["n"],
                "mean": r["mean"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
                "r_rb": r["r_rb"], "z": r["z"], "p": r["p"], "q": r["q"],
                "loo_worst_p": r["loo_worst_p"], "loo_max_shift": r["loo_max_shift"],
                "p_no_pat03": r["p_no_p03"],
                "patients": ",".join(r["patients"]),
            })
        ap("")
        ap(f"★ = q < {ALPHA} (FDR-BH within hypothesis).")
        ap("")

        # ── Cluster table
        ap("### Cluster-based permutation over k (sign-flip, 5,000 permutations)")
        ap("")
        ap("Primary: one-sample Wilcoxon z at each k, threshold |z| > 1.96.  "
           "Cluster mass = sum of z over contiguous supra-threshold k.  "
           "Cluster p = fraction of sign-flipped permutations with max cluster mass ≥ observed.")
        ap("")
        ap("| band | cluster k-range | len | mass | p (cluster-corr) |")
        ap("|------|:----------------|----:|-----:|----------------:|")
        any_cluster = False
        for band in BRAIN_BANDS_NAMES:
            cs = band_clusters.get(band, [])
            for c in cs:
                sig = "★" if c["p"] < ALPHA else ""
                ap(f"| {BRAIN_BAND_TEX_DICT[band]} | k={c['k_start']}–{c['k_end']} | "
                   f"{c['len']} | {c['mass']:.1f} | {c['p']:.4f}{sig} |")
                any_cluster = True
        if not any_cluster:
            ap("| — | no supra-threshold clusters | — | — | — |")
        ap("")
        ap(f"★ = p < {ALPHA} (cluster-corrected, two-sided-threshold one-sided test).")
        ap("")

    # H4: Kendall's W on per-patient band rankings of mean cross-phase VI.
    ap("## H4 — Frequency gradient (Kendall's W)")
    ap("")
    vi = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "vi_raw_profiles.csv")
    # mean cross-phase VI per (patient, band): average over cross pairs × all k
    # cross pairs classified in classify_pair
    from lrg_eegfc.config.const import classify_pair
    vi["ptype"] = [classify_pair(a, b) for a, b in zip(vi["phase_a"], vi["phase_b"])]
    cross = vi[vi["ptype"] == "cross"]
    table = cross.groupby(["patient", "band"])["vi"].mean().unstack("band")
    if table.shape[0] < 3:
        ap(f"n patients with cross-phase data = {table.shape[0]} (too few for Kendall).")
    else:
        # Rank bands per patient (smaller VI = rank 1)
        ranks = table[BRAIN_BANDS_NAMES].rank(axis=1, method="average")
        n, k = ranks.shape
        sum_ranks = ranks.sum(axis=0).to_numpy()
        S = np.sum((sum_ranks - sum_ranks.mean()) ** 2)
        W = 12 * S / (n ** 2 * (k ** 3 - k))
        chi2 = n * (k - 1) * W
        p_w = 1 - stats.chi2.cdf(chi2, df=k - 1)
        ap(f"n patients = {n}, n bands = {k}.")
        ap(f"Kendall's W = **{W:.3f}**  (χ²({k-1}) = {chi2:.2f}, p = {p_w:.4f})")
        ap("")
        ap("| band | mean rank | mean cross-VI |")
        ap("|------|----------:|--------------:|")
        for band in BRAIN_BANDS_NAMES:
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | {ranks[band].mean():.2f} | {table[band].mean():.3f} |")
    ap("")

    # Final write-out
    out_md = REPORTS_ROOT / "imcoh_vi" / "rigorous_tests.md"
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out_md}  ({out_md.stat().st_size:,} bytes)")

    out_csv = REPORTS_ROOT / "imcoh_vi" / "rigorous_tests.csv"
    pd.DataFrame(records).to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
