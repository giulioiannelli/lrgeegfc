#!/usr/bin/env python3
"""H2a′ — directed partition-level memory via conditional entropy.

Per (patient, band, k):
    δH(k) = H(P_rpost(k) | P_ttest(k)) - H(P_rpost(k) | P_rpre(k))

Trace hypothesis (one-sided): δH(k) < 0 — knowing the task partition
reduces uncertainty about rest_post more than knowing rest_pre does.

This is the asymmetric half of the VI(A,B) = H(A|B)+H(B|A) decomposition,
recovering a genuinely directed partition-level test after the symmetric
H2a failed at n=9. See `.agents/guides/02_methods/H2_METRICS.md` §1.

Pre-registered pass: ≥3 of 6 bands with q < 0.05 FDR → H2a′ goes to
main Results as "partition-level directed memory". <3 bands → moves to
supplementary as a robustness check.

Reads:  LRG imcoh_abs caches (partitions reconstructed from linkage).
Writes: data/reports/imcoh_vi/h2a_prime_conditional_entropy.{md,csv}
        data/reports/imcoh_vi/h2a_prime_conditional_entropy_raw.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics import conditional_entropy
from lrg_eegfc.workflow.lrg import load_lrg_result

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr
from h2d_coactivation_persistence import cluster_perm  # reused: sign-flip cluster-perm over k


K_RANGE = list(range(2, 50))


# ─────────────────────────── compute ───────────────────────────


def load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.linkage_matrix)


def delta_h(Z_rpre, Z_ttest, Z_rpost, k: int) -> float | None:
    """Return δH(k) = H(rpost | ttest) - H(rpost | rpre); None on mismatch."""
    if Z_rpre is None or Z_ttest is None or Z_rpost is None:
        return None
    n = Z_rpre.shape[0] + 1
    if Z_ttest.shape[0] + 1 != n or Z_rpost.shape[0] + 1 != n:
        return None
    lab_pre  = fcluster(Z_rpre,  k, criterion="maxclust")
    lab_tt   = fcluster(Z_ttest, k, criterion="maxclust")
    lab_post = fcluster(Z_rpost, k, criterion="maxclust")
    h_post_given_tt  = conditional_entropy(lab_post, lab_tt)
    h_post_given_pre = conditional_entropy(lab_post, lab_pre)
    return h_post_given_tt - h_post_given_pre


def collect_raw() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_rpre  = load_Z(pat, "rest_pre",  band)
            Z_ttest = load_Z(pat, "task_test", band)
            Z_rpost = load_Z(pat, "rest_post", band)
            for k in K_RANGE:
                dh = delta_h(Z_rpre, Z_ttest, Z_rpost, k)
                if dh is None:
                    continue
                rows.append({"patient": pat, "band": band, "k": k, "delta_h": dh})
    return pd.DataFrame(rows)


# ─────────────────────────── across-patient stats ───────────────────────────


def per_band_stats(rawdf: pd.DataFrame) -> list[dict]:
    """Per-band one-sample Wilcoxon on k-averaged δH (one-sided < 0).

    We flip sign so that "success" means positive values; this way we can
    reuse `wilcoxon_z` (alternative='greater') without extra plumbing.
    """
    pat_avg = (rawdf.groupby(["patient", "band"])["delta_h"].mean()
               .reset_index().rename(columns={"delta_h": "mean_dh"}))
    out = []
    for band in BRAIN_BANDS_NAMES:
        bdf = pat_avg[pat_avg["band"] == band]
        pats = bdf["patient"].tolist()
        vals = bdf["mean_dh"].to_numpy()
        if len(vals) < 3:
            out.append({"band": band, "n": len(vals)})
            continue
        # H1: δH < 0 → test -δH > 0.
        z, p = wilcoxon_z(-vals)
        r_rb = rank_biserial(-vals)
        mean, lo, hi = boot_ci_mean(vals)
        loo_ps = []
        for i in range(len(vals)):
            _, pi = wilcoxon_z(-np.delete(vals, i))
            if np.isfinite(pi):
                loo_ps.append(pi)
        loo_worst = max(loo_ps) if loo_ps else np.nan
        vals_no_p03 = np.array([v for pa, v in zip(pats, vals) if pa != "Pat_03"])
        _, p_no_p03 = (wilcoxon_z(-vals_no_p03)
                       if len(vals_no_p03) >= 3 else (np.nan, np.nan))
        out.append({
            "band": band, "n": len(vals),
            "mean_dh": mean, "ci_lo": lo, "ci_hi": hi,
            "r_rb": r_rb, "z": z, "p": p,
            "loo_worst_p": loo_worst, "p_no_p03": p_no_p03,
            "patients": pats, "values": vals.tolist(),
        })
    ps = [r.get("p", np.nan) for r in out]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(out, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return out


def per_band_cluster(rawdf: pd.DataFrame) -> list[dict]:
    """Cluster-perm over k on (-δH) so that positive-z = trace-consistent."""
    clusters_all = []
    for band in BRAIN_BANDS_NAMES:
        bdf = rawdf[rawdf["band"] == band]
        pats = sorted(bdf["patient"].unique())
        ks = sorted(bdf["k"].unique())
        if not pats or not ks:
            continue
        mat = np.full((len(pats), len(ks)), np.nan)
        for ip, pat in enumerate(pats):
            sub = bdf[bdf["patient"] == pat].set_index("k")
            for ik, k in enumerate(ks):
                if k in sub.index:
                    mat[ip, ik] = -sub.at[k, "delta_h"]  # flip sign
        keep_row = ~np.isnan(mat).all(axis=1)
        mat = mat[keep_row]
        ok_col = ~np.isnan(mat).any(axis=0)
        if not ok_col.any():
            continue
        for c in cluster_perm(mat[:, ok_col], np.array(ks)[ok_col]):
            clusters_all.append({"band": band, **c})
    return clusters_all


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    rawdf = collect_raw()
    raw_path = REPORTS_ROOT / "imcoh_vi" / "h2a_prime_conditional_entropy_raw.csv"
    rawdf.to_csv(raw_path, index=False)
    print(f"wrote {raw_path}  ({len(rawdf):,} rows)")

    band_rows = per_band_stats(rawdf)
    clusters = per_band_cluster(rawdf)

    lines: list[str] = []
    ap = lines.append
    ap("# H2a′ — directed partition-level memory (conditional entropy)")
    ap("")
    ap("**Metric per (patient, band, k):**")
    ap("")
    ap("$\\delta H(k) = H(\\mathcal{P}^{rpost}_k \\mid \\mathcal{P}^{ttest}_k) "
       "- H(\\mathcal{P}^{rpost}_k \\mid \\mathcal{P}^{rpre}_k)$")
    ap("")
    ap("**Hypothesis (one-sided):** δH(k) < 0 — knowing the task partition "
       "reduces uncertainty about rest_post *more* than knowing rest_pre does.")
    ap("")
    ap("**Test:** one-sample Wilcoxon signed-rank on -δH (so 'greater' = "
       "trace-consistent), FDR-BH across 6 bands, 10k bootstrap CI on the "
       "raw δH mean, LOO sensitivity, Pat_03-drop sensitivity. Scale "
       "localization: sign-flip cluster permutation over k (reusing the H2d "
       "implementation).")
    ap("")
    ap(f"k range: {K_RANGE[0]}..{K_RANGE[-1]}. Patients with node-count "
       "mismatch (Pat_10) or missing ttest (Pat_14) are excluded automatically.")
    ap("")
    ap("## Patient-level k-averaged test")
    ap("")
    ap("| band | n | mean δH | 95% CI | r_rb | z | p | q (BH) | LOO worst p | p no-Pat_03 |")
    ap("|------|--:|--------:|:------|-----:|--:|--:|-------:|------------:|------------:|")
    for r in band_rows:
        tex = BRAIN_BAND_TEX_DICT[r["band"]]
        if "mean_dh" not in r:
            ap(f"| {tex} | {r['n']} | — | — | — | — | — | — | — | — |")
            continue
        sig = "★" if r["q"] < 0.05 else ""
        ap(f"| {tex} | {r['n']} | {r['mean_dh']:+.4f} | "
           f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | {r['r_rb']:+.3f} | "
           f"{r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{sig} | "
           f"{r['loo_worst_p']:.4f} | {r['p_no_p03']:.4f} |")
    ap("")
    ap("★ = q < 0.05 (FDR-BH across 6 bands). Pre-registered pass: "
       "≥3 bands with ★.")
    ap("")

    ap("## Cluster-based permutation over k (sign-flipped, 5,000 perms)")
    ap("")
    ap("| band | k-range | len | mass | p |")
    ap("|------|:--------|----:|-----:|--:|")
    if not clusters:
        ap("| — | no supra-threshold clusters | — | — | — |")
    else:
        for c in clusters:
            sig = "★" if c["p"] < 0.05 else ""
            ap(f"| {BRAIN_BAND_TEX_DICT[c['band']]} | "
               f"k={c['k_start']}–{c['k_end']} | {c['len']} | "
               f"{c['mass']:.1f} | {c['p']:.4f}{sig} |")
    ap("")

    records = [{"hypothesis": "H2a_prime", **{k: r[k] for k in
                ("band", "n", "mean_dh", "ci_lo", "ci_hi", "r_rb",
                 "z", "p", "q", "loo_worst_p", "p_no_p03")}}
               for r in band_rows if "mean_dh" in r]
    csv_path = REPORTS_ROOT / "imcoh_vi" / "h2a_prime_conditional_entropy.csv"
    pd.DataFrame(records).to_csv(csv_path, index=False)

    md_path = REPORTS_ROOT / "imcoh_vi" / "h2a_prime_conditional_entropy.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md_path}")
    print(f"wrote {csv_path}")

    n_pass = sum(1 for r in band_rows if "q" in r and np.isfinite(r["q"]) and r["q"] < 0.05)
    print(f"\nPre-registered pass criterion: {n_pass}/6 bands with q<0.05.")
    if n_pass >= 3:
        print("  → FULL PASS: promote H2a′ to main Results (three-level convergence).")
    elif n_pass >= 1:
        print("  → PARTIAL: report in supplementary as 'subset-band partition memory'.")
    else:
        print("  → FAIL: keep as honest-limitation paragraph; ρ/Δρ remain primary.")


if __name__ == "__main__":
    main()
