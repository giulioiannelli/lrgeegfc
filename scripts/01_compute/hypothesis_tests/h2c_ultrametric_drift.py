#!/usr/bin/env python3
"""H2c — ultrametric drift direction test.

Per (patient, band):
  Δ_rest(i,j) = D_rpost(i,j) − D_rpre(i,j)
  Δ_task(i,j) = D_ttest(i,j) − D_rpre(i,j)
  ρ = Spearman correlation of Δ_rest and Δ_task over upper-triangle pairs.

If the task leaves a trace in post-rest, ρ > 0: pairs whose ultrametric
distance shrank during the task also shrank during the rest, and vice-versa.

Cross-patient: one-sided Wilcoxon signed-rank on ρ per band, FDR across 6
bands. Reports effect size (rank-biserial), bootstrap CI, LOO sensitivity,
with / without Pat_03, and a control test using Δ_learn = D_tlearn − D_rpre
instead of Δ_task (replicates in a second task phase).

Outputs:
    data/reports/imcoh_vi/h2c_ultrametric_drift.md
    data/reports/imcoh_vi/h2c_ultrametric_drift.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


N_BOOT = 10_000


def upper_tri(mat: np.ndarray) -> np.ndarray:
    """Return upper-triangle (k=1) as a 1-D vector.

    LRGResult.ultrametric_matrix is stored in condensed / 1-D form
    (length N*(N-1)/2). If given a square matrix, extract the upper tri.
    """
    m = np.asarray(mat)
    if m.ndim == 1:
        return m
    n = m.shape[0]
    iu = np.triu_indices(n, k=1)
    return m[iu]


from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr  # canonical


def load_D(patient: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(patient, phase, band,
                            fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.ultrametric_matrix)


def drift_rho(D_rpre, D_rpost, D_target) -> float:
    """Spearman ρ between Δ_rest and Δ_target over the upper triangle."""
    if D_rpre is None or D_rpost is None or D_target is None:
        return np.nan
    u_pre = upper_tri(D_rpre)
    u_post = upper_tri(D_rpost)
    u_tgt = upper_tri(D_target)
    if not (u_pre.shape == u_post.shape == u_tgt.shape):
        # Pat_10: rest 113-ch, task 116-ch → different N(N-1)/2.
        return np.nan
    d_rest = u_post - u_pre
    d_target = u_tgt - u_pre
    if d_rest.std() == 0 or d_target.std() == 0:
        return np.nan
    rho, _ = stats.spearmanr(d_rest, d_target)
    return float(rho) if np.isfinite(rho) else np.nan


def main() -> None:
    patients = list(PATIENTS_LIST)
    rows = []
    for pat in patients:
        for band in BRAIN_BANDS_NAMES:
            D_rpre = load_D(pat, "rest_pre", band)
            D_rpost = load_D(pat, "rest_post", band)
            D_tlearn = load_D(pat, "task_learn", band)
            D_ttest = load_D(pat, "task_test", band)
            rows.append({
                "patient": pat, "band": band,
                "rho_task":  drift_rho(D_rpre, D_rpost, D_ttest),
                "rho_learn": drift_rho(D_rpre, D_rpost, D_tlearn),
                "has_rpre":  D_rpre is not None,
                "has_rpost": D_rpost is not None,
                "has_ttest": D_ttest is not None,
                "has_tlearn": D_tlearn is not None,
            })
    df = pd.DataFrame(rows)
    df.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2c_ultrametric_drift_raw.csv", index=False)

    # Cross-patient test per band
    lines: list[str] = []
    ap = lines.append
    ap("# H2c — ultrametric drift direction (post-rest drifts toward task)")
    ap("")
    ap("**Metric per (patient, band):** Spearman ρ between the rest-period change "
       "vector Δ_rest(i,j) = D_rpost − D_rpre and the task change vector "
       "Δ_task(i,j) = D_ttest − D_rpre, over upper-triangle pairs of the LRG "
       "ultrametric matrix.")
    ap("")
    ap("**Hypothesis (one-sided):** ρ > 0 — pairs that got closer in task tend "
       "to have gotten closer in post-rest (and vice versa).")
    ap("")
    ap("**Test:** one-sample Wilcoxon signed-rank on per-patient ρ per band, "
       "FDR-BH across 6 bands, 10k bootstrap CI on mean ρ, LOO sensitivity, "
       "with/without Pat_03. A control using the first task phase (task_learn) "
       "replaces D_ttest → D_tlearn to check replication across task phases.")
    ap("")
    ap("Pat_10 cross-phase has n_nodes mismatch (rest 113, task 116) → "
       "excluded from this test (ρ undefined).")
    ap("")

    records = []
    for target in ("task", "learn"):
        ap(f"## Target = {'task_test' if target == 'task' else 'task_learn'}")
        ap("")
        rho_col = f"rho_{target}"
        # Per-band stats
        summary_rows = []
        for band in BRAIN_BANDS_NAMES:
            vals = df[(df["band"] == band)][rho_col].dropna().to_numpy()
            pats = df[(df["band"] == band)].dropna(subset=[rho_col])["patient"].tolist()
            if len(vals) < 3:
                summary_rows.append({"band": band, "n": len(vals)})
                continue
            z, p = wilcoxon_z(vals)
            r_rb = rank_biserial(vals)
            mean, lo, hi = boot_ci_mean(vals)
            # LOO
            loo_ps = []
            for i in range(len(vals)):
                _, pi = wilcoxon_z(np.delete(vals, i))
                if np.isfinite(pi):
                    loo_ps.append(pi)
            loo_worst = max(loo_ps) if loo_ps else np.nan
            # Without Pat_03
            vals_no_p03 = [v for pa, v in zip(pats, vals) if pa != "Pat_03"]
            if len(vals_no_p03) >= 3:
                _, p_no_p03 = wilcoxon_z(np.asarray(vals_no_p03))
            else:
                p_no_p03 = np.nan
            summary_rows.append({
                "band": band, "n": len(vals),
                "mean_rho": mean, "ci_lo": lo, "ci_hi": hi,
                "median_rho": float(np.median(vals)),
                "r_rb": r_rb, "z": z, "p": p,
                "loo_worst_p": loo_worst,
                "p_no_p03": p_no_p03,
                "patients": pats,
                "values": vals.tolist(),
            })
        ps = [r.get("p", np.nan) for r in summary_rows]
        finite_mask = [np.isfinite(p) for p in ps]
        finite_p = [p for p, ok in zip(ps, finite_mask) if ok]
        fdr_q = bh_fdr(finite_p)
        j = 0
        qs = []
        for ok in finite_mask:
            if ok:
                qs.append(fdr_q[j]); j += 1
            else:
                qs.append(np.nan)
        for r, q in zip(summary_rows, qs):
            r["q"] = q

        ap("| band | n | mean ρ | 95% CI | median ρ | r_rb | z | p | q (BH) | LOO worst p | p no-Pat_03 |")
        ap("|------|--:|-------:|:------|---------:|-----:|--:|--:|-------:|------------:|------------:|")
        for r in summary_rows:
            if "mean_rho" not in r:
                ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | — | — | — | — | — | — | — | — | — |")
                continue
            sig = "★" if r["q"] < 0.05 else ""
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | "
               f"{r['mean_rho']:+.3f} | [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] | "
               f"{r['median_rho']:+.3f} | {r['r_rb']:+.3f} | "
               f"{r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{sig} | "
               f"{r['loo_worst_p']:.4f} | {r['p_no_p03']:.4f} |")
            records.append({"target": target, "hypothesis": "H2c", **{k: r[k] for k in
                            ("band", "n", "mean_rho", "ci_lo", "ci_hi",
                             "median_rho", "r_rb", "z", "p", "q",
                             "loo_worst_p", "p_no_p03")}})
        ap("")
        ap(f"★ = q < 0.05 (FDR-BH across 6 bands).")
        ap("")

        # Per-patient dump for this target
        ap("**Per-patient ρ (raw)**")
        ap("")
        header = "| band | " + " | ".join(sorted({p for r in summary_rows
                                                   for p in r.get("patients", [])})) + " |"
        ap(header)
        pats_all = sorted({p for r in summary_rows for p in r.get("patients", [])})
        ap("|" + "|".join(["------"] + [":-----:"] * len(pats_all)) + "|")
        for r in summary_rows:
            if "values" not in r:
                continue
            pat_val = dict(zip(r["patients"], r["values"]))
            cells = [f"{pat_val.get(p, np.nan):+.2f}" if p in pat_val else "n/a" for p in pats_all]
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | " + " | ".join(cells) + " |")
        ap("")

    pd.DataFrame(records).to_csv(
        REPORTS_ROOT / "imcoh_vi" / "h2c_ultrametric_drift.csv", index=False
    )
    out = REPORTS_ROOT / "imcoh_vi" / "h2c_ultrametric_drift.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
