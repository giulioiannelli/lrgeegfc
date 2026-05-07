#!/usr/bin/env python3
"""Audit 43 -- Section-5 measure 09: cross-validation between LRG global probe
and per-leaf localization.

Bidirectional per-(patient, band) test:

(a) Enrichment direction. For pairs (i, j) where BOTH endpoints are
    calibrated trace-leaves (P_local_AND) under measure 12, compare the
    pair-level concordance score against pairs OUTSIDE that set:

        pair_score(i, j) = sign(dD_task) * sign(dD_rest)
                           * sqrt(|dD_task * dD_rest|)

        enrich(p, b)     = mean(pair_score | P_local) - mean(pair_score | ~P_local)

    Permutation null: shuffle the leaf is_trace flag 1000x and recompute.
    Two-sided p_perm = #{|enrich_perm| >= |enrich_real|} / N_perm.
    One-sided p_pos  = #{enrich_perm >= enrich_real} / N_perm.

(b) Concentration direction. The top-decile of pair_score (the strongest
    trace-direction pairs) should be enriched inside P_local:

        top10            = pairs with pair_score in top 10%
        conc(p, b)       = |top10 & P_local| / |top10|

    Same shuffle null. p_pos = #{conc_perm >= conc_real} / N_perm.

Both directions are also reported with P_local_OR (pair touches at least
one calibrated trace-leaf) as a sensitivity check.

Cohort verdict per band: paired Wilcoxon (real enrichment > 0) across the
10 patients, plus cohort fraction with permutation p < 0.05.

Outputs
-------
data/audit/lrg_cross_validation/per_patient_band.csv
data/audit/lrg_cross_validation/cohort_band_summary.csv
data/outputs/figures/section_5_lrg_trace/cross_validation/cohort_summary.pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

PAIR_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
LEAF_CAL = ROOT / "data" / "audit" / "per_leaf_rho_null" / "calibrated_trace_leaves.csv"
OUT_DIR = ROOT / "data" / "audit" / "lrg_cross_validation"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "cross_validation"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
N_PERM = 1000
RNG = np.random.default_rng(42)


def pair_score_vec(dT: np.ndarray, dR: np.ndarray) -> np.ndarray:
    """Magnitude-weighted concordance score per pair."""
    sigma = np.sign(dT) * np.sign(dR)
    mag = np.sqrt(np.abs(dT) * np.abs(dR))
    return sigma * mag


def cross_validate_one(pat: str, band: str, leaf_cal: pd.DataFrame) -> dict | None:
    npz_path = PAIR_DIR / f"{pat}_{band}.npz"
    if not npz_path.exists():
        return None
    d = np.load(npz_path)
    dT = d["dD_task"]; dR = d["dD_rest"]
    iu_i = d["iu_i"].astype(int); iu_j = d["iu_j"].astype(int)
    score = pair_score_vec(dT, dR)
    n_pairs = score.size

    sub = leaf_cal[(leaf_cal.patient == pat) & (leaf_cal.band == band)]
    if sub.empty:
        return None
    n_leaves = int(max(iu_i.max(), iu_j.max())) + 1
    in_trace = np.zeros(n_leaves, dtype=bool)
    for _, r in sub.iterrows():
        ell = int(r["leaf_id"])
        if ell < n_leaves and int(r["is_trace_p95"]) == 1:
            in_trace[ell] = True
    n_trace = int(in_trace.sum())
    if n_trace == 0:
        return {
            "patient": pat, "band": band, "n_pairs": n_pairs, "n_leaves": n_leaves,
            "n_trace_leaves": 0, "skip": True,
        }

    # Pair sets
    P_and = in_trace[iu_i] & in_trace[iu_j]
    P_or = in_trace[iu_i] | in_trace[iu_j]
    n_and = int(P_and.sum())
    n_or = int(P_or.sum())

    def enrich(P):
        if P.sum() < 2 or (~P).sum() < 2:
            return np.nan
        return float(score[P].mean() - score[~P].mean())

    enrich_and = enrich(P_and)
    enrich_or = enrich(P_or)

    # Top-decile concentration
    top10_thr = np.quantile(score, 0.90)
    top10 = score >= top10_thr
    n_top10 = int(top10.sum())

    def concentration(P):
        if n_top10 == 0:
            return np.nan
        return float((top10 & P).sum() / n_top10)

    conc_and = concentration(P_and)
    conc_or = concentration(P_or)

    # Permutation null: shuffle leaf is_trace flag, recompute
    e_and_null = np.empty(N_PERM)
    e_or_null = np.empty(N_PERM)
    c_and_null = np.empty(N_PERM)
    c_or_null = np.empty(N_PERM)
    for k in range(N_PERM):
        perm = RNG.permutation(n_leaves)
        in_trace_p = in_trace[perm]
        Pa = in_trace_p[iu_i] & in_trace_p[iu_j]
        Po = in_trace_p[iu_i] | in_trace_p[iu_j]
        e_and_null[k] = enrich(Pa) if not np.isnan(enrich(Pa)) else 0.0
        e_or_null[k] = enrich(Po) if not np.isnan(enrich(Po)) else 0.0
        c_and_null[k] = concentration(Pa)
        c_or_null[k] = concentration(Po)

    p_e_and = float((e_and_null >= enrich_and).mean()) if not np.isnan(enrich_and) else np.nan
    p_e_or = float((e_or_null >= enrich_or).mean()) if not np.isnan(enrich_or) else np.nan
    p_c_and = float((c_and_null >= conc_and).mean()) if not np.isnan(conc_and) else np.nan
    p_c_or = float((c_or_null >= conc_or).mean()) if not np.isnan(conc_or) else np.nan

    return {
        "patient": pat, "band": band,
        "n_pairs": n_pairs, "n_leaves": n_leaves, "n_trace_leaves": n_trace,
        "n_pairs_local_and": n_and, "n_pairs_local_or": n_or,
        "frac_local_and": n_and / n_pairs, "frac_local_or": n_or / n_pairs,
        "enrich_and": enrich_and, "enrich_or": enrich_or,
        "conc_top10_and": conc_and, "conc_top10_or": conc_or,
        "expected_conc_and": n_and / n_pairs,  # null expectation
        "expected_conc_or": n_or / n_pairs,
        "p_perm_enrich_and": p_e_and, "p_perm_enrich_or": p_e_or,
        "p_perm_conc_and": p_c_and, "p_perm_conc_or": p_c_or,
        "skip": False,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    leaf_cal = pd.read_csv(LEAF_CAL)
    rows = []
    for pat in PATIENTS:
        for band in BANDS:
            r = cross_validate_one(pat, band, leaf_cal)
            if r is None:
                continue
            rows.append(r)
            if r.get("skip"):
                print(f"[43] {pat} {band}: 0 trace-leaves, skipped")
                continue
            print(f"[43] {pat} {band}: n_trace={r['n_trace_leaves']:3d} | "
                  f"enrich_AND={r['enrich_and']:+.4f} (p={r['p_perm_enrich_and']:.3f}) | "
                  f"conc_AND={r['conc_top10_and']:.3f} (exp={r['expected_conc_and']:.3f}, "
                  f"p={r['p_perm_conc_and']:.3f})")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "per_patient_band.csv", index=False)

    # Cohort per-band summary
    band_rows = []
    df_v = df[~df.skip].copy()
    for band in BANDS:
        sub = df_v[df_v.band == band]
        if sub.empty:
            continue
        # Wilcoxon paired: real enrich vs zero
        e_and = sub["enrich_and"].dropna().to_numpy()
        e_or = sub["enrich_or"].dropna().to_numpy()
        # concentration: real conc vs expected (per-patient)
        c_and = sub["conc_top10_and"].dropna().to_numpy()
        e_c_and = sub.loc[~sub["conc_top10_and"].isna(), "expected_conc_and"].to_numpy()
        c_or = sub["conc_top10_or"].dropna().to_numpy()
        e_c_or = sub.loc[~sub["conc_top10_or"].isna(), "expected_conc_or"].to_numpy()

        def _wilcox(a, b=None, alt="greater"):
            if b is None:
                if a.size < 5:
                    return np.nan, np.nan
                try:
                    w = wilcoxon(a, alternative=alt, zero_method="wilcox")
                    return float(w.statistic), float(w.pvalue)
                except ValueError:
                    return np.nan, np.nan
            if a.size < 5:
                return np.nan, np.nan
            try:
                w = wilcoxon(a, b, alternative=alt, zero_method="wilcox")
                return float(w.statistic), float(w.pvalue)
            except ValueError:
                return np.nan, np.nan

        W_e_and, p_e_and = _wilcox(e_and, alt="greater")
        W_e_or, p_e_or = _wilcox(e_or, alt="greater")
        W_c_and, p_c_and = _wilcox(c_and, e_c_and, alt="greater")
        W_c_or, p_c_or = _wilcox(c_or, e_c_or, alt="greater")
        # cohort fraction with per-patient permutation p < 0.05
        n_pe_and = int((sub["p_perm_enrich_and"] < 0.05).sum())
        n_pe_or = int((sub["p_perm_enrich_or"] < 0.05).sum())
        n_pc_and = int((sub["p_perm_conc_and"] < 0.05).sum())
        n_pc_or = int((sub["p_perm_conc_or"] < 0.05).sum())
        band_rows.append({
            "band": band,
            "n_patients": int(len(sub)),
            "median_enrich_and": float(sub["enrich_and"].median(skipna=True)),
            "median_enrich_or": float(sub["enrich_or"].median(skipna=True)),
            "median_conc_top10_and": float(sub["conc_top10_and"].median(skipna=True)),
            "median_expected_conc_and": float(sub["expected_conc_and"].median()),
            "median_conc_top10_or": float(sub["conc_top10_or"].median(skipna=True)),
            "median_expected_conc_or": float(sub["expected_conc_or"].median()),
            "wilcoxon_p_enrich_and": p_e_and,
            "wilcoxon_p_enrich_or": p_e_or,
            "wilcoxon_p_conc_and_vs_exp": p_c_and,
            "wilcoxon_p_conc_or_vs_exp": p_c_or,
            "n_pat_perm_enrich_and_p05": n_pe_and,
            "n_pat_perm_enrich_or_p05": n_pe_or,
            "n_pat_perm_conc_and_p05": n_pc_and,
            "n_pat_perm_conc_or_p05": n_pc_or,
        })
    band_df = pd.DataFrame(band_rows)
    band_df.to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)
    print()
    print("Per-band cohort summary:")
    cols = ["band", "n_patients", "median_enrich_and", "median_conc_top10_and",
            "median_expected_conc_and", "wilcoxon_p_enrich_and",
            "wilcoxon_p_conc_and_vs_exp", "n_pat_perm_enrich_and_p05",
            "n_pat_perm_conc_and_p05"]
    print(band_df[cols].round(4).to_string(index=False))

    # ---------------- figure ------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(len(BANDS))

    # Left: per-band cohort enrichment box-plot (AND)
    ax = axes[0]
    data_and = [df_v[df_v.band == b]["enrich_and"].dropna().to_numpy() for b in BANDS]
    bp = ax.boxplot(data_and, positions=x, widths=0.5, patch_artist=True,
                    showfliers=False)
    for box in bp["boxes"]:
        box.set(facecolor="#9ec6e5", edgecolor="#1a4f73")
    for i, b in enumerate(BANDS):
        vals = data_and[i]
        ax.scatter(np.full_like(vals, i, dtype=float), vals, color="#1a4f73",
                   s=12, alpha=0.6, zorder=3)
    ax.axhline(0, color="0.4", lw=0.5, ls="--")
    ax.set_xticks(x); ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
    ax.set_ylabel(r"$\mathrm{enrich}(p,b)$ -- AND", fontsize=10)
    ax.set_title(r"pair-level concordance enrichment in trace-leaf pairs (AND)", fontsize=10)

    # Right: per-band concentration vs expected
    ax = axes[1]
    pos_real = x - 0.18
    pos_exp = x + 0.18
    data_c = [df_v[df_v.band == b]["conc_top10_and"].dropna().to_numpy() for b in BANDS]
    data_e = [df_v[df_v.band == b]["expected_conc_and"].dropna().to_numpy() for b in BANDS]
    bp_c = ax.boxplot(data_c, positions=pos_real, widths=0.32, patch_artist=True,
                      showfliers=False)
    bp_e = ax.boxplot(data_e, positions=pos_exp, widths=0.32, patch_artist=True,
                      showfliers=False)
    for box in bp_c["boxes"]:
        box.set(facecolor="#d62728", alpha=0.7, edgecolor="0.3")
    for box in bp_e["boxes"]:
        box.set(facecolor="#cccccc", edgecolor="0.3")
    ax.set_xticks(x); ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
    ax.set_ylabel(r"top-10\% pair concentration in $P_\mathrm{local}$ (AND)", fontsize=10)
    ax.set_title(r"red = real, grey = expected (chance) for AND-pair-set", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_summary.pdf")
    plt.close(fig)
    print(f"[43] wrote {FIG_DIR / 'cohort_summary.pdf'}")


if __name__ == "__main__":
    main()
