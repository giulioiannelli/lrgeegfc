#!/usr/bin/env python3
"""Emit the authoritative H1-H4 VI(k) results markdown for the post-reset
|ImCoh| pipeline.

Reads the two CSVs already produced by compute_imcoh_vi.py:
  data/reports/imcoh_vi/vi_raw_profiles.csv
  data/reports/imcoh_vi/hypothesis_contrasts.csv

Computes H4 (frequency-gradient concordance) in place, reusing the same
Kendall's W + pairwise Spearman rho formulas as
scripts/05_multiscale/analyze_h4_gradient.py so the two agree by construction.

H2 (H2a + H2b) is the central hypothesis of the work and is reported first.

Output:
  .agents/reports/H1_H4_VI_RESULTS_POST_RESET.md

Run: python scripts/01_compute/report_h1h4_vi.py
"""
from __future__ import annotations

import datetime as _dt
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT

BANDS = BRAIN_BANDS_NAMES
PATIENTS = list(PATIENTS_LIST)  # all 6 patients, no exclusions
VI_DIR = REPORTS_ROOT / "imcoh_vi"
OUT = ROOT / ".agents" / "reports" / "H1_H4_VI_RESULTS_POST_RESET.md"


# ── helpers ────────────────────────────────────────────────────────────
def _git_hash() -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(ROOT), capture_output=True, text=True, timeout=2,
        )
        return r.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _cache_mtime() -> str:
    try:
        # Use Pat_02 alpha_rsPre as a representative file
        f = IMCOH_LRG_CACHE / "Pat_02" / "alpha_rsPre_lrg_imcoh-abs.npz"
        if f.exists():
            ts = _dt.datetime.fromtimestamp(f.stat().st_mtime)
            return ts.strftime("%Y-%m-%d %H:%M")
    except Exception:
        pass
    return "unknown"


def _band_tex(band: str) -> str:
    return BRAIN_BAND_TEX_DICT.get(band, band)


def _fmt_pct(x: float) -> str:
    if np.isnan(x):
        return "  n/a"
    return f"{100 * x:+5.1f}%"


# ── H1-H3: direction & unanimity per (band, k) from hypothesis_contrasts.csv
def summarize_h1_h3(cdf: pd.DataFrame, hyp: str) -> pd.DataFrame:
    """Per-band summary of a hypothesis over k-range.

    Returns rows: band, n_patients_supporting_mean, frac_k_unanimous_pos,
    frac_k_unanimous_neg, mean_contrast, sd_contrast.
    """
    hdf = cdf[cdf["hypothesis"] == hyp].copy()
    rows = []
    for band in BANDS:
        bdf = hdf[hdf["band"] == band]
        if bdf.empty:
            rows.append({"band": band})
            continue
        k_cells = 0
        unanimous_pos = 0
        unanimous_neg = 0
        per_patient_frac_pos = {p: np.nan for p in PATIENTS}

        # Unanimity per k: all patients PRESENT at that cell agree on sign.
        # A patient may be absent from a specific (band, k) cell when the
        # required phase(s) are unavailable (e.g. Pat_13 no rest_pre,
        # Pat_14 no task_test, Pat_10 rest↔task skipped due to n_nodes
        # mismatch). This is data-driven, not an exclusion.
        for k, grp in bdf.groupby("k"):
            if len(grp) < 2:
                continue
            k_cells += 1
            signs = grp["sign"].values
            if (signs > 0).all():
                unanimous_pos += 1
            elif (signs < 0).all():
                unanimous_neg += 1

        # Per-patient: fraction of k where this patient supports the hypothesis
        for p in PATIENTS:
            pdf = bdf[bdf["patient"] == p]
            if pdf.empty:
                continue
            per_patient_frac_pos[p] = float((pdf["contrast"] > 0).mean())

        rows.append({
            "band": band,
            "k_cells": k_cells,
            "unanimous_pos": unanimous_pos,
            "unanimous_neg": unanimous_neg,
            "frac_unanimous_pos": unanimous_pos / k_cells if k_cells else np.nan,
            "frac_unanimous_neg": unanimous_neg / k_cells if k_cells else np.nan,
            "mean_contrast": float(bdf["contrast"].mean()),
            "sd_contrast": float(bdf["contrast"].std()),
            **{f"frac_pos_{p}": per_patient_frac_pos[p] for p in PATIENTS},
        })
    return pd.DataFrame(rows)


# ── H4: Kendall's W + pairwise Spearman rho across patients' band rankings
def summarize_h4(vi_df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Per-k Kendall's W across patients' band-rankings by mean cross-pair VI.

    Also returns the k-averaged analysis mirroring analyze_h4_gradient.py:
    W over bands ranked by patient-mean cross-VI averaged across all k.
    """
    cross = vi_df[vi_df["pair_type"] == "cross"]
    out_rows = []
    k_values = sorted(cross["k"].unique())
    for k in k_values:
        kdf = cross[cross["k"] == k]
        # per (patient, band) mean VI over cross pairs
        scores = (
            kdf.groupby(["patient", "band"])["vi"].mean().unstack("band")
        )
        band_cols = [b for b in BANDS if b in scores.columns]
        scores = scores[band_cols]
        pat_rows = [p for p in PATIENTS if p in scores.index]
        if len(pat_rows) < 2 or len(band_cols) < 3:
            continue
        rank_mat = np.array([
            scores.loc[p].rank(ascending=False).values for p in pat_rows
        ])  # (n_pat, n_bands)
        n_judges, n_items = rank_mat.shape
        rank_sums = rank_mat.sum(axis=0)
        S = np.sum((rank_sums - rank_sums.mean()) ** 2)
        W = 12 * S / (n_judges ** 2 * (n_items ** 3 - n_items))
        chi2 = n_judges * (n_items - 1) * W
        p_val = 1 - stats.chi2.cdf(chi2, df=n_items - 1)

        # Pairwise Spearman rho
        rhos = []
        for i in range(len(pat_rows)):
            for j in range(i + 1, len(pat_rows)):
                rho, _ = stats.spearmanr(rank_mat[i], rank_mat[j])
                if not np.isnan(rho):
                    rhos.append(rho)
        out_rows.append({
            "k": k,
            "W": W,
            "chi2": chi2,
            "p": p_val,
            "mean_rho": float(np.mean(rhos)) if rhos else np.nan,
            "min_rho": float(np.min(rhos)) if rhos else np.nan,
            "max_rho": float(np.max(rhos)) if rhos else np.nan,
        })

    # k-averaged (global) W: rank bands by mean cross-VI averaged over all k
    global_scores = (
        cross.groupby(["patient", "band"])["vi"].mean().unstack("band")
    )
    band_cols = [b for b in BANDS if b in global_scores.columns]
    global_scores = global_scores[band_cols]
    pat_rows = [p for p in PATIENTS if p in global_scores.index]
    rank_mat = np.array([
        global_scores.loc[p].rank(ascending=False).values for p in pat_rows
    ])
    n_judges, n_items = rank_mat.shape
    rank_sums = rank_mat.sum(axis=0)
    S = np.sum((rank_sums - rank_sums.mean()) ** 2)
    W_global = 12 * S / (n_judges ** 2 * (n_items ** 3 - n_items))
    chi2_global = n_judges * (n_items - 1) * W_global
    p_global = 1 - stats.chi2.cdf(chi2_global, df=n_items - 1)

    # Per-patient top-3 bands (most reorganized)
    top3 = {}
    for p in pat_rows:
        order = global_scores.loc[p].sort_values(ascending=False).index.tolist()
        top3[p] = order[:3]

    # Consensus ranking (mean rank across patients)
    consensus = pd.DataFrame({
        "band": band_cols,
        "mean_rank": rank_mat.mean(axis=0),
        "median_rank": np.median(rank_mat, axis=0),
        "mean_cross_VI": global_scores.mean(axis=0).values,
    }).sort_values("mean_rank").reset_index(drop=True)

    summary = {
        "W_global": W_global,
        "chi2_global": chi2_global,
        "p_global": p_global,
        "patients_used": pat_rows,
        "top3": top3,
        "consensus": consensus,
    }
    return pd.DataFrame(out_rows), summary


# ── markdown rendering ────────────────────────────────────────────────
def render(vi_df: pd.DataFrame, cdf: pd.DataFrame) -> str:
    lines: list[str] = []
    ap = lines.append

    ap(f"# H1-H4 on VI(k): post-ImCoh-reset authoritative results")
    ap("")
    ap("**FC method:** `imcoh_abs` = <|ImCoh|>_f (Nolte-2004 signed ImCoh with "
       "abs-then-band-average per frequency bin; Jensen's inequality ordering "
       "respected).")
    ap(f"**LRG cache:** `data/cache/imcoh_lrg/` (representative mtime "
       f"{_cache_mtime()}).  ")
    ap(f"**Patients (N={len(PATIENTS)}):** {', '.join(PATIENTS)}. "
       f"Per-patient caveats: Pat_03 recorded at 1024 Hz (all others 2048 Hz) — "
       f"kept as a documented outlier / negative control. "
       f"Pat_10 resting phases ship 113 channels vs 116 in task phases, so "
       f"cross-phase VI (rest↔task pairs) is skipped for Pat_10. "
       f"Pat_13 is missing rest_pre (corrupt .mat) and Pat_14 is missing "
       f"task_test (corrupt .mat) — any contrast requiring those phases is "
       f"n/a for the affected patient.")
    ap(f"**Bands:** {', '.join(BANDS)}.  ")
    ap(f"**k range:** 2 to N/2 per dendrogram (typically k ≤ 50).  ")
    ap(f"**Git:** `{_git_hash()}`  ")
    ap(f"**Report generated:** {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}  ")
    ap(f"**Source CSVs:** `data/reports/imcoh_vi/vi_raw_profiles.csv`, "
       f"`data/reports/imcoh_vi/hypothesis_contrasts.csv`.")
    ap("")
    ap("Signed contrast convention: **positive = hypothesis supported**. "
       "Unanimity at a (band, k) cell means every patient that contributes "
       "data at that cell has the same sign.  The number of contributing "
       "patients per cell is not constant across the cohort — the "
       "per-patient columns below show n/a where a patient lacks the "
       "phase(s) required by the contrast.")
    ap("")
    ap("---")
    ap("")

    # ── H2 — CENTRAL ──────────────────────────────────────────────────
    ap("## H2 — Task trace (central hypothesis)")
    ap("")
    ap("H2 captures whether a cognitive task leaves a persistent trace on "
       "the multiscale brain organization.  It is split into two "
       "complementary contrasts, both evaluated at every dendrogram cut k:")
    ap("")
    ap("- **H2a (task trace / reset asymmetry):** "
       "`VI(rest_pre, rest_post) > VI(task_test, rest_post)` — post-task rest is "
       "closer to task_test than to pre-task rest.")
    ap("- **H2b (approach vs exit asymmetry):** "
       "`VI(rest_pre, task_test) > VI(task_test, rest_post)` — entering the task "
       "reorganizes more than exiting it.")
    ap("")

    for hyp in ("H2a", "H2b"):
        ap(f"### {hyp}")
        ap("")
        s = summarize_h1_h3(cdf, hyp)
        ap("| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |"
           "  " + " | ".join(f"{p}" for p in PATIENTS) + " |")
        ap("|------|--------:|--------:|--------:|------:|------:|-----:|---:|"
           + "".join("------:|" for _ in PATIENTS))
        for _, r in s.iterrows():
            if "k_cells" not in r or pd.isna(r.get("k_cells", np.nan)):
                ap(f"| {_band_tex(r['band'])} | — | — | — | — | — | — | — | "
                   + " | ".join("—" for _ in PATIENTS) + " |")
                continue
            per_pat = [_fmt_pct(r[f"frac_pos_{p}"]) for p in PATIENTS]
            ap(f"| {_band_tex(r['band'])} | {int(r['k_cells'])} "
               f"| {int(r['unanimous_pos'])} | {int(r['unanimous_neg'])} "
               f"| {_fmt_pct(r['frac_unanimous_pos'])} "
               f"| {_fmt_pct(r['frac_unanimous_neg'])} "
               f"| {r['mean_contrast']:+.4f} | {r['sd_contrast']:.4f} | "
               + " | ".join(per_pat) + " |")
        ap("")
        # Short prose verdict
        total_cells = int(s["k_cells"].fillna(0).sum())
        total_pos = int(s["unanimous_pos"].fillna(0).sum())
        total_neg = int(s["unanimous_neg"].fillna(0).sum())
        ap(f"**Overall {hyp}:** {total_pos}/{total_cells} (band, k) cells "
           f"unanimous positive; {total_neg}/{total_cells} unanimous "
           f"negative; {total_cells - total_pos - total_neg} split.")
        ap("")

    ap("---")
    ap("")

    # ── H1 ────────────────────────────────────────────────────────────
    ap("## H1 — Task stability")
    ap("")
    ap("**Statement:** `VI(task_learn, task_test) < mean(VI over other phase "
       "pairs)`. Positive contrast = the task-learn/task-test hierarchy is "
       "closer than the average phase pair at that (band, k).")
    ap("")
    s = summarize_h1_h3(cdf, "H1")
    ap("| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |"
       "  " + " | ".join(f"{p}" for p in PATIENTS) + " |")
    ap("|------|--------:|--------:|--------:|------:|------:|-----:|---:|"
       + "".join("------:|" for _ in PATIENTS))
    for _, r in s.iterrows():
        if "k_cells" not in r or pd.isna(r.get("k_cells", np.nan)):
            ap(f"| {_band_tex(r['band'])} | — | — | — | — | — | — | — | "
               + " | ".join("—" for _ in PATIENTS) + " |")
            continue
        per_pat = [_fmt_pct(r[f"frac_pos_{p}"]) for p in PATIENTS]
        ap(f"| {_band_tex(r['band'])} | {int(r['k_cells'])} "
           f"| {int(r['unanimous_pos'])} | {int(r['unanimous_neg'])} "
           f"| {_fmt_pct(r['frac_unanimous_pos'])} "
           f"| {_fmt_pct(r['frac_unanimous_neg'])} "
           f"| {r['mean_contrast']:+.4f} | {r['sd_contrast']:.4f} | "
           + " | ".join(per_pat) + " |")
    ap("")
    total_cells = int(s["k_cells"].fillna(0).sum())
    total_pos = int(s["unanimous_pos"].fillna(0).sum())
    total_neg = int(s["unanimous_neg"].fillna(0).sum())
    ap(f"**Overall H1:** {total_pos}/{total_cells} unanimous positive, "
       f"{total_neg}/{total_cells} unanimous negative.")
    ap("")
    ap("---")
    ap("")

    # ── H3 ────────────────────────────────────────────────────────────
    ap("## H3 — Within-modality similarity exceeds cross-modality")
    ap("")
    ap("**Statement:** `mean_VI(within) < mean_VI(cross)`, with within = "
       "`{(rest_pre,rest_post), (task_learn,task_test)}` and cross = the four "
       "rest-to-task pairs.")
    ap("")
    s = summarize_h1_h3(cdf, "H3")
    ap("| band | k cells | unan. + | unan. − | frac+ | frac− | mean | sd |"
       "  " + " | ".join(f"{p}" for p in PATIENTS) + " |")
    ap("|------|--------:|--------:|--------:|------:|------:|-----:|---:|"
       + "".join("------:|" for _ in PATIENTS))
    for _, r in s.iterrows():
        if "k_cells" not in r or pd.isna(r.get("k_cells", np.nan)):
            ap(f"| {_band_tex(r['band'])} | — | — | — | — | — | — | — | "
               + " | ".join("—" for _ in PATIENTS) + " |")
            continue
        per_pat = [_fmt_pct(r[f"frac_pos_{p}"]) for p in PATIENTS]
        ap(f"| {_band_tex(r['band'])} | {int(r['k_cells'])} "
           f"| {int(r['unanimous_pos'])} | {int(r['unanimous_neg'])} "
           f"| {_fmt_pct(r['frac_unanimous_pos'])} "
           f"| {_fmt_pct(r['frac_unanimous_neg'])} "
           f"| {r['mean_contrast']:+.4f} | {r['sd_contrast']:.4f} | "
           + " | ".join(per_pat) + " |")
    ap("")
    total_cells = int(s["k_cells"].fillna(0).sum())
    total_pos = int(s["unanimous_pos"].fillna(0).sum())
    total_neg = int(s["unanimous_neg"].fillna(0).sum())
    ap(f"**Overall H3:** {total_pos}/{total_cells} unanimous positive, "
       f"{total_neg}/{total_cells} unanimous negative.")
    ap("")
    ap("---")
    ap("")

    # ── H4 ────────────────────────────────────────────────────────────
    ap("## H4 — Frequency gradient of reorganization")
    ap("")
    ap("**Statement:** across patients, is there a consistent band ordering "
       "by reorganization strength (mean cross-pair VI)?  Tested via "
       "Kendall's coefficient of concordance W (over patients as judges, "
       "bands as items) and pairwise Spearman ρ of per-patient band rankings.")
    ap("")

    h4_per_k, h4_global = summarize_h4(vi_df)

    ap("### Global (k-averaged) concordance")
    ap("")
    ap(f"- **Kendall's W = {h4_global['W_global']:.3f}** "
       f"(n_judges = {len(h4_global['patients_used'])}, n_items = {len(BANDS)})")
    ap(f"- Friedman χ²({len(BANDS) - 1}) = {h4_global['chi2_global']:.2f}, "
       f"p = {h4_global['p_global']:.4f}")
    _w = h4_global["W_global"]
    label = ("weak" if _w < 0.3 else "moderate" if _w < 0.5 else
             "good" if _w < 0.7 else "strong")
    ap(f"- Interpretation: **{label} agreement** (W < 0.3 weak, < 0.5 "
       f"moderate, < 0.7 good, ≥ 0.7 strong).")
    ap("")

    ap("**Consensus ranking (lower mean_rank = more reorganized):**")
    ap("")
    ap("| band | mean rank | median rank | mean cross-VI |")
    ap("|------|----------:|------------:|--------------:|")
    for _, r in h4_global["consensus"].iterrows():
        ap(f"| {_band_tex(r['band'])} | {r['mean_rank']:.2f} | "
           f"{r['median_rank']:.1f} | {r['mean_cross_VI']:.4f} |")
    ap("")

    ap("**Per-patient top-3 most-reorganized bands:**")
    ap("")
    for p in h4_global["patients_used"]:
        t3 = h4_global["top3"][p]
        ap(f"- **{p}:** {' > '.join(_band_tex(b) for b in t3)}")
    ap("")

    ap("### k-resolved concordance")
    ap("")
    ap("Kendall's W and pairwise Spearman ρ at each dendrogram cut k:")
    ap("")
    ap("| k | W | χ² | p | mean ρ | min ρ | max ρ |")
    ap("|---|------:|------:|------:|-------:|------:|------:|")
    step = max(1, len(h4_per_k) // 20)  # at most ~20 rows
    for _, r in h4_per_k.iloc[::step].iterrows():
        ap(f"| {int(r['k'])} | {r['W']:.3f} | {r['chi2']:.2f} | "
           f"{r['p']:.3f} | {r['mean_rho']:+.3f} | "
           f"{r['min_rho']:+.3f} | {r['max_rho']:+.3f} |")
    ap("")

    # ── Pat_03 outlier note ───────────────────────────────────────────
    ap("---")
    ap("")
    ap("## Notes")
    ap("")
    ap("- **Pat_03** was recorded at 1024 Hz (all others at 2048 Hz); its "
       "MSC-era FC is ~3× denser than the rest of the cohort and unstable. "
       "The ImCoh pipeline uses `nperseg_for_fs(fs)` (2048 for Pat_03, "
       "4096 for the others) so the 2-second segmentation is preserved. "
       "Pat_03 remains in the tables as a negative control: deviations "
       "from the cohort are expected there and do not invalidate a "
       "group-level result.")
    ap("- **Pat_06** was completed with task_learn + task_test on 2026-04-22 "
       "and now contributes to all 4-phase analyses.")
    ap("- **Pat_10** resting recordings carry 113 channels vs 116 in the task "
       "phases and in `channel_labels.csv` (3 channels dropped by the vendor "
       "at resting). Cross-phase VI is skipped for Pat_10 on any pair that "
       "mixes rest and task (n_nodes mismatch).")
    ap("- **Pat_13** rest_pre.mat is corrupted (neither scipy v5/7 nor valid "
       "HDF5); Pat_14 task_test.mat opens but carries no `Data` key. Both are "
       "pending vendor re-supply.")
    ap("- **Same-probe bias** (CLAUDE.md invariant 4) is *partially* "
       "mitigated by ImCoh relative to MSC but not fully — any community-"
       "level interpretation should still verify with probe-debiased FC. "
       "VI at coarse k inherits probe geometry; the k-resolved tables "
       "above expose this directly (compare low-k vs high-k behaviour).")
    ap("")
    ap("## Reproducing")
    ap("")
    ap("```bash")
    ap("conda activate lapbrain")
    ap("python scripts/01_compute/compute_imcoh_vi.py --dry-run       # "
       "coverage check")
    ap("python scripts/01_compute/compute_imcoh_vi.py                 # "
       "compute CSVs")
    ap("python scripts/01_compute/report_h1h4_vi.py                   # "
       "regenerate this markdown")
    ap("```")
    ap("")
    ap("For figures on the same data, the multiscale scripts now default "
       "to `--fc-method imcoh_abs` (output under the corresponding "
       "`.../<fc_method>/` subfolder):")
    ap("")
    ap("```bash")
    ap("python scripts/05_multiscale/multiscale_all_hypotheses.py")
    ap("python scripts/05_multiscale/continuous_multiscale_h2.py        "
       "# H2 central")
    ap("python scripts/05_multiscale/definitive_multiscale_h2.py        "
       "# H2 central")
    ap("python scripts/05_multiscale/continuous_all_pairs.py")
    ap("python scripts/05_multiscale/analyze_h4_gradient.py")
    ap("```")
    ap("")
    return "\n".join(lines)


def main() -> None:
    vi_df = pd.read_csv(VI_DIR / "vi_raw_profiles.csv")
    cdf = pd.read_csv(VI_DIR / "hypothesis_contrasts.csv")
    md = render(vi_df, cdf)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(md)
    print(f"Wrote {OUT}  ({len(md):,} chars)")


if __name__ == "__main__":
    main()
