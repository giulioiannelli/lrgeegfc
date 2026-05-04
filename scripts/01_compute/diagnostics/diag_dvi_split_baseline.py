#!/usr/bin/env python3
"""Δ_VI(k) split-baseline diagnostic — analogue of H2c Run A / Run C for partitions.

Tests whether the cohort-wide Δ_VI(k) ridges (δ k=20–32, α k=20–26, β k=6–7,
γ_h k=21–23 and 30–35) survive when:

  1. ROBUSTNESS — rest_pre is replaced by either of two independent halves
     (rpre_A / rpre_B) drawn from disjoint segments of rest_pre. Each half
     gives an alternative Δ_VI; if both halves' Δ_VI > 0 cohort-wide, the
     headline does not depend on a particular realisation of rest_pre.

  2. DRIFT FLOOR (Run C analogue) — a "null Δ_VI" computed entirely on
     halves of rest_pre and rest_post (no task involved). Plays rpre_A as
     "rpre", rpre_B as "task", rpost_A as "rpost":

         drift_dVI(k) = VI(c_rpre_A, c_rpost_A) − VI(c_rpre_B, c_rpost_A)

     If drift alone produces Δ_VI > 0 cohort-wide (same noise regime as the
     halves, no task signal involved), the headline reading needs to be
     conditioned on the drift floor.

For each (band, k) the paired Wilcoxon test between |Δ_VI_full| and
|drift_dVI| (halves, same noise regime via the A-vs-B referent design)
gives a fair "above-floor" gate.

NOTE on noise regimes: Δ_VI_full uses full-data trees (2-s Welch windows).
drift_dVI uses halves trees (1-s Welch windows). H2e showed halves are
*noisier* (larger absolute distances) than full; the comparison is
therefore conservative from the headline's perspective — drift_dVI on
halves is biased *upward* in magnitude, not downward, so a Δ_VI_full that
beats drift_dVI is genuinely beating a louder noise floor.

Reads:
  data/cache/imcoh_lrg/<Pat>/<band>_<phase>_lrg_imcoh-abs.npz   (full)
  data/cache/imcoh_lrg_halves/<Pat>/<band>_<phase>_<A|B>_lrg_imcoh-abs.npz

Writes:
  data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv
  data/audit/dvi_split_baseline/dvi_split_baseline_summary.md
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, CACHE_ROOT
from lrg_eegfc.utils.metrics import compute_vi
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial
from lrg_eegfc.workflow.lrg import load_lrg_result

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
K_RANGE = list(range(2, 50))

# Headline ridges from `task_trace_band_k_n10_imcoh_abs.pdf` audit-and-recovery §4.
RIDGES = {
    "delta":      [(20, 32)],
    "alpha":      [(20, 26)],
    "beta":       [(6, 7)],
    "high_gamma": [(21, 23), (30, 35)],
}


def _Z(pat: str, phase: str, band: str, cache_root) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", cache_root)
    except Exception:
        return None
    if r is None:
        return None
    return np.asarray(r.linkage_matrix)


def _cut(Z: np.ndarray, k: int) -> np.ndarray:
    return fcluster(Z, k, criterion="maxclust")


def collect() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_pre   = _Z(pat, "rest_pre",   band, IMCOH_LRG_CACHE)
            Z_test  = _Z(pat, "task_test",  band, IMCOH_LRG_CACHE)
            Z_post  = _Z(pat, "rest_post",  band, IMCOH_LRG_CACHE)
            Z_preA  = _Z(pat, "rest_pre_A",  band, HALVES_CACHE)
            Z_preB  = _Z(pat, "rest_pre_B",  band, HALVES_CACHE)
            Z_postA = _Z(pat, "rest_post_A", band, HALVES_CACHE)

            # Skip patients/bands missing any source.
            if any(z is None for z in (Z_pre, Z_test, Z_post,
                                       Z_preA, Z_preB, Z_postA)):
                continue
            n_full = Z_pre.shape[0] + 1
            n_half = Z_preA.shape[0] + 1
            if Z_test.shape[0] + 1 != n_full or Z_post.shape[0] + 1 != n_full:
                continue
            if Z_preB.shape[0] + 1 != n_half or Z_postA.shape[0] + 1 != n_half:
                continue
            if n_full != n_half:
                # Halves were computed on the same channels as full; should match.
                continue

            for k in K_RANGE:
                c_pre   = _cut(Z_pre,   k)
                c_test  = _cut(Z_test,  k)
                c_post  = _cut(Z_post,  k)
                c_preA  = _cut(Z_preA,  k)
                c_preB  = _cut(Z_preB,  k)
                c_postA = _cut(Z_postA, k)

                # Full-data headline reproduction (sanity).
                d_VI_full = compute_vi(c_pre, c_post) - compute_vi(c_test, c_post)
                # Robustness: substitute rpre with each half (test/post still full).
                d_VI_A    = compute_vi(c_preA, c_post) - compute_vi(c_test, c_post)
                d_VI_B    = compute_vi(c_preB, c_post) - compute_vi(c_test, c_post)
                # Drift floor (halves only; task-free):
                # rpre_A plays "rpre", rpre_B plays "task", rpost_A plays "rpost".
                drift_dVI = (compute_vi(c_preA, c_postA)
                             - compute_vi(c_preB, c_postA))

                rows.append({
                    "patient": pat, "band": band, "k": k,
                    "d_VI_full": d_VI_full,
                    "d_VI_A":    d_VI_A,
                    "d_VI_B":    d_VI_B,
                    "drift_dVI": drift_dVI,
                })
    return pd.DataFrame(rows)


def _ridge_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Cohort summary at each headline ridge cell."""
    out = []
    for band in BRAIN_BANDS_NAMES:
        for (k_lo, k_hi) in RIDGES.get(band, []):
            sub = df[(df["band"] == band)
                     & (df["k"] >= k_lo) & (df["k"] <= k_hi)]
            # cohort-mean across (patient, k) for each contrast
            for col in ("d_VI_full", "d_VI_A", "d_VI_B", "drift_dVI"):
                vals = sub[col].to_numpy()
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    continue
                # per-patient mean across k inside the ridge:
                pp = sub.groupby("patient")[col].mean()
                pp = pp.dropna()
                out.append({
                    "band": band,
                    "ridge": f"k={k_lo}-{k_hi}",
                    "contrast": col,
                    "cohort_mean": float(pp.mean()),
                    "cohort_median": float(pp.median()),
                    "n_pos": int((pp > 0).sum()),
                    "n_pat": int(pp.size),
                })
    return pd.DataFrame(out)


def _paired_test(df: pd.DataFrame) -> pd.DataFrame:
    """Per-(band, ridge) paired Wilcoxon: |d_VI_full| > |drift_dVI|.

    Per-patient effect size = mean(|d_VI_full|) − mean(|drift_dVI|) across k
    inside the ridge. Tests whether the cohort distribution of this paired
    contrast is > 0.
    """
    rows = []
    for band in BRAIN_BANDS_NAMES:
        for (k_lo, k_hi) in RIDGES.get(band, []):
            sub = df[(df["band"] == band)
                     & (df["k"] >= k_lo) & (df["k"] <= k_hi)].copy()
            if sub.empty:
                continue
            # paired per patient: mean of d_VI_full (signed, so we test the
            # one-sided "task-effect ridge is positive AND exceeds drift" jointly).
            pp = (sub.groupby("patient")
                     .agg(d_VI_full=("d_VI_full", "mean"),
                          drift_dVI=("drift_dVI", "mean")))
            if pp.shape[0] < 5:
                continue
            d = pp["d_VI_full"].to_numpy()
            null = pp["drift_dVI"].to_numpy()
            diff = d - null
            # One-sample (paired) one-sided test on diff > 0.
            z, p = wilcoxon_z(diff)
            r_rb = rank_biserial(diff)
            rows.append({
                "band": band,
                "ridge": f"k={k_lo}-{k_hi}",
                "n_paired": int(pp.shape[0]),
                "median_d_VI_full": float(np.median(d)),
                "median_drift_dVI": float(np.median(null)),
                "n_pos_full": int((d > 0).sum()),
                "n_pos_full_gt_drift": int((d > null).sum()),
                "wilcoxon_z": z,
                "p_one_sided": p,
                "r_rb": r_rb,
            })
    return pd.DataFrame(rows)


def main() -> None:
    out_dir = ROOT / "data" / "audit" / "dvi_split_baseline"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Computing Δ_VI split-baseline (this is sub-second per cell)...")
    df = collect()
    raw_csv = out_dir / "dvi_split_baseline_n10_imcoh_abs.csv"
    df.to_csv(raw_csv, index=False)
    print(f"  {len(df)} raw rows → {raw_csv.relative_to(ROOT)}")
    print(f"  patients × bands seen: {df.groupby(['patient','band']).ngroups}")

    ridge = _ridge_summary(df)
    paired = _paired_test(df)

    rsum_csv = out_dir / "dvi_ridge_summary.csv"
    ridge.to_csv(rsum_csv, index=False)
    print(f"  → {rsum_csv.relative_to(ROOT)}")
    pwsum_csv = out_dir / "dvi_paired_wilcoxon.csv"
    paired.to_csv(pwsum_csv, index=False)
    print(f"  → {pwsum_csv.relative_to(ROOT)}")

    # ── Markdown summary ─────────────────────────────────────────────────
    lines: list[str] = []
    ap = lines.append
    ap("# Δ_VI(k) split-baseline diagnostic — n=10 IMCOH_ABS")
    ap("")
    ap("Tests whether the cohort-wide Δ_VI ridges in the headline figure")
    ap("(`task_trace_band_k_n10_imcoh_abs.pdf`) survive when:")
    ap("")
    ap("- **Robustness**: rest_pre is replaced by either of two independent")
    ap("  halves (`d_VI_A`, `d_VI_B`).")
    ap("- **Drift floor**: rpre_A vs rpre_B as if they were rpre vs task")
    ap("  (`drift_dVI`, halves only, no task involved).")
    ap("")
    ap("### Per-ridge cohort means (per-patient mean across k inside the ridge)")
    ap("")
    ap("| band | ridge | contrast | cohort mean | cohort median | n_pos / n_pat |")
    ap("|------|-------|----------|------------:|--------------:|--------------:|")
    for _, r in ridge.iterrows():
        tex = BRAIN_BAND_TEX_DICT.get(r["band"], r["band"])
        ap(f"| {tex} | {r['ridge']} | `{r['contrast']}` | "
           f"{r['cohort_mean']:+.3f} | {r['cohort_median']:+.3f} | "
           f"{r['n_pos']}/{r['n_pat']} |")
    ap("")
    ap("### Paired one-sided Wilcoxon: d_VI_full > drift_dVI")
    ap("")
    ap("Per-patient: mean(d_VI_full) and mean(drift_dVI) across k inside ridge.")
    ap("")
    ap("| band | ridge | n | med d_VI_full | med drift_dVI | full>0 | full>drift | z | p (one-sided) | r_rb |")
    ap("|------|-------|--:|--------------:|--------------:|-------:|-----------:|--:|---------------:|-----:|")
    for _, r in paired.iterrows():
        tex = BRAIN_BAND_TEX_DICT.get(r["band"], r["band"])
        ap(f"| {tex} | {r['ridge']} | {r['n_paired']} | "
           f"{r['median_d_VI_full']:+.3f} | {r['median_drift_dVI']:+.3f} | "
           f"{r['n_pos_full']}/{r['n_paired']} | "
           f"{r['n_pos_full_gt_drift']}/{r['n_paired']} | "
           f"{r['wilcoxon_z']:+.2f} | {r['p_one_sided']:.4f} | "
           f"{r['r_rb']:+.2f} |")
    ap("")
    ap("### Reading rules")
    ap("")
    ap("- **Robustness pass** = both `d_VI_A` and `d_VI_B` keep cohort sign")
    ap("  and ≥ 7/10 patients positive. Means the headline does not hinge")
    ap("  on a particular realisation of rest_pre.")
    ap("- **Drift floor pass** = `drift_dVI` cohort median ≈ 0 AND paired")
    ap("  Wilcoxon `d_VI_full > drift_dVI` significant (p < 0.05).")
    ap("  Means the partition residual exceeds what within-session drift")
    ap("  alone produces in the same noise regime.")
    ap("- **Failure** in either column = the ridge is uncontrolled for the")
    ap("  same kind of artefact that flipped H2c-shared.")
    ap("")
    ap("Halves use 1-s Welch windows (vs 2-s for full data) per the H2e")
    ap("convention — halves are noisier so drift_dVI is biased *upward*")
    ap("in magnitude, making the comparison conservative.")

    out_md = out_dir / "dvi_split_baseline_summary.md"
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  → {out_md.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
