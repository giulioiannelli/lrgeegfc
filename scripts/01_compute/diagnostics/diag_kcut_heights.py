#!/usr/bin/env python3
"""Cut-height inventory — what does "k=28" actually mean across patients/phases?

For each (patient, band, phase, k), record:
  - h(k)      = absolute merge height at which fcluster(Z, k, maxclust)
                transitions from k+1 → k clusters
  - dmax(Z)   = root height
  - h_rel(k)  = h(k) / dmax(Z)
  - n_eff(k)  = Simpson effective number of clusters at k (Σ p_i²)⁻¹

Then summarise at the headline-ridge k-windows:
  - per-patient mean h_rel inside the ridge
  - cohort spread of h_rel (min, mean ± std, max)
  - per-phase comparison within patient

This is the data that lets us judge the user's concern — "dmax varies
across patients and across phases, so comparing on k is tricky." If
h_rel is tightly concentrated cohort-wide at a headline k-window, then
"k=28" maps to a similar fractional depth in every tree → the cohort
claim survives translation to h_rel. If h_rel is wide, integer-k is
the wrong axis.

Reads:  data/cache/imcoh_lrg/<Pat>/<band>_<phase>_lrg_imcoh-abs.npz
Writes: data/audit/dvi_split_baseline/kcut_heights_n10_imcoh_abs.csv
        data/audit/dvi_split_baseline/kcut_heights_summary.md
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.utils.metrics.tree import dmax_from_Z, simpson_neff
from lrg_eegfc.workflow.lrg import load_lrg_result

K_RANGE = list(range(2, 50))
PHASES = ("rest_pre", "task_test", "rest_post")

RIDGES = {
    "delta":      [(20, 32)],
    "alpha":      [(20, 26)],
    "beta":       [(6, 7)],
    "high_gamma": [(21, 23), (30, 35)],
}


def _Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return None if r is None else np.asarray(r.linkage_matrix)


def _h_at_k(Z: np.ndarray, k: int) -> float:
    """Cut height for `fcluster(Z, k, maxclust)`. Equals Z[-k, 2] for k ≥ 2."""
    n_leaves = Z.shape[0] + 1
    if k >= n_leaves:
        return 0.0
    if k < 1:
        return float("nan")
    return float(Z[-k, 2])


def collect() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            for phase in PHASES:
                Z = _Z(pat, phase, band)
                if Z is None:
                    continue
                n = Z.shape[0] + 1
                dmax = dmax_from_Z(Z)
                for k in K_RANGE:
                    h_k = _h_at_k(Z, k)
                    h_rel = h_k / dmax if dmax > 0 else float("nan")
                    labels = fcluster(Z, k, criterion="maxclust")
                    n_eff = simpson_neff(labels)
                    rows.append({
                        "patient": pat, "band": band, "phase": phase, "k": k,
                        "n_leaves": n, "h_k": h_k, "dmax": dmax,
                        "h_rel": h_rel, "n_eff": n_eff,
                    })
    return pd.DataFrame(rows)


def main() -> None:
    out_dir = ROOT / "data" / "audit" / "dvi_split_baseline"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = collect()
    csv = out_dir / "kcut_heights_n10_imcoh_abs.csv"
    df.to_csv(csv, index=False)
    print(f"{len(df)} rows → {csv.relative_to(ROOT)}")

    # ── dmax per patient × band × phase ─────────────────────────────
    dmax_df = (df.groupby(["patient", "band", "phase"])["dmax"]
                 .first().reset_index())
    dmax_pivot = dmax_df.pivot_table(
        index=["patient", "band"], columns="phase", values="dmax")

    # ── Per-ridge h_rel summary ─────────────────────────────────────
    lines: list[str] = []
    ap = lines.append
    ap("# Cut-height inventory — what does k mean across patients × phases?")
    ap("")
    ap("`h_k = Z[-k, 2]` is the merge height at which the partition transitions")
    ap("from k+1 → k clusters. `h_rel = h_k / dmax(Z)` is the fractional depth.")
    ap("`n_eff(k) = (Σ p_i²)⁻¹` is the Simpson effective number of clusters.")
    ap("")
    ap("## 1. dmax variation cohort-wide")
    ap("")
    ap("Per band, summary of dmax across (patient × phase). dmax is the")
    ap("absolute root height of each LRG tree under `imcoh_abs`.")
    ap("")
    ap("| band | min dmax | median dmax | max dmax | max/min |")
    ap("|------|---------:|------------:|---------:|--------:|")
    for band in BRAIN_BANDS_NAMES:
        bd = dmax_df[dmax_df["band"] == band]["dmax"].to_numpy()
        bd = bd[np.isfinite(bd)]
        if bd.size == 0:
            continue
        ap(f"| {BRAIN_BAND_TEX_DICT[band]} | {bd.min():.4f} | "
           f"{np.median(bd):.4f} | {bd.max():.4f} | {bd.max()/bd.min():.2f} |")
    ap("")

    ap("## 2. Within-patient phase variation of dmax")
    ap("")
    ap("For each (patient, band), how much does dmax differ across phases?")
    ap("Reported as `dmax_test / dmax_pre` and `dmax_post / dmax_pre`.")
    ap("Cohort percentiles per band.")
    ap("")
    ap("| band | n | median test/pre | IQR test/pre | median post/pre | IQR post/pre |")
    ap("|------|--:|----------------:|:-------------|----------------:|:-------------|")
    for band in BRAIN_BANDS_NAMES:
        sub = dmax_pivot.xs(band, level="band")
        if not all(p in sub.columns for p in PHASES):
            continue
        valid = sub.dropna()
        if valid.empty:
            continue
        r_test = valid["task_test"] / valid["rest_pre"]
        r_post = valid["rest_post"] / valid["rest_pre"]
        ap(f"| {BRAIN_BAND_TEX_DICT[band]} | {len(valid)} | "
           f"{r_test.median():.3f} | "
           f"[{r_test.quantile(0.25):.3f}, {r_test.quantile(0.75):.3f}] | "
           f"{r_post.median():.3f} | "
           f"[{r_post.quantile(0.25):.3f}, {r_post.quantile(0.75):.3f}] |")
    ap("")

    ap("## 3. h_rel at headline-ridge k cuts")
    ap("")
    ap("For each ridge cell, per-patient mean h_rel (across the k window")
    ap("AND across phases rpre/test/post). Cohort min/median/max gives the")
    ap("interpretation spread: a tight range means k tracks h_rel well; a")
    ap("wide range means 'k=N' interprets differently in different patients.")
    ap("")
    ap("| band | ridge | n_pat × n_phase | h_rel min | median | max | max/min |")
    ap("|------|-------|---------------:|----------:|-------:|----:|--------:|")
    for band in BRAIN_BANDS_NAMES:
        for (k_lo, k_hi) in RIDGES.get(band, []):
            sub = df[(df["band"] == band)
                     & (df["k"] >= k_lo) & (df["k"] <= k_hi)]
            if sub.empty:
                continue
            # Per (patient, phase) mean h_rel across the ridge k values:
            pp = (sub.groupby(["patient", "phase"])["h_rel"].mean()
                     .dropna())
            if pp.empty:
                continue
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | k={k_lo}-{k_hi} | "
               f"{len(pp)} | {pp.min():.3f} | {pp.median():.3f} | "
               f"{pp.max():.3f} | {pp.max()/max(pp.min(), 1e-9):.1f} |")
    ap("")

    ap("## 4. h_rel by phase at the δ k=20–32 headline ridge")
    ap("")
    ap("Per-patient, per-phase mean h_rel inside the δ k=20–32 ridge. If")
    ap("'k=20–32' translates to *similar* h_rel across phases within a")
    ap("patient, the integer-k cut is approximately phase-invariant.")
    ap("")
    sub = df[(df["band"] == "delta") & (df["k"] >= 20) & (df["k"] <= 32)]
    if not sub.empty:
        pp = (sub.groupby(["patient", "phase"])["h_rel"].mean()
                 .unstack("phase"))
        ap("| patient | rest_pre | task_test | rest_post | "
           "test/pre | post/pre |")
        ap("|---|---:|---:|---:|---:|---:|")
        for pat in pp.index:
            r = pp.loc[pat]
            r_t = r.get("task_test", float("nan")) / r.get("rest_pre", float("nan"))
            r_p = r.get("rest_post", float("nan")) / r.get("rest_pre", float("nan"))
            ap(f"| {pat} | {r.get('rest_pre', float('nan')):.3f} | "
               f"{r.get('task_test', float('nan')):.3f} | "
               f"{r.get('rest_post', float('nan')):.3f} | "
               f"{r_t:.3f} | {r_p:.3f} |")
        ap("")

    ap("## 5. n_eff at headline-ridge k cuts")
    ap("")
    ap("If n_eff(k) ≈ k cohort-wide, partitions are balanced. n_eff << k")
    ap("means a giant cluster + many singletons (k-artefact regime).")
    ap("")
    ap("| band | ridge | n_pat × n_phase | n_eff min | median | max |")
    ap("|------|-------|---------------:|----------:|-------:|----:|")
    for band in BRAIN_BANDS_NAMES:
        for (k_lo, k_hi) in RIDGES.get(band, []):
            sub = df[(df["band"] == band)
                     & (df["k"] >= k_lo) & (df["k"] <= k_hi)]
            if sub.empty:
                continue
            pp = (sub.groupby(["patient", "phase"])["n_eff"].mean()
                     .dropna())
            if pp.empty:
                continue
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | k={k_lo}-{k_hi} | "
               f"{len(pp)} | {pp.min():.2f} | {pp.median():.2f} | "
               f"{pp.max():.2f} |")
    ap("")

    out_md = out_dir / "kcut_heights_summary.md"
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  → {out_md.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
