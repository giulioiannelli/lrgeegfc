#!/usr/bin/env python3
"""H2e — session-drift noise floor for ρ (H2c) and Δρ (H2d).

Per patient × band, compare the cross-phase ρ (already in
`h2c_ultrametric_drift_raw.csv`) against a DRIFT-ONLY ρ computed on
split halves of the two resting recordings, and contextualize both
against a WITHIN-PHASE reliability ceiling:

    ρ_cross         = Spearman(D^rpost−D^rpre,  D^ttest−D^rpre)           [H2c]
    ρ_null_drift    = Spearman(D^rpre_B−D^rpre_A,  D^rpost_B−D^rpost_A)   [within-session drift only]
    ρ_within(rpre)  = Spearman(D^rpre_A,  D^rpre_B)                       [reliability ceiling]
    ρ_within(rpost) = Spearman(D^rpost_A, D^rpost_B)                      [reliability ceiling]

Primary test (per band, paired Wilcoxon, one-sided):
    ρ_cross  >  ρ_null_drift    — FDR-BH m=6
    Δρ_cross >  Δρ_null_drift   — FDR-BH m=6  (secondary; reuses H2d infrastructure)

Pre-registered pass: ≥4 of 6 bands with q<0.05 on the ρ paired test →
"task-driven correlation exceeds session-drift noise floor, not an
artefact of session order."

Halving T halves raw freq resolution → we use
`nperseg = nperseg_for_fs(fs) // 2` (1-s windows instead of 2-s). D
matrices on halved data are noisier than full-duration; H2e is a
statistical contrast, not a magnitude comparison.

Reads:  raw time series (via canonical loader + PATIENT_CHANNEL_DROP);
        `data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv` for ρ_cross;
        `data/reports/imcoh_vi/h2d_persistence_raw.csv` for Δρ_cross.
Writes: `data/reports/imcoh_vi/h2e_split_half.{md,csv}` + `_raw.csv`;
        split-half LRG caches under `data/cache/imcoh_lrg_halves/`.
"""
from __future__ import annotations

import argparse
import gc
import resource
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT,
    FS_OVERRIDES, PATIENTS_LIST, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, REPORTS_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.lrg import compute_lrg_analysis

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr
from _fc_split_half import compute_imcoh_abs_halves
from h2c_ultrametric_drift import upper_tri
from h2d_coactivation_persistence import compute_persistence, K_RANGE


HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
REST_PHASES = ("rest_pre", "rest_post")
DEFAULT_MEM_GB = 10.0  # hard cap; raises MemoryError if exceeded


def _set_mem_cap(gb: float) -> None:
    """Set RLIMIT_AS so the process fails fast instead of thrashing swap."""
    soft = int(gb * 1024 ** 3)
    try:
        resource.setrlimit(resource.RLIMIT_AS, (soft, soft))
    except (ValueError, OSError):
        pass  # best-effort; some kernels disallow lowering the cap


# ─────────────────────────── halves pipeline ───────────────────────────


def _load_X(pat: str, phase: str) -> tuple[np.ndarray, float] | None:
    try:
        X = load_timeseries(pat, phase, SEEG_DATAPATH)
    except (FileNotFoundError, OSError):
        return None
    if X is None:
        return None
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, 2048.0)
    return X, fs


def ensure_halves_lrg(pat: str, verbose: bool = False
                      ) -> dict[tuple[str, str, str], object]:
    """Compute (or load from cache) LRG on each (phase, half, band).

    Returns {(phase, half_tag, band): LRGResult}. Skips phases whose raw
    data is unavailable. Uses a 1-s Welch window (halved nperseg) for
    split-half spectral estimation. Aggressively frees FC matrices after
    the LRG result is cached to keep peak RSS bounded.
    """
    results: dict[tuple[str, str, str], object] = {}
    for phase in REST_PHASES:
        loaded = _load_X(pat, phase)
        if loaded is None:
            print(f"  {pat}/{phase}: timeseries missing, skip")
            continue
        X, fs = loaded
        nperseg_half = max(256, nperseg_for_fs(fs) // 2)
        if verbose:
            print(f"  {pat}/{phase}: X.shape={X.shape}, nperseg_half={nperseg_half}")
        halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
        del X
        gc.collect()
        for (band, tag), A in halves_fc.items():
            synth_phase = f"{phase}_{tag}"
            r = compute_lrg_analysis(
                A, pat, synth_phase, band, fc_method="imcoh_abs",
                cache_root=HALVES_CACHE, use_cache=True,
            )
            results[(phase, tag, band)] = r
        halves_fc.clear()
        del halves_fc
        gc.collect()
    return results


# ─────────────────────────── ρ/Δρ-null computation ───────────────────────────


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size == 0 or x.std() == 0 or y.std() == 0:
        return np.nan
    r, _ = stats.spearmanr(x, y)
    return float(r) if np.isfinite(r) else np.nan


def rho_quartet(halves: dict) -> dict[str, float]:
    """Compute ρ_null_drift, ρ_within_rpre, ρ_within_rpost per band.

    Keys of `halves`: (phase, tag, band) → LRGResult. Returns a flat
    dict keyed by (band, metric_name).
    """
    out: dict[tuple[str, str], float] = {}
    for band in BRAIN_BANDS_NAMES:
        try:
            d_pre_A = upper_tri(halves[("rest_pre",  "A", band)].ultrametric_matrix)
            d_pre_B = upper_tri(halves[("rest_pre",  "B", band)].ultrametric_matrix)
            d_post_A = upper_tri(halves[("rest_post", "A", band)].ultrametric_matrix)
            d_post_B = upper_tri(halves[("rest_post", "B", band)].ultrametric_matrix)
        except KeyError:
            continue
        if not (d_pre_A.shape == d_pre_B.shape == d_post_A.shape == d_post_B.shape):
            continue
        out[(band, "rho_null_drift")]   = _spearman(d_pre_B - d_pre_A,
                                                    d_post_B - d_post_A)
        out[(band, "rho_within_rpre")]  = _spearman(d_pre_A,  d_pre_B)
        out[(band, "rho_within_rpost")] = _spearman(d_post_A, d_post_B)
    return {f"{band}__{name}": v for (band, name), v in out.items()}


def delta_rho_null(halves: dict) -> dict[tuple[str, int], float]:
    """Δρ_null per (band, k) with triple (rpre_A, rpre_B, rpost_A).

    Substitution: play rpre_A as "rpre", rpre_B as "ttest", rpost_A as
    "rpost" — if drift alone produces Δρ > 0, this null captures it.
    Uses linkage_matrix from the halved LRG caches.
    """
    out: dict[tuple[str, int], float] = {}
    for band in BRAIN_BANDS_NAMES:
        try:
            Z_pre_A = halves[("rest_pre",  "A", band)].linkage_matrix
            Z_pre_B = halves[("rest_pre",  "B", band)].linkage_matrix
            Z_post_A = halves[("rest_post", "A", band)].linkage_matrix
        except KeyError:
            continue
        for k in K_RANGE:
            res = compute_persistence(Z_pre_A, Z_pre_B, Z_post_A, k)
            if res is not None:
                out[(band, k)] = res[0] - res[1]
    return out


# ─────────────────────────── raw assembly ───────────────────────────


def _load_cross_rho() -> pd.DataFrame:
    src = REPORTS_ROOT / "imcoh_vi" / "h2c_ultrametric_drift_raw.csv"
    return pd.read_csv(src)


def _load_cross_drho() -> pd.DataFrame:
    src = REPORTS_ROOT / "imcoh_vi" / "h2d_persistence_raw.csv"
    return pd.read_csv(src)


def collect_raw(verbose: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (rho_rows, drho_rows) with cross vs null aligned by (pat, band[, k])."""
    cross_rho  = _load_cross_rho()
    cross_drho = _load_cross_drho()
    rho_rows, drho_rows = [], []

    for pat in PATIENTS_LIST:
        halves = ensure_halves_lrg(pat, verbose=verbose)
        if not halves:
            print(f"  {pat}: no halves computed, skip")
            continue

        rho_null = rho_quartet(halves)
        drho_null = delta_rho_null(halves)

        # ρ rows
        sub_rho = cross_rho[cross_rho["patient"] == pat]
        for _, row in sub_rho.iterrows():
            band = row["band"]
            rho_rows.append({
                "patient": pat, "band": band,
                "rho_cross":       row["rho_task"],
                "rho_null_drift":  rho_null.get(f"{band}__rho_null_drift", np.nan),
                "rho_within_rpre": rho_null.get(f"{band}__rho_within_rpre", np.nan),
                "rho_within_rpost":rho_null.get(f"{band}__rho_within_rpost", np.nan),
            })

        # Δρ rows (per k)
        sub_drho = cross_drho[cross_drho["patient"] == pat]
        for _, row in sub_drho.iterrows():
            band, k = row["band"], int(row["k"])
            null_val = drho_null.get((band, k), np.nan)
            drho_rows.append({
                "patient": pat, "band": band, "k": k,
                "drho_cross": row["delta_rho"],
                "drho_null":  null_val,
            })

        # Release LRG result objects for this patient before moving on.
        halves.clear()
        del halves, rho_null, drho_null
        gc.collect()

    return pd.DataFrame(rho_rows), pd.DataFrame(drho_rows)


# ─────────────────────────── across-patient stats ───────────────────────────


def _paired_stats(diffs: np.ndarray, patients: list[str]) -> dict:
    """One-sample Wilcoxon on paired differences (cross - null), one-sided >0."""
    d = diffs[np.isfinite(diffs)]
    if len(d) < 3:
        return {"n": len(d)}
    z, p = wilcoxon_z(d)
    r_rb = rank_biserial(d)
    mean, lo, hi = boot_ci_mean(d)
    vals_no_p03 = np.array([v for pa, v in zip(patients, diffs)
                             if pa != "Pat_03" and np.isfinite(v)])
    _, p_no_p03 = (wilcoxon_z(vals_no_p03)
                   if len(vals_no_p03) >= 3 else (np.nan, np.nan))
    return {"n": len(d), "mean_diff": mean, "ci_lo": lo, "ci_hi": hi,
            "r_rb": r_rb, "z": z, "p": p, "p_no_p03": p_no_p03}


def rho_band_stats(rho_df: pd.DataFrame) -> list[dict]:
    rows = []
    for band in BRAIN_BANDS_NAMES:
        bdf = rho_df[rho_df["band"] == band]
        pats = bdf["patient"].tolist()
        diff = (bdf["rho_cross"] - bdf["rho_null_drift"]).to_numpy()
        mean_cross = np.nanmean(bdf["rho_cross"].to_numpy())
        mean_null  = np.nanmean(bdf["rho_null_drift"].to_numpy())
        mean_ceil  = np.nanmean(
            np.concatenate([bdf["rho_within_rpre"].to_numpy(),
                            bdf["rho_within_rpost"].to_numpy()])
        )
        s = _paired_stats(diff, pats)
        rows.append({"band": band, "mean_cross": mean_cross,
                     "mean_null": mean_null, "mean_ceiling": mean_ceil, **s})
    ps = [r.get("p", np.nan) for r in rows]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(rows, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return rows


def drho_band_stats(drho_df: pd.DataFrame) -> list[dict]:
    """K-averaged Δρ_cross − Δρ_null per patient×band, Wilcoxon across patients."""
    per = (drho_df.groupby(["patient", "band"])[["drho_cross", "drho_null"]]
           .mean().reset_index())
    rows = []
    for band in BRAIN_BANDS_NAMES:
        bdf = per[per["band"] == band]
        pats = bdf["patient"].tolist()
        diff = (bdf["drho_cross"] - bdf["drho_null"]).to_numpy()
        s = _paired_stats(diff, pats)
        rows.append({"band": band,
                     "mean_cross": float(np.nanmean(bdf["drho_cross"])),
                     "mean_null":  float(np.nanmean(bdf["drho_null"])),
                     **s})
    ps = [r.get("p", np.nan) for r in rows]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(rows, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return rows


# ─────────────────────────── report ───────────────────────────


def _rho_table(rho_rows, tex_dict):
    lines = ["| band | n | mean ρ_cross | mean ρ_null | mean ceiling | Δ cross−null | 95% CI | r_rb | z | p | q (BH) |",
             "|------|--:|-------------:|------------:|-------------:|-------------:|:------|-----:|--:|--:|-------:|"]
    for r in rho_rows:
        tex = tex_dict[r["band"]]
        if "mean_diff" not in r:
            lines.append(f"| {tex} | {r['n']} | — | — | — | — | — | — | — | — | — |")
            continue
        sig = "★" if r["q"] < 0.05 else ""
        lines.append(
            f"| {tex} | {r['n']} | {r['mean_cross']:+.3f} | {r['mean_null']:+.3f} | "
            f"{r['mean_ceiling']:+.3f} | {r['mean_diff']:+.3f} | "
            f"[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] | {r['r_rb']:+.3f} | "
            f"{r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{sig} |"
        )
    return lines


def _drho_table(rows, tex_dict):
    lines = ["| band | n | mean Δρ_cross | mean Δρ_null | Δ cross−null | 95% CI | r_rb | z | p | q (BH) |",
             "|------|--:|--------------:|-------------:|-------------:|:------|-----:|--:|--:|-------:|"]
    for r in rows:
        tex = tex_dict[r["band"]]
        if "mean_diff" not in r:
            lines.append(f"| {tex} | {r['n']} | — | — | — | — | — | — | — | — |")
            continue
        sig = "★" if r["q"] < 0.05 else ""
        lines.append(
            f"| {tex} | {r['n']} | {r['mean_cross']:+.4f} | {r['mean_null']:+.4f} | "
            f"{r['mean_diff']:+.4f} | [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | "
            f"{r['r_rb']:+.3f} | {r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{sig} |"
        )
    return lines


def write_report(rho_rows, drho_rows) -> None:
    out_dir = REPORTS_ROOT / "imcoh_vi"
    md = out_dir / "h2e_split_half.md"
    csv = out_dir / "h2e_split_half.csv"

    lines: list[str] = []
    ap = lines.append
    ap("# H2e — session-drift noise floor for ρ (and Δρ)")
    ap("")
    ap("Cross-phase ρ (H2c) and Δρ (H2d) are compared against a within-session "
       "drift null computed on split halves of the two resting recordings, and "
       "bracketed by a within-phase reliability ceiling. Split halves use a "
       "halved Welch `nperseg` (1-s windows) to preserve segment count; D "
       "matrices on half-duration data are noisier than full, so the null is "
       "modestly inflated — the test is a statistical contrast, not a "
       "magnitude comparison.")
    ap("")
    ap("Pat_14 can still be included (no task phase needed).")
    ap("")

    ap("## ρ paired test (primary: cross > null_drift)")
    ap("")
    ap("* `ρ_cross` = H2c Spearman(D^rpost−D^rpre, D^ttest−D^rpre).")
    ap("* `ρ_null_drift` = Spearman(D^rpre_B−D^rpre_A, D^rpost_B−D^rpost_A).")
    ap("* `ceiling` = mean of Spearman(D^phase_A, D^phase_B) over the two "
       "resting phases (split-half reliability of the ultrametric matrix).")
    ap("")
    lines.extend(_rho_table(rho_rows, BRAIN_BAND_TEX_DICT))
    ap("")
    ap("★ = q < 0.05 (FDR-BH across 6 bands). Pre-registered pass: ≥4 bands ★.")
    ap("")

    ap("## Δρ paired test (secondary: cross > null_drift)")
    ap("")
    ap("k-averaged Δρ_cross (H2d) vs Δρ_null computed on the triple "
       "(rpre_A, rpre_B, rpost_A) — drift-only analog of the H2d triple.")
    ap("")
    lines.extend(_drho_table(drho_rows, BRAIN_BAND_TEX_DICT))
    ap("")

    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")

    records = [{"metric": "rho", **{k: r.get(k) for k in
                ("band", "n", "mean_cross", "mean_null", "mean_ceiling",
                 "mean_diff", "ci_lo", "ci_hi", "r_rb", "z", "p", "q", "p_no_p03")}}
               for r in rho_rows]
    records += [{"metric": "drho", **{k: r.get(k) for k in
                 ("band", "n", "mean_cross", "mean_null",
                  "mean_diff", "ci_lo", "ci_hi", "r_rb", "z", "p", "q", "p_no_p03")}}
                for r in drho_rows]
    pd.DataFrame(records).to_csv(csv, index=False)
    print(f"wrote {csv}")


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mem-gb", type=float, default=DEFAULT_MEM_GB,
                   help="hard RSS cap (RLIMIT_AS); process fails fast on overrun.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    _set_mem_cap(args.mem_gb)
    HALVES_CACHE.mkdir(parents=True, exist_ok=True)
    rho_df, drho_df = collect_raw(verbose=args.verbose)
    rho_df.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2e_split_half_rho_raw.csv", index=False)
    drho_df.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2e_split_half_drho_raw.csv", index=False)

    rho_rows = rho_band_stats(rho_df)
    drho_rows = drho_band_stats(drho_df)
    write_report(rho_rows, drho_rows)

    n_pass_rho = sum(1 for r in rho_rows
                     if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05)
    print(f"\nρ test: {n_pass_rho}/6 bands pass (pre-reg: ≥4).")
    if n_pass_rho >= 4:
        print("  → PASS: H2c survives drift-floor control.")
    elif n_pass_rho >= 2:
        print("  → PARTIAL: drift is part of the signal — lean on H2d's ρ_inert baseline.")
    else:
        print("  → FAIL: ρ cannot be cleanly separated from session drift.")


if __name__ == "__main__":
    main()
