#!/usr/bin/env python3
"""H2-FROB — Frobenius-distance ratios between ultrametric matrices.

Magnitude-sensitive complement to the rank-based H2-RAW test. Per
(patient, band):

    r_a = ||D^rpost − D^ttest||_F / ||D^rpost − D^rpre||_F
    r_b = ||D^rpost − D^tlearn||_F / ||D^rpost − D^rpre||_F

If the trace claim holds at the magnitude level, r_a < 1 (rpost is
closer, in matrix-distance, to task than to rest_pre). r_b replicates
with the earlier task phase.

Tests (paired Wilcoxon across patients, FDR-BH m=6):
    H2a-FROB (target=ttest):  r_a < 1  →  log(r_a) < 0
    H2a-FROB (target=tlearn): r_b < 1  →  log(r_b) < 0

Reads:  LRG imcoh_abs caches.
Writes: data/reports/imcoh_vi/h2_frobenius_ratio.{md,csv,raw.csv}
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr
from h2c_ultrametric_drift import upper_tri


def load_D(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    return upper_tri(np.asarray(r.ultrametric_matrix)) if r is not None else None


def frob(d1, d2) -> float:
    if d1 is None or d2 is None or d1.shape != d2.shape:
        return np.nan
    return float(np.linalg.norm(d1 - d2))


def collect_raw() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            D_pre   = load_D(pat, "rest_pre",   band)
            D_post  = load_D(pat, "rest_post",  band)
            D_tt    = load_D(pat, "task_test",  band)
            D_tl    = load_D(pat, "task_learn", band)
            base = frob(D_post, D_pre)
            if not np.isfinite(base) or base == 0:
                continue
            rows.append({
                "patient": pat, "band": band,
                "r_ttest":  frob(D_post, D_tt) / base,
                "r_tlearn": frob(D_post, D_tl) / base,
                "base_pre_post": base,
            })
    return pd.DataFrame(rows)


def per_band(raw: pd.DataFrame, col: str) -> list[dict]:
    out = []
    for band in BRAIN_BANDS_NAMES:
        vals = raw[raw["band"] == band][col].dropna().to_numpy()
        if len(vals) < 3:
            out.append({"band": band, "n": len(vals)})
            continue
        # Test: r < 1, i.e. log(r) < 0. Flip sign so wilcoxon_z 'greater' works.
        lr = -np.log(vals)
        z, p = wilcoxon_z(lr)
        r_rb = rank_biserial(lr)
        mean, lo, hi = boot_ci_mean(vals)
        out.append({"band": band, "n": len(vals),
                    "mean_ratio": mean, "ci_lo": lo, "ci_hi": hi,
                    "r_rb": r_rb, "z": z, "p": p})
    ps = [r.get("p", np.nan) for r in out]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(out, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return out


def _table(rows):
    lines = ["| band | n | mean ratio | 95 % CI | r_rb | z | p | q (BH) |",
             "|------|--:|----------:|:-------|-----:|--:|--:|-------:|"]
    for r in rows:
        tx = BRAIN_BAND_TEX_DICT[r["band"]]
        if "mean_ratio" not in r:
            lines.append(f"| {tx} | {r['n']} | — | — | — | — | — | — |")
            continue
        star = "★" if r["q"] < 0.05 else ""
        lines.append(
            f"| {tx} | {r['n']} | {r['mean_ratio']:.3f} | "
            f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] | {r['r_rb']:+.3f} | "
            f"{r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{star} |"
        )
    return lines


def main() -> None:
    raw = collect_raw()
    raw.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2_frobenius_ratio_raw.csv", index=False)

    ttest  = per_band(raw, "r_ttest")
    tlearn = per_band(raw, "r_tlearn")

    lines: list[str] = []
    ap = lines.append
    ap("# H2-FROB — Frobenius-distance ratio between ultrametric matrices")
    ap("")
    ap("**Metric per (patient, band):** ratio r = ||D^rpost − D^target||_F / "
       "||D^rpost − D^rpre||_F over the upper-triangle cophenetic distances. "
       "If rpost is closer, in matrix distance, to target (task) than to "
       "rpre, then r < 1. **Magnitude-sensitive direct test; the complement "
       "of H2-RAW's rank-based version.**")
    ap("")
    ap("Test: log(r) < 0, one-sample Wilcoxon, FDR-BH m=6, 10k bootstrap "
       "CI on the raw ratio.")
    ap("")
    ap("## Target = task_test")
    ap("")
    lines.extend(_table(ttest))
    ap("")
    ap("## Target = task_learn (replication)")
    ap("")
    lines.extend(_table(tlearn))
    ap("")

    records = [{"target": "ttest",  **{k: r.get(k) for k in
                ("band", "n", "mean_ratio", "ci_lo", "ci_hi",
                 "r_rb", "z", "p", "q")}} for r in ttest]
    records += [{"target": "tlearn", **{k: r.get(k) for k in
                 ("band", "n", "mean_ratio", "ci_lo", "ci_hi",
                  "r_rb", "z", "p", "q")}} for r in tlearn]
    pd.DataFrame(records).to_csv(
        REPORTS_ROOT / "imcoh_vi" / "h2_frobenius_ratio.csv", index=False)

    md = REPORTS_ROOT / "imcoh_vi" / "h2_frobenius_ratio.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")

    n_tt = sum(1 for r in ttest  if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05)
    n_tl = sum(1 for r in tlearn if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05)
    print(f"  target=ttest:  {n_tt}/6 bands pass")
    print(f"  target=tlearn: {n_tl}/6 bands pass")


if __name__ == "__main__":
    main()
