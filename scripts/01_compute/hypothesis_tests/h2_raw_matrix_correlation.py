#!/usr/bin/env python3
"""H2-RAW — direct Spearman correlation between ultrametric matrices.

The most intuitive question, which we skipped:

    Is ρ(D^rpost, D^ttest)  >  ρ(D^rpost, D^rpre)?

i.e. does the post-task resting ultrametric look more like task than
like rest_pre, as a direct rank-correlation on the upper-triangle
matrices — no Δ-vector subtraction, no partition binarization. Pure
matrix similarity.

Compared to the current H2 family:
  * H2a (VI magnitude) — failed at n=9. Symmetric, categorical, lossy.
  * H2c (Spearman on Δ vectors) — passes but "directional" only in the
    signed-Δ sense. Does not directly answer "is rpost closer to task".
  * **H2-RAW (this script)** — directed at the role level, magnitude-
    sensitive, operates on continuous ultrametric matrices. Closest to
    what a reviewer's intuition asks for.

Per (patient, band, phase_pair) we compute Spearman ρ on the
upper-triangle cophenetic distance matrices. We then derive four
cross-patient tests:

  H1-raw:  ρ(tl, tt) > mean cross-type (H1 equivalent).
  H2a-raw: ρ(rpost, tt) > ρ(rpost, rpre)  (H2a equivalent, directed).
  H2b-raw: ρ(rpre, tt) > ρ(tt, rpost)    (H2b equivalent; Spearman is
           arg-swap symmetric, so this is a role-asymmetric contrast on
           the same matrix pairs).
  H3-raw:  mean ρ(within-type) > mean ρ(cross-type).

Reads:  LRG imcoh_abs caches via `load_lrg_result`.
Writes: data/reports/imcoh_vi/h2_raw_matrix_correlation.{md,csv}
        data/reports/imcoh_vi/h2_raw_matrix_correlation_raw.csv
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr
from h2c_ultrametric_drift import upper_tri


PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_TYPE = {"rest_pre": "rest", "rest_post": "rest",
              "task_learn": "task", "task_test": "task"}


# ─────────────────────────── per-pair ρ ───────────────────────────


def load_D(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.ultrametric_matrix) if r is not None else None


def spearman_matrices(D1, D2) -> float:
    if D1 is None or D2 is None:
        return np.nan
    u1, u2 = upper_tri(D1), upper_tri(D2)
    if u1.shape != u2.shape or u1.std() == 0 or u2.std() == 0:
        return np.nan
    r, _ = stats.spearmanr(u1, u2)
    return float(r) if np.isfinite(r) else np.nan


def collect_raw() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        Ds = {phase: load_D(pat, phase, band)
              for phase in PHASES
              for band in [None]}  # placeholder only — we iterate below
        for band in BRAIN_BANDS_NAMES:
            Dp = {phase: load_D(pat, phase, band) for phase in PHASES}
            for p1, p2 in combinations(PHASES, 2):
                rows.append({
                    "patient": pat, "band": band,
                    "p1": p1, "p2": p2,
                    "rho": spearman_matrices(Dp[p1], Dp[p2]),
                })
    return pd.DataFrame(rows)


# ─────────────────────────── cross-patient stats ───────────────────────────


def _rho_pair(df: pd.DataFrame, pat: str, band: str, p1: str, p2: str) -> float:
    """Unordered phase-pair lookup into the raw dataframe.

    `combinations(PHASES, 2)` uses tuple order, not alphabetical — so we
    try both orderings.
    """
    base = df[(df["patient"] == pat) & (df["band"] == band)]
    hit = base[(base["p1"] == p1) & (base["p2"] == p2)]
    if len(hit) == 0:
        hit = base[(base["p1"] == p2) & (base["p2"] == p1)]
    return float(hit["rho"].iloc[0]) if len(hit) else np.nan


def per_band_test(df: pd.DataFrame, test_fn, patients: list[str]
                  ) -> list[dict]:
    out = []
    for band in BRAIN_BANDS_NAMES:
        vals = np.array([test_fn(df, pat, band) for pat in patients], dtype=float)
        pats = [p for p, v in zip(patients, vals) if np.isfinite(v)]
        vals = vals[np.isfinite(vals)]
        if len(vals) < 3:
            out.append({"band": band, "n": len(vals)})
            continue
        z, p = wilcoxon_z(vals)
        r_rb = rank_biserial(vals)
        mean, lo, hi = boot_ci_mean(vals)
        out.append({"band": band, "n": len(vals), "mean": mean,
                    "ci_lo": lo, "ci_hi": hi, "r_rb": r_rb, "z": z, "p": p,
                    "patients": pats, "values": vals.tolist()})
    ps = [r.get("p", np.nan) for r in out]
    finite = [np.isfinite(pp) for pp in ps]
    qs = bh_fdr([pp for pp, ok in zip(ps, finite) if ok])
    j = 0
    for r, ok in zip(out, finite):
        r["q"] = qs[j] if ok else np.nan
        if ok:
            j += 1
    return out


def _h2a_raw(df, pat, band):
    """ρ(rpost, tt) − ρ(rpost, rpre). Positive ⇒ rpost closer to task."""
    return (_rho_pair(df, pat, band, "rest_post", "task_test")
            - _rho_pair(df, pat, band, "rest_post", "rest_pre"))


def _h2b_raw(df, pat, band):
    """ρ(rpre, tt) − ρ(tt, rpost). Positive ⇒ "approach" more similar than "exit"."""
    return (_rho_pair(df, pat, band, "rest_pre", "task_test")
            - _rho_pair(df, pat, band, "rest_post", "task_test"))


def _h1_raw(df, pat, band):
    """ρ(tl, tt) − mean cross-type ρ. Positive ⇒ task phases more self-similar."""
    r_task = _rho_pair(df, pat, band, "task_learn", "task_test")
    crosses = [_rho_pair(df, pat, band, p1, p2)
               for p1 in PHASES for p2 in PHASES
               if p1 < p2 and PHASE_TYPE[p1] != PHASE_TYPE[p2]]
    crosses = [c for c in crosses if np.isfinite(c)]
    if not crosses or not np.isfinite(r_task):
        return np.nan
    return r_task - float(np.mean(crosses))


def _h3_raw(df, pat, band):
    """mean ρ(within-type) − mean ρ(cross-type). Positive ⇒ phase-type is real."""
    within, cross = [], []
    for p1, p2 in combinations(PHASES, 2):
        v = _rho_pair(df, pat, band, p1, p2)
        if not np.isfinite(v):
            continue
        (within if PHASE_TYPE[p1] == PHASE_TYPE[p2] else cross).append(v)
    if not within or not cross:
        return np.nan
    return float(np.mean(within) - np.mean(cross))


# ─────────────────────────── report ───────────────────────────


def _table(rows, tex, metric_name):
    lines = [f"| band | n | mean {metric_name} | 95 % CI | r_rb | z | p | q (BH) |",
             "|------|--:|--------:|:-------|-----:|--:|--:|-------:|"]
    for r in rows:
        tx = tex[r["band"]]
        if "mean" not in r:
            lines.append(f"| {tx} | {r['n']} | — | — | — | — | — | — |")
            continue
        star = "★" if r["q"] < 0.05 else ""
        lines.append(
            f"| {tx} | {r['n']} | {r['mean']:+.4f} | "
            f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | {r['r_rb']:+.3f} | "
            f"{r['z']:+.2f} | {r['p']:.4f} | {r['q']:.4f}{star} |"
        )
    return lines


def write_report(raw_df, h1, h2a, h2b, h3) -> None:
    out_dir = REPORTS_ROOT / "imcoh_vi"
    md = out_dir / "h2_raw_matrix_correlation.md"
    csv = out_dir / "h2_raw_matrix_correlation.csv"

    lines: list[str] = []
    ap = lines.append
    ap("# H2-RAW — direct Spearman correlation between ultrametric matrices")
    ap("")
    ap("**Metric per (patient, band, phase_pair):** "
       "ρ = Spearman(upper-tri D^p1, upper-tri D^p2) on the cophenetic "
       "ultrametric matrices. No Δ-vector subtraction, no partition cut — "
       "the most direct continuous comparison between two LRG outputs.")
    ap("")
    ap("Four cross-patient tests (all paired one-sample Wilcoxon against 0, "
       "FDR-BH m=6, 10k bootstrap CI). The H2a-raw is the intuitive "
       "replacement of the failing H2a VI-magnitude test.")
    ap("")

    ap("## H1-raw — task phases self-similar vs cross-type")
    ap("")
    ap("`ρ(task_learn, task_test) − mean ρ(cross-type pair)`")
    ap("")
    lines.extend(_table(h1, BRAIN_BAND_TEX_DICT, "Δρ"))
    ap("")

    ap("## H2a-raw — rpost closer to task than to rest_pre")
    ap("")
    ap("`ρ(rest_post, task_test) − ρ(rest_post, rest_pre)`")
    ap("")
    ap("**This is the direct, magnitude-sensitive version of the failing "
       "H2a test. Positive ⇒ rpost's ultrametric is more similar to task's "
       "than to pre-task rest's.**")
    ap("")
    lines.extend(_table(h2a, BRAIN_BAND_TEX_DICT, "Δρ"))
    ap("")

    ap("## H2b-raw — approach (rpre↔tt) > exit (tt↔rpost)")
    ap("")
    ap("`ρ(rest_pre, task_test) − ρ(task_test, rest_post)`")
    ap("")
    lines.extend(_table(h2b, BRAIN_BAND_TEX_DICT, "Δρ"))
    ap("")

    ap("## H3-raw — within-type > cross-type")
    ap("")
    ap("`mean ρ(within-type pair) − mean ρ(cross-type pair)`")
    ap("")
    lines.extend(_table(h3, BRAIN_BAND_TEX_DICT, "Δρ"))
    ap("")

    records = []
    for tag, rows in [("H1_raw", h1), ("H2a_raw", h2a),
                       ("H2b_raw", h2b), ("H3_raw", h3)]:
        for r in rows:
            if "mean" not in r:
                continue
            records.append({"hypothesis": tag, **{k: r.get(k) for k in
                ("band", "n", "mean", "ci_lo", "ci_hi", "r_rb", "z", "p", "q")}})
    pd.DataFrame(records).to_csv(csv, index=False)
    raw_df.to_csv(out_dir / "h2_raw_matrix_correlation_raw.csv", index=False)

    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    print(f"wrote {csv}")

    # Pass counts per hypothesis
    for tag, rows in [("H1_raw", h1), ("H2a_raw", h2a),
                       ("H2b_raw", h2b), ("H3_raw", h3)]:
        n_pass = sum(1 for r in rows
                     if np.isfinite(r.get("q", np.nan)) and r["q"] < 0.05)
        print(f"  {tag}: {n_pass}/6 bands pass FDR")


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    raw = collect_raw()
    patients = list(PATIENTS_LIST)
    h1  = per_band_test(raw, _h1_raw,  patients)
    h2a = per_band_test(raw, _h2a_raw, patients)
    h2b = per_band_test(raw, _h2b_raw, patients)
    h3  = per_band_test(raw, _h3_raw,  patients)
    write_report(raw, h1, h2a, h2b, h3)


if __name__ == "__main__":
    main()
