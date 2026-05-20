#!/usr/bin/env python3
"""Audit 56 — Cohort KC trace, full tree vs epi-resected (V \\ E_p).

Direct testable follow-up from the 2026-05-08 Direction A finding
(`.agents/reports/2026-05-08_direction-a-induced-subtree.md`):
the epi-induced subtree at beta/lambda=1 shows RESET (heights revert
in rest_post) at the same band where the global cohort shows the
load-bearing TRACE (Result 2). If the epi nodes are confounders that
dampen the global signal, removing them should sharpen the trace.
If the epi nodes are carrying the trace, removing them should weaken
it.

Per (patient, band, lambda): load LRG linkage matrices for all three
phases. Compute Td_KC = d_KC(tt, post) - d_KC(pre, tt) on:

    full      = the cached full-tree linkage on N leaves
    resect    = the induced subtree on V \\ E_p (drop epi leaves)
    epi_only  = the induced subtree on E_p (audit_54's primary object,
                included here for direct comparison)

Cohort comparison: per-band, per-lambda median Td_KC + IQR + n_neg
across the three operationalisations. Lead with the cohort shape; BH
gating in a footnote.

Library reuse only. Same epi-mask building as audit_54.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.utils.io import load_channel_labels, load_epileptic_nodes
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, wilcoxon_z
from lrg_eegfc.utils.metrics.tree_distance import kc_distance
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14"]
BANDS = list(BRAIN_BANDS_NAMES)
PHASES = ("rest_pre", "task_test", "rest_post")
KC_LAMBDAS = (0.0, 0.5, 1.0)
OUT_DIR = ROOT / "data" / "audit" / "trace_minus_epi"


def _probe_of(label: str) -> str:
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


@dataclass
class PatientMasks:
    patient: str
    channels: list[str]
    epi_mask: np.ndarray
    probes: np.ndarray


def _build_masks(patient: str) -> PatientMasks:
    ch = load_channel_labels(patient)
    epi_set = set(load_epileptic_nodes(patient))
    epi_mask = np.array([c in epi_set for c in ch])
    probes = np.array([_probe_of(c) for c in ch], dtype=object)
    return PatientMasks(patient, ch, epi_mask, probes)


def _induced_linkage_from_coph(coph_full: np.ndarray,
                                leaf_indices: np.ndarray) -> np.ndarray:
    sub = coph_full[np.ix_(leaf_indices, leaf_indices)]
    return linkage(squareform(sub, checks=False), method="average")


def _td_kc(Z_pre: np.ndarray, Z_tt: np.ndarray, Z_post: np.ndarray,
           lam: float) -> tuple[float, float, float]:
    d_pre_tt = kc_distance(Z_pre, Z_tt, lam=lam, normalize=True)
    d_tt_post = kc_distance(Z_tt, Z_post, lam=lam, normalize=True)
    d_pre_post = kc_distance(Z_pre, Z_post, lam=lam, normalize=True)
    return d_pre_tt, d_tt_post, d_pre_post


def compute_for_patient(masks: PatientMasks) -> list[dict]:
    pat = masks.patient
    epi_idx = np.flatnonzero(masks.epi_mask)
    non_epi_idx = np.flatnonzero(~masks.epi_mask)
    rows: list[dict] = []

    if epi_idx.size < 3 or non_epi_idx.size < 3:
        return rows

    for band in BANDS:
        Zs: dict[str, np.ndarray] = {}
        coph: dict[str, np.ndarray] = {}
        for phase in PHASES:
            r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs")
            if r is None or r.linkage_matrix is None:
                Zs = {}
                break
            Z = np.asarray(r.linkage_matrix)
            Zs[phase] = Z
            coph[phase] = squareform(cophenet(Z))
        if not Zs:
            continue

        Z_resect = {p: _induced_linkage_from_coph(coph[p], non_epi_idx)
                    for p in PHASES}
        Z_epi = {p: _induced_linkage_from_coph(coph[p], epi_idx)
                 for p in PHASES}

        for lam in KC_LAMBDAS:
            for variant, Zd in (("full", Zs), ("resect", Z_resect),
                                  ("epi_only", Z_epi)):
                d_a, d_b, d_c = _td_kc(
                    Zd["rest_pre"], Zd["task_test"], Zd["rest_post"], lam)
                rows.append(dict(
                    patient=pat, band=band, lam=float(lam),
                    variant=variant,
                    n_leaves=int(masks.epi_mask.size if variant == "full"
                                  else (non_epi_idx.size if variant == "resect"
                                        else epi_idx.size)),
                    d_pre_tt=float(d_a),
                    d_tt_post=float(d_b),
                    d_pre_post=float(d_c),
                    Td_KC=float(d_b - d_a),
                ))
    return rows


def cohort_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per (variant, band, lam): median, IQR, n_neg, raw Wilcoxon p (less)."""
    rows: list[dict] = []
    for variant in ("full", "resect", "epi_only"):
        for lam in KC_LAMBDAS:
            for band in BANDS:
                v = df[(df.variant == variant) & (df.band == band) &
                       (df.lam == lam)].Td_KC.dropna().to_numpy()
                if v.size < 3:
                    continue
                med = float(np.median(v))
                q1 = float(np.percentile(v, 25))
                q3 = float(np.percentile(v, 75))
                iqr = q3 - q1
                pct_iqr = (abs(med) / iqr * 100.0) if iqr > 1e-9 else float("nan")
                n_neg = int(np.sum(v < 0))
                _, p = wilcoxon_z(-v)  # alternative='less' on raw v
                rows.append(dict(
                    variant=variant, band=band, lam=float(lam),
                    n=int(v.size),
                    median=med, q1=q1, q3=q3, iqr=iqr,
                    pct_iqr=pct_iqr, n_neg=n_neg,
                    p_one_sided=float(p),
                ))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_bh"] = np.nan
    for variant in ("full", "resect", "epi_only"):
        for lam in KC_LAMBDAS:
            sel = (out.variant == variant) & (out.lam == lam)
            ps = out.loc[sel, "p_one_sided"].tolist()
            valid = [i for i, x in enumerate(ps) if np.isfinite(x)]
            if not valid:
                continue
            qs_v = bh_fdr([ps[i] for i in valid])
            qs_full: list[float] = [float("nan")] * len(ps)
            for i, q in zip(valid, qs_v):
                qs_full[i] = q
            out.loc[sel, "q_bh"] = qs_full
    return out


def comparison_table(cohort: pd.DataFrame) -> pd.DataFrame:
    """Side-by-side full vs resect per (band, lam) with Δ-shifts."""
    rows: list[dict] = []
    for band in BANDS:
        for lam in KC_LAMBDAS:
            f = cohort[(cohort.variant == "full") & (cohort.band == band) &
                       (cohort.lam == lam)]
            r = cohort[(cohort.variant == "resect") & (cohort.band == band) &
                       (cohort.lam == lam)]
            e = cohort[(cohort.variant == "epi_only") & (cohort.band == band) &
                       (cohort.lam == lam)]
            if f.empty or r.empty:
                continue
            row = dict(
                band=band, lam=float(lam),
                full_med=float(f["median"].iloc[0]),
                full_pctIQR=float(f["pct_iqr"].iloc[0]),
                full_n_neg=int(f["n_neg"].iloc[0]),
                resect_med=float(r["median"].iloc[0]),
                resect_pctIQR=float(r["pct_iqr"].iloc[0]),
                resect_n_neg=int(r["n_neg"].iloc[0]),
                delta_med=float(r["median"].iloc[0]) - float(f["median"].iloc[0]),
            )
            if not e.empty:
                row.update(epi_med=float(e["median"].iloc[0]),
                            epi_pctIQR=float(e["pct_iqr"].iloc[0]),
                            epi_n_neg=int(e["n_neg"].iloc[0]))
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for pat in COHORT:
        try:
            masks = _build_masks(pat)
        except Exception as e:
            print(f"[audit_56] skip {pat}: cannot build masks ({e})")
            continue
        if int(masks.epi_mask.sum()) < 3:
            print(f"[audit_56] skip {pat}: |E_p|<3")
            continue
        n_non = int((~masks.epi_mask).sum())
        print(f"[audit_56] {pat}: |E_p|={int(masks.epi_mask.sum())}, "
              f"|V \\ E_p|={n_non}, N={masks.epi_mask.size}")
        rows.extend(compute_for_patient(masks))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "Td_KC_per_patient.csv", index=False)
    cohort = cohort_summary(df)
    cohort.to_csv(OUT_DIR / "Td_KC_cohort_summary.csv", index=False)
    cmp = comparison_table(cohort)
    cmp.to_csv(OUT_DIR / "Td_KC_full_vs_resect.csv", index=False)

    print(f"\n[audit_56] CSVs at {OUT_DIR}")

    # Lead with the comparison: per-band, per-lambda full vs resect
    print("\n=== Per-(band, lam) cohort Td_KC: FULL vs RESECT (V\\E_p) ===")
    print(f"{'band':>10s} {'lam':>4s} | "
          f"{'full med (%IQR, n_neg)':>26s} | "
          f"{'resect med (%IQR, n_neg)':>28s} | "
          f"{'Δmed':>7s} | direction")
    print("-" * 95)
    for _, row in cmp.iterrows():
        f_str = f"{row.full_med:+.3f} ({row.full_pctIQR:.0f}%, {row.full_n_neg}/9)"
        rs_str = f"{row.resect_med:+.3f} ({row.resect_pctIQR:.0f}%, {row.resect_n_neg}/9)"
        d_str = f"{row.delta_med:+.3f}"
        if row.full_med < 0 and row.delta_med < 0:
            tag = "trace strengthens"
        elif row.full_med < 0 and row.delta_med > 0:
            tag = "trace weakens"
        elif row.full_med > 0 and row.delta_med > 0:
            tag = "reset strengthens"
        elif row.full_med > 0 and row.delta_med < 0:
            tag = "reset weakens / flip"
        else:
            tag = ""
        print(f"{row.band:>10s} {row.lam:>4.1f} | {f_str:>26s} | {rs_str:>28s} | "
              f"{d_str:>7s} | {tag}")

    # Beta deep-dive (the prediction)
    print("\n=== Prediction check: BETA trace (full vs resect) ===")
    for lam in KC_LAMBDAS:
        sub = cmp[(cmp.band == "beta") & (cmp.lam == lam)]
        if sub.empty:
            continue
        row = sub.iloc[0]
        print(f"  beta  lam={lam}: full Td_med={row.full_med:+.3f} "
              f"({row.full_n_neg}/9 trace) → resect Td_med={row.resect_med:+.3f} "
              f"({row.resect_n_neg}/9 trace).  Δ = {row.delta_med:+.3f}")


if __name__ == "__main__":
    main()
