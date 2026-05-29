#!/usr/bin/env python3
"""Audit 54 — Epi cross-phase rigidity (TARR taxonomy on the fixed leaf set E_p).

Direction A of the LRG-epilepsy research programme. Scope at
.agents/guides/task-persistence-investigation/2026-05-08_epi-cross-phase-rigidity.md.

Per (patient, band, lambda, variant): build the induced subtree on the
epi node set E_p in each of {rest_pre, task_test, rest_post}; compute the
three pairwise Kendall-Colijn distances; classify the per-patient call
into one of {anchor, trace, reset, rearrange}; build a size-matched
random-leaf-subset null (R=100); emit per-patient z-scores and the
continuous trace score Td_E = d(pre, tt) - d(tt, post). T_d > 0 = trace.

Library reuse only:
    load_channel_labels, load_epileptic_nodes  -> E_p
    load_lrg_result                            -> linkage matrices
    induced_linkage                            -> Z|_S (UPGMA on cophenetic sub-matrix)
    kc_distance                                -> three triangle distances
    wilcoxon_z, bh_fdr                         -> cohort tests

Outputs (data/audit/epi_rigidity/):
    M_class_per_patient.csv        per (patient, band, lambda, variant)
    M_class_cohort.csv             per (band, lambda, variant, class)
    M_trace_score_cohort.csv       per (band, lambda, variant)
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

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
VARIANTS = ("all", "cp")
CLASS_NAMES = ("anchor", "trace", "reset", "rearrange")
EPSILON_KC = 0.10
N_NULL = 100
RNG_SEED = 20260508
MIN_LEAF_SIZE = 3

OUT_DIR = ROOT / "data" / "audit" / "epi_rigidity"


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


def _epi_indices(masks: PatientMasks, variant: str) -> np.ndarray:
    """E_p (variant='all') or E_p^cp (variant='cp', drop nodes whose only
    epi peers are on the same probe)."""
    epi_idx = np.flatnonzero(masks.epi_mask)
    if variant == "all":
        return epi_idx
    if variant == "cp":
        epi_probes = masks.probes[epi_idx]
        keep: list[int] = []
        for k, idx in enumerate(epi_idx):
            other_probes = np.delete(epi_probes, k)
            if np.any(other_probes != epi_probes[k]):
                keep.append(int(idx))
        return np.asarray(keep, dtype=int)
    raise ValueError(f"unknown variant: {variant}")


def _induced_linkage_from_coph(coph_full: np.ndarray,
                                leaf_indices: np.ndarray) -> np.ndarray:
    """Linkage matrix of the induced subtree using a precomputed cophenetic
    full matrix. Faster than rebuilding cophenetic every call."""
    sub = coph_full[np.ix_(leaf_indices, leaf_indices)]
    return linkage(squareform(sub, checks=False), method="average")


def _classify(d_pre_tt: float, d_tt_post: float, d_pre_post: float,
              eps: float) -> str:
    pre_tt_small = d_pre_tt <= eps
    tt_post_small = d_tt_post <= eps
    pre_post_small = d_pre_post <= eps
    if pre_tt_small and tt_post_small and pre_post_small:
        return "anchor"
    if (not pre_tt_small) and tt_post_small and (not pre_post_small):
        return "trace"
    if (not pre_tt_small) and (not tt_post_small) and pre_post_small:
        return "reset"
    return "rearrange"


def _triangle(Z_pre: np.ndarray, Z_tt: np.ndarray, Z_post: np.ndarray,
              lam: float) -> tuple[float, float, float]:
    return (
        kc_distance(Z_pre, Z_tt, lam=lam, normalize=True),
        kc_distance(Z_tt, Z_post, lam=lam, normalize=True),
        kc_distance(Z_pre, Z_post, lam=lam, normalize=True),
    )


def compute_for_patient(masks: PatientMasks,
                         rng: np.random.Generator) -> list[dict]:
    pat = masks.patient
    n_chan = len(masks.channels)
    rows: list[dict] = []

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

        for variant in VARIANTS:
            S = _epi_indices(masks, variant)
            if S.size < MIN_LEAF_SIZE:
                continue
            Z_sub_phase = {p: _induced_linkage_from_coph(coph[p], S) for p in PHASES}

            null_subsets = [
                rng.choice(n_chan, size=int(S.size), replace=False)
                for _ in range(N_NULL)
            ]
            null_Z = [
                {p: _induced_linkage_from_coph(coph[p], S_r) for p in PHASES}
                for S_r in null_subsets
            ]

            for lam in KC_LAMBDAS:
                d_a, d_b, d_c = _triangle(
                    Z_sub_phase["rest_pre"], Z_sub_phase["task_test"],
                    Z_sub_phase["rest_post"], lam)
                cls_E = _classify(d_a, d_b, d_c, EPSILON_KC)
                # T_d > 0 = trace (rsPost closer to task than rsPre).
                Td_E = d_a - d_b

                class_count = {c: 0 for c in CLASS_NAMES}
                Td_null: list[float] = []
                for Zn in null_Z:
                    da, db, dc = _triangle(
                        Zn["rest_pre"], Zn["task_test"], Zn["rest_post"], lam)
                    class_count[_classify(da, db, dc, EPSILON_KC)] += 1
                    Td_null.append(da - db)
                Td_null_arr = np.asarray(Td_null)

                P_null = {c: class_count[c] / N_NULL for c in CLASS_NAMES}
                z_class: dict[str, float] = {}
                for c in CLASS_NAMES:
                    p_c = P_null[c]
                    sigma = max(np.sqrt(p_c * (1.0 - p_c) / N_NULL), 1e-3)
                    z_class[c] = ((1.0 if cls_E == c else 0.0) - p_c) / sigma
                z_trace = (Td_E - float(Td_null_arr.mean())) / max(
                    float(Td_null_arr.std(ddof=1)), 1e-9)

                row: dict[str, object] = dict(
                    patient=pat, band=band, lam=float(lam), variant=variant,
                    n_S=int(S.size),
                    d_pre_tt=float(d_a), d_tt_post=float(d_b),
                    d_pre_post=float(d_c),
                    cls_E=cls_E, Td_E=float(Td_E),
                    z_trace=float(z_trace),
                    Td_null_mean=float(Td_null_arr.mean()),
                    Td_null_std=float(Td_null_arr.std(ddof=1)),
                )
                for c in CLASS_NAMES:
                    row[f"P_null_{c}"] = float(P_null[c])
                    row[f"z_class_{c}"] = float(z_class[c])
                rows.append(row)
    return rows


def cohort_class_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for variant in VARIANTS:
        for lam in KC_LAMBDAS:
            for band in BANDS:
                sub = df[(df.variant == variant) & (df.lam == lam) &
                         (df.band == band)]
                if sub.empty:
                    continue
                for c in CLASS_NAMES:
                    z = sub[f"z_class_{c}"].to_numpy()
                    z_finite = z[np.isfinite(z)]
                    if z_finite.size < 3:
                        rows.append(dict(
                            band=band, lam=float(lam), variant=variant, cls=c,
                            n=int(z_finite.size), n_pos=0,
                            z_mean=float("nan"), p_one_sided=float("nan")))
                        continue
                    n_pos = int(np.sum(z_finite > 1.96))
                    z_score, p = wilcoxon_z(z_finite)
                    rows.append(dict(
                        band=band, lam=float(lam), variant=variant, cls=c,
                        n=int(z_finite.size), n_pos=n_pos,
                        z_mean=float(z_finite.mean()),
                        p_one_sided=float(p)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_bh"] = np.nan
    for variant in VARIANTS:
        for lam in KC_LAMBDAS:
            sel = (out.variant == variant) & (out.lam == lam)
            ps = out.loc[sel, "p_one_sided"].tolist()
            valid = [i for i, x in enumerate(ps) if np.isfinite(x)]
            if not valid:
                continue
            ps_v = [ps[i] for i in valid]
            qs_v = bh_fdr(ps_v)
            qs_full: list[float] = [float("nan")] * len(ps)
            for i, q in zip(valid, qs_v):
                qs_full[i] = q
            out.loc[sel, "q_bh"] = qs_full
    return out


def cohort_trace_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Cohort tests for the trace score.

    Two tests per (band, lambda, variant):
      - z_trace_vs_null: per-patient Td_E z-scored against random subset null,
                         then cohort Wilcoxon (less-than). Direction A's
                         primary scope test (controls for global tree changes).
      - Td_vs_zero:      per-patient raw Td_E cohort Wilcoxon (less-than).
                         Direct comparable to audit_48 fig_05 / Section 5 KC.
    """
    rows: list[dict] = []
    for variant in VARIANTS:
        for lam in KC_LAMBDAS:
            for band in BANDS:
                sub = df[(df.variant == variant) & (df.lam == lam) &
                         (df.band == band)]
                if sub.empty:
                    continue
                zt = sub["z_trace"].to_numpy()
                td = sub["Td_E"].to_numpy()
                zt = zt[np.isfinite(zt)]
                td = td[np.isfinite(td)]
                if zt.size < 3 or td.size < 3:
                    rows.append(dict(
                        band=band, lam=float(lam), variant=variant,
                        n=int(min(zt.size, td.size)), n_pos=0, n_pos_raw=0,
                        z_trace_mean=float("nan"),
                        Td_mean=float("nan"), Td_median=float("nan"),
                        p_z_trace=float("nan"), p_Td_vs_zero=float("nan")))
                    continue
                # z_trace > 0 = trace; wilcoxon_z tests 'greater' by default.
                _, p_zt = wilcoxon_z(zt)
                _, p_td = wilcoxon_z(td)
                rows.append(dict(
                    band=band, lam=float(lam), variant=variant,
                    n=int(zt.size),
                    n_pos=int(np.sum(zt > 1.96)),
                    n_pos_raw=int(np.sum(td > 0.0)),
                    z_trace_mean=float(zt.mean()),
                    Td_mean=float(td.mean()),
                    Td_median=float(np.median(td)),
                    p_z_trace=float(p_zt),
                    p_Td_vs_zero=float(p_td)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_z_trace"] = np.nan
    out["q_Td_vs_zero"] = np.nan
    for col_p, col_q in [("p_z_trace", "q_z_trace"),
                         ("p_Td_vs_zero", "q_Td_vs_zero")]:
        for variant in VARIANTS:
            for lam in KC_LAMBDAS:
                sel = (out.variant == variant) & (out.lam == lam)
                ps = out.loc[sel, col_p].tolist()
                valid = [i for i, x in enumerate(ps) if np.isfinite(x)]
                if not valid:
                    continue
                ps_v = [ps[i] for i in valid]
                qs_v = bh_fdr(ps_v)
                qs_full: list[float] = [float("nan")] * len(ps)
                for i, q in zip(valid, qs_v):
                    qs_full[i] = q
                out.loc[sel, col_q] = qs_full
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    rows: list[dict] = []
    for pat in COHORT:
        try:
            masks = _build_masks(pat)
        except Exception as e:
            print(f"[audit_54] skip {pat}: cannot build masks ({e})")
            continue
        if int(masks.epi_mask.sum()) < MIN_LEAF_SIZE:
            print(f"[audit_54] skip {pat}: |E_p|={int(masks.epi_mask.sum())} < {MIN_LEAF_SIZE}")
            continue
        n_probes = len(set(masks.probes[masks.epi_mask]))
        n_cp = int(_epi_indices(masks, "cp").size)
        print(f"[audit_54] {pat}: |E_p|={int(masks.epi_mask.sum())}, "
              f"|E_p^cp|={n_cp}, n_probes={n_probes}")
        rows.extend(compute_for_patient(masks, rng))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "M_class_per_patient.csv", index=False)
    cls_df = cohort_class_stats(df)
    cls_df.to_csv(OUT_DIR / "M_class_cohort.csv", index=False)
    trace_df = cohort_trace_stats(df)
    trace_df.to_csv(OUT_DIR / "M_trace_score_cohort.csv", index=False)

    print(f"\n[audit_54] CSVs written to {OUT_DIR}")
    print(f"  rows in M_class_per_patient: {len(df)}")
    print(f"  rows in M_class_cohort:      {len(cls_df)}")
    print(f"  rows in M_trace_score_cohort:{len(trace_df)}")

    if not cls_df.empty:
        sig = cls_df[(cls_df.q_bh < 0.05) & (cls_df.n_pos >= 7)].sort_values(
            ["variant", "lam", "band", "cls"])
        if not sig.empty:
            print("\n[audit_54] Cohort-positive class cells (BH q<0.05, n_pos>=7):")
            for _, r in sig.iterrows():
                print(f"  {r.variant}/lam={r.lam:.1f}/{r.band}/{r.cls}: "
                      f"n_pos={r.n_pos}, p={r.p_one_sided:.4g}, q={r.q_bh:.4g}")
        else:
            print("\n[audit_54] No cohort-positive class cells "
                  "(BH q<0.05 + n_pos>=7) under any (band, lambda, variant).")

    if not trace_df.empty:
        # z_trace test (random-null)
        sig_zt = trace_df[(trace_df.q_z_trace < 0.05) &
                          (trace_df.n_pos >= 7)].sort_values(
            ["variant", "lam", "band"])
        if not sig_zt.empty:
            print("\n[audit_54] Cohort-positive z_trace cells (BH q<0.05, n_pos>=7):")
            for _, r in sig_zt.iterrows():
                print(f"  z_trace {r.variant}/lam={r.lam:.1f}/{r.band}: "
                      f"n_pos={r.n_pos}, p={r.p_z_trace:.4g}, q={r.q_z_trace:.4g}")
        else:
            print("\n[audit_54] z_trace test: no cohort-positive cells.")

        # Td_E vs 0 (raw, comparable to audit_48 fig_05)
        sig_td = trace_df[(trace_df.q_Td_vs_zero < 0.05) &
                          (trace_df.n_pos_raw >= 7)].sort_values(
            ["variant", "lam", "band"])
        if not sig_td.empty:
            print("\n[audit_54] Cohort-positive Td_E vs 0 cells "
                  "(BH q<0.05, n_pos_raw>=7):")
            for _, r in sig_td.iterrows():
                print(f"  Td_vs_0 {r.variant}/lam={r.lam:.1f}/{r.band}: "
                      f"n_pos={r.n_pos_raw}/{r.n}, "
                      f"median={r.Td_median:+.3f}, "
                      f"p={r.p_Td_vs_zero:.4g}, q={r.q_Td_vs_zero:.4g}")
        else:
            print("\n[audit_54] Td_E vs 0 test: no cohort-positive cells.")


if __name__ == "__main__":
    main()
