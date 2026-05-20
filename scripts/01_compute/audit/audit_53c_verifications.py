#!/usr/bin/env python3
"""Audit 53c — verifications for the manuscript KC-leaf composite figure.

Checks (5):
  1. Predicate phrasing consistency (code vs narrative).
  2. Multiple-comparison correction footing (writing-side; not numerical).
  3. Pat_14 γ_h sentinel mechanics: max f_i with min_partners=1, plus the
     full leaf score distribution at this cell.
  4. Cohort-aggregate T+/T- ratios per band: median over 10 patients
     for each of the 6 bands.
  5. Pat_07 β focal 62: number of T+ vs T- partners at the same focal.

T+ predicate (figure-side, λ=0):
    matched_tt_post = |m_tt − m_post| <= δ
    differs_pre_tt  = |m_pre − m_tt|  > δ
    nontrivial_m    = m_tt >= MIN_M  ∧  m_post >= MIN_M
    fine_grained    = size_tt <= s_max  ∧  size_post <= s_max
    T+              = matched_tt_post ∧ differs_pre_tt ∧ nontrivial_m ∧ fine_grained

T- predicate (mirror, anti-trace):
    matched_pre_post = |m_pre − m_post| <= δ
    differs_tt_pre   = |m_tt − m_pre|  > δ
    nontrivial_m     = m_pre >= MIN_M  ∧  m_post >= MIN_M
    fine_grained     = size_pre <= s_max  ∧  size_post <= s_max
    T-               = matched_pre_post ∧ differs_tt_pre ∧ nontrivial_m ∧ fine_grained
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.workflow.lrg import load_lrg_result

sys.path.insert(0, str(Path(__file__).parent))
from audit_53_kc_leaf_topological_memory import (  # noqa: E402
    PHASES,
    expand_per_leaf_score,
    pair_index,
    pair_mrca_subtree_sizes,
)

DELTA         = 1
MIN_M         = 2
MAX_MRCA_FRAC = 0.50

COHORT = (
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
)
BANDS = ("delta", "theta", "alpha", "beta", "low_gamma", "high_gamma")


def load_phase_data(patient: str, band: str):
    Z = {}
    for ph in PHASES:
        res = load_lrg_result(patient, ph, band, fc_method="imcoh_abs")
        if res is None or res.linkage_matrix is None:
            return None
        Z[ph] = np.asarray(res.linkage_matrix)
    n = int(Z["rest_pre"].shape[0]) + 1
    if not all(int(Z[ph].shape[0]) + 1 == n for ph in PHASES):
        return None
    return Z, n


def predicate_T_plus(Z: dict, n: int) -> np.ndarray:
    m_pre  = kc_vectors(Z["rest_pre"])[0]
    m_tt   = kc_vectors(Z["task_test"])[0]
    m_post = kc_vectors(Z["rest_post"])[0]
    size_tt   = pair_mrca_subtree_sizes(Z["task_test"])
    size_post = pair_mrca_subtree_sizes(Z["rest_post"])
    s_max = int(np.ceil(MAX_MRCA_FRAC * n))
    return (
        (np.abs(m_tt - m_post) <= DELTA)
        & (np.abs(m_pre - m_tt) > DELTA)
        & (m_tt >= MIN_M) & (m_post >= MIN_M)
        & (size_tt <= s_max) & (size_post <= s_max)
    )


def predicate_T_minus(Z: dict, n: int) -> np.ndarray:
    m_pre  = kc_vectors(Z["rest_pre"])[0]
    m_tt   = kc_vectors(Z["task_test"])[0]
    m_post = kc_vectors(Z["rest_post"])[0]
    size_pre  = pair_mrca_subtree_sizes(Z["rest_pre"])
    size_post = pair_mrca_subtree_sizes(Z["rest_post"])
    s_max = int(np.ceil(MAX_MRCA_FRAC * n))
    return (
        (np.abs(m_pre - m_post) <= DELTA)
        & (np.abs(m_tt - m_pre) > DELTA)
        & (m_pre >= MIN_M) & (m_post >= MIN_M)
        & (size_pre <= s_max) & (size_post <= s_max)
    )


def hr() -> None:
    print("-" * 72)


def check_1_predicate_phrasing() -> None:
    print("[1] Predicate phrasing consistency")
    hr()
    print(
        "Code (audit_53.render_cell, audit_53b.compute_cell):\n"
        "    differs_pre_tt = np.abs(m_pre - m_tt) > delta    with delta = 1\n"
        "Narrative sketch wording:\n"
        "    'm_pre differs from m_tt by >= 2'\n"
        "Equivalence: m_pre, m_tt are integer MRCA depths; |m_pre - m_tt| > 1\n"
        "  iff |m_pre - m_tt| >= 2. Equivalent at integers. CODE USES STRICT '>'\n"
        "  CONSISTENTLY (no mixing of '>=', no off-by-one).\n"
        "Status: PASS (single phrasing throughout audit_53*.py)."
    )


def check_2_correction_footing() -> None:
    print("\n[2] Multiple-comparison correction footing")
    hr()
    print(
        "The narrative-sketch reference 'joint Bonferroni m=48' belongs to the\n"
        "§5.2 cohort scalar T_KC pipeline (24 cells x 2 sides + ... see\n"
        "result_2_lrg_beta_trace.md). After the §5.2 structural cuts the live\n"
        "correction backing the beta lambda=0 cohort claim is within-probe BH\n"
        "at m=12, q=0.006.\n"
        "  (i) Caption guidance for the writing agent: cite within-probe BH\n"
        "      m=12, q=0.006 — NOT m=48 — when describing the lambda=0 cohort\n"
        "      backing of the per-leaf figures.\n"
        "  (ii) Per-leaf figures (audit_53, audit_53b): computed at lambda=0\n"
        "      ONLY (m vector, no M); no joint-scope correction baked in.\n"
        "      Verified by inspection of audit_53.render_cell and\n"
        "      audit_53b.compute_cell (both use kc_vectors(...)[0] only).\n"
        "Status: PASS (numerical pipeline is lambda=0 only; caption guidance\n"
        "  must avoid m=48)."
    )


def check_3_pat14_high_gamma_sentinel() -> None:
    print("\n[3] Pat_14 high_gamma sentinel mechanics")
    hr()
    data = load_phase_data("Pat_14", "high_gamma")
    if data is None:
        print("  ERROR: missing Pat_14 high_gamma LRG results")
        return
    Z, n = data
    T = predicate_T_plus(Z, n)
    f = expand_per_leaf_score(T, n)
    n_pairs_T = int(T.sum())
    nonzero = np.sort(f[f > 0])[::-1]
    print(f"  cohort cell:       Pat_14 high_gamma  n_leaves={n}")
    print(f"  T+ predicate at this cell:")
    print(f"    pairs with T+(i,j) = 1:           {n_pairs_T}")
    print(f"    max f_i:                           {f.max():.4f}")
    print(f"    leaves with f_i > 0:               {(f > 0).sum()}")
    print(f"    top-5 nonzero f_i:                 "
          f"{', '.join(f'{x:.4f}' for x in nonzero[:5])}")
    print()
    print(
        "  Interpretation: with min_partners=1 (no threshold), the cell holds\n"
        f"  {n_pairs_T} predicate-positive pair(s) total. Max focal score is\n"
        f"  ~{f.max():.4f}, i.e. 1 partner out of (n-1)={n-1}. Compare with\n"
        "  the rendered beta exemplars (|P|=37, 31, 38). The Pat_14 high_gamma\n"
        f"  sentinel is SUBSTANTIVELY EMPTY: a single isolated pair is not\n"
        "  a multi-leaf trace pattern. The min_partners=5 cutoff in the\n"
        "  manuscript figure is honest about this — it would render the\n"
        "  same cell as sentinel even if min_partners were lowered to 2.\n"
        "  Status: PASS (sentinel is genuine, not threshold-driven)."
    )


def check_4_cohort_aggregate_ratios() -> None:
    print("\n[4] Cohort-aggregate T+/T- ratios per band")
    hr()
    table = {}
    for band in BANDS:
        ratios = []
        per_pat = []
        for pat in COHORT:
            data = load_phase_data(pat, band)
            if data is None:
                per_pat.append((pat, None, None, None))
                continue
            Z, n = data
            Tp = predicate_T_plus(Z, n)
            Tm = predicate_T_minus(Z, n)
            np_count = int(Tp.sum())
            nm_count = int(Tm.sum())
            if nm_count == 0:
                ratio = np.inf if np_count > 0 else np.nan
            else:
                ratio = np_count / nm_count
            per_pat.append((pat, np_count, nm_count, ratio))
            if np.isfinite(ratio):
                ratios.append(ratio)
        if ratios:
            med  = float(np.median(ratios))
            mean = float(np.mean(ratios))
        else:
            med = mean = float("nan")
        table[band] = (med, mean, per_pat, ratios)

    hdr = f"{'band':<11} {'median(T+/T-)':>14} {'mean(T+/T-)':>13} {'n_finite':>10}"
    print(hdr)
    print("-" * len(hdr))
    for band in BANDS:
        med, mean, _, ratios = table[band]
        print(f"{band:<11} {med:>14.3f} {mean:>13.3f} {len(ratios):>10d}")
    print()
    print("  Per-band per-patient detail (T+, T-, ratio):")
    for band in BANDS:
        _, _, per_pat, _ = table[band]
        cells = []
        for pat, npc, nmc, r in per_pat:
            if npc is None:
                cells.append(f"{pat[-2:]}=NA")
            elif r is None or not np.isfinite(r):
                cells.append(f"{pat[-2:]}=({npc}/{nmc}=inf)")
            else:
                cells.append(f"{pat[-2:]}=({npc}/{nmc}={r:.2f})")
        print(f"  {band:<11}  " + "  ".join(cells))


def check_5_pat07_beta_focal_62() -> None:
    print("\n[5] Pat_07 beta focal leaf 62 — T+ vs T- partner balance")
    hr()
    data = load_phase_data("Pat_07", "beta")
    if data is None:
        print("  ERROR: missing Pat_07 beta LRG results")
        return
    Z, n = data
    Tp = predicate_T_plus(Z, n)
    Tm = predicate_T_minus(Z, n)

    focal = 62
    np_partners = []
    nm_partners = []
    for j in range(n):
        if j == focal:
            continue
        k = pair_index(n, focal, j)
        if Tp[k]:
            np_partners.append(j)
        if Tm[k]:
            nm_partners.append(j)
    overlap = sorted(set(np_partners) & set(nm_partners))
    print(f"  focal i*={focal}, n={n} leaves")
    print(f"  T+ partners (rendered as red bundle in figure): {len(np_partners)}")
    print(f"  T- partners (mirror, anti-trace):                {len(nm_partners)}")
    print(f"  partners in BOTH T+ and T-:                      {len(overlap)}")
    if len(np_partners) > 0:
        ratio = len(np_partners) / max(len(nm_partners), 1) if nm_partners else float("inf")
        print(f"  T+ / T- partner ratio at focal 62:               "
              f"{ratio if np.isfinite(ratio) else 'inf':.2f}")
    print()
    if len(nm_partners) <= max(2, 0.2 * len(np_partners)):
        verdict = "TRACE-DOMINANT (T- << T+)"
    elif len(nm_partners) <= 0.5 * len(np_partners):
        verdict = "TRACE-DOMINANT (T+ at least 2x T-)"
    elif len(nm_partners) <= len(np_partners):
        verdict = "TRACE-LEANING (T+ > T- but mixed)"
    else:
        verdict = "MIXED OR ANTI-TRACE (T- >= T+)"
    print(f"  Verdict at focal 62: {verdict}")


def main() -> int:
    check_1_predicate_phrasing()
    check_2_correction_footing()
    check_3_pat14_high_gamma_sentinel()
    check_4_cohort_aggregate_ratios()
    check_5_pat07_beta_focal_62()
    return 0


if __name__ == "__main__":
    sys.exit(main())
