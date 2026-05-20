#!/usr/bin/env python3
"""Audit 71 — C4 cross-probe Wilcoxon gate (replaces n_trace_xprobe ≥ 6/10).

Implements CONTROLS.md §C4 (locked 2026-05-19 Decision 9, post writing-agent
feedback): the cross-probe restriction gate is the **paired one-sided
Wilcoxon** on `(rho_split − rho_xprobe)` per patient under
`H_1: rho_split > rho_xprobe`, **failing to reject** (p ≥ 0.05) ⇒ no
significant degradation, **plus** `sign(rho_xprobe_median) ==
sign(rho_split_median)`. No patient-count threshold anywhere
(per `feedback_no_hardcoded_test_thresholds.md`).

LOO diagnostic per `feedback_no_single_patient_p_driven.md`:
`wilcoxon_loo_max_p_split_gt_xprobe` = worst-case Wilcoxon p after
dropping each patient once. Descriptive only, never a gate.

Inputs
------
- `data/audit/ctm_triangle/Td_per_patient_per_band.csv` columns
  `patient, band, rho_split, rho_split_cross_probe`.

Outputs
-------
- `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv` — per-band cohort
  Wilcoxon results + LOO max + C4 verdict.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
IN_CSV = PROJECT_ROOT / "data/audit/ctm_triangle/Td_per_patient_per_band.csv"
OUT_CSV = PROJECT_ROOT / "data/audit/ctm_triangle/c4_wilcoxon_cohort.csv"


def wilcoxon_one_sided_greater(d: np.ndarray) -> float:
    """Return paired Wilcoxon one-sided p under H_1: d > 0."""
    d = d[np.isfinite(d)]
    if d.size < 3 or np.all(d == 0):
        return 1.0
    try:
        _, p = wilcoxon(d, alternative="greater", zero_method="wilcox", correction=False)
        return float(p)
    except Exception:
        return 1.0


def loo_max_p(d_full: np.ndarray, patients: list[str]) -> tuple[float, str]:
    """Return (max LOO p, argmax patient) for paired Wilcoxon one-sided greater."""
    n = len(d_full)
    p_loo = np.zeros(n)
    for i in range(n):
        keep = np.ones(n, dtype=bool); keep[i] = False
        p_loo[i] = wilcoxon_one_sided_greater(d_full[keep])
    argmax = int(np.argmax(p_loo))
    return float(p_loo[argmax]), patients[argmax]


def main():
    df = pd.read_csv(IN_CSV)
    rows = []
    for band in BANDS:
        sub = df[df["band"] == band].copy()
        sub = sub.sort_values("patient").reset_index(drop=True)
        patients = sub["patient"].tolist()
        rho_split = sub["rho_split"].to_numpy()
        rho_xprobe = sub["rho_split_cross_probe"].to_numpy()
        d = rho_split - rho_xprobe

        # Primary test: paired Wilcoxon, H1: rho_split > rho_xprobe
        p_split_gt_xprobe = wilcoxon_one_sided_greater(d)

        # LOO max p
        loo_max, loo_patient = loo_max_p(d, patients)

        # Sign agreement on cohort medians
        med_split = float(np.median(rho_split))
        med_xprobe = float(np.median(rho_xprobe))
        sign_agree = bool(np.sign(med_split) == np.sign(med_xprobe))

        # C4 verdict: pass iff FAILS to reject (p >= 0.05) AND sign agreement
        c4_pass = bool((p_split_gt_xprobe >= 0.05) and sign_agree)

        rows.append({
            "band": band,
            "n_patients": len(patients),
            "rho_split_median": med_split,
            "rho_xprobe_median": med_xprobe,
            "paired_wilcoxon_p_split_gt_xprobe": p_split_gt_xprobe,
            "wilcoxon_loo_max_p_split_gt_xprobe": loo_max,
            "wilcoxon_loo_argmax_patient": loo_patient,
            "sign_agreement": sign_agree,
            "c4_pass": c4_pass,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_CSV, index=False)
    print(out_df.to_string(index=False))
    print(f"\nWrote: {OUT_CSV}")


if __name__ == "__main__":
    main()
