"""Q2: absolute distance scale alongside Z.

Reads the existing per-patient audit.csv + null.npz produced by
audit_25_raw_fc_phase_distance.py, and adds:

- Within-rest_pre null median + IQR (Q3-Q1) of the 50-split distribution.
- Observed d for each phase pair (already in audit.csv as `d_obs`).
- Ratio = d_obs / null_median.

Writes:
- ``data/audit/raw_fc_phase_distance/<Patient>/audit_with_scale.csv``
- ``data/audit/raw_fc_phase_distance/cohort_scale_summary.csv``
  with cohort-median d_obs across all n=10 patients per band x distance x pair.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import DATA_ROOT

N_COHORT = len(PATIENTS_4PHASE)  # single source of truth — never hardcode

DISTANCES = ("P", "S", "F")
PHASE_PAIRS = (
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
)

OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"

per_patient_rows = []
for p in PATIENTS_4PHASE:
    audit_csv = OUT_BASE / p / "audit.csv"
    null_npz = OUT_BASE / p / "null.npz"
    if not audit_csv.exists() or not null_npz.exists():
        print(f"[{p}] missing artefacts, skipping", file=sys.stderr)
        continue
    df = pd.read_csv(audit_csv)
    nulls = np.load(null_npz)

    for band in BRAIN_BANDS_NAMES:
        for d_label in DISTANCES:
            key = f"{band}_{d_label}"
            if key not in nulls.files:
                continue
            arr = np.asarray(nulls[key])
            arr_finite = arr[np.isfinite(arr)]
            if arr_finite.size < 3:
                continue
            null_median = float(np.median(arr_finite))
            q1 = float(np.percentile(arr_finite, 25))
            q3 = float(np.percentile(arr_finite, 75))
            null_iqr = q3 - q1
            null_mad = float(np.median(np.abs(arr_finite - null_median)))

            for phi_A, phi_B in PHASE_PAIRS:
                row = df[
                    (df.band == band)
                    & (df.distance == d_label)
                    & (df.phase_A == phi_A)
                    & (df.phase_B == phi_B)
                ]
                if row.empty:
                    continue
                d_obs = float(row.iloc[0].d_obs)
                z = float(row.iloc[0].z)
                ratio = d_obs / null_median if null_median > 0 else float("nan")
                per_patient_rows.append({
                    "patient": p,
                    "band": band,
                    "distance": d_label,
                    "phase_A": phi_A,
                    "phase_B": phi_B,
                    "d_obs": d_obs,
                    "null_median": null_median,
                    "null_q1": q1,
                    "null_q3": q3,
                    "null_iqr": null_iqr,
                    "null_mad": null_mad,
                    "z": z,
                    "ratio_obs_over_nullmed": ratio,
                })

per_patient_df = pd.DataFrame(per_patient_rows)
per_patient_df.to_csv(OUT_BASE / "per_patient_with_scale.csv", index=False)

# Per-patient extended audit is intentionally NOT written —
# `per_patient_with_scale.csv` (long-format) at the cohort level
# already contains every per-patient row. Splitting it per-patient
# duplicates data and clutters the per-patient subfolders.

# Cohort summary: median across all 10 patients (Pat_03 included; Z is
# dimensionless w.r.t. fs and ratios are unitless)
in_pool = per_patient_df
cohort_rows = []
for band in BRAIN_BANDS_NAMES:
    for d_label in DISTANCES:
        for phi_A, phi_B in PHASE_PAIRS:
            sub = in_pool[
                (in_pool.band == band)
                & (in_pool.distance == d_label)
                & (in_pool.phase_A == phi_A)
                & (in_pool.phase_B == phi_B)
            ]
            if sub.empty:
                continue
            cohort_rows.append({
                "band": band,
                "distance": d_label,
                "phase_A": phi_A,
                "phase_B": phi_B,
                "n_patients": int(sub.shape[0]),
                "median_d_obs": float(sub.d_obs.median()),
                "median_null_median": float(sub.null_median.median()),
                "median_null_iqr": float(sub.null_iqr.median()),
                "median_ratio": float(sub.ratio_obs_over_nullmed.median()),
                "median_z": float(sub.z.median()),
            })
cohort_df = pd.DataFrame(cohort_rows)
cohort_df.to_csv(OUT_BASE / "cohort_scale_summary.csv", index=False)

print("=== Pat_06 — absolute scale per band × distance × pair ===")
pat06 = per_patient_df[per_patient_df.patient == "Pat_06"].copy()
pat06["d_obs"] = pat06["d_obs"].round(4)
pat06["null_median"] = pat06["null_median"].round(4)
pat06["null_iqr"] = pat06["null_iqr"].round(5)
pat06["ratio_obs_over_nullmed"] = pat06["ratio_obs_over_nullmed"].round(2)
pat06["z"] = pat06["z"].round(1)
print(pat06[["band", "distance", "phase_A", "phase_B",
             "d_obs", "null_median", "null_iqr",
             "ratio_obs_over_nullmed", "z"]].to_string(index=False))

print(f"\n=== Cohort scale summary (n={N_COHORT}; Pat_03 included) ===")
cohort_short = cohort_df.copy()
cohort_short["median_d_obs"] = cohort_short["median_d_obs"].round(4)
cohort_short["median_null_median"] = cohort_short["median_null_median"].round(4)
cohort_short["median_null_iqr"] = cohort_short["median_null_iqr"].round(5)
cohort_short["median_ratio"] = cohort_short["median_ratio"].round(2)
cohort_short["median_z"] = cohort_short["median_z"].round(1)
print(cohort_short.to_string(index=False))
print(f"\nwrote {OUT_BASE / 'per_patient_with_scale.csv'}")
print(f"wrote {OUT_BASE / 'cohort_scale_summary.csv'}")
print(f"wrote per-patient audit_with_scale.csv files in {OUT_BASE}/Pat_*/")
