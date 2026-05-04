"""Q4: d_S vs d_F triangle-persistence contrast per band.

Reads per-patient audit.csv, extracts T_d under each distance, and:

- Tabulates T_d^(d_S) vs T_d^(d_F) per patient × band (n=10 patients; Pat_03 included).
- Computes Spearman rho and Pearson r across patients within each band.
- Flags bands where |T_d^(d_S)| ≪ |T_d^(d_F)| (amplitude-drift candidates)
  vs bands where T_d^(d_S) and T_d^(d_F) agree in sign and rough magnitude
  (structural-shift candidates).

Decision rule for the per-band flag:
- "amplitude-drift" if median |T_d^F| > 2 × median |T_d^S| AND median |T_d^S| < 0.05
- "structural-shift" if (sign agreement >= 7/10 patients) AND (Spearman rho > 0.5)
- "ambiguous" otherwise
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import DATA_ROOT

N_COHORT = len(PATIENTS_4PHASE)  # single source of truth — never hardcode

OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"

rows = []
for p in PATIENTS_4PHASE:
    csv = OUT_BASE / p / "audit.csv"
    if not csv.exists():
        continue
    df = pd.read_csv(csv)
    tri = df[df.flag == "T_d"]
    for band in BRAIN_BANDS_NAMES:
        for d_label in ("S", "F", "P"):
            sub = tri[(tri.band == band) & (tri.distance == d_label)]
            if sub.empty:
                continue
            rows.append({"patient": p, "band": band, "distance": d_label, "T_d": float(sub.iloc[0].d_obs)})

td_long = pd.DataFrame(rows)
td_wide = td_long.pivot_table(index=["patient", "band"], columns="distance", values="T_d").reset_index()
td_wide.to_csv(OUT_BASE / "Td_per_patient_per_band.csv", index=False)

in_pool = td_wide.copy()  # n=10; Pat_03 included (Z dimensionless w.r.t. fs)

# Per-band contrast & flags
band_rows = []
for band in BRAIN_BANDS_NAMES:
    sub = in_pool[in_pool.band == band].dropna(subset=["S", "F"])
    if sub.shape[0] < 5:
        continue
    Td_S = sub["S"].to_numpy()
    Td_F = sub["F"].to_numpy()
    Td_P = sub["P"].to_numpy() if "P" in sub.columns else np.full_like(Td_S, np.nan)

    rho_SF, _ = spearmanr(Td_S, Td_F)
    r_SF, _ = pearsonr(Td_S, Td_F)
    sign_agree = int((np.sign(Td_S) == np.sign(Td_F)).sum())
    med_abs_S = float(np.median(np.abs(Td_S)))
    med_abs_F = float(np.median(np.abs(Td_F)))
    med_S = float(np.median(Td_S))
    med_F = float(np.median(Td_F))

    if med_abs_F > 2 * med_abs_S and med_abs_S < 0.05:
        flag = "amplitude-drift"
    elif sign_agree >= 7 and rho_SF > 0.5:
        flag = "structural-shift"
    else:
        flag = "ambiguous"

    band_rows.append({
        "band": band,
        "n_patients": int(sub.shape[0]),
        "median_T_d_S": round(med_S, 4),
        "median_T_d_F": round(med_F, 4),
        "median_abs_T_d_S": round(med_abs_S, 4),
        "median_abs_T_d_F": round(med_abs_F, 4),
        "sign_agree_S_F": f"{sign_agree}/{sub.shape[0]}",
        "spearman_rho_S_F": round(rho_SF, 3),
        "pearson_r_S_F": round(r_SF, 3),
        "flag": flag,
    })

contrast_df = pd.DataFrame(band_rows)
contrast_df.to_csv(OUT_BASE / "Td_dS_vs_dF_band_contrast.csv", index=False)

print(f"=== T_d per patient per band (n={N_COHORT}; columns: P, S, F) ===")
in_pool_show = in_pool.copy()
for c in ("P", "S", "F"):
    if c in in_pool_show.columns:
        in_pool_show[c] = in_pool_show[c].round(4)
print(in_pool_show.sort_values(["band", "patient"]).to_string(index=False))

print(f"\n=== Per-band contrast: d_S vs d_F (n={N_COHORT}) ===")
print(contrast_df.to_string(index=False))

print(f"\nwrote {OUT_BASE / 'Td_per_patient_per_band.csv'}")
print(f"wrote {OUT_BASE / 'Td_dS_vs_dF_band_contrast.csv'}")
