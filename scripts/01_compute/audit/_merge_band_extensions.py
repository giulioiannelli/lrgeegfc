#!/usr/bin/env python3
"""Merge two-band δ/θ outputs from audit_63 / audit_65 into the
canonical 5-band CSVs by concatenating with the .bak_alpha_beta_lowgamma
backups taken before the new runs.

Usage:
    python _merge_band_extensions.py audit_63
    python _merge_band_extensions.py audit_65
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

AUDITS = {
    "audit_63": {
        "dir": ROOT / "data" / "audit"
               / "matched_strength_surrogate_split_baseline",
        "files": ["per_patient_per_band.csv", "cohort_summary.csv"],
    },
    "audit_65": {
        "dir": ROOT / "data" / "audit"
               / "kc_matched_strength_surrogate",
        "files": ["per_patient_per_band.csv", "cohort_summary.csv",
                  "joint_signature.csv"],
    },
}


def merge_one(name: str) -> None:
    spec = AUDITS[name]
    out_dir = spec["dir"]
    for fn in spec["files"]:
        new_path = out_dir / fn
        bak_path = out_dir / f"{fn}.bak_alpha_beta_lowgamma"
        if not new_path.exists():
            print(f"[merge] SKIP missing new file: {new_path}")
            continue
        if not bak_path.exists():
            print(f"[merge] SKIP missing backup: {bak_path}")
            continue
        new_df = pd.read_csv(new_path)
        bak_df = pd.read_csv(bak_path)
        if "band" not in new_df.columns or "band" not in bak_df.columns:
            print(f"[merge] {fn}: no band column, concat as-is")
            merged = pd.concat([bak_df, new_df], ignore_index=True)
        else:
            new_bands = set(new_df.band.unique())
            bak_bands = set(bak_df.band.unique())
            overlap = new_bands & bak_bands
            if overlap:
                print(f"[merge] {fn}: WARNING bands overlap: {overlap} — "
                      f"new run will override backup rows")
                bak_df = bak_df[~bak_df.band.isin(overlap)]
            merged = pd.concat([bak_df, new_df], ignore_index=True)
            # Reorder bands canonically when possible
            cat = pd.CategoricalDtype(
                [b for b in BAND_ORDER if b in merged.band.unique()],
                ordered=True)
            merged["band"] = merged["band"].astype(cat)
            sort_cols = ["band"]
            if "patient" in merged.columns:
                sort_cols.append("patient")
            if "lam" in merged.columns:
                sort_cols.append("lam")
            elif "lambda" in merged.columns:
                sort_cols.append("lambda")
            merged = merged.sort_values(sort_cols).reset_index(drop=True)
            merged["band"] = merged["band"].astype(str)
        merged_path = out_dir / fn
        merged.to_csv(merged_path, index=False)
        print(f"[merge] {merged_path}: {len(merged)} rows "
              f"({sorted(merged.band.unique()) if 'band' in merged.columns else 'no band'})")


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else None
    if target not in AUDITS:
        print(f"usage: {sys.argv[0]} {{{'|'.join(AUDITS)}}}")
        sys.exit(1)
    merge_one(target)
