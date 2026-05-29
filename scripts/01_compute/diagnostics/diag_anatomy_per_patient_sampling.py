#!/usr/bin/env python
"""Per-patient sampling check across ALL locked anatomy claims.

Question (raised 2026-05-29 after β coverage diagnostic): the locked
ANATOMY_LEDGER.md flags up to 11 named DK regions as "strong localized"
per band, but the per-band coverage diagnostic on β revealed that several
regions are sampled by ≤2 patients (some by only LRG-anti patients).

This script applies the same coverage check to *every* locked anatomy
network and flags:
- SINGLE-PATIENT  — region sampled by 1 patient (anatomy enrichment is
  effectively a single-patient finding).
- THIN            — region sampled by 2 patients.
- ROBUST          — region sampled by ≥3 patients.

For each band's network we also report which patients sample each region
so that a follow-up audit can cross-reference with the per-band per-patient
trace direction (e.g., β: was a region sampled only by LRG-anti patients?).

Outputs:
- ``data/audit/anatomy_per_patient_sampling/per_band_per_region.csv`` —
  long table (band, probe, region, n_patients_sampling, patient_list).
- ``data/audit/anatomy_per_patient_sampling/per_band_summary.csv`` —
  one row per (band, probe) with counts of SINGLE / THIN / ROBUST regions.
- ``data/audit/anatomy_per_patient_sampling/README.md`` — verdict + table
  identifying regions at risk of overselling.

Sources for region lists: `.agents/preprint/locked/ANATOMY_LEDGER.md`
(2026-05-19 lock + 2026-05-28 audit edits).

Runtime: ~10 seconds.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.io.regions import load_channel_regions


OUT_ROOT = Path("data/audit/anatomy_per_patient_sampling")


ANATOMY_NETWORKS = {
    ("beta", "cophenet"): [
        "ctx-lh-isthmuscingulate",
        "ctx-lh-superiorfrontal",
        "ctx-rh-insula",
        "ctx-lh-parahippocampal",
        "ctx-rh-postcentral",
        "ctx-lh-entorhinal",
        "ctx-rh-rostralanteriorcingulate",
    ],
    ("beta", "grassmann"): [
        "Hip",
        "ctx-lh-insula",
        "ctx-lh-middletemporal",
        "ctx-lh-lateralorbitofrontal",
        "ctx-rh-medialorbitofrontal",
        "ctx-lh-superiortemporal",
        "ctx-rh-rostralmiddlefrontal",
    ],
    ("alpha", "cophenet"): [
        "ctx-lh-caudalanteriorcingulate",
        "ctx-lh-parahippocampal",
        "ctx-lh-rostralanteriorcingulate",
        "ctx-rh-medialorbitofrontal",
        "ctx-rh-caudalmiddlefrontal",
        "ctx-rh-postcentral",
        "ctx-lh-caudalmiddlefrontal",
        "ctx-rh-caudalanteriorcingulate",
        "ctx-rh-precuneus",
        "ctx-rh-superiorparietal",
        "ctx-rh-posteriorcingulate",
    ],
    ("low_gamma", "grassmann"): [
        "ctx-lh-lateraloccipital",
        "ctx-lh-middletemporal",
        "ctx-lh-rostralmiddlefrontal",
        "ctx-lh-superiortemporal",
        "ctx-rh-medialorbitofrontal",
        "ctx-rh-parstriangularis",
        "ctx-lh-cuneus",
    ],
    ("delta", "grassmann_full"): [
        "ctx-lh-inferiortemporal",
        "ctx-lh-inferiorparietal",
        "ctx-rh-parstriangularis",
        "ctx-lh-superiortemporal",
    ],
    ("delta", "grassmann_epiX"): [
        "ctx-lh-superiorparietal",
        "ctx-rh-rostralmiddlefrontal",
        "ctx-lh-superiorfrontal",
    ],
}


def _flag(n_patients: int) -> str:
    if n_patients == 0:
        return "UNSAMPLED"
    if n_patients == 1:
        return "SINGLE-PATIENT"
    if n_patients == 2:
        return "THIN"
    return "ROBUST"


def main() -> None:
    pat_regions = {pat: load_channel_regions(pat) for pat in PATIENTS_4PHASE}

    long_rows = []
    summary_rows = []
    for (band, probe), regions in ANATOMY_NETWORKS.items():
        flag_counts = {"SINGLE-PATIENT": 0, "THIN": 0, "ROBUST": 0,
                       "UNSAMPLED": 0}
        flagged_regions: list[str] = []
        for region in regions:
            sampling = []
            for pat, df in pat_regions.items():
                n = int((df["region"] == region).sum())
                if n > 0:
                    sampling.append((pat, n))
            n_patients = len(sampling)
            flag = _flag(n_patients)
            flag_counts[flag] += 1
            long_rows.append({
                "band": band,
                "probe": probe,
                "region": region,
                "n_patients_sampling": n_patients,
                "patient_list": "; ".join(
                    f"{p}({n})" for p, n in sampling) or "—",
                "flag": flag,
            })
            if flag in {"SINGLE-PATIENT", "THIN", "UNSAMPLED"}:
                flagged_regions.append(f"{region} [{flag}]")

        summary_rows.append({
            "band": band,
            "probe": probe,
            "n_regions": len(regions),
            "n_robust": flag_counts["ROBUST"],
            "n_thin": flag_counts["THIN"],
            "n_single_patient": flag_counts["SINGLE-PATIENT"],
            "n_unsampled": flag_counts["UNSAMPLED"],
            "frac_robust": flag_counts["ROBUST"] / len(regions),
            "flagged_regions": "; ".join(flagged_regions) or "—",
        })

    long_df = pd.DataFrame(long_rows)
    summary_df = pd.DataFrame(summary_rows)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(OUT_ROOT / "per_band_per_region.csv", index=False)
    summary_df.to_csv(OUT_ROOT / "per_band_summary.csv", index=False)

    print(f"{'band':<10s} {'probe':<16s} {'n':>3s} {'robust':>7s} "
          f"{'thin':>5s} {'single':>7s} {'unsamp':>7s}")
    print("-" * 70)
    for _, r in summary_df.iterrows():
        print(f"{r['band']:<10s} {r['probe']:<16s} {r['n_regions']:>3d} "
              f"{r['n_robust']:>7d} {r['n_thin']:>5d} "
              f"{r['n_single_patient']:>7d} {r['n_unsampled']:>7d}")

    print("\nFlagged regions per (band, probe):")
    for _, r in summary_df.iterrows():
        if r["flagged_regions"] != "—":
            print(f"  {r['band']:<10s} {r['probe']:<16s}: {r['flagged_regions']}")

    print(f"\nOutputs → {OUT_ROOT}/")


if __name__ == "__main__":
    main()
