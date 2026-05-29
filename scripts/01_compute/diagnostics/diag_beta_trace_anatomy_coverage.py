#!/usr/bin/env python
"""β-trace anatomy coverage per patient — explains Pat_10 / Pat_14 / Pat_15
non-trace at the cophenet probe via implant sampling, or rules it out.

The β cophenet trace is anatomically localized to 7 Desikan-Killiany regions
(`ANATOMY_LEDGER.md` §β cophenet anatomy, audit_71 2026-05-19):

  - ctx-lh-isthmuscingulate
  - ctx-lh-superiorfrontal
  - ctx-rh-insula
  - ctx-lh-parahippocampal
  - ctx-rh-postcentral
  - ctx-lh-entorhinal
  - ctx-rh-rostralanteriorcingulate

Hypothesis (mechanistic): patients whose stereoEEG implant samples this
network sparsely or asymmetrically should fail to carry the cohort β
cophenet trace at the per-patient level. Pat_10 (`ρ_split^coph` = −0.091),
Pat_14 (−0.049), and Pat_15 (+0.083, n.s.) are the three LRG-anti or
non-significant patients at β cophenet at full cohort.

This diagnostic counts, per patient:
- total contacts in the 7-region β-trace network
- per-region count (sparse coverage reveals missing nodes)
- hemispheric balance (4 LH regions + 3 RH regions in the network →
  balanced sampling should show similar LH/RH counts; asymmetry can
  obscure the trace)
- total contacts overall (denominator for normalized coverage)

Comparison groups:
- LRG-anti / non-sig at β cophenet: Pat_10, Pat_14, Pat_15
- LRG-trace at β cophenet (n_above = 7/10): Pat_02, 03, 05, 06, 07, 08, 13

Outputs:
- ``data/audit/beta_trace_anatomy_coverage/per_patient.csv`` — wide table.
- ``data/audit/beta_trace_anatomy_coverage/per_patient_per_region.csv`` —
  long table (patient × region).
- ``data/audit/beta_trace_anatomy_coverage/README.md`` — verdict + reading.

Runtime: < 30 seconds.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.io.regions import load_channel_regions


OUT_ROOT = Path("data/audit/beta_trace_anatomy_coverage")

BETA_COPHENET_REGIONS = [
    "ctx-lh-isthmuscingulate",
    "ctx-lh-superiorfrontal",
    "ctx-rh-insula",
    "ctx-lh-parahippocampal",
    "ctx-rh-postcentral",
    "ctx-lh-entorhinal",
    "ctx-rh-rostralanteriorcingulate",
]

LRG_ANTI_OR_NS = {"Pat_10", "Pat_14", "Pat_15"}
LRG_TRACE = {"Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08", "Pat_13"}


def _patient_coverage(patient: str) -> tuple[dict, pd.DataFrame]:
    regs = load_channel_regions(patient)
    n_total = len(regs)
    in_net = regs["region"].isin(BETA_COPHENET_REGIONS)
    n_in_net = int(in_net.sum())

    per_region = {}
    long_rows = []
    for r in BETA_COPHENET_REGIONS:
        c = int((regs["region"] == r).sum())
        per_region[r] = c
        long_rows.append({"patient": patient, "region": r, "n_contacts": c})

    lh_regs = [r for r in BETA_COPHENET_REGIONS if "-lh-" in r]
    rh_regs = [r for r in BETA_COPHENET_REGIONS if "-rh-" in r]
    n_lh = int(regs["region"].isin(lh_regs).sum())
    n_rh = int(regs["region"].isin(rh_regs).sum())

    summary = {
        "patient": patient,
        "group": "LRG-anti/ns" if patient in LRG_ANTI_OR_NS else "LRG-trace",
        "n_contacts_total": n_total,
        "n_in_beta_net": n_in_net,
        "frac_in_beta_net": n_in_net / n_total if n_total else np.nan,
        "n_LH_in_net": n_lh,
        "n_RH_in_net": n_rh,
        "n_regions_covered": int(sum(1 for c in per_region.values() if c > 0)),
        **{f"n_{r.split('-')[-1]}_{r.split('-')[1]}": c
            for r, c in per_region.items()},
    }
    return summary, pd.DataFrame(long_rows)


def _write_summary(df: pd.DataFrame, out_path: Path) -> None:
    by_group = df.groupby("group")[
        ["n_in_beta_net", "frac_in_beta_net", "n_LH_in_net",
         "n_RH_in_net", "n_regions_covered"]
    ].agg(["median", "min", "max"])

    md = f"""---
name: beta-trace-anatomy-coverage
era: IMCOH_ABS_COHORT_N10
status: current
kind: diagnostic-report
date: 2026-05-29
band: beta
probe: cophenet
network: 7 DK regions from ANATOMY_LEDGER β cophenet
---

# β-trace anatomy coverage per patient

Tests whether Pat_10 / Pat_14 / Pat_15 — the three patients not in the
β cophenet trace direction at full cohort — sample the β-trace 7-region
network differently from the 7 trace-carrying patients.

## Per-patient table

| Patient | Group | Total contacts | In β-net | Fraction | LH/RH in net | Regions covered (of 7) |
|---|---|---|---|---|---|---|
"""
    for _, row in df.iterrows():
        md += (f"| {row['patient']} | {row['group']} | "
               f"{row['n_contacts_total']} | {row['n_in_beta_net']} | "
               f"{row['frac_in_beta_net']:.3f} | "
               f"{row['n_LH_in_net']}/{row['n_RH_in_net']} | "
               f"{row['n_regions_covered']} |\n")

    md += "\n## Group comparison (median / min / max)\n\n"
    md += "| Group | n_in_net (median) | frac_in_net (median) | n_LH (median) | n_RH (median) | n_regions covered (median) |\n"
    md += "|---|---|---|---|---|---|\n"
    for grp in ["LRG-trace", "LRG-anti/ns"]:
        sub = df[df["group"] == grp]
        md += (f"| {grp} (n={len(sub)}) | "
               f"{int(sub['n_in_beta_net'].median())} "
               f"[{int(sub['n_in_beta_net'].min())}–"
               f"{int(sub['n_in_beta_net'].max())}] | "
               f"{sub['frac_in_beta_net'].median():.3f} | "
               f"{int(sub['n_LH_in_net'].median())} | "
               f"{int(sub['n_RH_in_net'].median())} | "
               f"{int(sub['n_regions_covered'].median())} |\n")

    md += """
## How to read

- **n_in_beta_net** is the number of recording contacts a patient has
  inside the 7-region β-trace network. Higher = better sampled.
- **frac_in_beta_net** normalizes by total contacts.
- **LH/RH** balance — the network has 4 left-hemi regions and 3 right-hemi
  regions. Severe asymmetry could mute the trace if the patient samples
  only one side.
- **Regions covered** — how many of the 7 regions have at least one
  contact. Sparse coverage of the network = sparse sampling of the trace.

## Verdict criteria

- **Mechanism confirmed** if LRG-anti/ns group median `n_in_beta_net`
  is < 50% of LRG-trace group median, OR median `n_regions_covered` is
  ≤ 3 (vs LRG-trace ≥ 5), OR severe LH/RH asymmetry (single hemi 0
  contacts).
- **Mechanism rejected** if both groups are comparable across all four
  metrics. Move to direction (2) within-baseline cophenet stability.

## Caveats

This is a coverage diagnostic, not a causal test. A patient with low
coverage in this network might still carry the trace if the actual β
reorganization extends beyond the 7-region network (e.g., to subcortical
nodes not in the DK parcellation). Conversely, a patient with high
coverage but no trace might have biological reasons unrelated to
sampling.

Pat_15 has 0 epi contacts and right-hemisphere-only implant — this is
already documented biology (`feedback_no_patient_dropout.md`) and the
coverage diagnostic is expected to confirm the LH-region sparse sampling.
"""
    out_path.write_text(md)


def main() -> None:
    rows = []
    long_dfs = []
    for pat in PATIENTS_4PHASE:
        print(f"[diag] {pat} ...", end="", flush=True)
        try:
            summary, long_df = _patient_coverage(pat)
        except Exception as e:
            print(f" FAIL ({e})")
            continue
        rows.append(summary)
        long_dfs.append(long_df)
        print(f" n_in_net={summary['n_in_beta_net']} "
              f"LH/RH={summary['n_LH_in_net']}/{summary['n_RH_in_net']} "
              f"regions_covered={summary['n_regions_covered']}/7")

    df = pd.DataFrame(rows)
    long_df = pd.concat(long_dfs, ignore_index=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_ROOT / "per_patient.csv", index=False)
    long_df.to_csv(OUT_ROOT / "per_patient_per_region.csv", index=False)
    _write_summary(df, OUT_ROOT / "README.md")

    print("\n[diag] group comparison:")
    for grp in ["LRG-trace", "LRG-anti/ns"]:
        sub = df[df["group"] == grp]
        print(f"  {grp:<12s} (n={len(sub)}): "
              f"n_in_net median={int(sub['n_in_beta_net'].median())} "
              f"[{int(sub['n_in_beta_net'].min())}-"
              f"{int(sub['n_in_beta_net'].max())}], "
              f"regions covered median={int(sub['n_regions_covered'].median())}/7, "
              f"LH/RH median={int(sub['n_LH_in_net'].median())}/"
              f"{int(sub['n_RH_in_net'].median())}")
    print(f"\n[diag] outputs → {OUT_ROOT}/")


if __name__ == "__main__":
    main()
