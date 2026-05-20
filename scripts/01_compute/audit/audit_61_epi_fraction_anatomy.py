#!/usr/bin/env python3
"""Audit 61 — descriptive epileptic-zone fractions for §5.5 region pools.

C4 of the 2026-05-10 coding-agent batch. Five region pools (one per band
cell flagged in §5.5):

    band       region               cohort pool
    low_gamma  ctx-lh-fusiform      38   (Bonferroni-survived at m=48)
    beta       Hip                  27   (uncorrected; named 'Hippocampus' in C4 spec)
    beta       ctx-lh-fusiform      38   (uncorrected)
    beta       ctx-lh-superiortemporal  53 (uncorrected)
    alpha      ctx-lh-parsopercularis   26 (uncorrected suggestive)

For each region, compute the fraction of contacts in the cohort pool
flagged epileptic (red font in the patient implant XLSX, surfaced via
`load_epileptic_nodes`). Descriptive only — no hypergeometric, no
Bonferroni expansion, no fold into TARR (per
`feedback_epilepsy_not_trace_locked.md` and `16_epilepsy_pointer.md`).

Cohort pool definition matches §5.5 (audit_56_anatomy_distribution.py
load_pool): per-patient `implant_pat_NN.csv`, first comma-separated
Desikan-Killany token per row, drop Wm/Unk.

Output:
    data/audit/epi_fraction_anatomy/region_epi_fractions.csv
    data/audit/epi_fraction_anatomy/per_contact_table.csv
    data/audit/epi_fraction_anatomy/README.md
"""
from __future__ import annotations

import pandas as pd

from lrg_eegfc.utils.io import load_channel_labels, load_epileptic_nodes
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

TARGETS = [
    ("low_gamma", "ctx-lh-fusiform"),
    ("beta",      "Hip"),
    ("beta",      "ctx-lh-fusiform"),
    ("beta",      "ctx-lh-superiortemporal"),
    ("alpha",     "ctx-lh-parsopercularis"),
]

OUT_DIR = ROOT / "data" / "audit" / "epi_fraction_anatomy"


def load_pool() -> pd.DataFrame:
    rows = []
    for pat in COHORT:
        impl_path = (ROOT / "data" / "raw" / "stereoeeg_patients" / pat
                     / f"implant_pat_{pat[-2:]}.csv")
        impl = pd.read_csv(impl_path)
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        epi_set = set(load_epileptic_nodes(pat))
        labels = [str(x).strip() for x in impl["label"]]
        for label, dk in zip(labels, impl[dk_col]):
            region = str(dk).split(",")[0].strip()
            rows.append({
                "patient": pat,
                "label": label,
                "region": region,
                "is_epi": label in epi_set,
            })
    pool = pd.DataFrame(rows)
    pool = pool[~pool.region.isin(["Wm", "Unk"])].reset_index(drop=True)
    return pool


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pool = load_pool()
    pool.to_csv(OUT_DIR / "per_contact_table.csv", index=False)

    out_rows = []
    for band, region in TARGETS:
        sub = pool[pool.region == region]
        n_contacts = int(len(sub))
        n_epi = int(sub.is_epi.sum())
        n_pat_total = int(sub.patient.nunique())
        n_pat_epi = int(sub.loc[sub.is_epi, "patient"].nunique())
        frac = n_epi / n_contacts if n_contacts else float("nan")
        out_rows.append({
            "band": band,
            "region": region,
            "n_contacts": n_contacts,
            "n_epi": n_epi,
            "frac_epi": round(frac, 4) if n_contacts else None,
            "n_patients_in_pool": n_pat_total,
            "n_patients_with_epi_contact": n_pat_epi,
        })

    df = pd.DataFrame(out_rows)
    df.to_csv(OUT_DIR / "region_epi_fractions.csv", index=False)

    print(df.to_string(index=False))

    readme_lines = [
        "---",
        "name: epi_fraction_anatomy",
        "scope: descriptive_disclosure_for_section_5_5",
        "date: 2026-05-10",
        "status: descriptive",
        "---",
        "",
        "# Epileptic-zone fractions for §5.5 region pools",
        "",
        "Descriptive disclosure only. The cohort pool of each region is",
        "tallied against the patient implant epileptic-node flag (red font",
        "in the implant XLSX). No hypergeometric test, no Bonferroni",
        "expansion, no integration with the §5.5 anatomy enrichment.",
        "",
        "## Numbers",
        "",
        "| band | region | n_contacts | n_epi | frac_epi | patients_in_pool | patients_with_epi |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in out_rows:
        readme_lines.append(
            f"| {r['band']} | {r['region']} | {r['n_contacts']} | {r['n_epi']} "
            f"| {r['frac_epi']:.3f} | {r['n_patients_in_pool']} "
            f"| {r['n_patients_with_epi_contact']} |"
        )
    readme_lines.extend([
        "",
        "## §5.5 disclosure-sentence drafts",
        "",
    ])
    for r in out_rows:
        readme_lines.append(
            f"- **{r['band']} / {r['region']}**: of the cohort's "
            f"{r['n_contacts']} contacts at this region, {r['n_epi']} "
            f"({r['frac_epi']*100:.1f}%) are flagged epileptic-zone in "
            f"the patient-level implant catalog. The trace-leaf "
            f"enrichment reported above is computed against the full "
            f"cohort cortical pool without epileptic-zone exclusion; "
            f"an exclusion-sensitivity analysis is owed."
        )
    readme_lines.extend([
        "",
        "## Provenance",
        "",
        "- Cohort: " + ", ".join(COHORT),
        "- Pool: per-patient `implant_pat_NN.csv` first Desikan-Killany "
        "token, dropping `Wm` / `Unk` (matches §5.5 / audit_56).",
        "- Epi flag: `load_epileptic_nodes` (red-font label in implant XLSX).",
        "- Build script: `scripts/01_compute/audit/audit_61_epi_fraction_anatomy.py`",
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(readme_lines))


if __name__ == "__main__":
    main()
