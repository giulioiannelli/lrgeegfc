"""Channel → Desikan-Killiany region mapping from per-patient implant CSVs.

The implant CSV stores DK labels in a single packed column with format
`{region1, w1, region2, w2, ..., PTD, value}` where weights sum to ~100.
The dominant region is the one with the highest non-PTD weight; PTD is
the proximal-tissue-density score and is excluded from the argmax.

Channel labels in `channel_labels.csv` use the form `"F 1,G2"` (with a
space between contact prefix and number, and a `G2` reference suffix).
Implant labels use `"F1"`. Normalisation strips whitespace and the
`,Gn` suffix.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import SEEG_DATAPATH

_LABEL_SUFFIX = re.compile(r",G\d+$")


def _normalise_label(label: str) -> str:
    """`"F 1,G2"` → `"F1"`. Idempotent on already-normalised labels."""
    s = str(label).strip()
    s = _LABEL_SUFFIX.sub("", s)
    s = s.replace(" ", "")
    return s


def _parse_dk_cell(cell: str) -> tuple[str, float]:
    """Parse `" ctx-rh-medialorbitofrontal,68.0,Unk,31.0,Wm,1.0,PTD, 0.98"`
    into `(dominant_region, weight)` where weight excludes PTD."""
    if not isinstance(cell, str) or not cell.strip():
        return "unknown", float("nan")
    parts = [p.strip() for p in cell.split(",")]
    pairs: list[tuple[str, float]] = []
    i = 0
    while i + 1 < len(parts):
        name = parts[i]
        try:
            w = float(parts[i + 1])
        except ValueError:
            i += 1
            continue
        if name.lower() != "ptd":
            pairs.append((name, w))
        i += 2
    if not pairs:
        return "unknown", float("nan")
    pairs.sort(key=lambda kv: -kv[1])
    return pairs[0]


_REGION_GROUPS = {
    # Coarse anatomical lobes/networks; falls back to white matter / unknown.
    "frontal": ("ctx-lh-superiorfrontal", "ctx-lh-rostralmiddlefrontal",
                "ctx-lh-caudalmiddlefrontal", "ctx-lh-medialorbitofrontal",
                "ctx-lh-lateralorbitofrontal", "ctx-lh-parsopercularis",
                "ctx-lh-parstriangularis", "ctx-lh-parsorbitalis",
                "ctx-lh-frontalpole", "ctx-lh-precentral",
                "ctx-rh-superiorfrontal", "ctx-rh-rostralmiddlefrontal",
                "ctx-rh-caudalmiddlefrontal", "ctx-rh-medialorbitofrontal",
                "ctx-rh-lateralorbitofrontal", "ctx-rh-parsopercularis",
                "ctx-rh-parstriangularis", "ctx-rh-parsorbitalis",
                "ctx-rh-frontalpole", "ctx-rh-precentral"),
    "parietal": ("ctx-lh-superiorparietal", "ctx-lh-inferiorparietal",
                  "ctx-lh-supramarginal", "ctx-lh-precuneus",
                  "ctx-lh-postcentral", "ctx-lh-paracentral",
                  "ctx-rh-superiorparietal", "ctx-rh-inferiorparietal",
                  "ctx-rh-supramarginal", "ctx-rh-precuneus",
                  "ctx-rh-postcentral", "ctx-rh-paracentral"),
    "temporal": ("ctx-lh-superiortemporal", "ctx-lh-middletemporal",
                  "ctx-lh-inferiortemporal", "ctx-lh-temporalpole",
                  "ctx-lh-transversetemporal", "ctx-lh-bankssts",
                  "ctx-lh-fusiform", "ctx-lh-entorhinal",
                  "ctx-lh-parahippocampal",
                  "ctx-rh-superiortemporal", "ctx-rh-middletemporal",
                  "ctx-rh-inferiortemporal", "ctx-rh-temporalpole",
                  "ctx-rh-transversetemporal", "ctx-rh-bankssts",
                  "ctx-rh-fusiform", "ctx-rh-entorhinal",
                  "ctx-rh-parahippocampal"),
    "occipital": ("ctx-lh-lateraloccipital", "ctx-lh-cuneus",
                   "ctx-lh-pericalcarine", "ctx-lh-lingual",
                   "ctx-rh-lateraloccipital", "ctx-rh-cuneus",
                   "ctx-rh-pericalcarine", "ctx-rh-lingual"),
    "cingulate": ("ctx-lh-rostralanteriorcingulate",
                   "ctx-lh-caudalanteriorcingulate",
                   "ctx-lh-posteriorcingulate", "ctx-lh-isthmuscingulate",
                   "ctx-rh-rostralanteriorcingulate",
                   "ctx-rh-caudalanteriorcingulate",
                   "ctx-rh-posteriorcingulate", "ctx-rh-isthmuscingulate"),
    "insula": ("ctx-lh-insula", "ctx-rh-insula"),
    "subcortical": ("Left-Hippocampus", "Right-Hippocampus",
                     "Left-Amygdala", "Right-Amygdala",
                     "Left-Thalamus-Proper", "Right-Thalamus-Proper",
                     "Left-Caudate", "Right-Caudate",
                     "Left-Putamen", "Right-Putamen",
                     "Left-Pallidum", "Right-Pallidum",
                     "Left-Accumbens-area", "Right-Accumbens-area"),
    "white_matter": ("Wm", "Left-Cerebral-White-Matter",
                      "Right-Cerebral-White-Matter"),
}


def _coarse_lobe(region: str) -> str:
    if not region or region == "unknown":
        return "unknown"
    for lobe, members in _REGION_GROUPS.items():
        if region in members:
            return lobe
    if region.startswith("ctx-"):
        return "other_cortex"
    return "other"


def _hemisphere(region: str) -> str:
    if "-lh-" in region or region.startswith("Left-"):
        return "L"
    if "-rh-" in region or region.startswith("Right-"):
        return "R"
    return "?"


def load_channel_regions(patient: str,
                            data_root: Path = SEEG_DATAPATH) -> pd.DataFrame:
    """Return a DataFrame indexed by FC-matrix row order with columns:
        ``label_raw``, ``label`` (normalised), ``region`` (DK name),
        ``region_weight``, ``lobe`` (coarse), ``hemisphere``, ``x``, ``y``, ``z``.

    Row order matches `channel_labels.csv`, which is the order the FC
    matrix and LRG dendrogram leaves use.
    """
    pat_dir = data_root / patient
    patnum = int(patient.split("_")[-1])
    ch_csv = pat_dir / "channel_labels.csv"
    impl_csv = pat_dir / f"implant_pat_{patnum:02d}.csv"
    if not impl_csv.exists():
        impl_csv = pat_dir / f"Implant_pat_{patnum:02d}.csv"

    ch = pd.read_csv(ch_csv)
    impl = pd.read_csv(impl_csv)

    # find the DK column (some patients have "Desikan-K" truncated)
    dk_col = next((c for c in impl.columns if "esikan" in c), None)
    if dk_col is None:
        raise ValueError(f"{patient}: no Desikan column found in {impl_csv.name}")

    impl = impl.copy()
    impl["_norm"] = impl["label"].map(_normalise_label)
    parsed = impl[dk_col].map(_parse_dk_cell)
    impl["region"] = [p[0] for p in parsed]
    impl["region_weight"] = [p[1] for p in parsed]

    ch = ch.copy()
    ch["label_raw"] = ch["label"]
    ch["_norm"] = ch["label"].map(_normalise_label)

    merged = ch.merge(impl[["_norm", "region", "region_weight", "x", "y", "z"]],
                      on="_norm", how="left")
    merged["region"] = merged["region"].fillna("unknown")
    merged["lobe"] = merged["region"].map(_coarse_lobe)
    merged["hemisphere"] = merged["region"].map(_hemisphere)
    merged = merged.rename(columns={"_norm": "label"})
    cols = ["label_raw", "label", "region", "region_weight",
             "lobe", "hemisphere", "x", "y", "z"]
    return merged[cols].reset_index(drop=True)
