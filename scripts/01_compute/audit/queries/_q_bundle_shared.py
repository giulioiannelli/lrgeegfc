"""Shared header for the q_bundle figure family.

Split out of q_bundle_figures.py on 2026-05-29 (Phase 4-B split 3/7).
Holds the import block, output paths, phase-pair palette, the
`lookup_pair` helper, and the data loaders consumed by each
`q_fig{N}_*.py` script.

Replaces the previous `sys.path.insert(scripts/archive/2026-02_imcoh-dev-notes)`
injection with proper library imports.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: F401  (re-exported for fig modules)
import numpy as np  # noqa: F401
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    PATIENTS_4PHASE,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    NODE_SIZE_DEFAULT,
)
from lrg_eegfc.config.paths import DATA_ROOT, SEEG_DATAPATH
from lrg_eegfc.workflow.fc import load_fc_matrix  # noqa: F401
from lrg_eegfc.utils.probe import extract_probe_labels  # noqa: F401
from lrgsglib.plotlib import imshow_colorbar_caxdivider  # noqa: F401

# Layout + drawing helpers (promoted from figures_for_notes/_shared.py on
# 2026-05-29 to lrg_eegfc.visuals.network_layouts / network_drawing —
# Phase 4-B split 1/7).
from lrg_eegfc.visuals.network_layouts import compute_network_layout  # noqa: F401
from lrg_eegfc.visuals.network_drawing import (  # noqa: F401
    draw_network_edges,
    _probe_color_map,
)


def probe_sort_indices(channel_labels):
    """Return indices that sort channels by probe, then by contact number."""
    import numpy as np
    from lrg_eegfc.utils.probe import extract_probe_labels
    probes = extract_probe_labels(channel_labels)
    unique = sorted(set(probes))
    return np.argsort([unique.index(p) * 1000 + i for i, p in enumerate(probes)])


def probe_boundaries(sorted_probes):
    """Return boundary positions (first index of each new probe)."""
    return [i for i in range(1, len(sorted_probes))
            if sorted_probes[i] != sorted_probes[i - 1]]


# draw_probe_outlines lives in the library — re-export for split-fig callers.
from lrg_eegfc.visuals import draw_probe_outlines  # noqa: F401, E402


def load_channel_labels(patient: str) -> list[str]:
    """Load cleaned monopolar channel labels for *patient*.

    Inline copy of the helper that used to live in the archived
    `_shared.py` -- re-implementing it here lets us drop the
    `sys.path.insert(scripts/archive/...)` hack.
    """
    for ext in ("csv", "txt"):
        fpath = SEEG_DATAPATH / patient / f"channel_labels.{ext}"
        if not fpath.exists():
            continue
        labels: list[str] = []
        with open(fpath) as f:
            for i, line in enumerate(f):
                line = line.strip().strip('"')
                if not line:
                    continue
                if i == 0 and line.lower() == "label":
                    continue
                label = line.split(",")[0].strip().strip('"').replace(" ", "")
                if label:
                    labels.append(label)
        if labels:
            return labels
    raise FileNotFoundError(f"No channel_labels for {patient}")


PHASE_SHORT = {
    "rest_pre":   "RPre",
    "task_learn": "TL",
    "task_test":  "TT",
    "rest_post":  "RPost",
}

N_COHORT = len(PATIENTS_4PHASE)


def resolve_substrate() -> tuple[str, str, Path, Path]:
    """Parse CLI substrate arg + return (SUBSTRATE, SUFFIX, SRC, OUT)."""
    substrate = sys.argv[1] if len(sys.argv) > 1 else "imcoh_abs"
    assert substrate in ("imcoh_abs", "imcoh_sq"), \
        f"unknown substrate {substrate}"
    if substrate == "imcoh_abs":
        src = DATA_ROOT / "audit" / "raw_fc_phase_distance"
    else:
        src = DATA_ROOT / "audit" / "raw_fc_phase_distance_imcoh_sq"
    suffix = f"_{substrate}"
    out = ROOT / ".agents" / "writing-bundles" / "raw-fc"
    out.mkdir(parents=True, exist_ok=True)
    return substrate, suffix, src, out


PAIR_COLORS = {
    ("task_learn", "task_test"): "#1f77b4",  # blue
    ("task_test", "rest_post"):  "#d62728",  # red - the trace pair
    ("rest_pre",  "task_test"):  "#ff7f0e",  # orange
    ("rest_pre",  "rest_post"):  "#7f7f7f",  # gray - drift floor
}
PAIR_LABELS = {
    ("task_learn", "task_test"): "TL$\\leftrightarrow$TT",
    ("task_test", "rest_post"):  "TT$\\leftrightarrow$RPost",
    ("rest_pre",  "task_test"):  "RPre$\\leftrightarrow$TT",
    ("rest_pre",  "rest_post"):  "RPre$\\leftrightarrow$RPost",
}
# Distance keys + labels used across multiple bundle figures (was at
# module-top in fig6 section of the original q_bundle_figures.py and
# reused by figS1, figS2, figS3).
DISTANCE_KEYS = ["S", "P", "F"]
DISTANCE_LABEL = {
    "S": r"$d_S$  Spearman",
    "P": r"$d_P$  Pearson",
    "F": r"$d_F$  Frobenius",
}
PHASES_ORDER = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_X = {ph: i for i, ph in enumerate(PHASES_ORDER)}
PHASE_LABEL_SHORT = {"rest_pre": "RPre", "task_learn": "TL",
                     "task_test": "TT", "rest_post": "RPost"}
REF_LINES = [0.0, 0.25, 0.5, 0.75, 1.0]

PAIR_ORDER = [
    ("task_learn", "task_test"),
    ("task_test",  "rest_post"),
    ("rest_pre",   "task_test"),
    ("rest_pre",   "rest_post"),
]


def lookup_pair(df, band, A, B, distance="S"):
    """Return per-patient `d_obs` Series indexed by patient for a phase pair."""
    sub = df[
        (df.distance == distance)
        & (df.band == band)
        & (((df.phase_A == A) & (df.phase_B == B))
           | ((df.phase_A == B) & (df.phase_B == A)))
    ]
    return sub.groupby("patient").d_obs.first()


def load_bundle_data(src: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load (distances_long, td, contrast) DataFrames from the substrate src dir."""
    rows = []
    for p in PATIENTS_4PHASE:
        f = src / p / "distance_4phase.csv"
        if not f.exists():
            print(f"[{p}] missing distance_4phase.csv, skipping")
            continue
        rows.append(pd.read_csv(f))
    distances_long = pd.concat(rows, ignore_index=True)
    td = pd.read_csv(src / "Td_per_patient_per_band.csv")
    contrast = pd.read_csv(src / "Td_dS_vs_dF_band_contrast.csv").set_index("band")
    return distances_long, td, contrast
