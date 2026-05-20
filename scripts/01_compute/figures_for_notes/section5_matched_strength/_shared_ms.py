"""Shared helpers for the §5 matched-strength refinement figures.

Centralizes band order, cohort, paths, colours so the four per-probe
scripts stay compact and behave identically.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

# Walk upward to find project root (contains pyproject.toml) — mirrors
# lrg_eegfc.utils.scripting.setup_script_env without configuring matplotlib.
def _find_root() -> Path:
    here = Path(__file__).resolve().parent
    for _ in range(10):
        if (here / "pyproject.toml").exists():
            return here
        here = here.parent
    raise RuntimeError("project root with pyproject.toml not found")


ROOT = _find_root()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

# Canonical low-frequency-first band order used in §5.7.
ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

# Visual constants ----------------------------------------------------------
CLR_OBS = "#1f3d6e"        # dark blue — observed
CLR_SURR = "#6e6e6e"        # neutral grey — strength surrogate
CLR_SURR_FILL = "#cccccc"
CLR_EPIX = "#7f3b2c"        # rust brown — epi-zone-excluded
CLR_ZERO = "#888888"        # zero-line / median tick
CLR_HIGHLIGHT = "#1f7a1f"   # green: matched-strength surviving window

BAND_LABELS = {b: BRAIN_BAND_TEX_DICT[b] for b in ALL_BANDS}


# Paths --------------------------------------------------------------------
REPORT_DIR = ROOT / "data" / "reports" / "section_5_matched_strength_refinement"
FIG_DIR = REPORT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Source CSV roots (one per probe layer)
RAW_FC_DIR = ROOT / "data" / "audit" / "raw_fc_matched_strength"
LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
ALPHA_EPIX_DIR = ROOT / "data" / "audit" / "alpha_epi_exclusion"
GRASSMANN_DIR = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
GRASSMANN_EPIX_DIR = ROOT / "data" / "audit" / "grassmann_epi_exclusion"
KC_DIR = ROOT / "data" / "audit" / "kc_matched_strength_surrogate"


def cohort_ratio(obs_med: float, surr_med: float) -> float:
    """Cohort effect-size ratio |obs|/|surr| with safe handling of zero."""
    if surr_med == 0 or not np.isfinite(surr_med):
        return float("nan")
    return abs(obs_med) / abs(surr_med)


def patient_short(pat: str) -> str:
    return pat.replace("Pat_", "P")


def load_strip(path: Path, bands: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df.band.isin(bands)].copy()
    df["band"] = pd.Categorical(df.band, categories=bands, ordered=True)
    df = df.sort_values(["band", "patient"]).reset_index(drop=True)
    return df


def longest_contiguous(mask: np.ndarray) -> tuple[int, int, int]:
    """Length, start, end (inclusive) of longest True run in *mask*."""
    best = 0
    best_start = -1
    best_end = -1
    cur = 0
    cur_start = -1
    for i, m in enumerate(mask):
        if m:
            if cur == 0:
                cur_start = i
            cur += 1
            if cur > best:
                best = cur
                best_start = cur_start
                best_end = i
        else:
            cur = 0
    return best, best_start, best_end
