"""Refresh the halves LRG cache so every NPZ carries `eigenvalues` and
`eigenvectors` (added 2026-04-29 to support the eigenvector-direct
pivot E1/E2/E3 rungs).

Walks ``data/cache/imcoh_lrg_halves/<pat>/<band>_<phase>_<A|B>_lrg_imcoh-abs.npz``
for the n=10 cohort and recomputes each cell with ``overwrite_cache=True``.
Mirrors the FC step from ``h2e_split_half.ensure_halves_lrg`` but skips the
downstream stats so this is a pure cache-refresh utility.

Usage:
    python scripts/01_compute/audit/refresh_lrg_halves_eigvec.py [-v]
"""
from __future__ import annotations

import argparse
import gc
import sys

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    FS_OVERRIDES,
    PATIENTS_LIST,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.lrg import compute_lrg_analysis

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # noqa: E402

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
REST_PHASES = ("rest_pre", "rest_post")


def _load_X(pat: str, phase: str):
    try:
        X = load_timeseries(pat, phase, SEEG_DATAPATH)
    except (FileNotFoundError, OSError):
        return None
    if X is None:
        return None
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, 2048.0)
    return X, fs


def refresh_patient(pat: str, verbose: bool = False) -> int:
    n = 0
    for phase in REST_PHASES:
        loaded = _load_X(pat, phase)
        if loaded is None:
            print(f"  {pat}/{phase}: timeseries missing, skip")
            continue
        X, fs = loaded
        nperseg_half = max(256, nperseg_for_fs(fs) // 2)
        if verbose:
            print(f"  {pat}/{phase}: X.shape={X.shape}, nperseg_half={nperseg_half}")
        halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
        del X
        gc.collect()
        for (band, tag), A in halves_fc.items():
            synth_phase = f"{phase}_{tag}"
            compute_lrg_analysis(
                A, pat, synth_phase, band, fc_method="imcoh_abs",
                cache_root=HALVES_CACHE,
                use_cache=True, overwrite_cache=True,
                verbose=False,
            )
            n += 1
        halves_fc.clear()
        del halves_fc
        gc.collect()
    return n


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    total = 0
    for pat in PATIENTS_LIST:
        print(f"=== {pat} ===")
        total += refresh_patient(pat, verbose=args.verbose)
    print(f"\nDone: refreshed {total} halves NPZs.")


if __name__ == "__main__":
    main()
