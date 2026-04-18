#!/usr/bin/env python3
"""Regenerate imcoh_abs LRG cache from FC matrices (one patient per subprocess)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import compute_lrg_analysis

BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patient")
    args = ap.parse_args()
    pat = args.patient

    n = 0
    for phase in PHASES:
        for band in BANDS:
            A = load_fc_matrix(pat, phase, band, "imcoh_abs")
            if A is None:
                continue
            # use_cache=True writes to cache (populate-on-miss); path is empty
            compute_lrg_analysis(A, pat, phase, band, "imcoh_abs",
                                 use_cache=True, overwrite_cache=False,
                                 verbose=False)
            n += 1
    print(f"{pat}: regenerated {n} LRG cache files", flush=True)


if __name__ == "__main__":
    main()
