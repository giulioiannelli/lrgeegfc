#!/usr/bin/env python3
"""Data-based check: the tensor re-run IS the E1 pipeline, not a new one.

``lane_e_verify_rankcorr.py`` proves the algebra on random graphs. This proves
the *pipeline*: it derives the four functionals from the stored rank-correlation
tensor of ``rankcorr_grid`` and compares them, cell by cell, with the numbers
``lane_e_01_encinf_knob_grid.py`` wrote independently -- observed values, the
whole matched-strength surrogate ensemble (the seeds are shared, so realization
``r`` must match realization ``r``), and the role swap.

If these agree, the re-run introduced no change other than what is stored, and
every result derived from the tensor inherits E1's provenance. Expected
agreement is ~1e-7: the tensor is stored as float32.
"""
from __future__ import annotations

import argparse

import numpy as np

from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_from_corr_stack as funcs_of,
)
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = ROOT / "data" / "paper_final" / "lane_e_encinf"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-cells", type=int, default=0, help="0 = all")
    ap.add_argument("--tol", type=float, default=1e-5)
    a = ap.parse_args()

    cells = sorted((BASE / "rankcorr_grid" / "cells").glob("*.npz"))
    if a.max_cells:
        cells = cells[:a.max_cells]
    worst = dict(obs=0.0, surr=0.0, swap=0.0)
    n = 0
    for f in cells:
        old = BASE / "grid" / "cells" / f.name
        if not old.exists():
            continue
        zn, zo = np.load(f), np.load(old)
        got = dict(obs=funcs_of(zn["Robs"].astype(float)),
                   surr=funcs_of(zn["Rsurr"].astype(float)),
                   swap=funcs_of(zn["Robs"].astype(float), swap=True))
        for i, k in enumerate(list(zo["funcs"])):
            for key in worst:
                worst[key] = max(worst[key],
                                 float(np.nanmax(np.abs(got[key][k] - zo[key][..., i]))))
        n += 1
    print(f"cells compared: {n} of {len(cells)} tensor cells")
    for k, v in worst.items():
        print(f"  max |tensor {k:4s} - E1 {k:4s}| = {v:.3e}")
    bad = [k for k, v in worst.items() if v > a.tol]
    print("\nVERDICT:", "PASS" if (n and not bad) else f"FAIL {bad or 'no overlap'}")
    raise SystemExit(0 if (n and not bad) else 1)


if __name__ == "__main__":
    main()
