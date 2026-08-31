#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/src:$ROOT/scripts/01_compute/paper_final${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export W0S_OUT_ANALYSIS="$ROOT/scratch/fake_out"
exec /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python - "$@" <<'PY'
import sys, os
sys.argv = ["fake"]
os.environ.setdefault("W0S_FAKE", "1")
import numpy as np, pathlib
import w0s_00_artifacts as A
_orig = A.load_cells
A.load_cells = lambda *a, **k: A.CellSet(
    pathlib.Path("scratch/fake_grid"),
    ["Pat_02","Pat_03","Pat_05","Pat_06","Pat_07","Pat_08","Pat_10","Pat_13","Pat_14","Pat_15"],
    ["theta","beta"])
import w0s_02_dilution_diagnostics as D
D.load_cells = A.load_cells
D.main()
import w0s_03_scale_local_verdict as V
V.load_cells = A.load_cells
V.N_PERM = 400
V.N_CAL_DRAWS = 15
V.N_NEFF_DRAWS = 15
V.main()
PY
