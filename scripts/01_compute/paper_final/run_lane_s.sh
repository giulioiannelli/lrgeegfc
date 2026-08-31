#!/usr/bin/env bash
# Lane W0-S end to end. Run from the repository root (the data paths are
# relative to it). Safe to re-run: the master grid reuses any per-cell surrogate
# ensemble already on disk (W0S_RESUME=1 by default), and the draws are seeded
# per cell so a resumed run is identical to an unbroken one.
#
#   ./scripts/01_compute/paper_final/run_lane_s.sh              # everything
#   ./scripts/01_compute/paper_final/run_lane_s.sh w0s_03_scale_local_verdict.py
#
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
export W0S_WORKERS="${W0S_WORKERS:-8}"
export W0S_R="${W0S_R:-100}"

if [ $# -gt 0 ]; then
    exec "$HERE/run_py.sh" "$HERE/$1" "${@:2}"
fi

"$HERE/run_py.sh" "$HERE/w0s_verify_readouts.py"          # equivalence, ~10 s
"$HERE/run_py.sh" "$HERE/w0s_01_scale_locality_grid.py"   # master grid
"$HERE/run_py.sh" "$HERE/w0s_02_dilution_diagnostics.py"  # S1
"$HERE/run_py.sh" "$HERE/w0s_03_scale_local_verdict.py"   # S0 criteria, S3, S4
"$HERE/run_py.sh" "$HERE/w0s_04_carrier_overlap.py"       # carriers + yardstick
"$HERE/run_py.sh" "$HERE/w0s_05_figures.py"               # figures
