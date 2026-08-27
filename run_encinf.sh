#!/usr/bin/env bash
# Run the W0-C five-phase encoding/inference scale grid. Resumable.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export W0C_WORKERS="${W0C_WORKERS:-10}"
exec "$HERE/run_py.sh" "$HERE/scripts/01_compute/paper_final/w0c_06_encinf_scale_grid.py" "$@"
