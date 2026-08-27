#!/usr/bin/env bash
# Resume/complete the W0-C master grid. Safe to re-run: cells that already have
# a surrogate ensemble on disk are reused (W0C_RESUME=1 by default).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export W0C_WORKERS="${W0C_WORKERS:-4}"
exec "$HERE/run_py.sh" "$HERE/w0c_01_gate_and_tau_grid.py" "$@"
