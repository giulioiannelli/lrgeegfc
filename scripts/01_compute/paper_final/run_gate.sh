#!/usr/bin/env bash
# Run the W0-C analysis stage. W0C_BASE selects the output/base directory.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$HERE/run_py.sh" "$HERE/$1" "${@:2}"
