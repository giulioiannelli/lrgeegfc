#!/usr/bin/env bash
# Canonical runner for the W0-B null ladder.
# PYTHONPATH must point at THIS worktree's src/ (the editable install points at
# the main checkout, and PYTHONPATH wins over the .pth).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export PYTHONPATH="$HERE/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python "$@"
