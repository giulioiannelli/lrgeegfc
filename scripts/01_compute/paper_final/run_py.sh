#!/usr/bin/env bash
# Worktree python runner: forces the worktree's src/ ahead of the editable
# install that points at the main checkout, and pins BLAS to 1 thread so the
# multiprocessing pools in scripts/ do not oversubscribe.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python "$@"
