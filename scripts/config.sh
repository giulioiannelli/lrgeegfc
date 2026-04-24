#!/bin/bash
# =============================================================================
# Pipeline Configuration
# =============================================================================
# Default parameters for testing and batch processing

# Test mode defaults (single case for quick validation)
TEST_PATIENT="Pat_03"
TEST_PHASE="rest_pre"
TEST_BAND="beta"

# Batch mode defaults.
# BATCH_PATIENTS is intentionally a small smoke-test subset; the full
# canonical roster (n=10 as of 2026-04-22) is auto-detected by
# list_patients() in config/const.py and used by the CLI when
# --patients is omitted.
BATCH_PATIENTS=(Pat_02 Pat_03)
BATCH_PHASES=(rest_pre task_learn task_test rest_post)
BATCH_BANDS=(delta theta alpha beta low_gamma high_gamma)

# FC method choices
FC_METHODS=(corr msc)

# Visualization types
PLOT_TYPES=(summary full)

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Export for use in other scripts
export TEST_PATIENT TEST_PHASE TEST_BAND
export BATCH_PATIENTS BATCH_PHASES BATCH_BANDS
export FC_METHODS PLOT_TYPES
export RED GREEN YELLOW BLUE CYAN NC
