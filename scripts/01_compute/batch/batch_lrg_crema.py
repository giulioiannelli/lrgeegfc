#!/usr/bin/env python
"""Batch LRG analysis on CReMa-validated MSC matrices.

Two phases:
  1. Compute: run LRG on CReMa MSC for all patients (subprocess per patient)
  2. Plot:   generate all LRG figures (full panel, entropy, dendrogram, ultrametric)

Results go to:
  - Cache:   data/lrg_cache_crema/
  - Figures: data/figures/lrg_crema/
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
import time

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, DEFAULT_NPERSEG
from lrg_eegfc.config.paths import FIGURES_ROOT, LRG_CREMA_CACHE

PATIENTS_PHASES = {
    "Pat_02": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_03": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_05": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "rest_post", "task_learn"],
    "Pat_08": ["rest_pre", "rest_post", "task_learn", "task_test"],
}

BANDS = BRAIN_BANDS_NAMES

# CReMa loading defaults (must match batch_ecm_adaptive.py)
NPERSEG = DEFAULT_NPERSEG
ECM_ALPHA_MIN = 0.01
ECM_ALPHA_MAX = 0.50
ECM_N_ENSEMBLE = 100
ECM_WEIGHT_SCALE = 1000

LRG_CACHE = str(LRG_CREMA_CACHE)
FIG_DIR = str(FIGURES_ROOT / "lrg_crema")

# -----------------------------------------------------------------------
# Worker: compute LRG for one patient (all bands × phases)
# -----------------------------------------------------------------------
COMPUTE_WORKER = textwrap.dedent("""\
import sys, gc, time
import numpy as np
sys.path.insert(0, "src")

from pathlib import Path
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import compute_lrg_analysis, load_lrg_result

patient = "{patient}"
phases = {phases}
bands = {bands}
lrg_cache = Path("{lrg_cache}")

total = len(phases) * len(bands)
done = 0

for phase in phases:
    for band in bands:
        done += 1
        tag = f"{{patient}}/{{band}}/{{phase}}"

        # Check if already cached
        existing = load_lrg_result(patient, phase, band, "msc", lrg_cache)
        if existing is not None:
            print(f"  [{{done}}/{{total}}] {{tag}} — cached ({{existing.n_nodes}} nodes)")
            continue

        # Load CReMa-validated MSC
        W = load_msc_matrix(
            patient, phase, band,
            sparsify="ecm_adaptive", nperseg={nperseg},
            ecm_alpha_min={ecm_alpha_min}, ecm_alpha_max={ecm_alpha_max},
            ecm_n_ensemble={ecm_n_ensemble}, ecm_weight_scale={ecm_weight_scale},
        )
        if W is None:
            print(f"  [{{done}}/{{total}}] {{tag}} — SKIP (no CReMa cache)")
            continue

        t0 = time.time()
        try:
            result = compute_lrg_analysis(
                adjacency_matrix=W,
                patient=patient,
                phase=phase,
                band=band,
                fc_method="msc",
                cache_root=lrg_cache,
                use_cache=True,
                overwrite_cache=False,
                verbose=False,
            )
            elapsed = time.time() - t0
            print(f"  [{{done}}/{{total}}] {{tag}} — {{result.n_nodes}} nodes, {{elapsed:.1f}}s")
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  [{{done}}/{{total}}] {{tag}} — FAILED: {{e}} ({{elapsed:.0f}}s)")

        del W
        gc.collect()

print(f"DONE {{patient}}")
""")

# -----------------------------------------------------------------------
# Worker: generate plots for one patient (all bands × phases)
# -----------------------------------------------------------------------
PLOT_WORKER = textwrap.dedent("""\
import sys, gc, warnings
import numpy as np
sys.path.insert(0, "src")
warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib
matplotlib.use("Agg")

from pathlib import Path
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.lrg import (
    plot_lrg_entropy_curves,
    plot_lrg_dendrogram,
    plot_ultrametric_heatmap,
)
from lrg_eegfc.visuals.lrg_panels import plot_lrg_full_panel

patient = "{patient}"
phases = {phases}
bands = {bands}
lrg_cache = Path("{lrg_cache}")
fig_dir = Path("{fig_dir}")

total = len(phases) * len(bands)
done = 0

for phase in phases:
    for band in bands:
        done += 1
        tag = f"{{patient}}/{{band}}/{{phase}}"

        # Check LRG result exists
        result = load_lrg_result(patient, phase, band, "msc", lrg_cache)
        if result is None:
            print(f"  [{{done}}/{{total}}] {{tag}} — SKIP (no LRG cache)")
            continue

        # Load CReMa MSC for full panel heatmap
        W = load_msc_matrix(
            patient, phase, band,
            sparsify="ecm_adaptive", nperseg={nperseg},
            ecm_alpha_min={ecm_alpha_min}, ecm_alpha_max={ecm_alpha_max},
            ecm_n_ensemble={ecm_n_ensemble}, ecm_weight_scale={ecm_weight_scale},
        )

        out_dir = fig_dir / patient
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"{{band}}_{{phase}}_lrg_msc"

        try:
            # 1. Full 5-panel
            plot_lrg_full_panel(
                patient, phase, band, "msc",
                cache_root=lrg_cache,
                output_path=out_dir / f"{{prefix}}_full.png",
                fc_matrix=W,
                verbose=False,
            )

            # 2. Entropy curves
            plot_lrg_entropy_curves(
                patient, phase, band, "msc",
                cache_root=lrg_cache,
                output_path=out_dir / f"{{prefix}}_entropy.png",
            )

            # 3. Dendrogram
            plot_lrg_dendrogram(
                patient, phase, band, "msc",
                cache_root=lrg_cache,
                output_path=out_dir / f"{{prefix}}_dendrogram.png",
            )

            # 4. Ultrametric heatmap
            plot_ultrametric_heatmap(
                patient, phase, band, "msc",
                cache_root=lrg_cache,
                output_path=out_dir / f"{{prefix}}_ultrametric.png",
            )

            print(f"  [{{done}}/{{total}}] {{tag}} — 4 figures OK")

        except Exception as e:
            print(f"  [{{done}}/{{total}}] {{tag}} — FAILED: {{e}}")

        del W
        gc.collect()

print(f"DONE {{patient}}")
""")


def run_phase(phase_name, worker_template):
    """Run a worker template for all patients."""
    total_combos = sum(len(ph) * len(BANDS) for ph in PATIENTS_PHASES.values())
    print(f"\n{'='*60}")
    print(f"  {phase_name}: {total_combos} combos across {len(PATIENTS_PHASES)} patients")
    print(f"{'='*60}\n")

    for patient, phases in PATIENTS_PHASES.items():
        n = len(phases) * len(BANDS)
        print(f"--- {patient}: {n} combos ({len(phases)} phases x {len(BANDS)} bands) ---")

        code = worker_template.format(
            patient=patient,
            phases=repr(phases),
            bands=repr(BANDS),
            lrg_cache=LRG_CACHE,
            fig_dir=FIG_DIR,
            nperseg=NPERSEG,
            ecm_alpha_min=ECM_ALPHA_MIN,
            ecm_alpha_max=ECM_ALPHA_MAX,
            ecm_n_ensemble=ECM_N_ENSEMBLE,
            ecm_weight_scale=ECM_WEIGHT_SCALE,
        )

        t0 = time.time()
        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, timeout=3600,
        )
        elapsed = time.time() - t0

        for line in proc.stdout.strip().splitlines():
            print(line)

        if proc.returncode != 0:
            print(f"  FAILED (rc={proc.returncode}) in {elapsed:.0f}s")
            for line in proc.stderr.strip().splitlines()[-10:]:
                print(f"  ERR: {line}")
        else:
            print(f"  Completed in {elapsed:.0f}s\n")


def main():
    total_start = time.time()

    # Phase 1: Compute LRG
    run_phase("PHASE 1: COMPUTE LRG on CReMa MSC", COMPUTE_WORKER)

    # Phase 2: Generate plots
    run_phase("PHASE 2: GENERATE LRG PLOTS", PLOT_WORKER)

    total_elapsed = time.time() - total_start
    print(f"\nTotal time: {total_elapsed/60:.1f} minutes")
    print(f"LRG cache: {LRG_CACHE}")
    print(f"Figures:   {FIG_DIR}")


if __name__ == "__main__":
    main()
