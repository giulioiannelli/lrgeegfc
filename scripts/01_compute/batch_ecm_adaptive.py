#!/usr/bin/env python
"""Batch-compute CReMa adaptive sparsification for all patients.

Runs each patient in a separate subprocess to free memory between patients.
Dense MSC must already be cached (nperseg=4096).
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
import time

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, DEFAULT_NPERSEG

PATIENTS_PHASES = {
    "Pat_02": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_03": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_05": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "rest_post", "task_learn"],
    "Pat_08": ["rest_pre", "rest_post", "task_learn", "task_test"],
}

BANDS = BRAIN_BANDS_NAMES
NPERSEG = DEFAULT_NPERSEG

# Worker code: process all bands/phases for ONE patient
WORKER = textwrap.dedent("""\
import sys, gc, time
import numpy as np
sys.path.insert(0, "src")

from lrg_eegfc.utils.fc.msc.sparsify import ecm_sparsify_adaptive
from lrg_eegfc.workflow.msc import load_msc_matrix, get_msc_cache_path, DEFAULT_MSC_CACHE_ROOT

patient = "{patient}"
phases = {phases}
bands = {bands}
nperseg = {nperseg}

total = len(phases) * len(bands)
done = 0

for phase in phases:
    for band in bands:
        done += 1
        tag = f"{{patient}}/{{band}}/{{phase}}"

        # Check if already cached
        cache_path = get_msc_cache_path(
            patient, phase, band, DEFAULT_MSC_CACHE_ROOT,
            sparsify="ecm_adaptive", nperseg=nperseg,
            ecm_alpha_min=0.01, ecm_alpha_max=0.50,
            ecm_n_ensemble=100, ecm_weight_scale=1000,
        )
        if cache_path.exists():
            A = np.load(cache_path)
            nz = (A > 0).sum() // 2
            print(f"  [{{done}}/{{total}}] {{tag}} — cached ({{nz}} edges)")
            continue

        # Load dense
        W = load_msc_matrix(patient, phase, band, sparsify="none", nperseg=nperseg)
        if W is None:
            print(f"  [{{done}}/{{total}}] {{tag}} — SKIP (no dense cache)")
            continue

        t0 = time.time()
        try:
            A, alpha_used = ecm_sparsify_adaptive(W, alpha_min=0.01, alpha_max=0.50,
                                                   tol=0.01, n_ensemble=100, weight_scale=1000)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(cache_path, A)
            nz = (A > 0).sum() // 2
            elapsed = time.time() - t0
            print(f"  [{{done}}/{{total}}] {{tag}} — alpha={{alpha_used:.3f}}, edges={{nz}}, {{elapsed:.0f}}s")
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  [{{done}}/{{total}}] {{tag}} — FAILED: {{e}} ({{elapsed:.0f}}s)")

        del W
        gc.collect()

print(f"DONE {{patient}}")
""")


def main():
    total_start = time.time()
    total_combos = sum(len(phases) * len(BANDS) for phases in PATIENTS_PHASES.values())
    print(f"Batch CReMa adaptive: {total_combos} combos across {len(PATIENTS_PHASES)} patients")
    print(f"Each patient runs in its own subprocess\n")

    for patient, phases in PATIENTS_PHASES.items():
        n = len(phases) * len(BANDS)
        print(f"{'='*60}")
        print(f"{patient}: {n} combos ({len(phases)} phases × {len(BANDS)} bands)")
        print(f"{'='*60}")

        code = WORKER.format(
            patient=patient,
            phases=repr(phases),
            bands=repr(BANDS),
            nperseg=NPERSEG,
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

    total_elapsed = time.time() - total_start
    print(f"\nTotal time: {total_elapsed/60:.1f} minutes")


if __name__ == "__main__":
    main()
