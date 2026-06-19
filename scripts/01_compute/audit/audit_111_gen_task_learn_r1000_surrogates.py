#!/usr/bin/env python3
"""audit_111 — generate the R=1000 matched-strength surrogate ensemble for the
β ``task_learn`` phase, completing the R=1000 cache for the inference-mark
localization (audit_110).

audit_92 generated R=1000 for β phases {rest_pre_A, rest_pre_B, task_test,
rest_post} — but NOT ``task_learn`` (it pre-dates the consolidation arc). The
inference-mark localization needs ``task_learn`` (encoding / inference vectors),
so this adds the one missing β phase at R=1000, cached at the canonical path so
audit_110 reads it via ``--R 1000``. R=1000 is β-only (audit_92 only did β), so
the inference-mark R=1000 verdict is β-only by construction — which matches the
arc's β-specific result and the OFC localization's β-only multiplicity family.

Phase FC loaded IDENTICALLY to audit_92.load_W (``load_fc_matrix`` for the full
``task_learn`` phase) so the R=1000 task_learn surrogate is the strength-matched
null of exactly the graph the observed inference vectors are built from. The
R=200 task_learn ensemble (audit_103 --null) is untouched (different filename).

Cost: ~1000 × 1 phase × 10 patients; idempotent (cache-skips existing files).
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import load_or_compute_surrogate_eigs
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND = "beta"
PHASE = "task_learn"
R, SWAP, SEED = 1000, 20, 20260511


def main():
    rng = np.random.default_rng(SEED)
    for pat in COHORT:
        try:
            W = np.asarray(load_fc_matrix(pat, PHASE, BAND, "imcoh_abs"),
                           dtype=np.float64)
        except FileNotFoundError:
            print(f"  [skip] {pat} {PHASE}: FC missing")
            continue
        ev, _ = load_or_compute_surrogate_eigs(
            pat, BAND, PHASE, W, R, SWAP, SEED, rng,
            fc_method="imcoh_abs", verbose=True)
        ok = int(np.isfinite(ev).all(axis=1).sum())
        print(f"  {pat} {PHASE}: N={W.shape[0]}  ok={ok}/{R}")
    print("done -> R=1000 beta task_learn surrogate ensemble cached")


if __name__ == "__main__":
    main()
