#!/usr/bin/env python3
"""audit_92 — generate an R=1000 matched-strength surrogate ensemble for beta.

Purpose: remove the R=200 permutation floor (min p = 1/201 ≈ 0.005) so the
BH-FDR-corrected q for the beta OFC localization (audit_83) is no longer pinned
at the resolution limit. A fresh, independent ensemble of 1000 strength-
preserving surrogates per (patient, phase) is generated and cached at the
canonical surrogate path with R=1000 in the filename, so it does NOT clobber the
locked R=200 ensemble. audit_83 then reads it via ``--R 1000``.

The phase FCs are loaded IDENTICALLY to audit_83.obs_trace (the canonical
observed): full phases via ``load_fc_matrix``; rest_pre halves via the audit_63
half-FC ``.npy`` cache — so the R=1000 surrogate is the strength-matched null of
exactly the graphs the observed trace is built from.

Cost: ~1000 × 4 phases × 10 patients surrogates; ~110 MB/file, ~4.3 GB total.
Read-only on FCs; writes only surrogate eig caches. Idempotent (cache-skips).
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import load_or_compute_surrogate_eigs
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = setup_script_env()
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND = "beta"
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
R, SWAP, SEED = 1000, 20, 20260511


def load_W(pat, phase):
    if phase in ("rest_pre_A", "rest_pre_B"):
        return np.load(HALVES_FC_CACHE / pat
                       / f"{BAND}_{phase}_imcoh_abs.npy").astype(np.float64)
    return np.asarray(load_fc_matrix(pat, phase, BAND, "imcoh_abs"),
                      dtype=np.float64)


def main():
    rng = np.random.default_rng(SEED)
    for pat in COHORT:
        for ph in PHASES:
            try:
                W = load_W(pat, ph)
            except FileNotFoundError:
                print(f"  [skip] {pat} {ph}: FC missing")
                continue
            ev, _ = load_or_compute_surrogate_eigs(
                pat, BAND, ph, W, R, SWAP, SEED, rng,
                fc_method="imcoh_abs", verbose=True)
            ok = int(np.isfinite(ev).all(axis=1).sum())
            print(f"  {pat} {ph}: N={W.shape[0]}  ok={ok}/{R}")
    print("done -> R=1000 beta surrogate ensemble cached")


if __name__ == "__main__":
    main()
