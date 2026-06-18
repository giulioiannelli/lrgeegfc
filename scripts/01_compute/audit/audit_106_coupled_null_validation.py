#!/usr/bin/env python3
"""audit_106 — coupled-cross-phase null: VALIDATION GATE.

STEP 2 of the cross-phase taxonomy (scope 2026-06-12). Before the coupled null
can be used to test the *similarity* channels (anchor, reset), it must pass a
validation gate: the coupled surrogate (shared 4-cycle swap sequence across
phases) must PRESERVE the observed cross-phase cophenetic similarity (ρ^coph),
while the independent-per-phase null (audit_63 default) DEGRADES it. If the
coupled null does not bracket observed, the coupling is mis-specified and the
similarity channels fall back to the geometry baseline only.

We measure four cross-phase quantities on the LRG cophenetic vectors:
  AB        = Spearman(D_preA, D_preB)   within-baseline ceiling
  pre_post  = Spearman(D_preA, D_post)   shared-backbone similarity  <- GATE
  task_post = Spearman(D_task, D_post)
  rho_split = Spearman(D_task-D_preA, D_post-D_preB)   the trace (difference)

Expectation: coupled preserves AB / pre_post / task_post near observed; independent
drops them (esp. rho_split -> ~0, reproducing the audit_63 independent-null result).

Reuses lrg_eegfc.utils.surrogate.coupled_surrogate_cophenet + the audit_63
FC->cophenet pipeline. NO tree-cutting.

Usage:
    python audit_106_coupled_null_validation.py --band beta --patients Pat_02 Pat_08 Pat_15 --R 20
Outputs data/audit/cross_phase_taxonomy/coupled_null_validation.csv
"""
from __future__ import annotations

import argparse
import importlib.util

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.utils.surrogate import coupled_surrogate_cophenet

ROOT = setup_script_env()
_spec = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
SEED = 20260612


def xphase_metrics(DA, DB, DT, DP) -> dict:
    return dict(
        AB=float(spearmanr(DA, DB).statistic),
        pre_post=float(spearmanr(DA, DP).statistic),
        task_post=float(spearmanr(DT, DP).statistic),
        rho_split=float(spearmanr(DT - DA, DP - DB).statistic),
    )


def surrogate_metric_dist(surco: dict, n_surr: int) -> dict:
    """Per-metric array over realizations from a coupled/independent ensemble."""
    acc = {k: [] for k in ("AB", "pre_post", "task_post", "rho_split")}
    for r in range(n_surr):
        vecs = [surco[ph][r] for ph in PHASES]
        if any(v is None for v in vecs):
            continue
        m = xphase_metrics(*vecs)
        for k in acc:
            acc[k].append(m[k])
    return {k: np.array(v) for k, v in acc.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="beta")
    ap.add_argument("--patients", nargs="+",
                    default=["Pat_02", "Pat_08", "Pat_15"])
    ap.add_argument("--R", type=int, default=20)
    ap.add_argument("--swap-factor", type=int, default=20)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    for pat in args.patients:
        a63.ensure_half_fcs(pat, [args.band])

    rows = []
    print(f"\n==== coupled-null validation gate : {args.band}  R={args.R} ====")
    print(f"{'patient':<8}{'metric':<10}{'observed':>10}{'coupled med[IQR]':>22}"
          f"{'independent med[IQR]':>24}")
    for pat in args.patients:
        Ws = {ph: a63.load_phase_fc(pat, ph, args.band) for ph in PHASES}
        Dobs = {ph: a63.lrg_ultrametric_condensed(Ws[ph]) for ph in PHASES}
        m_obs = xphase_metrics(Dobs["rest_pre_A"], Dobs["rest_pre_B"],
                               Dobs["task_test"], Dobs["rest_post"])

        cpl = coupled_surrogate_cophenet(Ws, args.R, args.swap_factor, SEED,
                                         coupled=True, verbose=True)
        ind = coupled_surrogate_cophenet(Ws, args.R, args.swap_factor, SEED,
                                         coupled=False, verbose=True)
        d_cpl = surrogate_metric_dist(cpl, args.R)
        d_ind = surrogate_metric_dist(ind, args.R)

        for k in ("AB", "pre_post", "task_post", "rho_split"):
            cm, ci = np.median(d_cpl[k]), np.subtract(*np.quantile(d_cpl[k], [.75, .25]))
            im, ii = np.median(d_ind[k]), np.subtract(*np.quantile(d_ind[k], [.75, .25]))
            print(f"{pat:<8}{k:<10}{m_obs[k]:>+10.3f}"
                  f"{f'{cm:+.3f}[{ci:.3f}]':>22}{f'{im:+.3f}[{ii:.3f}]':>24}")
            rows.append(dict(patient=pat, band=args.band, metric=k,
                             observed=m_obs[k], coupled_med=float(cm),
                             coupled_iqr=float(ci), independent_med=float(im),
                             independent_iqr=float(ii)))
        print()

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "coupled_null_validation.csv", index=False)

    # gate verdict on the shared-backbone similarity (pre_post)
    pp = df[df.metric == "pre_post"]
    ratio = (pp.coupled_med / pp.observed).median()
    margin = (pp.coupled_med - pp.independent_med).median()
    print(f"[GATE] pre_post: coupled/observed median ratio = {ratio:.2f} "
          f"(want ~1, >0.7 ok); coupled-independent margin = {margin:+.3f} "
          f"(want >0)")
    rs = df[df.metric == "rho_split"]
    print(f"[sanity] rho_split independent median = "
          f"{rs.independent_med.median():+.3f} (want ~0, matches audit_63); "
          f"coupled median = {rs.coupled_med.median():+.3f}")
    gate_ok = ratio > 0.7 and margin > 0
    print(f"[GATE] {'PASS' if gate_ok else 'FAIL'} — coupled null "
          f"{'preserves' if gate_ok else 'does NOT preserve'} the cross-phase "
          f"backbone above the independent null.")
    print(f"\nOutputs -> {OUT}/coupled_null_validation.csv")


if __name__ == "__main__":
    main()
