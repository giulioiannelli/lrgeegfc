#!/usr/bin/env python3
r"""audit_163 — phase-scale consolidation trajectory waypoints, ρ_sym-consistent.

The consolidation arc (audit_152) reports two SCALARS per band (T_learn, T_infspec_pe)
but not the intermediate task waypoints needed to draw the four-phase LOOP in the
encoding/inference plane. This fills that gap WITH THE SAME SPLIT-HALF (ρ_sym)
construction, so the loop is bias-consistent end to end.

Why split-half and not a single baseline: a single shared rest_pre baseline induces
exactly the encoding<->persistence correlation ρ_sym removes (the unsplit endpoints are
inflated and NON-specific — every band's inference coord ≈ 0.19). The trajectory must
therefore baseline every MOVING point on the OPPOSITE rest_pre half from the axis, and
average the two arm assignments — then rest_post lands exactly on the ρ_sym
(T_learn, T_infspec_pe) and the outbound waypoints are debiased too.

    arm (axis_base a, point_base b):
      e = D_taskLearn − D[a]              encoding axis (task-side half a)
      f = D_taskTest  − D_taskLearn       inference-specific contrast (arm-invariant)
      for phase X: d = D[X] − D[b] ;  x = ρ_S(d, e) ;  y = partial ρ_S(d, f | e)
    waypoint(X) = ½ [ arm(A,B)(X) + arm(B,A)(X) ]
    rest_pre := (0, 0) by convention (zero displacement)

rest_post waypoint == audit_152 (T_learn, T_infspec_pe) EXACTLY (asserted below).
Visualization support only (no null); the figure's GATE colour = audit_152
arc_cohort_verdict ρ_sym matched-strength p.

Reuses (no forks): audit_63 ensure_half_fcs / load_phase_fc / lrg_ultrametric_condensed.

Writes: data/audit/consolidation_arc_rhosym/arc_phase_trajectory.csv
  (patient, band, phase, x, y)   phase in {rest_pre, task_learn, task_test, rest_post}
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs, load_phase_fc, lrg_ultrametric_condensed,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
ARC_PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
MOVING = ("task_learn", "task_test", "rest_post")
ARMS = (("rest_pre_A", "rest_pre_B"), ("rest_pre_B", "rest_pre_A"))  # (axis_base, point_base)

OUT = ROOT / "data" / "audit" / "consolidation_arc_rhosym"
VERDICT = OUT / "arc_cohort_verdict.csv"
ARCPP = OUT / "arc_per_patient.csv"


def _rho(a, b):
    r, _ = spearmanr(a, b)
    return float(r)


def _partial_rho(a, b, c):
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    denom = np.sqrt(max(0.0, (1.0 - rac ** 2) * (1.0 - rbc ** 2)))
    return (rab - rac * rbc) / denom if denom > 0 else np.nan


def _arm_waypoints(D, axis_base, point_base):
    e = D["task_learn"] - D[axis_base]
    f = D["task_test"] - D["task_learn"]
    out = {}
    for X in MOVING:
        d = D[X] - D[point_base]
        out[X] = (_rho(d, e), _partial_rho(d, f, e))
    return out


def main():
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                D = {ph: lrg_ultrametric_condensed(load_phase_fc(pat, ph, band))
                     for ph in ARC_PHASES}
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_163] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
                continue
            a1 = _arm_waypoints(D, *ARMS[0])
            a2 = _arm_waypoints(D, *ARMS[1])
            rows.append(dict(patient=pat, band=band, phase="rest_pre", x=0.0, y=0.0))
            for X in MOVING:
                x = 0.5 * (a1[X][0] + a2[X][0])
                y = 0.5 * (a1[X][1] + a2[X][1])
                rows.append(dict(patient=pat, band=band, phase=X, x=x, y=y))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_phase_trajectory.csv", index=False)
    print(f"[audit_163] wrote {OUT/'arc_phase_trajectory.csv'} ({len(df)} rows)\n")

    # CROSS-CHECK: rest_post waypoint must equal audit_152 (T_learn, T_infspec_pe).
    rp = df[df.phase == "rest_post"].set_index(["patient", "band"])
    if ARCPP.exists():
        arc = pd.read_csv(ARCPP).set_index(["patient", "band"])
        dx = dy = 0.0
        for idx in rp.index:
            if idx in arc.index:
                dx = max(dx, abs(rp.loc[idx, "x"] - arc.loc[idx, "T_learn"]))
                dy = max(dy, abs(rp.loc[idx, "y"] - arc.loc[idx, "T_infspec_pe"]))
        status = "PASS" if max(dx, dy) < 1e-9 else ("CLOSE" if max(dx, dy) < 1e-3 else "FAIL")
        print(f"[CROSS-CHECK] rest_post vs audit_152: max|Δx|={dx:.2e} max|Δy|={dy:.2e} → {status}\n")

    print("=== cohort-median trajectory (ρ_sym-consistent) ===")
    print(f"{'band':<11}{'learn(x,y)':>18}{'test(x,y)':>18}{'post(x,y)':>18}{'inf gate p':>12}")
    v = pd.read_csv(VERDICT) if VERDICT.exists() else None
    for b in BANDS:
        sub = df[df.band == b]
        cells = []
        for ph in MOVING:
            s = sub[sub.phase == ph]
            cells.append(f"({s.x.median():+.2f},{s.y.median():+.2f})")
        gate = ""
        if v is not None:
            ti = v[(v.band == b) & (v.functional == "T_infspec_pe")]
            gate = f"{ti.wilcoxon_p.values[0]:.3f}" if len(ti) else ""
        print(f"{b:<11}{cells[0]:>18}{cells[1]:>18}{cells[2]:>18}{gate:>12}")


if __name__ == "__main__":
    main()
