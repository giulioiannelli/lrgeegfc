#!/usr/bin/env python3
"""Audit 103 — four-phase cophenetic consolidation arc (observed pass).

Extends the audit_83 split-baseline full-graph trace from a binary
(does rest_post echo task_test?) to the full arc
``rest_pre → task_learn → task_test → rest_post``, and decomposes the persistent
offline reorganization into an ENCODING component and an INFERENCE-SPECIFIC
component. Scope + 5-point preamble:
``.agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md``.

Pure reuse of verified primitives (audit_63 loaders + LRG ultrametric). NO new
connectivity estimator, NO behavioral data. This is the OBSERVED pass: it
computes the arc functionals, AUDITS phase-length asymmetry, and CROSS-CHECKS
``T_test`` against the cached audit_83 full-graph ``obs_stat`` (anti-hallucination
gate). Matched-strength surrogate null + figures are the next stage, run only
after the observed pattern is read with the user (no pre-registered gate).

Notation (see scope doc)
------------------------
D_x  = LRG cophenetic condensed vector for phase x (τ=1/λ_max, UPGMA).
e = D_TL − D_preA   (encoding)         g = D_TT − D_preA   (inference-online)
f = D_TT − D_TL     (inference-specific) p = D_RP − D_preB  (persistent)
T_test    = ρ(g, p)   ← reproduces audit_83 full-graph obs_stat (CROSS-CHECK)
T_learn   = ρ(e, p)   T_infspec = ρ(f, p)   C_LT = ρ(e, g)
partials  : ρ(f,p·e), ρ(e,p·f)   (first-order partial Spearman, rank world)

Positive = TRACE (locked sign convention).

Output
------
data/audit/consolidation_arc/arc_per_patient.csv   (tidy: one row per pat×band)
+ printed cohort summary table.
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS

# Reuse the EXACT proven audit_63 loaders + LRG ultrametric (import, no fork).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)

OUT = ROOT / "data" / "audit" / "consolidation_arc"
OUT.mkdir(parents=True, exist_ok=True)

# audit_83 full-graph cophenetic obs, for the cross-check gate.
AUDIT83_CSV = ROOT / "data" / "audit" / "wm_stratified" / "cophenetic_raw_per_patient.csv"


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman, matching audit_83._rho exactly."""
    r, _ = spearmanr(a, b)
    return float(r)


def _partial_rho(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """First-order partial Spearman of (a, b) controlling c (rank world)."""
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    denom = np.sqrt(max(0.0, (1.0 - rac**2) * (1.0 - rbc**2)))
    return (rab - rac * rbc) / denom if denom > 0 else np.nan


def _arc_row(pat: str, band: str) -> dict:
    # Cophenetic condensed vectors for the five phase views.
    D = {x: lrg_ultrametric_condensed(load_phase_fc(pat, x, band))
         for x in ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")}
    e = D["task_learn"] - D["rest_pre_A"]   # encoding
    g = D["task_test"] - D["rest_pre_A"]    # inference-online
    f = D["task_test"] - D["task_learn"]    # inference-specific
    p = D["rest_post"] - D["rest_pre_B"]    # persistent
    return {
        "patient": pat, "band": band, "n_pairs": int(p.size),
        "T_test":    _rho(g, p),   # ← cross-check vs audit_83 full-graph obs_stat
        "T_learn":   _rho(e, p),
        "T_infspec": _rho(f, p),
        "C_LT":      _rho(e, g),
        "T_infspec_pe": _partial_rho(f, p, e),  # inference-specific | encoding
        "T_learn_pf":   _partial_rho(e, p, f),  # encoding | inference-specific
    }


def main() -> None:
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)  # audit_63 half-FC cache (no-op if present)
        for band in BANDS:
            try:
                rows.append(_arc_row(pat, band))
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_103] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
                rows.append({"patient": pat, "band": band, "error": str(exc)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_per_patient.csv", index=False)
    print(f"\n[audit_103] wrote {OUT/'arc_per_patient.csv'}  ({len(df)} rows)\n")

    # --- CROSS-CHECK GATE: T_test must reproduce audit_83 full-graph cophenetic ---
    if AUDIT83_CSV.exists():
        ref = pd.read_csv(AUDIT83_CSV)
        ref = ref[(ref["substrate"] == "cophenetic") & (ref["config"] == "full")]
        ref = ref.set_index(["patient", "band"])["obs_stat"]
        merged = df.dropna(subset=["T_test"]).set_index(["patient", "band"])
        diffs = []
        for idx, r in merged.iterrows():
            if idx in ref.index:
                diffs.append(abs(r["T_test"] - float(ref.loc[idx])))
        if diffs:
            md = max(diffs)
            status = "PASS" if md < 1e-6 else ("CLOSE" if md < 0.02 else "FAIL")
            print(f"[CROSS-CHECK] T_test vs audit_83 full-graph: n={len(diffs)} "
                  f"cells, max|Δ|={md:.2e}  → {status}")
            if status == "FAIL":
                print("  !! T_test does NOT reproduce audit_83 — new numbers NOT trusted.")
        else:
            print("[CROSS-CHECK] no overlapping (patient,band) cells found.")
    else:
        print(f"[CROSS-CHECK] audit_83 csv not found at {AUDIT83_CSV} — skipped.")

    # --- Cohort summary per band: median, k/10 positive, IQR ---
    print("\n=== COHORT ARC SUMMARY (median [IQR]; k/10 positive) ===")
    cols = ["T_test", "T_learn", "T_infspec", "T_infspec_pe", "T_learn_pf"]
    hdr = f"{'band':<11}" + "".join(f"{c:>22}" for c in cols)
    print(hdr)
    for band in BANDS:
        sub = df[df["band"] == band]
        cells = []
        for c in cols:
            v = sub[c].dropna().values if c in sub else np.array([])
            if v.size:
                med = np.median(v)
                q1, q3 = np.percentile(v, [25, 75])
                kpos = int((v > 0).sum())
                cells.append(f"{med:+.2f}[{q1:+.2f},{q3:+.2f}] {kpos}/{v.size}")
            else:
                cells.append("--")
        print(f"{band:<11}" + "".join(f"{c:>22}" for c in cells))

    print("\n[read] T_test = established trace (audit_83). "
          "T_infspec vs T_learn = does rest_post consolidate INFERENCE (f) or "
          "ENCODING (e)? Partials T_infspec_pe / T_learn_pf = unique contribution.")


if __name__ == "__main__":
    main()
