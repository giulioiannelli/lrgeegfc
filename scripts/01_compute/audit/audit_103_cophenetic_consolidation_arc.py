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

import argparse
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
# Matched-strength surrogate machinery — reuse audit_83's per-phase stack builder
# and the _wm_stratify cache-path / per-cell-RNG vocabulary verbatim, so the
# canonical full-graph ensembles (preA/preB/task_test/rest_post, seed 20260511)
# are CACHE HITS and only task_learn is generated fresh.
from audit_83_wm_stratified_cophenetic import _surr_stacks  # type: ignore
import _wm_stratify as ws  # type: ignore

ARC_PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
FUNCTIONALS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")

# ws.cell_rng indexes _epi_stratify._PHASE_IDX = {preA,preB,task_test,rest_post}
# which has NO task_learn → KeyError. Mirror cell_rng's SeedSequence EXACTLY
# (so the 4 canonical phases are byte-identical — and cache hits regardless),
# extending the phase index with task_learn=4 for the one fresh ensemble.
_ARC_PHASE_IDX = {"rest_pre_A": 0, "rest_pre_B": 1, "task_test": 2,
                  "rest_post": 3, "task_learn": 4}


def _arc_cell_rng(pat: str, band: str, phase: str) -> np.random.Generator:
    ss = np.random.SeedSequence([
        ws.config_seed("full"), 600 + ws.ALL_CONFIGS.index("full"),
        int(pat.split("_")[-1]), ws.ALL_BANDS.index(band), _ARC_PHASE_IDX[phase],
    ])
    return np.random.default_rng(ss)

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


# ---------------------------------------------------------------------------
# Stage 2 — matched-strength surrogate null on the new arc functionals.
# Reuses audit_83's _surr_stacks + ws cache paths/RNG verbatim: preA/preB/
# task_test/rest_post are cache hits (seed 20260511); only task_learn generates.
# ---------------------------------------------------------------------------
def _obs_functionals(D: dict) -> dict:
    e = D["task_learn"] - D["rest_pre_A"]
    g = D["task_test"] - D["rest_pre_A"]
    f = D["task_test"] - D["task_learn"]
    p = D["rest_post"] - D["rest_pre_B"]
    return {"T_test": _rho(g, p), "T_learn": _rho(e, p),
            "T_infspec": _rho(f, p), "T_infspec_pe": _partial_rho(f, p, e)}


def _surr_functionals(sc: dict, R: int) -> dict:
    """Per-surrogate null realizations of each functional, from cophenetic
    surrogate stacks ``sc[phase]`` (shape (R, n_pairs)); pairs surrogate index r
    across phases exactly as audit_83 does. NaN (strength-violating) rows skipped."""
    e_s = sc["task_learn"] - sc["rest_pre_A"]
    g_s = sc["task_test"] - sc["rest_pre_A"]
    f_s = sc["task_test"] - sc["task_learn"]
    p_s = sc["rest_post"] - sc["rest_pre_B"]
    out = {k: [] for k in FUNCTIONALS}
    for r in range(R):
        ev, gv, fv, pv = e_s[r], g_s[r], f_s[r], p_s[r]
        if not (np.isfinite(ev).all() and np.isfinite(gv).all()
                and np.isfinite(fv).all() and np.isfinite(pv).all()):
            continue
        out["T_test"].append(_rho(gv, pv))
        out["T_learn"].append(_rho(ev, pv))
        out["T_infspec"].append(_rho(fv, pv))
        out["T_infspec_pe"].append(_partial_rho(fv, pv, ev))
    return {k: np.asarray(v, float) for k, v in out.items()}


def run_null(verbose: bool = True) -> None:
    R, sf = ws.N_SURROGATES, ws.SWAP_FACTOR
    print(f"[audit_103 null] R={R} swap_factor={sf}  "
          f"(preA/preB/task_test/rest_post cached; task_learn generated)\n")
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
                D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in ARC_PHASES}
                obs = _obs_functionals(D)
                sc = {}
                for ph in ARC_PHASES:
                    coph, _raw = _surr_stacks(
                        Ws[ph], ws.surr_eig_path("full", pat, band, ph, R, sf),
                        _arc_cell_rng(pat, band, ph), R, sf, False)
                    sc[ph] = coph
                surr = _surr_functionals(sc, R)
                row = {"patient": pat, "band": band}
                for k in FUNCTIONALS:
                    s = surr[k][np.isfinite(surr[k])]
                    row[f"{k}_obs"] = obs[k]
                    row[f"{k}_surr_p50"] = float(np.median(s)) if s.size else np.nan
                    row[f"{k}_p"] = float(np.mean(s >= obs[k])) if s.size else np.nan
                    row[f"{k}_nsurr"] = int(s.size)
                rows.append(row)
                if verbose:
                    print(f"[null] {pat:<7} {band:<10} " + "  ".join(
                        f"{k}={obs[k]:+.2f}(p={row[f'{k}_p']:.3f})"
                        for k in ("T_learn", "T_infspec", "T_infspec_pe")))
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_103 null] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_null_per_patient.csv", index=False)
    print(f"\n[audit_103 null] wrote {OUT/'arc_null_per_patient.csv'} ({len(df)} rows)\n")

    # Cohort verdict per band per functional (paired Wilcoxon obs vs surr-median,
    # reusing ws.cohort_verdict — the exact audit_83 cohort logic). LO-Pat_15 too.
    print("=== COHORT MATCHED-STRENGTH VERDICT (positive = trace; obs vs surr-median) ===")
    print(f"{'band':<11}{'functional':<15}{'obs_med':>9}{'surr_med':>9}"
          f"{'n>surr':>8}{'wilcox_p':>10}{'LOp15_p':>9}")
    for band in BANDS:
        sub = df[df["band"] == band]
        for k in FUNCTIONALS:
            d = sub.dropna(subset=[f"{k}_obs", f"{k}_surr_p50"])
            if d.empty:
                continue
            n_above = int((d[f"{k}_p"] < 0.05).sum())
            v = ws.cohort_verdict(d[f"{k}_obs"].values, d[f"{k}_surr_p50"].values, n_above)
            d2 = d[d["patient"] != "Pat_15"]
            n_above2 = int((d2[f"{k}_p"] < 0.05).sum())
            v2 = ws.cohort_verdict(d2[f"{k}_obs"].values, d2[f"{k}_surr_p50"].values, n_above2)
            print(f"{band:<11}{k:<15}{v['med_obs']:>+9.3f}{v['med_surr']:>+9.3f}"
                  f"{n_above:>5}/{v['n_defined']:<2}{v['wilcoxon_p']:>10.4f}"
                  f"{v2['wilcoxon_p']:>9.4f}")
    print("\n[read] wilcox_p = one-sided paired Wilcoxon (obs > strength-matched "
          "surrogate median) across cohort. T_infspec_pe = inference-specific "
          "persistence controlling encoding — the claim under test.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", action="store_true",
                    help="run matched-strength surrogate null (slow: generates task_learn)")
    args = ap.parse_args()
    if args.null:
        run_null()
    else:
        main()
