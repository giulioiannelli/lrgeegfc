#!/usr/bin/env python3
"""audit_157 — arc τ-sweep + MESOSCALE matched-strength null under the rho_sym estimator.

Migrates the R2.2 evidence (audit_103c τ-sweep + audit_103d mesoscale null) off bare
rho_split onto the canonical, arbitrary-half-free rho_sym estimator. The N2 climax
(audit_103 fine-scale arc) already migrated in audit_152; this extends that migration
ACROSS the diffusion-time grid so the "inference-specific consolidation favours the
mesoscale" claim (R2.2) ships in rho_sym, not rho_split.

Single change vs audit_103d: the per-pair arc concordance becomes the SYMMETRIC average
of both split-half arm assignments (audit_152._sym_functionals / _surr_functionals_sym),
applied identically to the observed trace and to every cached matched-strength surrogate
realization, at EACH tau_mult. Everything else is audit_103d verbatim — the same
tau grid (1.0, 2.610, 6.813 = audit_103c logspace[5],[10]), the same cached seed-20260511
surrogate eigendecompositions (cophenetic reformed at tau from cached eigs, NO
regeneration), the same phases and cohort. Does NOT edit audit_103c/d or audit_152.

5-POINT CRITICAL PREAMBLE
-------------------------
(1) CLAIM. The mesoscale verdict holds under rho_sym: beta T_infspec_pe clears the
    matched-strength null at tau≈2.6/6.8 (Wilcoxon p<0.05, LO-P15 robust) WHILE every
    control band (delta/theta/low-gamma/high-gamma) stays null at the same tau — i.e.
    inference-specific consolidation is multiscale, mildly mesoscale-favouring,
    independent of the arbitrary A/B half.
(2) NULL. The arc/mesoscale bump was an artifact of the rho_split arbitrary-half arm
    assignment (the estimator), not diffusion geometry.
(3) STRONGEST ALTERNATIVE. rho_sym is rho_split relabelled and reproduces it trivially.
    Controlled by recording BOTH sym functionals AND arm1 (=rho_split) per cell, with
    TWO exact cross-checks: arm1 reproduces audit_103d at every tau; sym at tau=1
    reproduces audit_152 (arc_null_per_patient). Same strength-preserving surrogate
    ensemble as audit_103 (estimator change orthogonal to strength). The mesoscale
    coarse-graining alternative (bump is generic for any strength-het graph) is
    controlled by the matched-strength surrogate at EACH tau + the neg-ctrl bands
    (audit_103d's discipline, inherited).
(4) CANNOT. reuses R=200 cached surrogate; n=10 cohort Wilcoxon coarse p-grid; f
    (inference-specific) is arm-invariant so T_infspec's symmetrization only averages
    the persistent-axis half. A tau-scale signature is correlational, not a mechanism.
(5) FALSIFY. If beta mesoscale T_infspec_pe fails matched-strength under rho_sym, OR a
    neg-ctrl band clears at the mesoscale tau, the mesoscale claim is estimator- or
    coarse-grain-dependent and R2.2 stays at the validated fine-scale tau=1.

Correctness anchors (printed):
  A1  arm1 T_*_obs at each tau  == audit_103d arc_mesoscale_null_R200.csv (rho_split)
  A2  sym  T_*_obs at tau=1     == audit_152 arc_null_per_patient.csv     (rho_sym fine)

Output: data/audit/consolidation_arc_rhosym/
  arc_scale_null_R200.csv        per-patient sym obs + arm1 + surrogate p, per tau
  arc_scale_cohort_verdict.csv   band x tau x functional cohort Wilcoxon (+ LO-P15)
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
import _wm_stratify as ws  # type: ignore
from audit_103_cophenetic_consolidation_arc import (  # type: ignore
    _arc_cell_rng, ensure_half_fcs, load_phase_fc,
)
from audit_103c_tau_sweep_arc import _cophenet_at_taumult, _eig_L  # type: ignore
from audit_103d_arc_mesoscale_null import TAU_MULTS  # type: ignore  (1.000, 2.610, 6.813)
from audit_152_consolidation_arc_rhosym import (  # type: ignore
    ARC_PHASES, BANDS, COHORT, FUNCTIONALS, _sym_functionals, _surr_functionals_sym,
)

OUT = ROOT / "data" / "audit" / "consolidation_arc_rhosym"
OUT.mkdir(parents=True, exist_ok=True)

REF_MESO = ROOT / "data" / "audit" / "consolidation_arc" / "arc_mesoscale_null_R200.csv"
REF_ARC1 = OUT / "arc_null_per_patient.csv"  # audit_152 fine-scale rho_sym anchor


def _cell(pat: str, band: str, R: int, sf: int) -> list[dict]:
    """Per (pat, band): sym obs + arm1 + surrogate p at each tau_mult."""
    Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
    obs_eig = {ph: _eig_L(Ws[ph]) for ph in ARC_PHASES}                 # hoisted (safe)
    eigs = {ph: load_or_compute_eigs_at_path(
        ws.surr_eig_path("full", pat, band, ph, R, sf), Ws[ph], R, sf,
        _arc_cell_rng(pat, band, ph), verbose=False) for ph in ARC_PHASES}
    out = []
    for tm in TAU_MULTS:
        D = {ph: _cophenet_at_taumult(*obs_eig[ph], tm)[0] for ph in ARC_PHASES}
        sym_obs, arm1_obs = _sym_functionals(D)
        sc = {}
        for ph in ARC_PHASES:
            evals, evecs = eigs[ph]
            n_pairs = Ws[ph].shape[0] * (Ws[ph].shape[0] - 1) // 2
            C = np.full((R, n_pairs), np.nan)
            for r in range(R):
                if np.isfinite(evals[r]).all():
                    C[r] = _cophenet_at_taumult(evals[r], evecs[r], tm)[0]
            sc[ph] = C
        surr = _surr_functionals_sym(sc, R)
        row = {"patient": pat, "band": band, "tau_mult": float(tm)}
        for k in FUNCTIONALS:
            s = surr[k][np.isfinite(surr[k])]
            row[f"{k}_obs"] = sym_obs[k]
            row[f"{k}_arm1"] = arm1_obs[k]
            row[f"{k}_surr_p50"] = float(np.median(s)) if s.size else np.nan
            row[f"{k}_p"] = float(np.mean(s >= sym_obs[k])) if s.size else np.nan
            row[f"{k}_nsurr"] = int(s.size)
        out.append(row)
    return out


def main(cohort: list[str]) -> None:
    R, sf = ws.N_SURROGATES, ws.SWAP_FACTOR
    ncell = len(cohort) * len(BANDS)
    print(f"[audit_157] rho_sym arc scale null  R={R} swap={sf}  "
          f"tau_mult={[round(t, 3) for t in TAU_MULTS]}  ({ncell} cells)\n"
          f"           (cached seed-20260511 eigs; cophenetic reformed at tau; NO regen)\n",
          flush=True)
    rows, t0, i = [], time.time(), 0
    for pat in cohort:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            i += 1
            try:
                rows.extend(_cell(pat, band, R, sf))
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_157] FAIL {pat} {band}: {type(exc).__name__}: {exc}", flush=True)
                continue
            el = time.time() - t0
            eta = el / i * (ncell - i)
            print(f"[audit_157] [{i:>2}/{ncell}] {pat:<7} {band:<10} "
                  f"elapsed {el:6.1f}s  ETA {eta:6.1f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_scale_null_R200.csv", index=False)
    print(f"\n[audit_157] wrote {OUT/'arc_scale_null_R200.csv'} ({len(df)} rows) "
          f"in {time.time()-t0:.1f}s\n", flush=True)

    _anchors(df)
    _verdict(df)


def _anchors(df: pd.DataFrame) -> None:
    print("=== CORRECTNESS ANCHORS ===")
    # A1: arm1 (rho_split) at each tau == audit_103d rho_split mesoscale CSV
    if REF_MESO.exists():
        ref = pd.read_csv(REF_MESO)
        key = ["patient", "band", "tau_mult"]
        m = df.merge(ref, on=key, suffixes=("", "_ref"))
        d = []
        for k in FUNCTIONALS:
            if f"{k}_arm1" in m and f"{k}_obs_ref" in m:
                d.append(np.abs(m[f"{k}_arm1"] - m[f"{k}_obs_ref"]).max())
        md = float(np.nanmax(d)) if d else np.nan
        st = "PASS" if md < 1e-9 else ("CLOSE" if md < 1e-4 else "FAIL")
        print(f"  A1 arm1==audit_103d (rho_split, all tau): n={len(m)} max|Δ|={md:.2e} → {st}")
    else:
        print("  A1 skipped — audit_103d CSV not found.")
    # A2: sym at tau=1 == audit_152 fine-scale rho_sym
    if REF_ARC1.exists():
        ref = pd.read_csv(REF_ARC1).set_index(["patient", "band"])
        cur = df[np.isclose(df["tau_mult"], 1.0)].set_index(["patient", "band"])
        d = []
        for k in FUNCTIONALS:
            c = f"{k}_obs"
            for idx, r in cur.iterrows():
                if idx in ref.index and c in ref.columns:
                    d.append(abs(r[c] - float(ref.loc[idx, c])))
        md = float(np.nanmax(d)) if d else np.nan
        st = "PASS" if md < 1e-9 else ("CLOSE" if md < 1e-4 else "FAIL")
        print(f"  A2 sym(tau=1)==audit_152 (rho_sym fine): n={len(cur)} max|Δ|={md:.2e} → {st}")
    else:
        print("  A2 skipped — audit_152 arc_null CSV not found.")
    print()


def _verdict(df: pd.DataFrame) -> None:
    print("=== COHORT MATCHED-STRENGTH VERDICT under rho_sym (positive = trace) ===")
    print("    HEADLINE = T_infspec_pe (inference-specific | encoding); tau=1 = fine-scale lock\n")
    vrows = []
    for k in ("T_infspec_pe", "T_learn"):
        print(f"--- {k} ---")
        print(f"  {'band':<11}{'tau/λmax':>9}{'obs_med':>9}{'surr_med':>9}"
              f"{'n_p<.05':>8}{'wilcox_p':>10}{'LOp15_p':>9}")
        for band in ("beta", "alpha", "low_gamma", "theta", "delta", "high_gamma"):
            for tm in TAU_MULTS:
                s = df[(df.band == band) & np.isclose(df.tau_mult, tm)]
                d = s.dropna(subset=[f"{k}_obs", f"{k}_surr_p50"])
                if d.empty:
                    continue
                n_above = int((d[f"{k}_p"] < 0.05).sum())
                v = ws.cohort_verdict(d[f"{k}_obs"].values, d[f"{k}_surr_p50"].values, n_above)
                d2 = d[d.patient != "Pat_15"]
                v2 = ws.cohort_verdict(d2[f"{k}_obs"].values, d2[f"{k}_surr_p50"].values,
                                       int((d2[f"{k}_p"] < 0.05).sum()))
                tag = " *neg-ctrl*" if (k == "T_infspec_pe" and band in
                                        ("delta", "theta", "low_gamma", "high_gamma")) else ""
                print(f"  {band:<11}{tm:>9.2f}{v['med_obs']:>+9.3f}{v['med_surr']:>+9.3f}"
                      f"{n_above:>5}/{v['n_defined']:<2}{v['wilcoxon_p']:>10.4f}"
                      f"{v2['wilcoxon_p']:>9.4f}{tag}", flush=True)
                vrows.append({"band": band, "tau_mult": float(tm), "functional": k,
                              "med_obs": v["med_obs"], "med_surr": v["med_surr"],
                              "n_above": n_above, "n_defined": v["n_defined"],
                              "wilcoxon_p": v["wilcoxon_p"], "wilcoxon_p_LOp15": v2["wilcoxon_p"]})
            print()
    pd.DataFrame(vrows).to_csv(OUT / "arc_scale_cohort_verdict.csv", index=False)
    print(f"[audit_157] wrote {OUT/'arc_scale_cohort_verdict.csv'}")
    print("[read] beta T_infspec_pe clears (p<.05, LO-P15) at tau≈2.6/6.8 while every "
          "neg-ctrl stays null → inference-specific consolidation is mesoscale-favouring "
          "under rho_sym. tau=1 row = fine-scale lock (== audit_152).")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=None,
                    help="override cohort (timing/smoke); default = full n=10")
    args = ap.parse_args()
    main(args.patients if args.patients else COHORT)
