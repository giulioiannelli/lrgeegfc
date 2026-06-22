#!/usr/bin/env python3
"""audit_103d — matched-strength null on the arc functionals at the MESOSCALE τ.

Scope: extends `.agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-
consolidation-arc.md` (P2 "scale-of-inference") — the τ-sweep (`audit_103c`,
observed-only) surfaced a mild β inference-specific bump at τ≈2.6 and τ≈6.8,
slightly above the fine-scale τ=1. The sweep's own rule: "Asserting a τ≠1 peak as
a result requires the matched-strength surrogate at that τ — deferred." This is
that deferred follow-up: the MANDATORY referee at the peak τ.

5-POINT CRITICAL PREAMBLE
-------------------------
(1) CLAIM under test. "The β inference-specific consolidation (T_infspec·e) is
    stronger / still significant at a mesoscale diffusion time (τ≈2.6–6.8 / λ_max)
    than at the fine scale τ=1 — i.e. inference is (partly) a multi-step,
    coarse-grained phenomenon." NOT yet a result; the τ-sweep only saw an
    observed bump.
(2) NULL. The canonical matched-strength surrogate ensemble (4-cycle ±δ rewiring,
    seed 20260511, R=200), the SAME ensemble that validated the fine-scale arc
    (`audit_103 --null`) and every cohort cophenetic claim — just with its
    cophenetic reformed at the mesoscale τ from the cached eigendecomposition
    (zero regeneration). Per-patient p = mean(surr ≥ obs); cohort = one-sided
    paired Wilcoxon(obs > surr-median), exactly audit_103's logic.
(3) STRONGEST ALTERNATIVE the null must control for. That a mesoscale "bump" is a
    generic property of coarse-graining ANY strength-heterogeneous graph (the
    heat kernel homogenizes, cophenetic distances compress, concordance can drift)
    — not anything inference-specific. The matched-strength null reproduces the
    observed node strengths at EVERY τ, so if a band's bump is pure
    strength/coarse-graining geometry the surrogate bumps identically and the
    obs−surr gap stays flat. The β-specificity is the discriminator.
(4) NEGATIVE-CONTROL DISCIPLINE (the audit_103b lesson, locked). A new τ is a new
    test surface; run a KNOWN-NULL band through it. δ is dead-null for inference
    at τ=1 (p≈0.9) — if δ's T_infspec·e "passes" at the mesoscale τ, the mesoscale
    surface is manufacturing significance and the bump is retired, exactly as the
    truncation null was. All 6 bands are run; δ/θ/low-γ/high-γ are the controls.
(5) FALSIFICATION / LIMITS. PASS = β T_infspec·e clears matched-strength at τ≈2.6
    AND/OR 6.8 (Wilcoxon p<0.05, LO-P15 robust) WHILE the control bands stay null
    → a real, reportable "inference is mesoscale" claim. FAIL = β does not clear,
    or a control band also clears → the bump is observed-only noise / coarse-grain
    geometry, and the headline stays at the validated fine-scale τ=1. CANNOT do:
    prove a causal multi-step-integration mechanism (a τ-scale signature is
    correlational); rescue the τ=1 verdict if it somehow weakened (it is already
    locked by audit_103 --null).

Correctness anchor: at τ_mult=1 the per-patient T_infspec_pe obs + p reproduce
`arc_null_per_patient.csv` (audit_103 --null) — same ensemble, same τ.

Output: data/audit/consolidation_arc/arc_mesoscale_null_R200.csv
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
import _wm_stratify as ws  # type: ignore
from audit_103_cophenetic_consolidation_arc import (  # type: ignore
    ARC_PHASES, BANDS, COHORT, FUNCTIONALS, OUT, _arc_cell_rng,
    _obs_functionals, _surr_functionals, ensure_half_fcs, load_phase_fc,
)
from audit_103c_tau_sweep_arc import _cophenet_at_taumult, _eig_L  # type: ignore

# τ multiples: 1.0 = fine-scale anchor; 2.610 / 6.813 = the two β bumps, taken
# EXACTLY from the audit_103c grid (np.logspace(0,1,13)[5] and [10]).
_GRID = np.logspace(0.0, 1.0, 13)
TAU_MULTS = [1.0, float(_GRID[5]), float(_GRID[10])]   # 1.000, 2.610, 6.813
ARC_NULL_CSV = OUT / "arc_null_per_patient.csv"        # τ=1 correctness anchor


def main() -> None:
    R, sf = ws.N_SURROGATES, ws.SWAP_FACTOR
    print(f"[audit_103d] mesoscale matched-strength null  R={R} swap={sf}  "
          f"τ_mult={[round(t,3) for t in TAU_MULTS]}\n"
          f"            (cached ensemble seed 20260511 — cophenetic reformed at τ, "
          f"no regeneration)\n")
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
                eigs = {ph: load_or_compute_eigs_at_path(
                    ws.surr_eig_path("full", pat, band, ph, R, sf), Ws[ph], R, sf,
                    _arc_cell_rng(pat, band, ph), verbose=False) for ph in ARC_PHASES}
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_103d] FAIL load {pat} {band}: {type(exc).__name__}: {exc}")
                continue
            for tm in TAU_MULTS:
                D = {ph: _cophenet_at_taumult(*_eig_L(Ws[ph]), tm)[0] for ph in ARC_PHASES}
                obs = _obs_functionals(D)
                sc = {}
                for ph in ARC_PHASES:
                    evals, evecs = eigs[ph]
                    n_pairs = Ws[ph].shape[0] * (Ws[ph].shape[0] - 1) // 2
                    C = np.full((R, n_pairs), np.nan)
                    for r in range(R):
                        if np.isfinite(evals[r]).all():
                            C[r] = _cophenet_at_taumult(evals[r], evecs[r], tm)[0]
                    sc[ph] = C
                surr = _surr_functionals(sc, R)
                row = {"patient": pat, "band": band, "tau_mult": float(tm)}
                for k in FUNCTIONALS:
                    s = surr[k][np.isfinite(surr[k])]
                    row[f"{k}_obs"] = obs[k]
                    row[f"{k}_surr_p50"] = float(np.median(s)) if s.size else np.nan
                    row[f"{k}_p"] = float(np.mean(s >= obs[k])) if s.size else np.nan
                    row[f"{k}_nsurr"] = int(s.size)
                rows.append(row)
            meso = next((r for r in rows[-len(TAU_MULTS):]
                         if np.isclose(r["tau_mult"], TAU_MULTS[1])), None)
            extra = (f"  T_infspec_pe@τ{TAU_MULTS[1]:.1f}: "
                     f"obs={meso['T_infspec_pe_obs']:+.2f} p={meso['T_infspec_pe_p']:.3f}"
                     if meso else "")
            print(f"[audit_103d] {pat:<7} {band:<10} done{extra}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_mesoscale_null_R200.csv", index=False)
    print(f"\n[audit_103d] wrote {OUT/'arc_mesoscale_null_R200.csv'} ({len(df)} rows)\n")

    # --- Correctness anchor: τ_mult=1 reproduces audit_103 --null ---
    if ARC_NULL_CSV.exists():
        ref = pd.read_csv(ARC_NULL_CSV).set_index(["patient", "band"])
        cur = df[np.isclose(df["tau_mult"], 1.0)].set_index(["patient", "band"])
        diffs = []
        for idx, r in cur.iterrows():
            if idx in ref.index and "T_infspec_pe_obs" in ref.columns:
                diffs.append(abs(r["T_infspec_pe_obs"] - float(ref.loc[idx, "T_infspec_pe_obs"])))
        md = max(diffs) if diffs else np.nan
        status = "PASS" if md < 1e-9 else ("CLOSE" if md < 1e-4 else "FAIL")
        print(f"[ANCHOR] τ_mult=1 T_infspec_pe_obs vs audit_103 --null: max|Δ|={md:.2e} → {status}")

    # --- Cohort verdict per band per τ (the headline = T_infspec_pe) ---
    def _wilcox_gt(a, b):
        a = np.asarray(a, float); b = np.asarray(b, float)
        m = np.isfinite(a) & np.isfinite(b); a, b = a[m], b[m]
        if a.size < 3 or np.allclose(a, b):
            return np.nan
        try:
            return float(wilcoxon(a, b, alternative="greater").pvalue)
        except ValueError:
            return np.nan

    print("\n=== COHORT MATCHED-STRENGTH VERDICT  (positive = trace; obs vs surr-median) ===")
    print("    HEADLINE functional = T_infspec_pe (inference-specific, controlling encoding)\n")
    for k in ("T_infspec_pe", "T_learn"):
        print(f"--- {k} ---")
        print(f"  {'band':<11}{'τ/λmax':>8}{'obs_med':>9}{'surr_med':>9}"
              f"{'n_p<.05':>8}{'wilcox_p':>10}{'LOp15_p':>9}")
        for band in ("beta", "alpha", "low_gamma", "theta", "delta", "high_gamma"):
            for tm in TAU_MULTS:
                s = df[(df.band == band) & np.isclose(df.tau_mult, tm)]
                if s.empty:
                    continue
                o, sm = s[f"{k}_obs"].values, s[f"{k}_surr_p50"].values
                npos = int((s[f"{k}_p"] < 0.05).sum())
                s15 = s[s.patient != "Pat_15"]
                tag = " *neg-ctrl*" if (k == "T_infspec_pe" and band in
                                       ("delta", "theta", "low_gamma", "high_gamma")) else ""
                print(f"  {band:<11}{tm:>8.2f}{np.nanmedian(o):>+9.3f}{np.nanmedian(sm):>+9.3f}"
                      f"{npos:>5}/{len(s):<2}{_wilcox_gt(o, sm):>10.4f}"
                      f"{_wilcox_gt(s15[f'{k}_obs'], s15[f'{k}_surr_p50']):>9.4f}{tag}")
            print()
    print("[read] β T_infspec_pe clears (p<0.05, LO-P15) at τ≈2.6/6.8 WHILE every "
          "neg-ctrl band stays null → 'inference is mesoscale' is real. If a "
          "neg-ctrl band also clears at the mesoscale τ → bump is coarse-grain "
          "geometry, retire it (the audit_103b lesson). Compare τ=1 row = the "
          "locked fine-scale verdict.")


if __name__ == "__main__":
    main()
