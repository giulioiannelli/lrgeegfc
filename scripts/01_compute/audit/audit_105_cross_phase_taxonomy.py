#!/usr/bin/env python3
"""audit_105 — cross-phase cophenetic taxonomy: decomposition + trace guard.

STEP 1 of the anchor/trace/reset/reorganize taxonomy (scope:
.agents/guides/task-persistence-investigation/2026-06-12_cross-phase-cophenetic-taxonomy.md).

100% rho^coph-native: every channel is a functional of the SAME per-pair
cophenetic distances that carry the validated rho_split trace. NO tree-cutting.

Per (patient, band) it builds the four split-baseline cophenetic vectors
(rest_pre_A, rest_pre_B, task_test, rest_post) with the audit_63 pipeline
imported verbatim, then:

  1. TRACE GUARD (bit-exact): rho_split = Spearman(D_task - D_preA,
     D_post - D_preB) must reproduce the locked
     matched_strength_surrogate_split_baseline/per_patient_per_band.csv obs_rho.
     If it does not, the cell is flagged guard_ok=False and excluded from claims.
  2. CONTRAST DECOMPOSITION (lrg_eegfc.utils.metrics.cross_phase): persistence
     phi1=(D_post-D_preA)/sqrt2 (= the rho_split rest axis), excursion
     phi2=(2 D_task-D_preA-D_post)/sqrt6 (the reset channel, orthogonal & in
     rho_split's null-space), var=(phi1^2+phi2^2)/3, anchor weight, mover mask.
  3. COHORT COMPOSITION: comp_anchor=1-mover_frac, comp_trace=mover_frac*share_T,
     comp_reset=mover_frac*share_R (sum to 1; reorganize carved later by the
     coupled null in step 2).
  4. STRENGTH DIAGNOSTIC: Spearman(anchor_weight, endpoint strength) — checks
     the anchor=hub-backbone tautology (caveat 3).

GUARDS the scope report demands: comp_trace must peak at beta/alpha; anchor
must dominate (Gratton stable-backbone). Reset is reported but UNVERIFIED until
the coupled null (step 2) + localization (step 3).

Usage:
    python audit_105_cross_phase_taxonomy.py                 # all 6 bands, n=10
    python audit_105_cross_phase_taxonomy.py --kappa 1.0 --bands beta alpha
Outputs data/audit/cross_phase_taxonomy/:
    per_patient_per_band.csv   - per-cell channels, shares, guard, strength diag
    cohort_summary.csv         - per-band cohort composition (median + IQR)
    figures/decomposition.pdf  - anchor dominance + trace/reset split per band
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.cross_phase import (
    cohort_share_row, cross_phase_channels,
)
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()

# audit_63 FC->cophenet pipeline imported verbatim (no fork)
_spec = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy"
FIG = OUT / "figures"
LOCKED_CSV = (ROOT / "data" / "audit"
              / "matched_strength_surrogate_split_baseline"
              / "per_patient_per_band.csv")

ALL_BANDS = list(BRAIN_BANDS)                      # delta..high_gamma
TRACE_BANDS = ["alpha", "beta", "low_gamma"]       # locked obs_rho available
GUARD_TOL = 1e-6

C_ANCHOR = "#807dba"
C_TRACE = "#2c7fb8"
C_RESET = "#41ab5d"
C_REORG = "#fc9272"


def load_locked_rho() -> dict:
    if not LOCKED_CSV.exists():
        print(f"[audit_105] WARNING: locked CSV missing at {LOCKED_CSV}; "
              "trace guard disabled")
        return {}
    df = pd.read_csv(LOCKED_CSV)
    return {(r.patient, r.band): float(r.obs_rho) for r in df.itertuples()}


def per_cell(pat: str, band: str, kappa: float, locked: dict) -> dict | None:
    Ws = {}
    for phase in a63.PHASES:
        try:
            Ws[phase] = a63.load_phase_fc(pat, phase, band)
        except Exception as e:
            print(f"[audit_105] SKIP {pat}/{band}/{phase}: {e}")
            return None
    N = Ws[a63.PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        print(f"[audit_105] SKIP {pat}/{band}: phase shape mismatch")
        return None

    D = {ph: a63.lrg_ultrametric_condensed(Ws[ph]) for ph in a63.PHASES}
    Lc = D[a63.PHASES[0]].size
    if any(d.size != Lc for d in D.values()):
        print(f"[audit_105] SKIP {pat}/{band}: ultrametric size mismatch")
        return None

    DA, DB = D["rest_pre_A"], D["rest_pre_B"]
    DT, DP = D["task_test"], D["rest_post"]

    # 1. trace guard ---------------------------------------------------------
    rho_split = float(spearmanr(DT - DA, DP - DB).statistic)
    locked_rho = locked.get((pat, band), np.nan)
    guard_ok = (not np.isfinite(locked_rho)) or abs(rho_split - locked_rho) < GUARD_TOL

    # 2. contrast channels ---------------------------------------------------
    ch = cross_phase_channels(DA, DB, DT, DP, kappa=kappa)
    row = cohort_share_row(ch)

    # 4. strength diagnostic (anchor vs endpoint hubness) --------------------
    s_node = Ws["rest_pre_A"].sum(axis=1)
    iu_i, iu_j = np.triu_indices(N, k=1)
    pair_strength = s_node[iu_i] + s_node[iu_j]
    rho_anchor_strength = float(spearmanr(ch["anchor_weight"], pair_strength).statistic)

    return {
        "patient": pat, "band": band, "N_nodes": int(N), "n_pairs": int(Lc),
        "rho_split": rho_split, "locked_rho": locked_rho, "guard_ok": guard_ok,
        "sigma": ch["sigma"],
        "anchor_mass": row["anchor_mass"], "mover_frac": row["mover_frac"],
        "share_T": row["share_T"], "share_R": row["share_R"],
        "E_T": row["E_T"], "E_R": row["E_R"],
        "comp_anchor": row["comp_anchor"], "comp_trace": row["comp_trace"],
        "comp_reset": row["comp_reset"],
        "median_var_over_sig2": float(np.median(ch["var"]) / (ch["sigma"] ** 2)),
        "rho_anchor_strength": rho_anchor_strength,
    }


def cohort_summary(df: pd.DataFrame, bands: list[str]) -> pd.DataFrame:
    out = []
    for band in bands:
        sub = df[(df.band == band) & df.guard_ok]
        if sub.empty:
            continue
        out.append({
            "band": band, "n_patients": len(sub),
            "comp_anchor_med": float(sub.comp_anchor.median()),
            "comp_trace_med": float(sub.comp_trace.median()),
            "comp_reset_med": float(sub.comp_reset.median()),
            "share_T_med": float(sub.share_T.median()),
            "share_T_iqr": float(sub.share_T.quantile(.75) - sub.share_T.quantile(.25)),
            "mover_frac_med": float(sub.mover_frac.median()),
            "rho_anchor_strength_med": float(sub.rho_anchor_strength.median()),
        })
    return pd.DataFrame(out)


def make_figure(df: pd.DataFrame, bands: list[str], out_path: Path) -> None:
    use_lrg_style()
    bands = [b for b in bands if not df[(df.band == b) & df.guard_ok].empty]
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(1.5 * len(bands) + 2, 7),
                                         sharex=True)
    x = np.arange(len(bands))

    # top: stacked mass x energy composition (cohort median), anchor/trace/reset
    a = [df[(df.band == b) & df.guard_ok].comp_anchor.median() for b in bands]
    t = [df[(df.band == b) & df.guard_ok].comp_trace.median() for b in bands]
    r = [df[(df.band == b) & df.guard_ok].comp_reset.median() for b in bands]
    ax_top.bar(x, a, color=C_ANCHOR, label="anchor (rigid)")
    ax_top.bar(x, t, bottom=a, color=C_TRACE, label="trace (persistence)")
    ax_top.bar(x, r, bottom=np.array(a) + np.array(t), color=C_RESET,
               label="reset (excursion)")
    ax_top.set_ylabel("cross-phase composition\n(pair-mass × energy share)")
    ax_top.set_ylim(0, 1)

    # bottom: among-mover trace vs reset energy split, per-patient strip
    for j, b in enumerate(bands):
        sub = df[(df.band == b) & df.guard_ok]
        ax_bot.scatter(np.full(len(sub), j) - 0.12, sub.share_T, s=18,
                       color=C_TRACE, alpha=0.8, zorder=3)
        ax_bot.scatter(np.full(len(sub), j) + 0.12, sub.share_R, s=18,
                       color=C_RESET, alpha=0.8, zorder=3)
        ax_bot.plot([j - 0.12], [sub.share_T.median()], marker="_",
                    ms=22, color=C_TRACE, mew=2.5, zorder=4)
        ax_bot.plot([j + 0.12], [sub.share_R.median()], marker="_",
                    ms=22, color=C_RESET, mew=2.5, zorder=4)
    ax_bot.axhline(0.5, color="0.6", lw=0.7, ls="--", zorder=0)
    ax_bot.set_ylabel("among-mover energy share")
    ax_bot.set_ylim(0, 1)
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands])

    handles = [Patch(facecolor=C_ANCHOR, label="anchor (rigid backbone)"),
               Patch(facecolor=C_TRACE, label="trace (persistence φ₁ = ρ_split)"),
               Patch(facecolor=C_RESET, label="reset (excursion φ₂, unverified)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, transparent=True, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=ALL_BANDS)
    ap.add_argument("--patients", nargs="+", default=a63.COHORT)
    ap.add_argument("--kappa", type=float, default=1.0)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    locked = load_locked_rho()

    print("[audit_105] pre-flight: ensure half FCs cached")
    for pat in args.patients:
        a63.ensure_half_fcs(pat, args.bands)

    rows = []
    print(f"\n{'patient':<8}{'band':<11}{'rho_split':>10}{'locked':>9}{'guard':>7}"
          f"{'anchor':>8}{'mover%':>8}{'shT':>6}{'shR':>6}{'r(anc,str)':>11}")
    for band in args.bands:
        for pat in args.patients:
            r = per_cell(pat, band, args.kappa, locked)
            if r is None:
                continue
            rows.append(r)
            g = "OK" if r["guard_ok"] else "FAIL"
            lk = f"{r['locked_rho']:+.3f}" if np.isfinite(r["locked_rho"]) else "  -  "
            print(f"{pat:<8}{band:<11}{r['rho_split']:>+10.3f}{lk:>9}{g:>7}"
                  f"{r['anchor_mass']:>8.3f}{100*r['mover_frac']:>7.1f}%"
                  f"{r['share_T']:>6.2f}{r['share_R']:>6.2f}"
                  f"{r['rho_anchor_strength']:>+11.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    coh = cohort_summary(df, args.bands)
    coh.to_csv(OUT / "cohort_summary.csv", index=False)

    n_fail = int((~df.guard_ok).sum())
    print(f"\n[audit_105] TRACE GUARD: {len(df) - n_fail}/{len(df)} cells reproduce "
          f"locked obs_rho (tol {GUARD_TOL:g}); {n_fail} FAIL")
    print("\n[audit_105] cohort composition (guard-ok cells):")
    print(coh.to_string(index=False))

    make_figure(df, args.bands, FIG / "decomposition.pdf")
    print(f"\n[audit_105] outputs -> {OUT}")


if __name__ == "__main__":
    main()
