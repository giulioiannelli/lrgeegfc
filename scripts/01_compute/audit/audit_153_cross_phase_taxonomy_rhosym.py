#!/usr/bin/env python3
"""audit_153 — cross-phase cophenetic taxonomy (anchor/trace/reset) under rho_sym.

Migrates the taxonomy decomposition (audit_105) off bare rho_split onto the canonical,
arbitrary-half-free rho_sym estimator (feedback_rho_sym_canonical_estimator). The
taxonomy channels are functionals of the split-baseline cophenetic distances, so which
rest_pre half is the "primary" (task-side) baseline and which is the rest-side baseline
is the same arbitrary choice rho_sym removes. Both equally-valid arm assignments are
computed and their compositions averaged:

    arm1 (A primary / task-side, B rest-side):  cross_phase_channels(DA, DB, DT, DP)
    arm2 (B primary / task-side, A rest-side):  cross_phase_channels(DB, DA, DT, DP)
    composition = 1/2 [ share(arm1) + share(arm2) ]

The primary reported trace scalar is rho_sym = 1/2[rho(DT-DA,DP-DB)+rho(DT-DB,DP-DA)].
The trace GUARD is kept bit-exact on ARM1 rho_split vs the locked audit_63 CSV (the
pipeline-integrity anchor audit_105 already validates); rho_sym is derived on top of
that verified pipeline and cross-checked (loose) against the audit_150 gate output.
Does NOT edit audit_105; new script, new output dir.

5-point critical preamble
1. Claim: the composition (anchor dominates; trace = phi1 persistence peaks at
   beta/alpha) holds under rho_sym as under rho_split.
2. Null: the composition depended on the arbitrary A/B primary choice.
3. Strongest alternative: rho_sym trivially relabels rho_split. Controlled by the
   bit-exact arm1 guard vs locked audit_63 and by reporting both arms' compositions.
4. Cannot: the reset channel remains UNVERIFIED here (its coupled null is a separate
   step, retracted for the 3-way dissociation — only trace->OFC survives, audit_151);
   this script re-states the decomposition + anchor dominance only.
5. Falsify: if arm1 rho_split stops reproducing the locked audit_63 obs_rho, the
   pipeline drifted and the numbers are untrusted.

Outputs data/audit/cross_phase_taxonomy_rhosym/:
  per_patient_per_band.csv   per-cell sym composition, rho_sym, arm1 guard, strength diag
  cohort_summary.csv         per-band cohort composition (median + IQR)
  figures/decomposition.pdf  anchor dominance + trace/reset split per band
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

# audit_63 FC->cophenet pipeline imported verbatim (no fork) — bit-exact guard anchor
_spec = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "cross_phase_taxonomy_rhosym"
FIG = OUT / "figures"
LOCKED_CSV = (ROOT / "data" / "audit"
              / "matched_strength_surrogate_split_baseline"
              / "per_patient_per_band.csv")
# rho_sym reference (audit_150) — loose cross-check only (different observed pipeline)
RHOSYM_CSV = ROOT / "data" / "audit" / "rho_sym_gate" / "per_patient_per_band.csv"

ALL_BANDS = list(BRAIN_BANDS)
GUARD_TOL = 1e-6

C_ANCHOR = "#807dba"
C_TRACE = "#2c7fb8"
C_RESET = "#41ab5d"

_SHARE_KEYS = ("anchor_mass", "mover_frac", "share_T", "share_R", "E_T", "E_R",
               "comp_anchor", "comp_trace", "comp_reset")


def load_locked_rho() -> dict:
    if not LOCKED_CSV.exists():
        print(f"[audit_153] WARNING: locked CSV missing at {LOCKED_CSV}; guard disabled")
        return {}
    df = pd.read_csv(LOCKED_CSV)
    return {(r.patient, r.band): float(r.obs_rho) for r in df.itertuples()}


def load_rhosym_ref() -> dict:
    if not RHOSYM_CSV.exists():
        return {}
    df = pd.read_csv(RHOSYM_CSV)
    return {(r.patient, r.band): float(r.obs_rho) for r in df.itertuples()}


def per_cell(pat: str, band: str, kappa: float, locked: dict, rhosym_ref: dict) -> dict | None:
    Ws = {}
    for phase in a63.PHASES:
        try:
            Ws[phase] = a63.load_phase_fc(pat, phase, band)
        except Exception as e:
            print(f"[audit_153] SKIP {pat}/{band}/{phase}: {e}")
            return None
    N = Ws[a63.PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        print(f"[audit_153] SKIP {pat}/{band}: phase shape mismatch")
        return None

    D = {ph: a63.lrg_ultrametric_condensed(Ws[ph]) for ph in a63.PHASES}
    Lc = D[a63.PHASES[0]].size
    if any(d.size != Lc for d in D.values()):
        print(f"[audit_153] SKIP {pat}/{band}: ultrametric size mismatch")
        return None

    DA, DB = D["rest_pre_A"], D["rest_pre_B"]
    DT, DP = D["task_test"], D["rest_post"]

    # 1. trace scalars: arm1 = rho_split (bit-exact guard); rho_sym = mean of both arms
    rho_split = float(spearmanr(DT - DA, DP - DB).statistic)   # arm1
    rho_ba = float(spearmanr(DT - DB, DP - DA).statistic)      # arm2
    rho_sym = 0.5 * (rho_split + rho_ba)
    locked_rho = locked.get((pat, band), np.nan)
    guard_ok = (not np.isfinite(locked_rho)) or abs(rho_split - locked_rho) < GUARD_TOL
    rho_sym_ref = rhosym_ref.get((pat, band), np.nan)

    # 2. symmetric composition: average both arm assignments
    ch_ab = cross_phase_channels(DA, DB, DT, DP, kappa=kappa)   # A primary
    ch_ba = cross_phase_channels(DB, DA, DT, DP, kappa=kappa)   # B primary
    row_ab, row_ba = cohort_share_row(ch_ab), cohort_share_row(ch_ba)
    row = {k: 0.5 * (row_ab[k] + row_ba[k]) for k in _SHARE_KEYS}
    sigma = ch_ab["sigma"]   # RMS|DA-DB| is arm-invariant

    # 4. strength diagnostic (anchor vs endpoint hubness), averaged over arms
    s_node = Ws["rest_pre_A"].sum(axis=1)
    iu_i, iu_j = np.triu_indices(N, k=1)
    pair_strength = s_node[iu_i] + s_node[iu_j]
    ras = 0.5 * (float(spearmanr(ch_ab["anchor_weight"], pair_strength).statistic)
                 + float(spearmanr(ch_ba["anchor_weight"], pair_strength).statistic))

    return {
        "patient": pat, "band": band, "N_nodes": int(N), "n_pairs": int(Lc),
        "rho_sym": rho_sym, "rho_split": rho_split, "rho_ba": rho_ba,
        "locked_rho": locked_rho, "guard_ok": guard_ok, "rho_sym_ref": rho_sym_ref,
        "sigma": sigma,
        "anchor_mass": row["anchor_mass"], "mover_frac": row["mover_frac"],
        "share_T": row["share_T"], "share_R": row["share_R"],
        "E_T": row["E_T"], "E_R": row["E_R"],
        "comp_anchor": row["comp_anchor"], "comp_trace": row["comp_trace"],
        "comp_reset": row["comp_reset"],
        "rho_anchor_strength": ras,
    }


def cohort_summary(df: pd.DataFrame, bands: list[str]) -> pd.DataFrame:
    out = []
    for band in bands:
        sub = df[(df.band == band) & df.guard_ok]
        if sub.empty:
            continue
        out.append({
            "band": band, "n_patients": len(sub),
            "rho_sym_med": float(sub.rho_sym.median()),
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
    a = [df[(df.band == b) & df.guard_ok].comp_anchor.median() for b in bands]
    t = [df[(df.band == b) & df.guard_ok].comp_trace.median() for b in bands]
    r = [df[(df.band == b) & df.guard_ok].comp_reset.median() for b in bands]
    ax_top.bar(x, a, color=C_ANCHOR, label="anchor (rigid)")
    ax_top.bar(x, t, bottom=a, color=C_TRACE, label="trace (persistence)")
    ax_top.bar(x, r, bottom=np.array(a) + np.array(t), color=C_RESET, label="reset")
    ax_top.set_ylabel("cross-phase composition\n(pair-mass × energy share, ρ_sym)")
    ax_top.set_ylim(0, 1)
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
               Patch(facecolor=C_TRACE, label="trace (persistence φ₁, ρ_sym)"),
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
    rhosym_ref = load_rhosym_ref()

    print("[audit_153] pre-flight: ensure half FCs cached")
    for pat in args.patients:
        a63.ensure_half_fcs(pat, args.bands)

    rows = []
    print(f"\n{'patient':<8}{'band':<11}{'rho_sym':>9}{'rho_split':>10}{'guard':>7}"
          f"{'anchor':>8}{'mover%':>8}{'shT':>6}{'shR':>6}{'r(anc,str)':>11}")
    for band in args.bands:
        for pat in args.patients:
            r = per_cell(pat, band, args.kappa, locked, rhosym_ref)
            if r is None:
                continue
            rows.append(r)
            g = "OK" if r["guard_ok"] else "FAIL"
            print(f"{pat:<8}{band:<11}{r['rho_sym']:>+9.3f}{r['rho_split']:>+10.3f}{g:>7}"
                  f"{r['anchor_mass']:>8.3f}{100*r['mover_frac']:>7.1f}%"
                  f"{r['share_T']:>6.2f}{r['share_R']:>6.2f}"
                  f"{r['rho_anchor_strength']:>+11.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    coh = cohort_summary(df, args.bands)
    coh.to_csv(OUT / "cohort_summary.csv", index=False)

    n_fail = int((~df.guard_ok).sum())
    print(f"\n[audit_153] TRACE GUARD (arm1 rho_split vs locked audit_63): "
          f"{len(df) - n_fail}/{len(df)} cells bit-exact (tol {GUARD_TOL:g}); {n_fail} FAIL")
    # loose cross-check rho_sym vs audit_150
    ok = df.dropna(subset=["rho_sym_ref"])
    if not ok.empty:
        md = float((ok.rho_sym - ok.rho_sym_ref).abs().max())
        print(f"[audit_153] rho_sym vs audit_150 gate: n={len(ok)} cells, "
              f"max|Δ|={md:.2e} (loose; different observed pipeline)")
    print("\n[audit_153] cohort composition under rho_sym (guard-ok cells):")
    print(coh.to_string(index=False))

    make_figure(df, args.bands, FIG / "decomposition.pdf")
    print(f"\n[audit_153] outputs -> {OUT}")


if __name__ == "__main__":
    main()
