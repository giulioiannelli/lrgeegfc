#!/usr/bin/env python3
"""audit_152 — four-phase cophenetic consolidation arc under the rho_sym estimator.

Migrates the N2 climax (audit_103) off bare rho_split onto the canonical, arbitrary-
half-free rho_sym estimator (feedback_rho_sym_canonical_estimator). Every arc
functional is a Spearman between a task-side contrast (baselined on one rest_pre
half) and the persistent rest-side contrast (baselined on the OTHER half). Bare
rho_split hard-codes which half feeds which arm; rho_sym averages the two equally-
valid assignments:

    for the split-baseline pair (A, B):
      arm1 (A->task-side, B->rest-side):  e =D_TL-D_A  g =D_TT-D_A  p =D_RP-D_B
      arm2 (B->task-side, A->rest-side):  e'=D_TL-D_B  g'=D_TT-D_B  p'=D_RP-D_A
      f = D_TT - D_TL   (inference-specific; NO rest_pre half -> arm-invariant)
    T_x^sym = 1/2 [ rho(arm1) + rho(arm2) ]

applied IDENTICALLY to observed and to every matched-strength surrogate row. This
does NOT edit audit_103; new script, new output dir. Reuses audit_83's cached
surrogate cophenetic stacks verbatim (preA/preB/task_test/rest_post = canonical
seed-20260511 cache hits; task_learn = audit_103 --null cache hit) so NO surrogate
is regenerated — only the arm recombination changes.

Notation (see scope .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md)
  e = encoding    g = inference-online    f = inference-specific    p = persistent
  T_test    = rho(g,p)   <- reproduces audit_83 full-graph obs_stat (arm1 = CROSS-CHECK)
  T_learn   = rho(e,p)   T_infspec = rho(f,p)
  T_infspec_pe = partial rho(f,p | e)   (inference-specific | encoding)
Positive = TRACE (locked sign convention).

5-point critical preamble
1. Claim: the arc verdicts (beta T_infspec_pe clears matched-strength = inference-
   specific consolidation; alpha/beta T_learn = learning phase leaves its own trace)
   hold under rho_sym exactly as under rho_split.
2. Null: the arc depended on the arbitrary A/B half assignment (the estimator, not
   the biology).
3. Strongest alternative: rho_sym is rho_split relabelled and reproduces it trivially.
   Controlled by recording BOTH the sym functionals AND arm1 (=rho_split) per cell,
   the arm1 T_test cross-check against the locked audit_83 full-graph obs_stat, and
   the SAME strength-preserving surrogate audit_103 uses (estimator change is
   orthogonal to strength).
4. Cannot: reuses R=200 cached surrogate; n=10 cohort Wilcoxon coarse p-grid; f
   (inference-specific contrast) is arm-invariant so T_infspec's symmetrization only
   averages the persistent-axis half — a weaker de-biasing than T_test/T_learn (both
   arms differ there).
5. Falsify: if beta T_infspec_pe or alpha/beta T_learn fails matched-strength under
   rho_sym, the arc is estimator-dependent.

Output: data/audit/consolidation_arc_rhosym/
  arc_per_patient.csv        observed sym functionals + arm1 (rho_split) cross-check
  arc_null_per_patient.csv   sym obs + surrogate p per functional (--null)
  README.md
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
from audit_83_wm_stratified_cophenetic import _surr_stacks  # type: ignore
import _wm_stratify as ws  # type: ignore

ARC_PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
FUNCTIONALS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")

# Mirror audit_103's SeedSequence EXACTLY so the 4 canonical phases are byte-identical
# cache hits and task_learn (index 4) reuses the audit_103 --null ensemble.
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

OUT = ROOT / "data" / "audit" / "consolidation_arc_rhosym"
OUT.mkdir(parents=True, exist_ok=True)

AUDIT83_CSV = ROOT / "data" / "audit" / "wm_stratified" / "cophenetic_raw_per_patient.csv"

# arm assignments: (task-side baseline, rest-side baseline)
_ARMS = (("rest_pre_A", "rest_pre_B"), ("rest_pre_B", "rest_pre_A"))


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    r, _ = spearmanr(a, b)
    return float(r)


def _partial_rho(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """First-order partial Spearman of (a, b) controlling c (rank world)."""
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    denom = np.sqrt(max(0.0, (1.0 - rac**2) * (1.0 - rbc**2)))
    return (rab - rac * rbc) / denom if denom > 0 else np.nan


def _arm_functionals(D_TL, D_TT, D_RP, D_bt, D_br) -> dict:
    """Functionals for ONE arm: D_bt = task-side baseline, D_br = rest-side baseline."""
    e = D_TL - D_bt          # encoding
    g = D_TT - D_bt          # inference-online
    f = D_TT - D_TL          # inference-specific (arm-invariant)
    p = D_RP - D_br          # persistent
    return {"T_test": _rho(g, p), "T_learn": _rho(e, p),
            "T_infspec": _rho(f, p), "T_infspec_pe": _partial_rho(f, p, e)}


def _sym_functionals(D: dict) -> tuple[dict, dict]:
    """(sym, arm1). sym = mean over the two arm assignments; arm1 = rho_split arm."""
    arm1 = _arm_functionals(D["task_learn"], D["task_test"], D["rest_post"],
                            D["rest_pre_A"], D["rest_pre_B"])
    arm2 = _arm_functionals(D["task_learn"], D["task_test"], D["rest_post"],
                            D["rest_pre_B"], D["rest_pre_A"])
    sym = {k: 0.5 * (arm1[k] + arm2[k]) for k in FUNCTIONALS}
    return sym, arm1


def _arc_row(pat: str, band: str) -> dict:
    D = {x: lrg_ultrametric_condensed(load_phase_fc(pat, x, band)) for x in ARC_PHASES}
    sym, arm1 = _sym_functionals(D)
    row = {"patient": pat, "band": band, "n_pairs": int(D["rest_pre_A"].size)}
    for k in FUNCTIONALS:
        row[k] = sym[k]                 # primary = rho_sym
    row["T_test_split"] = arm1["T_test"]  # arm1 = rho_split cross-check vs audit_83
    return row


def main() -> None:
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                rows.append(_arc_row(pat, band))
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_152] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
                rows.append({"patient": pat, "band": band, "error": str(exc)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_per_patient.csv", index=False)
    print(f"\n[audit_152] wrote {OUT/'arc_per_patient.csv'}  ({len(df)} rows)\n")

    # CROSS-CHECK: arm1 T_test_split must reproduce audit_83 full-graph cophenetic obs.
    if AUDIT83_CSV.exists() and "T_test_split" in df:
        ref = pd.read_csv(AUDIT83_CSV)
        ref = ref[(ref["substrate"] == "cophenetic") & (ref["config"] == "full")]
        ref = ref.set_index(["patient", "band"])["obs_stat"]
        merged = df.dropna(subset=["T_test_split"]).set_index(["patient", "band"])
        diffs = [abs(r["T_test_split"] - float(ref.loc[idx]))
                 for idx, r in merged.iterrows() if idx in ref.index]
        if diffs:
            md = max(diffs)
            status = "PASS" if md < 1e-6 else ("CLOSE" if md < 0.02 else "FAIL")
            print(f"[CROSS-CHECK] arm1 T_test vs audit_83 full-graph: n={len(diffs)} "
                  f"cells, max|Δ|={md:.2e}  → {status}")
            if status == "FAIL":
                print("  !! arm1 does NOT reproduce audit_83 — pipeline drift, numbers suspect.")
    else:
        print("[CROSS-CHECK] audit_83 csv not found — skipped.")

    print("\n=== COHORT ARC SUMMARY under rho_sym (median [IQR]; k/10 positive) ===")
    cols = ["T_test", "T_learn", "T_infspec", "T_infspec_pe"]
    print(f"{'band':<11}" + "".join(f"{c:>24}" for c in cols))
    for band in BANDS:
        sub = df[df["band"] == band]
        cells = []
        for c in cols:
            v = sub[c].dropna().values if c in sub else np.array([])
            if v.size:
                q1, q3 = np.percentile(v, [25, 75])
                cells.append(f"{np.median(v):+.2f}[{q1:+.2f},{q3:+.2f}] {int((v>0).sum())}/{v.size}")
            else:
                cells.append("--")
        print(f"{band:<11}" + "".join(f"{c:>24}" for c in cells))


# ---------------------------------------------------------------------------
# Stage 2 — matched-strength surrogate null on the SYM arc functionals.
# ---------------------------------------------------------------------------
def _surr_functionals_sym(sc: dict, R: int) -> dict:
    """Per-surrogate SYM realizations of each functional from cophenetic stacks
    ``sc[phase]`` (R, n_pairs). Pairs surrogate index r across phases exactly as
    audit_103; strength-violating (NaN) rows skipped. Both arms recombine the SAME
    5 phase stacks."""
    out = {k: [] for k in FUNCTIONALS}
    for r in range(R):
        cols = {ph: sc[ph][r] for ph in ARC_PHASES}
        if not all(np.isfinite(v).all() for v in cols.values()):
            continue
        arm1 = _arm_functionals(cols["task_learn"], cols["task_test"], cols["rest_post"],
                                cols["rest_pre_A"], cols["rest_pre_B"])
        arm2 = _arm_functionals(cols["task_learn"], cols["task_test"], cols["rest_post"],
                                cols["rest_pre_B"], cols["rest_pre_A"])
        for k in FUNCTIONALS:
            out[k].append(0.5 * (arm1[k] + arm2[k]))
    return {k: np.asarray(v, float) for k, v in out.items()}


def run_null(verbose: bool = True) -> None:
    R, sf = ws.N_SURROGATES, ws.SWAP_FACTOR
    print(f"[audit_152 null] R={R} swap_factor={sf}  rho_sym  "
          f"(preA/preB/task_test/rest_post + task_learn cached)\n")
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
                D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in ARC_PHASES}
                obs, _arm1 = _sym_functionals(D)
                sc = {}
                for ph in ARC_PHASES:
                    coph, _raw = _surr_stacks(
                        Ws[ph], ws.surr_eig_path("full", pat, band, ph, R, sf),
                        _arc_cell_rng(pat, band, ph), R, sf, False)
                    sc[ph] = coph
                surr = _surr_functionals_sym(sc, R)
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
                print(f"[audit_152 null] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_null_per_patient.csv", index=False)
    print(f"\n[audit_152 null] wrote {OUT/'arc_null_per_patient.csv'} ({len(df)} rows)\n")

    print("=== COHORT MATCHED-STRENGTH VERDICT under rho_sym (positive = trace) ===")
    print(f"{'band':<11}{'functional':<15}{'obs_med':>9}{'surr_med':>9}"
          f"{'n>surr':>8}{'wilcox_p':>10}{'LOp15_p':>9}")
    verdict_rows = []
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
            verdict_rows.append({"band": band, "functional": k,
                                 "med_obs": v["med_obs"], "med_surr": v["med_surr"],
                                 "n_above": n_above, "n_defined": v["n_defined"],
                                 "wilcoxon_p": v["wilcoxon_p"], "wilcoxon_p_LOp15": v2["wilcoxon_p"]})
    pd.DataFrame(verdict_rows).to_csv(OUT / "arc_cohort_verdict.csv", index=False)
    print("\n[read] T_infspec_pe = inference-specific persistence controlling encoding "
          "(the beta-only claim); T_learn = learning-phase own trace (alpha/beta).")
    write_readme(pd.DataFrame(verdict_rows))


def write_readme(verdict: pd.DataFrame) -> None:
    L = ["---", "name: consolidation_arc_rhosym",
         "scope: four_phase_cophenetic_consolidation_arc_under_rho_sym",
         "date: 2026-07-06", "status: current", "---", "",
         "# Consolidation arc under rho_sym (canonical estimator)", "",
         "**Head.** N2 climax migrated off bare rho_split onto rho_sym (mean of both",
         "split-half arm assignments; feedback_rho_sym_canonical_estimator). Same",
         "audit_103 matched-strength surrogate (R=200, cached; NO regeneration), same",
         "cross-check to audit_83 full-graph. Verdict per band per functional below.", ""]
    if not verdict.empty:
        L += ["| band | functional | obs_med | surr_med | n>surr | wilcox_p | LO-P15_p |",
              "|---|---|---|---|---|---|---|"]
        for _, r in verdict.iterrows():
            L.append(f"| {r['band']} | {r['functional']} | {r['med_obs']:+.3f} | "
                     f"{r['med_surr']:+.3f} | {int(r['n_above'])}/{int(r['n_defined'])} | "
                     f"{r['wilcoxon_p']:.4f} | {r['wilcoxon_p_LOp15']:.4f} |")
    L += ["", "## Provenance",
          "- estimator: rho_sym (mean of both A/B arm assignments); arm1 = rho_split cross-check.",
          "- f (inference-specific) has NO rest_pre half -> arm-invariant; only its p axis symmetrizes.",
          f"- surrogate: audit_83 _surr_stacks, R={ws.N_SURROGATES} swap_factor={ws.SWAP_FACTOR}, cached.",
          "- Build: scripts/01_compute/audit/audit_152_consolidation_arc_rhosym.py",
          "- Does NOT modify audit_103; supersedes it as estimator of record."]
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", action="store_true",
                    help="run matched-strength surrogate null on the sym functionals")
    args = ap.parse_args()
    main()          # observed pass + cross-check always runs (fast)
    if args.null:
        run_null()
