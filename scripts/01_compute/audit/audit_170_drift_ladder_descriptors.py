#!/usr/bin/env python3
"""audit_170 — within-recording DRIFT null for the FULL pairwise-descriptor ladder.

Extends audit_167 (raw + cophenetic only) to the whole classical descriptor bank
(node strength, weighted clustering, eigenvector centrality, PageRank, closeness) so
the drift control covers every rung audit_166 scored under matched-strength.
audit_166 found raw/eigcent/closeness "detect" the beta inference-specific component
under matched-strength; audit_167 then showed the RAW detection is drift (raw real
sits BELOW its within-rest drift floor while cophenetic sits above it; coph-vs-raw
drift-filtering advantage beta p=0.032). This asks the same of the node scalars: does
their matched-strength "detection" survive the drift control, or is it (like raw)
reproduced by within-rest drift?

Method — IDENTICAL null to audit_167, only phi(W) widens across rungs:
  drift arc = a single rest_pre recording windowed into 5 ordered segments
              (preA=w0 preB=w1 TL=w2 TT=w3 RP=w4) -> pure pre-task drift.
  real arc  = the 5 cached task phases.
  per rung  : T_test = rho_sym(A,B,TT,RP), T_infspec_pe = partial_rho(f,p|e), on
              drift and on real, per (patient, band).
  reading   : real_ispe <= drift_ispe => the rung's "detection" IS drift (disappears).
              coph advantage = coph(real-drift) - rung(real-drift) > 0 => the
              multiscale rung separates inference-from-drift better than the rung
              (bias-robust differential; the ABSOLUTE gate is biased +ve for a
              conditional trace, so the differential is the load-bearing statistic).

Optimisation vs audit_167: the imcoh cross-spectrum is nperseg-fixed and band-
independent, so it is computed ONCE per window (full + 5 windows = 6 Welch/patient)
and every band is sliced from it -- bit-identical to audit_167's per-band recompute,
~6x fewer Welch calls. Patients run in parallel (independent).

5-point preamble (inherits audit_167; same null, wider representation menu):
1. Claim: no classical low-level descriptor carries the beta inference-specific
   component once within-rest drift is controlled; only cophenetic sits above drift.
2. Null (drift): a same-length arc from pure pre-task rest reproduces the rung's
   trace. Matched-strength (audit_166) does NOT control drift -> orthogonal.
3. Strongest alt: monotonic session nonstationarity makes any later-vs-earlier
   contrast correlate for a low-level descriptor (as for raw FC).
4. Reach: rest_pre precedes any task -> a positive T on the drift arc CANNOT be
   consolidation. 1/5-length windows -> noisier FC -> drift biased LOW (conservative).
   The conditional T_infspec_pe is +ve-biased on a no-signal arc, so the ABSOLUTE
   gate is a caution only; the load-bearing statistic is the coph-vs-rung DIFFERENTIAL.
5. Falsify (of "only multiscale"): a classical rung sitting significantly ABOVE its
   drift floor in beta where raw does not -> reported honestly.

Reads  : rest_pre timeseries + cached phase FCs (audit_63 loaders).
Writes : data/audit/drift_ladder_descriptors/{per_patient_per_band,cohort_summary}.csv + README
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.metrics.graph_descriptors import (
    closeness_centrality,
    eigenvector_centrality,
    node_strength,
    pagerank,
    raw_edges,
    weighted_clustering_onnela,
)
from lrg_eegfc.utils.surrogate.matched_strength import (
    cophenetic_condensed_from_adjacency as coph,
)

# reuse audit_167's band-extract, condensed-triu, rho_sym arc verbatim (no fork)
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_167_drift_null_ladder import (  # type: ignore
    BANDS,
    COHORT,
    NWIN,
    PH,
    _arc,
    _band_abs,
    _iu,
    _rho,
)
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore

RUNGS = [
    ("raw_fc", raw_edges),
    ("strength", node_strength),
    ("clustering", weighted_clustering_onnela),
    ("eigcent", eigenvector_centrality),
    ("pagerank", pagerank),
    ("closeness", closeness_centrality),
    ("cophenetic", coph),
]
OUT = ROOT / "data" / "audit" / "drift_ladder_descriptors"


def _phi_arc(phi, mats):
    """Apply descriptor phi to the 5 phase FC matrices -> (T_test, T_infspec_pe)."""
    try:
        D = {ph: np.asarray(phi(mats[ph]), float) for ph in PH}
        if any(not np.all(np.isfinite(D[ph])) for ph in PH):
            return np.nan, np.nan
        return _arc(D)
    except Exception:  # noqa: BLE001
        return np.nan, np.nan


def per_patient(pat):
    """All bands for one patient. Welch computed ONCE per window, sliced per band."""
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nper = nperseg_for_fs(fs)
    T = X.shape[1]
    # one cross-spectrum per recording: full + 5 ordered windows
    fr_full, C_full = compute_msc_welch(X, fs, nperseg=nper, metric="imcoh")
    win_specs = []
    for i in range(NWIN):
        Xi = X[:, i * T // NWIN:(i + 1) * T // NWIN]
        win_specs.append(compute_msc_welch(Xi, fs, nperseg=min(nper, Xi.shape[1]),
                                            metric="imcoh"))
    rows = []
    for band in BANDS:
        flo, fhi = BRAIN_BANDS[band]
        Ws = [_band_abs(C, fr, flo, fhi) for fr, C in win_specs]
        full = _band_abs(C_full, fr_full, flo, fhi)
        cached_pre = load_phase_fc(pat, "rest_pre", band)
        val = (_rho(_iu(full), _iu(cached_pre))
               if cached_pre.shape == full.shape else np.nan)
        drift_mats = {ph: Ws[i] for i, ph in enumerate(PH)}
        real_mats = {ph: load_phase_fc(pat, ph, band) for ph in PH}
        for name, phi in RUNGS:
            tt_dr, ip_dr = _phi_arc(phi, drift_mats)
            tt_re, ip_re = _phi_arc(phi, real_mats)
            rows.append(dict(patient=pat, band=band, descriptor=name, val=val,
                             drift_Ttest=tt_dr, real_Ttest=tt_re,
                             drift_ispe=ip_dr, real_ispe=ip_re))
    return rows


def _wp(x, alt="greater"):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    if np.nanmax(np.abs(x)) < 1e-9:
        return 1.0
    try:
        return float(wilcoxon(x, alternative=alt).pvalue)
    except Exception:  # noqa: BLE001
        return float("nan")


def cohort_summary(df):
    out = []
    coph_df = df[df.descriptor == "cophenetic"].set_index(["patient", "band"])
    for name, _ in RUNGS:
        d = df[df.descriptor == name]
        for band in BANDS:
            s = d[d.band == band]
            if s.empty:
                continue
            row = {"descriptor": name, "band": band, "n_pat": len(s)}
            for tgt in ("Ttest", "ispe"):
                real = s[f"real_{tgt}"].values
                drift = s[f"drift_{tgt}"].values
                row[f"{tgt}_real_med"] = float(np.nanmedian(real))
                row[f"{tgt}_drift_med"] = float(np.nanmedian(drift))
                row[f"{tgt}_real_gt_drift_p"] = _wp(real - drift)
                row[f"{tgt}_n_real_gt_drift"] = int(np.nansum(real > drift))
                row[f"{tgt}_drift_alone_p"] = _wp(drift)
            if name != "cophenetic":
                m = s.set_index(["patient", "band"]).join(
                    coph_df[["real_ispe", "drift_ispe"]], rsuffix="_coph")
                adv = ((m["real_ispe_coph"] - m["drift_ispe_coph"])
                       - (m["real_ispe"] - m["drift_ispe"])).values
                row["coph_adv_ispe_p"] = _wp(adv)
                row["coph_adv_ispe_n"] = int(np.nansum(adv > 0))
            else:
                row["coph_adv_ispe_p"] = np.nan
                row["coph_adv_ispe_n"] = 0
            out.append(row)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="first N patients only")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    cohort = COHORT[:args.limit] if args.limit else COHORT
    ncpu = min(len(cohort), max(1, (os.cpu_count() or 4) - 2))
    print(f"[audit_170] within-rest DRIFT ladder, {len(RUNGS)} rungs, "
          f"{len(cohort)} patients x {len(BANDS)} bands, {ncpu} workers\n", flush=True)

    rows, t0 = [], time.time()
    with Pool(ncpu) as pool:
        for i, pr in enumerate(pool.imap_unordered(per_patient, cohort), 1):
            rows.extend(pr)
            el = time.time() - t0
            print(f"[{i}/{len(cohort)}] {pr[0]['patient']} done  {el:.0f}s "
                  f"ETA {el/i*(len(cohort)-i):.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    summ = cohort_summary(df)
    summ.to_csv(OUT / "cohort_summary.csv", index=False)

    print("\n=== BETA inference-specific: real vs within-rest drift, per rung ===")
    print(f"{'rung':<12}{'real':>8}{'drift':>8}{'r>d p':>8}{'n>':>4}"
          f"{'drift0 p':>9}   {'coph_adv p':>11}{'n':>4}")
    b = summ[summ.band == "beta"].set_index("descriptor")
    for name, _ in RUNGS:
        if name not in b.index:
            continue
        r = b.loc[name]
        adv = "" if name == "cophenetic" else f"{r.coph_adv_ispe_p:>11.3f}{int(r.coph_adv_ispe_n):>4}"
        print(f"{name:<12}{r.ispe_real_med:>+8.3f}{r.ispe_drift_med:>+8.3f}"
              f"{r.ispe_real_gt_drift_p:>8.3f}{int(r.ispe_n_real_gt_drift):>4}"
              f"{r.ispe_drift_alone_p:>9.3f}   {adv}")
    print("\nREAD: real<=drift or r>d p large => rung's inference 'detection' is DRIFT (disappears).")
    print("      coph_adv p<0.05 => cophenetic separates inference-from-drift better than the rung.")
    _write_readme(summ, time.time() - t0, len(cohort))
    print(f"\n[audit_170] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}")


def _write_readme(summ, rt, npat):
    b = summ[summ.band == "beta"].set_index("descriptor")
    L = ["---", "name: drift_ladder_descriptors",
         "scope: within-rest drift null for the full pairwise-descriptor ladder (extends audit_167)",
         "date: 2026-07-10", "status: current", "---", "",
         "# Drift ladder — does any low-level descriptor's inference-specific trace",
         "# survive the within-rest drift control? (extends audit_167 raw+coph)", "",
         "**Head.** Same within-rest 5-window drift null as audit_167, applied to every",
         "audit_166 rung. A rung whose real inference-specific value sits at/below its",
         "drift floor is drift, not consolidation. Load-bearing statistic = the",
         "coph-vs-rung differential (subtracts the shared conditional-functional bias).", "",
         "## Beta inference-specific (T_infspec_pe)", "",
         "| rung | real | drift | real>drift p | n>/10 | drift-alone p | coph adv p |",
         "|---|---|---|---|---|---|---|"]
    for name, _ in RUNGS:
        if name not in b.index:
            continue
        r = b.loc[name]
        adv = "--" if name == "cophenetic" else f"{r.coph_adv_ispe_p:.3f}"
        L.append(f"| {name} | {r.ispe_real_med:+.3f} | {r.ispe_drift_med:+.3f} | "
                 f"{r.ispe_real_gt_drift_p:.3f} | {int(r.ispe_n_real_gt_drift)} | "
                 f"{r.ispe_drift_alone_p:.3f} | {adv} |")
    L += ["", "## Reading",
          "- real <= drift (or real>drift p large) => the rung's inference detection is DRIFT.",
          "- coph adv p < 0.05 => cophenetic separates inference-from-drift better than the rung",
          "  (bias-robust differential; the ABSOLUTE drift gate is biased +ve for a conditional trace).",
          "- matched-strength (audit_166) is ORTHOGONAL to drift; a rung can clear MS yet be drift.", "",
          "## Provenance",
          "- null: within-rest_pre 5-window drift arc (audit_167 machinery, verbatim).",
          "- estimator: rho_sym arc (audit_152); positive = trace.",
          "- descriptors: src/lrg_eegfc/utils/metrics/graph_descriptors.py.",
          f"- cohort n={npat}, Wilcoxon one-sided; runtime {rt:.0f}s.",
          "- Build: scripts/01_compute/audit/audit_170_drift_ladder_descriptors.py."]
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
