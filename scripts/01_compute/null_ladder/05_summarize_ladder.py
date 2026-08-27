#!/usr/bin/env python3
"""W0-B: assemble the null ladder into the calibration tables.

Reads every rung written under ``data/paper_final/w0b_nulls/`` plus the incumbent
matched-strength gate, and emits per-band x per-scale x per-rung tables.

Reporting rules obeyed here:
  - per-scale, never best-scale: the full 16-point s-grid is written out, and
    the headline counts "scales cleared / 16", not "the best scale".
  - the cohort gate is a one-sided Wilcoxon on the MARGIN (obs - surr_p50), the
    same statistic script 13 uses, so rungs are directly comparable.
  - BH-FDR is applied across the whole (band x scale) family WITHIN a rung, since
    that is the coordinated family a reader would scan.
  - per-patient counts are reported descriptively and are never the gate.
  - leave-one-patient-out of the VERDICT: for each scale, is the gate still
    cleared with any single patient removed?

The incumbent matched-strength numbers are read from the existing
``sparsified_arc/ms_mst020`` run for CALIBRATION ONLY (to place the new rungs on
a common axis). Every number attributed to a new rung is computed here from that
rung's own per-patient CSV.
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr

OUT = ROOT / "data" / "paper_final" / "w0b_nulls"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
RUNGS = ["ms", "n1", "n2", "n1b", "n3_order", "n3_free", "n4"]
LABEL = {
    "ms": "matched-strength (incumbent, matrix-level)",
    "n1": "N1 lag-randomized coherency (circular shift, freq domain)",
    "n2": "N2 per-(chan,freq) phase randomization",
    "n1b": "N1b segment-lattice shift (shift predictor / independence)",
    "n3_order": "N3 order-preserving block rotation (EXACT, drift-aware)",
    "n3_free": "N3 free block permutation",
    "n4": "N4 within-rest placebo (duration-matched)",
}


def load_rung(rung):
    if rung == "ms":
        p = ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv"
        if not p.exists():
            return None
        d = pd.read_csv(p)
        d["rung"] = "ms"
        return d
    p = OUT / rung / "per_patient_scale.csv"
    return pd.read_csv(p) if p.exists() else None


def gate_table(d, rung):
    """Per (band, scale): cohort Wilcoxon on the margin, LOO of the verdict."""
    rows = []
    for band in sorted(d.band.unique()):
        for s in np.sort(d.s.unique()):
            x = d[(d.band == band) & np.isclose(d.s, s)].dropna(subset=["obs_rho", "surr_p50"])
            if len(x) < 5:
                continue
            m = x.obs_rho.values - x.surr_p50.values
            try:
                p = float(wilcoxon(m, alternative="greater")[1]) if np.any(m != 0) else np.nan
            except Exception:
                p = np.nan
            loo = []
            for i in range(len(m)):
                mm = np.delete(m, i)
                try:
                    loo.append(float(wilcoxon(mm, alternative="greater")[1])
                               if np.any(mm != 0) else np.nan)
                except Exception:
                    loo.append(np.nan)
            loo = np.array(loo, float)
            rows.append(dict(
                rung=rung, band=band, s=float(s), n_pat=len(x), gate_p=p,
                loo_p_max=float(np.nanmax(loo)) if loo.size else np.nan,
                loo_all_clear=bool(np.all(loo < 0.05)) if np.isfinite(loo).all() else False,
                obs_med=float(np.median(x.obs_rho)), surr_med=float(np.median(x.surr_p50)),
                margin_med=float(np.median(m)),
                n_pat_above=int((x.p < 0.05).sum()),
                obs_iqr=float(np.percentile(x.obs_rho, 75) - np.percentile(x.obs_rho, 25)),
            ))
    g = pd.DataFrame(rows)
    if g.empty:
        return g
    ok = g.gate_p.notna()
    g.loc[ok, "gate_q"] = bh_fdr(g.loc[ok, "gate_p"].values)      # whole band x scale family
    return g


def main():
    tabs = {}
    for rung in RUNGS:
        d = load_rung(rung)
        if d is None:
            print(f"[skip] {rung}: no per_patient_scale.csv", flush=True)
            continue
        tabs[rung] = gate_table(d, rung)
    if not tabs:
        print("nothing to summarize"); return

    allg = pd.concat(tabs.values(), ignore_index=True)
    allg.to_csv(OUT / "ladder_gate_all.csv", index=False)

    nS = int(allg.groupby(["rung", "band"]).size().max())
    print(f"\n{'='*100}\nNULL LADDER — scales cleared out of {nS} "
          f"(cohort one-sided Wilcoxon on obs - surr_p50)\n{'='*100}", flush=True)
    hdr = f"{'band':12s}" + "".join(f"{r:>12s}" for r in tabs)
    for metric, fn, note in [
        ("p<0.05", lambda g: (g.gate_p < 0.05).sum(), "raw p"),
        ("q<0.05", lambda g: (g.gate_q < 0.05).sum(), "BH across band x scale, within rung"),
        ("LOO-robust", lambda g: g.loo_all_clear.sum(), "gate holds dropping ANY single patient"),
    ]:
        print(f"\n--- {metric}  ({note}) ---\n{hdr}", flush=True)
        for band in BANDS:
            line = f"{band:12s}"
            for r in tabs:
                gb = tabs[r][tabs[r].band == band]
                line += f"{('-' if gb.empty else str(int(fn(gb)))):>12s}"
            print(line, flush=True)

    print(f"\n{'='*100}\nEFFECT SIZE — cohort median margin (obs - surr_p50), "
          f"min / median / max over scales\n{'='*100}", flush=True)
    print(hdr, flush=True)
    for band in BANDS:
        line = f"{band:12s}"
        for r in tabs:
            gb = tabs[r][tabs[r].band == band]
            line += ("-" if gb.empty else
                     f"{gb.margin_med.median():+.3f}").rjust(12)
        print(line, flush=True)

    print(f"\n{'='*100}\nNULL HEIGHT — cohort median surrogate rho_sym at s=1 and s=180\n{'='*100}",
          flush=True)
    print(f"{'band':12s}" + "".join(f"{r:>12s}" for r in tabs), flush=True)
    for band in BANDS:
        for lbl, pick in (("s=1", 0), ("s=180", -1)):
            line = f"{band[:9]+' '+lbl:12s}"
            for r in tabs:
                gb = tabs[r][tabs[r].band == band].sort_values("s")
                line += ("-" if gb.empty else f"{gb.surr_med.values[pick]:+.3f}").rjust(12)
            print(line, flush=True)

    # per-scale detail for the two claim bands
    for band in ("alpha", "beta"):
        print(f"\n{'='*100}\nPER-SCALE DETAIL — {band}\n{'='*100}", flush=True)
        print(f"{'s':>7s}" + "".join(f"{r+' p':>14s}" for r in tabs), flush=True)
        ss = np.sort(allg.s.unique())
        for s in ss:
            line = f"{s:7.1f}"
            for r in tabs:
                gb = tabs[r][(tabs[r].band == band) & np.isclose(tabs[r].s, s)]
                line += ("-" if gb.empty else f"{gb.gate_p.values[0]:.4f}").rjust(14)
            print(line, flush=True)

    print(f"\nwrote {OUT/'ladder_gate_all.csv'}", flush=True)
    for r in tabs:
        print(f"  {r:9s} {LABEL[r]}", flush=True)


if __name__ == "__main__":
    main()
