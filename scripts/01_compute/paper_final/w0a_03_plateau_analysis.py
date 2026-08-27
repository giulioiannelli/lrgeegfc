#!/usr/bin/env python3
"""W0-A / A2-A3 — read the stability surface and apply the frozen plateau rule.

Pure analysis of the CSVs written by w0a_02_sparsification_stability.py. No new
compute, no new choices: the plateau definition applied here is verbatim R2 as
frozen in w0a_00_preregistration.md before any number was seen.

  R2.2 PLATEAU (frozen). For band b and functional T, the maximal contiguous
  set of f over which the PER-SCALE cohort gate verdict (gate_p < 0.05 at each
  of the 16 scales -- no best-scale collapse) is CONSTANT, requiring
  f_max / f_min >= 2 AND requiring the verdict to be reproduced by at least two
  different density MECHANISMS at matched density.

  R2.1 A parameter-free filter is adopted iff its (density, verdict) point lies
  INSIDE that plateau.
  R2.3 Otherwise the readout is knob-integrated: the median over the plateau,
  with the plateau width stated and no single f reported.
  R2.4 If no octave-wide invariant region exists for a band, that band's
  verdict is KNOB-DEPENDENT and is reported as such, not rescued by picking f.

The strict rule demands an exact 16-bit verdict vector match, which is a harsh
criterion -- one scale flipping breaks a plateau. The strict answer is reported
FIRST and is the gate. A graded stability profile (per-scale agreement with a
reference config, Hamming distance across the knob) is reported after it as
description, explicitly labelled as post-hoc and not part of the rule.

Note on multiplicity: the plateau is an INVARIANCE statement evaluated under a
fixed decision rule applied identically to every config. It is not itself a
significance claim, and the alpha=0.05 per-scale threshold here is a fixed
yardstick, not a corrected test. The cohort gate specification (BH across the
band x scale family, LOO, bootstrap CI) belongs to lane W0-C.

Outputs (data/paper_final/w0a_substrate/a2_stability/{transform}/):
  plateaus.csv           functional, band, f_lo, f_hi, octaves, verdict_vector,
                         n_scales_cleared, cross_mechanism_ok
  parameter_free.csv     functional, band, filter, density, inside_plateau, ...
  knob_integrated.csv    functional, band, plateau median of the statistic
  stability_profile.csv  graded agreement of every config with the reference
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

FUNCTIONALS = ("T_probe", "T_encode", "T_probespec", "T_probespec_pe")
PARAM_FREE = ("tmfg", "perc")
ALPHA = 0.05
MIN_OCTAVES = 1.0                      # f_max / f_min >= 2


def verdict_vector(g, config, functional, band):
    """16-bit per-scale cleared/not vector for one (config, functional, band)."""
    x = g[(g.config == config) & (g.functional == functional) & (g.band == band)]
    if x.empty:
        return None
    x = x.sort_values("s")
    return tuple((x.gate_p.values < ALPHA).astype(int))


def vstr(v):
    return "".join(str(b) for b in v) if v is not None else "-" * 16


def find_plateaus(g, functional, band, family="mst"):
    """Maximal contiguous runs of identical verdict vector along the f grid."""
    fam = g[(g.method == family) & (g.functional == functional) & (g.band == band)]
    if fam.empty:
        return []
    fracs = sorted(fam.param.unique())
    vecs = [verdict_vector(g, f"{family}@{f:g}", functional, band) for f in fracs]
    runs, i = [], 0
    while i < len(fracs):
        j = i
        while j + 1 < len(fracs) and vecs[j + 1] == vecs[i]:
            j += 1
        runs.append(dict(f_lo=fracs[i], f_hi=fracs[j], n_points=j - i + 1,
                         verdict=vecs[i],
                         octaves=float(np.log2(fracs[j] / fracs[i])) if fracs[i] > 0 else 0.0))
        i = j + 1
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transform", default=os.environ.get("W0A_TRANSFORM", "abs"))
    a = ap.parse_args()
    OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a2_stability" / a.transform
    gate = pd.read_csv(OUT / "cohort_gate.csv")
    st = pd.read_csv(OUT / "structure.csv")
    per = pd.read_csv(OUT / "per_patient_scale.csv")
    bands = list(dict.fromkeys(gate.band))
    dens = st.groupby("config").density.median()
    ncomp = st.groupby("config").n_components.median()
    print(f"=== W0-A / A2 stability surface — transform = imcoh_{a.transform} ===", flush=True)
    print(f"configs={gate.config.nunique()} bands={len(bands)} "
          f"functionals={gate.functional.nunique()} scales={gate.s.nunique()} "
          f"patients={per.patient.nunique()}\n", flush=True)

    # ------------------------------------------------------------------ #
    # 1. the surface: scales cleared, per config x band x functional
    # ------------------------------------------------------------------ #
    for func in FUNCTIONALS:
        gF = gate[gate.functional == func]
        if gF.empty:
            continue
        print(f"--- [{func}] scales cleared (gate_p<{ALPHA}) of 16 | verdict bitstring ---",
              flush=True)
        print("  " + "config".ljust(15) + "dens  ncmp  "
              + "".join(b[:5].rjust(7) for b in bands), flush=True)
        for cfg in gate.config.unique():
            row = gF[gF.config == cfg]
            if row.empty:
                continue
            cells = "".join(f"{int((row[row.band == b].gate_p < ALPHA).sum()):7d}" for b in bands)
            print(f"  {cfg.ljust(15)}{dens.get(cfg, np.nan):.3f} {ncomp.get(cfg, np.nan):5.1f}{cells}",
                  flush=True)
        print("", flush=True)

    # ------------------------------------------------------------------ #
    # 2. STRICT plateau rule (R2.2)
    # ------------------------------------------------------------------ #
    print("=== R2.2 STRICT plateau: maximal contiguous f with an IDENTICAL 16-scale "
          "verdict vector (mst-union family) ===", flush=True)
    prows = []
    for func in FUNCTIONALS:
        for band in bands:
            runs = find_plateaus(gate, func, band)
            if not runs:
                continue
            best = max(runs, key=lambda r: (r["octaves"], r["n_points"]))
            qualifies = best["octaves"] >= MIN_OCTAVES
            prows.append(dict(functional=func, band=band, f_lo=best["f_lo"],
                              f_hi=best["f_hi"], n_points=best["n_points"],
                              octaves=best["octaves"],
                              n_scales_cleared=int(sum(best["verdict"])),
                              verdict=vstr(best["verdict"]),
                              qualifies_octave=bool(qualifies),
                              n_runs=len(runs)))
            print(f"  {func:15s} {band:11s} widest run f=[{best['f_lo']:g},{best['f_hi']:g}] "
                  f"({best['octaves']:.2f} oct, {best['n_points']} pts) "
                  f"cleared={sum(best['verdict']):2d}/16  {vstr(best['verdict'])}  "
                  f"[{'PLATEAU' if qualifies else 'no octave-wide plateau'}] "
                  f"({len(runs)} distinct verdicts across the grid)", flush=True)
    pd.DataFrame(prows).to_csv(OUT / "plateaus.csv", index=False)

    # ------------------------------------------------------------------ #
    # 3. GRADED stability profile (post-hoc description, NOT the rule)
    # ------------------------------------------------------------------ #
    print("\n=== graded stability (POST-HOC description, not the frozen rule): "
          "per-scale agreement with mst@0.2 ===", flush=True)
    srows = []
    for func in FUNCTIONALS:
        for band in bands:
            ref = verdict_vector(gate, "mst@0.2", func, band)
            if ref is None:
                continue
            for cfg in gate.config.unique():
                v = verdict_vector(gate, cfg, func, band)
                if v is None:
                    continue
                agree = float(np.mean(np.array(v) == np.array(ref)))
                srows.append(dict(functional=func, band=band, config=cfg,
                                  method=gate[gate.config == cfg].method.iloc[0],
                                  param=gate[gate.config == cfg].param.iloc[0],
                                  density=dens.get(cfg, np.nan),
                                  n_cleared=int(sum(v)), agree_with_ref=agree,
                                  hamming=int(np.sum(np.array(v) != np.array(ref))),
                                  verdict=vstr(v)))
    sp = pd.DataFrame(srows)
    sp.to_csv(OUT / "stability_profile.csv", index=False)
    for func in FUNCTIONALS:
        x = sp[(sp.functional == func) & (sp.method == "mst")]
        if x.empty:
            continue
        print(f"  [{func}] mean per-scale agreement with mst@0.2 across the mst f-grid, "
              f"by band:", flush=True)
        for band in bands:
            y = x[x.band == band].sort_values("param")
            prof = "  ".join(f"{p:g}:{v:.2f}" for p, v in zip(y.param, y.agree_with_ref))
            print(f"    {band:11s} {prof}", flush=True)

    # ------------------------------------------------------------------ #
    # 3b. CONTINUOUS margin surface (post-hoc description, NOT the rule)
    # ------------------------------------------------------------------ #
    # A binary clear/not verdict at n=10 is a coarse readout: the Wilcoxon
    # p-grid on 10 pairs is discrete and a cell sitting near alpha flips on
    # nothing. The cohort MARGIN (median over patients of obs - surr_p50) is the
    # underlying effect size and is far better behaved. Two things are reported:
    # its level, and the SHAPE of its 16-scale profile -- because the project's
    # scientific claims are about tau-dependence (beta scale-invariant vs alpha
    # scale-tuned), so a substrate that preserves the profile shape preserves the
    # claim even where a marginal cell flips its binary verdict.
    print("\n=== continuous margin surface (POST-HOC description, not the frozen rule) ===",
          flush=True)
    print("    per (functional, band): median margin over scales, and Spearman of the "
          "16-scale margin PROFILE against mst@0.2", flush=True)
    mrows2 = []
    for func in FUNCTIONALS:
        for band in bands:
            ref = gate[(gate.config == "mst@0.2") & (gate.functional == func)
                       & (gate.band == band)].sort_values("s")
            if ref.empty:
                continue
            rv = ref.margin_med.values
            for cfg in gate.config.unique():
                x = gate[(gate.config == cfg) & (gate.functional == func)
                         & (gate.band == band)].sort_values("s")
                if len(x) != len(rv):
                    continue
                from scipy.stats import spearmanr
                sp = float(spearmanr(x.margin_med.values, rv).statistic)
                mrows2.append(dict(functional=func, band=band, config=cfg,
                                   method=x.method.iloc[0], param=x.param.iloc[0],
                                   density=dens.get(cfg, np.nan),
                                   margin_med_over_scales=float(np.median(x.margin_med)),
                                   margin_max=float(np.max(x.margin_med)),
                                   profile_spearman_vs_ref=sp))
    ms = pd.DataFrame(mrows2)
    ms.to_csv(OUT / "margin_surface.csv", index=False)
    for func in FUNCTIONALS:
        x = ms[(ms.functional == func) & (ms.method == "mst")]
        if x.empty:
            continue
        print(f"  [{func}] median margin across the mst f-grid, by band "
              f"(profile Spearman vs mst@0.2 in brackets):", flush=True)
        for band in bands:
            y = x[x.band == band].sort_values("param")
            prof = " ".join(f"{p:g}:{v:+.3f}[{c:+.2f}]" for p, v, c
                            in zip(y.param, y.margin_med_over_scales,
                                   y.profile_spearman_vs_ref))
            print(f"    {band:11s} {prof}", flush=True)

    # ------------------------------------------------------------------ #
    # 4. cross-MECHANISM reproduction at matched density (R2.2 second clause)
    # ------------------------------------------------------------------ #
    print("\n=== R2.2 cross-mechanism check: mst-union vs plain threshold at MATCHED "
          "density (and disparity nearby) ===", flush=True)
    mrows = []
    for func in FUNCTIONALS:
        for band in bands:
            for f in (0.05, 0.1, 0.2, 0.4):
                vm = verdict_vector(gate, f"mst@{f:g}", func, band)
                vt = verdict_vector(gate, f"thresh@{f:g}", func, band)
                if vm is None or vt is None:
                    continue
                agree = float(np.mean(np.array(vm) == np.array(vt)))
                mrows.append(dict(functional=func, band=band, frac=f,
                                  mst_cleared=int(sum(vm)), thresh_cleared=int(sum(vt)),
                                  agree=agree,
                                  mst_ncomp=float(ncomp.get(f"mst@{f:g}", np.nan)),
                                  thresh_ncomp=float(ncomp.get(f"thresh@{f:g}", np.nan))))
    mm = pd.DataFrame(mrows)
    mm.to_csv(OUT / "cross_mechanism.csv", index=False)
    for func in FUNCTIONALS:
        x = mm[mm.functional == func]
        if x.empty:
            continue
        print(f"  [{func}] mst vs thresh per-scale agreement, by band x frac:", flush=True)
        for band in bands:
            y = x[x.band == band]
            prof = "  ".join(f"f={r.frac:g}:{r.agree:.2f}({r.mst_cleared}v{r.thresh_cleared})"
                             for r in y.itertuples())
            print(f"    {band:11s} {prof}", flush=True)

    # ------------------------------------------------------------------ #
    # 5. R2.1 parameter-free candidates
    # ------------------------------------------------------------------ #
    print("\n=== R2.1 parameter-free filters: density, verdict, and whether they sit "
          "inside the plateau ===", flush=True)
    pf = pd.DataFrame(prows).set_index(["functional", "band"]) if prows else None
    frows = []
    for filt in PARAM_FREE:
        d = dens.get(filt, np.nan)
        print(f"  --- {filt} (median density {d:.4f}, "
              f"median n_components {ncomp.get(filt, np.nan):.1f}) ---", flush=True)
        for func in FUNCTIONALS:
            for band in bands:
                v = verdict_vector(gate, filt, func, band)
                if v is None:
                    continue
                inside = np.nan
                pl_lo = pl_hi = np.nan
                if pf is not None and (func, band) in pf.index:
                    r = pf.loc[(func, band)]
                    pl_lo, pl_hi = float(r.f_lo), float(r.f_hi)
                    d_lo = dens.get(f"mst@{pl_lo:g}", np.nan)
                    d_hi = dens.get(f"mst@{pl_hi:g}", np.nan)
                    inside = bool(vstr(v) == r.verdict and bool(r.qualifies_octave)
                                  and d_lo <= d <= d_hi)
                frows.append(dict(functional=func, band=band, filter=filt, density=d,
                                  n_cleared=int(sum(v)), verdict=vstr(v),
                                  plateau_f_lo=pl_lo, plateau_f_hi=pl_hi,
                                  inside_plateau=inside))
        x = pd.DataFrame(frows)
        x = x[x["filter"] == filt]
        for func in FUNCTIONALS:
            y = x[x.functional == func]
            if y.empty:
                continue
            print(f"    [{func}] " + "  ".join(
                f"{r.band[:5]}:{r.n_cleared}/16{'*' if r.inside_plateau is True else ''}"
                for r in y.itertuples()), flush=True)
    pd.DataFrame(frows).to_csv(OUT / "parameter_free.csv", index=False)

    # ------------------------------------------------------------------ #
    # 6. R2.3 knob-integrated readout over the mst family
    # ------------------------------------------------------------------ #
    print("\n=== R2.3 knob-integrated readout: statistic MEDIANED over the plateau "
          "(no single f reported) ===", flush=True)
    krows = []
    for func in FUNCTIONALS:
        for band in bands:
            if pf is None or (func, band) not in pf.index:
                continue
            r = pf.loc[(func, band)]
            lo, hi = float(r.f_lo), float(r.f_hi)
            sel = gate[(gate.method == "mst") & (gate.functional == func)
                       & (gate.band == band) & (gate.param >= lo) & (gate.param <= hi)]
            if sel.empty:
                continue
            # per scale: median over the plateau of the cohort margin, and the
            # fraction of plateau fractions at which the scale clears
            byscale = sel.groupby("s").agg(margin=("margin_med", "median"),
                                           frac_clear=("gate_p", lambda z: float(np.mean(z < ALPHA))),
                                           gate_p_med=("gate_p", "median"))
            krows.append(dict(functional=func, band=band, f_lo=lo, f_hi=hi,
                              octaves=float(r.octaves),
                              n_scales_always_clear=int((byscale.frac_clear == 1.0).sum()),
                              n_scales_never_clear=int((byscale.frac_clear == 0.0).sum()),
                              n_scales_ambiguous=int(((byscale.frac_clear > 0)
                                                      & (byscale.frac_clear < 1)).sum()),
                              margin_med=float(byscale.margin.median()),
                              gate_p_med=float(byscale.gate_p_med.median())))
            print(f"  {func:15s} {band:11s} f in [{lo:g},{hi:g}] ({r.octaves:.2f} oct): "
                  f"always-clear {int((byscale.frac_clear == 1.0).sum()):2d}/16, "
                  f"never {int((byscale.frac_clear == 0.0).sum()):2d}/16, "
                  f"ambiguous {int(((byscale.frac_clear > 0) & (byscale.frac_clear < 1)).sum()):2d}/16, "
                  f"median margin {byscale.margin.median():+.4f}", flush=True)
    pd.DataFrame(krows).to_csv(OUT / "knob_integrated.csv", index=False)
    print(f"\n[a2-analysis] wrote plateaus.csv, stability_profile.csv, "
          f"cross_mechanism.csv, parameter_free.csv, knob_integrated.csv -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
