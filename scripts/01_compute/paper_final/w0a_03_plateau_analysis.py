#!/usr/bin/env python3
"""W0-A / A2-A3 — read the stability surface and apply the frozen plateau rule.

Pure analysis of the CSVs written by w0a_02_sparsification_stability.py. No new
compute on the graphs; the only new computation is the sign-flip cluster null on
the already-computed per-patient margin profiles.

FROZEN RULE R2 (verbatim from w0a_00_preregistration.md, written before any
number was seen):
  R2.2 PLATEAU. For band b and functional T, the maximal contiguous set of f over
  which the per-scale cohort gate verdict is CONSTANT, requiring f_max/f_min >= 2
  AND requiring the verdict to be reproduced by at least two different density
  MECHANISMS at matched density.
  R2.1 A parameter-free filter is adopted iff its (density, verdict) point lies
  INSIDE that plateau.
  R2.3 Otherwise the readout is knob-integrated: the median over the plateau.
  R2.4 If no octave-wide invariant region exists, the band is KNOB-DEPENDENT.

THREE-WAY PLATEAU LABEL (added 2026-08-27). A plateau is only informative once
its CONTENT is named, because an all-zeros run is invariance too:
  invariantly-positive  the verdict is constant AND positive across the region
  invariantly-null      constant AND null -- a real result for a null band, and
                        a FAILURE of that substrate for a signal band
  knob-dependent        no octave-wide constant region
Only an invariantly-positive plateau on a signal band, with the null bands
invariantly-null over the same region, licenses a substrate choice.

MULTIPLICITY (added 2026-08-27, after sibling lane W0-C). A 16-point scale sweep
is not 16 independent tests -- neighbouring scales are near-duplicates -- so a
per-scale cleared-COUNT inherits whatever multiplicity family it was computed
under, and a plateau that exists under one family and vanishes under another is
not a plateau. Three families are therefore reported side by side:
  (P) PRIMARY, family-free: ONE sign-flip cluster-mass test per (config, band,
      functional) over the whole 16-scale margin profile (Maris-Oostenveld;
      whole-patient-profile flips preserve the along-axis correlation). This
      collapses the axis to a single test and does not inherit the per-scale
      multiplicity problem at all. The plateau verdict is defined on THIS.
  (S1) per-scale raw alpha=0.05 -- the family used while the sweep ran, kept so
      earlier rows stay comparable.
  (S2) per-scale BH across the whole (band x scale) grid within a config.
The effective number of independent tests on the scale axis is reported so the
reader can see how conservative S2 is.

CALIBRATION WARNING, carried from W0-C: a conditional (partial-correlation)
statistic can return systematically positive values on input containing no
signal -- a sham arc built inside pre-task rest produced p = 0.007. T_probespec_pe
is exactly such a statistic. Its cells are printed but are marked UNCALIBRATED
and must not be gated until a data-based placebo has been run.

Every gate here is on the MARGIN (obs - surr_p50), never a raw observed value:
per-patient observed values co-vary with the height of that patient's own null.

Outputs (data/paper_final/w0a_substrate/a2_stability/{transform}/):
  cluster_gate.csv       PRIMARY: one cluster p per config x band x functional
  plateaus.csv           three-way plateau label per band x functional
  multiplicity_check.csv the same verdict under P / S1 / S2
  stability_profile.csv  graded per-scale agreement with the reference config
  margin_surface.csv     continuous margin level + 16-scale profile shape
  cross_mechanism.csv    mst-union vs plain threshold at matched density
  parameter_free.csv     TMFG / percolation vs the plateau
  knob_integrated.csv    plateau-median readout
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, cluster_stats

FUNCTIONALS = ("T_probe", "T_encode", "T_probespec", "T_probespec_pe")
UNCALIBRATED = ("T_probespec_pe",)          # conditional statistic; see docstring
PARAM_FREE = ("tmfg", "perc")
ALPHA = 0.05
MIN_OCTAVES = 1.0
N_PERM = 10_000
Z_THRESH = 1.0
REF_CONFIG = "mst@0.2"


# --------------------------------------------------------------------------- #
# PRIMARY gate: one sign-flip cluster test per (config, band, functional)
# --------------------------------------------------------------------------- #
def margin_matrix(per, config, functional, band):
    """(K patients, n_scales) matrix of per-patient margins obs - surr_p50."""
    x = per[(per.config == config) & (per.functional == functional)
            & (per.band == band)]
    if x.empty:
        return None
    piv = x.pivot_table(index="patient", columns="s",
                        values="obs_rho", aggfunc="first")
    piv50 = x.pivot_table(index="patient", columns="s",
                          values="surr_p50", aggfunc="first")
    M = (piv - piv50).values
    return M if M.size and np.isfinite(M).any() else None


def _cluster_mass(Y, z_thresh):
    mu = Y.mean(axis=0)
    sd = Y.std(axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd > 1e-12, mu / (sd / np.sqrt(Y.shape[0])), 0.0)
    cl = cluster_stats(z, z_thresh)
    return max((c[2] for c in cl), default=0.0), z


def axis_cluster_p(M, rng, n_perm=N_PERM, z_thresh=Z_THRESH):
    """Sign-flip cluster-mass p for a whole swept axis (Maris-Oostenveld).

    Whole-patient-profile sign flips, so the along-axis correlation is preserved
    and the 16 scales cost ONE test rather than 16.
    """
    ok = np.isfinite(M).all(axis=0)
    X = M[:, ok]
    if X.shape[0] < 3 or X.shape[1] < 1:
        return np.nan, np.nan
    obs, _ = _cluster_mass(X, z_thresh)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = _cluster_mass(X * rng.choice((-1.0, 1.0), size=(X.shape[0], 1)),
                                z_thresh)[0]
    return float((1 + int(np.sum(null >= obs))) / (n_perm + 1)), float(obs)


def effective_tests(M):
    """How many independent tests the correlated scale axis is really worth."""
    ok = np.isfinite(M).all(axis=0)
    X = M[:, ok]
    if X.shape[1] < 2 or X.shape[0] < 3:
        return np.nan, np.nan
    C = np.corrcoef(X, rowvar=False)
    C = np.where(np.isfinite(C), C, 0.0)
    lam = np.clip(np.linalg.eigvalsh(C), 0.0, None)
    n_eff = float(lam.sum() ** 2 / np.sum(lam ** 2)) if lam.sum() > 0 else np.nan
    iu = np.triu_indices(C.shape[0], 1)
    return n_eff, float(np.mean(C[iu]))


# --------------------------------------------------------------------------- #
# plateau over the f grid, with the three-way label
# --------------------------------------------------------------------------- #
def label_run(verdicts):
    """invariantly-positive / invariantly-null for a constant verdict run."""
    return "invariantly-positive" if verdicts[0] else "invariantly-null"


def find_plateaus(bits_by_frac):
    """Maximal contiguous runs of identical verdict along the sorted f grid."""
    fracs = sorted(bits_by_frac)
    runs, i = [], 0
    while i < len(fracs):
        j = i
        while j + 1 < len(fracs) and bits_by_frac[fracs[j + 1]] == bits_by_frac[fracs[i]]:
            j += 1
        runs.append(dict(f_lo=fracs[i], f_hi=fracs[j], n_points=j - i + 1,
                         verdict=bits_by_frac[fracs[i]],
                         octaves=float(np.log2(fracs[j] / fracs[i])) if fracs[i] > 0 else 0.0))
        i = j + 1
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transform", default=os.environ.get("W0A_TRANSFORM", "abs"))
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--null-bands", default="theta",
                    help="comma list of bands required to stay NULL for a fraction to be "
                         "admissible. The project's band-selectivity claim is "
                         "'theta,low_gamma'; 'theta' alone is the weaker criterion.")
    a = ap.parse_args()
    OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a2_stability" / a.transform
    per = pd.read_csv(OUT / "per_patient_scale.csv")
    st = pd.read_csv(OUT / "structure.csv")
    gate = pd.read_csv(OUT / "cohort_gate.csv")
    bands = [b for b in ("delta", "theta", "alpha", "beta", "low_gamma", "high_gamma")
             if b in set(per.band)]
    configs = list(dict.fromkeys(per.config))
    dens = st.groupby("config").density.median()
    ncomp = st.groupby("config").n_components.median()
    n_scales = per.s.nunique()
    rng = np.random.default_rng(20260827)

    print(f"=== W0-A / A2 stability surface — transform = imcoh_{a.transform} ===", flush=True)
    print(f"configs={len(configs)} bands={bands} functionals={per.functional.nunique()} "
          f"scales={n_scales} patients={per.patient.nunique()}", flush=True)

    # ---------------------------------------------------------------- #
    # PRIMARY: one cluster test per (config, band, functional)
    # ---------------------------------------------------------------- #
    print(f"\n=== PRIMARY gate: sign-flip cluster mass over the whole {n_scales}-scale "
          f"axis ({a.n_perm} perms) — ONE test per config x band x functional ===",
          flush=True)
    crows = []
    for cfg in configs:
        for func in FUNCTIONALS:
            for band in bands:
                M = margin_matrix(per, cfg, func, band)
                if M is None:
                    continue
                p, mass = axis_cluster_p(M, rng, n_perm=a.n_perm)
                n_eff, moff = effective_tests(M)
                crows.append(dict(config=cfg, functional=func, band=band,
                                  cluster_p=p, cluster_mass=mass,
                                  n_eff_scales=n_eff, mean_scale_corr=moff,
                                  density=float(dens.get(cfg, np.nan)),
                                  margin_med=float(np.nanmedian(M)),
                                  method=per[per.config == cfg].method.iloc[0],
                                  param=per[per.config == cfg].param.iloc[0]))
    cg = pd.DataFrame(crows)
    cg.to_csv(OUT / "cluster_gate.csv", index=False)
    ne = cg.n_eff_scales.median()
    print(f"  effective independent tests on the {n_scales}-scale axis: median "
          f"n_eff = {ne:.2f} (mean cross-scale corr {cg.mean_scale_corr.median():+.2f}) "
          f"-> per-scale BH over-charges by ~{n_scales/ne:.0f}x", flush=True)
    for func in FUNCTIONALS:
        tag = "  [UNCALIBRATED — conditional statistic, do not gate]" if func in UNCALIBRATED else ""
        print(f"\n  --- [{func}] cluster p per config x band{tag} ---", flush=True)
        print("    " + "config".ljust(15) + "dens   "
              + "".join(b[:5].rjust(9) for b in bands), flush=True)
        for cfg in configs:
            row = cg[(cg.config == cfg) & (cg.functional == func)]
            if row.empty:
                continue
            cells = ""
            for b in bands:
                y = row[row.band == b]
                cells += (f"{y.cluster_p.iloc[0]:9.4f}" if len(y) else " " * 9)
            print(f"    {cfg.ljust(15)}{dens.get(cfg, np.nan):.3f}{cells}", flush=True)

    # ---------------------------------------------------------------- #
    # secondary families S1 / S2, for the robustness check
    # ---------------------------------------------------------------- #
    mrows = []
    for cfg in configs:
        for func in FUNCTIONALS:
            g = gate[(gate.config == cfg) & (gate.functional == func)]
            if g.empty:
                continue
            q = dict(zip(zip(g.band, g.s), bh_fdr(list(g.gate_p.fillna(1.0)))))
            for band in bands:
                gb = g[g.band == band]
                if gb.empty:
                    continue
                cp = cg[(cg.config == cfg) & (cg.functional == func) & (cg.band == band)]
                mrows.append(dict(
                    config=cfg, functional=func, band=band,
                    P_cluster_p=float(cp.cluster_p.iloc[0]) if len(cp) else np.nan,
                    P_verdict=int(float(cp.cluster_p.iloc[0]) < ALPHA) if len(cp) else -1,
                    S1_n_cleared=int((gb.gate_p < ALPHA).sum()),
                    S1_verdict=int((gb.gate_p < ALPHA).any()),
                    S2_n_cleared=int(sum(q[(band, s)] < ALPHA for s in gb.s)),
                    S2_verdict=int(any(q[(band, s)] < ALPHA for s in gb.s))))
    mc = pd.DataFrame(mrows)
    mc.to_csv(OUT / "multiplicity_check.csv", index=False)
    agree_P_S1 = float((mc.P_verdict == mc.S1_verdict).mean())
    agree_P_S2 = float((mc.P_verdict == mc.S2_verdict).mean())
    print(f"\n=== multiplicity-family robustness ===", flush=True)
    print(f"  verdict agreement: PRIMARY(cluster) vs S1(raw per-scale) = {agree_P_S1:.3f}; "
          f"vs S2(whole-grid BH) = {agree_P_S2:.3f}  over {len(mc)} config x band x functional cells",
          flush=True)

    # ---------------------------------------------------------------- #
    # the plateau, on the PRIMARY verdict, with the three-way label
    # ---------------------------------------------------------------- #
    print(f"\n=== R2.2 plateau on the PRIMARY (cluster) verdict, mst-union family, "
          f"three-way label ===", flush=True)
    prows = []
    for func in FUNCTIONALS:
        for band in bands:
            fam = cg[(cg.method == "mst") & (cg.functional == func) & (cg.band == band)]
            if fam.empty:
                continue
            bits = {float(r.param): int(r.cluster_p < ALPHA) for r in fam.itertuples()}
            runs = find_plateaus(bits)
            pos = [r for r in runs if r["verdict"] == 1]
            best = max(runs, key=lambda r: (r["octaves"], r["n_points"]))
            best_pos = max(pos, key=lambda r: (r["octaves"], r["n_points"])) if pos else None
            if best["octaves"] >= MIN_OCTAVES:
                lab = label_run([best["verdict"]])
            else:
                lab = "knob-dependent"
            prows.append(dict(functional=func, band=band, label=lab,
                              f_lo=best["f_lo"], f_hi=best["f_hi"],
                              octaves=best["octaves"], verdict=best["verdict"],
                              n_distinct_runs=len(runs),
                              n_f_positive=int(sum(bits.values())), n_f=len(bits),
                              pos_f_lo=best_pos["f_lo"] if best_pos else np.nan,
                              pos_f_hi=best_pos["f_hi"] if best_pos else np.nan,
                              pos_octaves=best_pos["octaves"] if best_pos else np.nan))
            bitstr = "".join(str(bits[f]) for f in sorted(bits))
            tag = "  [UNCALIBRATED]" if func in UNCALIBRATED else ""
            print(f"  {func:15s} {band:11s} {lab:22s} widest run f=[{best['f_lo']:g},"
                  f"{best['f_hi']:g}] ({best['octaves']:.2f} oct)  positive at "
                  f"{sum(bits.values())}/{len(bits)} fractions  bits(f asc)={bitstr}"
                  f"{tag}", flush=True)
    pd.DataFrame(prows).to_csv(OUT / "plateaus.csv", index=False)

    # ---------------------------------------------------------------- #
    # continuous margin surface (level + 16-scale profile shape)
    # ---------------------------------------------------------------- #
    print(f"\n=== continuous margin surface — level, and Spearman of the {n_scales}-scale "
          f"profile vs {REF_CONFIG} ===", flush=True)
    srows = []
    for func in FUNCTIONALS:
        for band in bands:
            ref = gate[(gate.config == REF_CONFIG) & (gate.functional == func)
                       & (gate.band == band)].sort_values("s")
            if ref.empty:
                continue
            rv = ref.margin_med.values
            for cfg in configs:
                x = gate[(gate.config == cfg) & (gate.functional == func)
                         & (gate.band == band)].sort_values("s")
                if len(x) != len(rv):
                    continue
                srows.append(dict(functional=func, band=band, config=cfg,
                                  method=x.method.iloc[0], param=x.param.iloc[0],
                                  density=float(dens.get(cfg, np.nan)),
                                  margin_med=float(np.median(x.margin_med)),
                                  margin_max=float(np.max(x.margin_med)),
                                  profile_spearman=float(spearmanr(x.margin_med.values, rv).statistic)))
    ms = pd.DataFrame(srows)
    ms.to_csv(OUT / "margin_surface.csv", index=False)
    for func in FUNCTIONALS:
        x = ms[(ms.functional == func) & (ms.method == "mst")]
        if x.empty:
            continue
        print(f"  [{func}] median margin across the mst f-grid "
              f"(profile Spearman vs {REF_CONFIG} in brackets):", flush=True)
        for band in bands:
            y = x[x.band == band].sort_values("param")
            print(f"    {band:11s} " + " ".join(
                f"{p:g}:{v:+.3f}[{c:+.2f}]" for p, v, c
                in zip(y.param, y.margin_med, y.profile_spearman)), flush=True)

    # ---------------------------------------------------------------- #
    # cross-mechanism at matched density (R2.2 second clause)
    # ---------------------------------------------------------------- #
    print(f"\n=== R2.2 cross-mechanism: mst-union vs plain threshold at MATCHED density "
          f"(PRIMARY verdict; ncomp = components) ===", flush=True)
    xrows = []
    for func in FUNCTIONALS:
        for band in bands:
            for f in sorted(set(cg[cg.method == "thresh"].param.dropna())):
                a_ = cg[(cg.config == f"mst@{f:g}") & (cg.functional == func) & (cg.band == band)]
                b_ = cg[(cg.config == f"thresh@{f:g}") & (cg.functional == func) & (cg.band == band)]
                if a_.empty or b_.empty:
                    continue
                xrows.append(dict(functional=func, band=band, frac=f,
                                  mst_p=float(a_.cluster_p.iloc[0]),
                                  thresh_p=float(b_.cluster_p.iloc[0]),
                                  agree=int((a_.cluster_p.iloc[0] < ALPHA)
                                            == (b_.cluster_p.iloc[0] < ALPHA)),
                                  thresh_ncomp=float(ncomp.get(f"thresh@{f:g}", np.nan))))
    xm = pd.DataFrame(xrows)
    xm.to_csv(OUT / "cross_mechanism.csv", index=False)
    for func in FUNCTIONALS:
        x = xm[xm.functional == func]
        if x.empty:
            continue
        print(f"  [{func}] " + " | ".join(
            f"{b[:5]} " + " ".join(f"f={r.frac:g}:{'=' if r.agree else 'X'}"
                                   f"({r.mst_p:.3f}/{r.thresh_p:.3f})"
                                   for r in x[x.band == b].itertuples())
            for b in bands if not x[x.band == b].empty), flush=True)

    # ---------------------------------------------------------------- #
    # R2.1 parameter-free filters
    # ---------------------------------------------------------------- #
    print(f"\n=== R2.1 parameter-free filters vs the plateau (PRIMARY verdict) ===", flush=True)
    pf = pd.DataFrame(prows).set_index(["functional", "band"]) if prows else None
    frows = []
    for filt in PARAM_FREE:
        d = float(dens.get(filt, np.nan))
        print(f"  --- {filt} (median density {d:.4f}, components {ncomp.get(filt, np.nan):.1f}) ---",
              flush=True)
        for func in FUNCTIONALS:
            cells = []
            for band in bands:
                row = cg[(cg.config == filt) & (cg.functional == func) & (cg.band == band)]
                if row.empty:
                    continue
                p = float(row.cluster_p.iloc[0])
                v = int(p < ALPHA)
                inside = np.nan
                if pf is not None and (func, band) in pf.index:
                    r = pf.loc[(func, band)]
                    d_lo = float(dens.get(f"mst@{r.f_lo:g}", np.nan))
                    d_hi = float(dens.get(f"mst@{r.f_hi:g}", np.nan))
                    inside = bool(v == int(r.verdict) and r.label != "knob-dependent"
                                  and d_lo <= d <= d_hi)
                frows.append(dict(functional=func, band=band, filter=filt, density=d,
                                  cluster_p=p, verdict=v, inside_plateau=inside))
                cells.append(f"{band[:5]}:{p:.3f}{'*' if v else ''}")
            print(f"    [{func:15s}] " + "  ".join(cells), flush=True)
    pd.DataFrame(frows).to_csv(OUT / "parameter_free.csv", index=False)

    # ---------------------------------------------------------------- #
    # R2.3 knob-integrated readout
    # ---------------------------------------------------------------- #
    print(f"\n=== R2.3 knob-integrated readout over the mst f-grid (no single f) ===",
          flush=True)
    krows = []
    for func in FUNCTIONALS:
        for band in bands:
            fam = cg[(cg.method == "mst") & (cg.functional == func) & (cg.band == band)]
            if fam.empty:
                continue
            frac_pos = float((fam.cluster_p < ALPHA).mean())
            krows.append(dict(functional=func, band=band,
                              n_f=len(fam), frac_f_positive=frac_pos,
                              cluster_p_med=float(fam.cluster_p.median()),
                              margin_med=float(fam.margin_med.median())))
            print(f"  {func:15s} {band:11s} positive at {frac_pos*100:5.1f}% of fractions; "
                  f"median cluster p = {fam.cluster_p.median():.4f}; "
                  f"median margin = {fam.margin_med.median():+.4f}", flush=True)
    pd.DataFrame(krows).to_csv(OUT / "knob_integrated.csv", index=False)

    # ---------------------------------------------------------------- #
    # ADMISSIBLE WINDOW — the actual substrate decision
    # ---------------------------------------------------------------- #
    # A substrate is admissible only where the SIGNAL bands are positive AND the
    # NULL bands are null, on the same graph. Per R2's three-way label, an
    # all-zeros plateau on a signal band is a failure, and a positive verdict on
    # a null band is a leak. Reported on T_probe, the standard trace, because
    # that is the functional the contract is chosen on; the other functionals are
    # then read ON the chosen window rather than used to choose it.
    SIGNAL = ("alpha", "beta")
    NULLB = tuple(b.strip() for b in a.null_bands.split(",") if b.strip())
    print(f"\n=== ADMISSIBLE WINDOW on T_probe: signal {SIGNAL} positive AND "
          f"null {NULLB} null, at the same f ===", flush=True)
    fam = cg[(cg.method == "mst") & (cg.functional == "T_probe")]
    fracs = sorted(set(fam.param.dropna()))
    arows = []
    for f in fracs:
        row = {"frac": f, "density": float(dens.get(f"mst@{f:g}", np.nan))}
        ok_sig, ok_null = True, True
        for b in bands:
            y = fam[(fam.band == b) & np.isclose(fam.param, f)]
            p = float(y.cluster_p.iloc[0]) if len(y) else np.nan
            row[b] = p
            if b in SIGNAL and not (p < ALPHA):
                ok_sig = False
            if b in NULLB and (p < ALPHA):
                ok_null = False
        row["admissible"] = bool(ok_sig and ok_null)
        arows.append(row)
    ad = pd.DataFrame(arows)
    ad.to_csv(OUT / f"admissible_window_null-{'-'.join(NULLB)}.csv", index=False)
    for r in ad.itertuples():
        marks = " ".join(f"{b[:5]}={getattr(r, b):.3f}{'*' if getattr(r, b) < ALPHA else ' '}"
                         for b in bands)
        print(f"  f={r.frac:<5g} d={r.density:.3f}  {marks}   "
              f"{'ADMISSIBLE' if r.admissible else ''}", flush=True)
    good = ad[ad.admissible]
    if len(good):
        lo, hi = float(good.frac.min()), float(good.frac.max())
        contiguous = list(good.frac) == [f for f in fracs if lo <= f <= hi]
        centre = float(np.sqrt(lo * hi))
        nearest = min(fracs, key=lambda f: abs(np.log(f) - np.log(centre)))
        print(f"\n  --> admissible f in [{lo:g}, {hi:g}] "
              f"({np.log2(hi/lo):.2f} octaves, {len(good)} of {len(fracs)} fractions, "
              f"{'contiguous' if contiguous else 'NOT contiguous'}); "
              f"density [{good.density.min():.3f}, {good.density.max():.3f}]", flush=True)
        print(f"  --> log-centre f = {centre:.3f}; nearest swept fraction = {nearest:g}",
              flush=True)
        print(f"  --> incumbent f = 0.20 is "
              f"{'INSIDE' if 0.20 in list(good.frac) else 'OUTSIDE'} the window, at its "
              f"{'upper edge' if abs(np.log(0.20) - np.log(hi)) < 1e-9 else 'interior'}",
              flush=True)
        # knob-integrated readout over the admissible window, all functionals
        print(f"\n  knob-integrated readout over the admissible window "
              f"(median over f in [{lo:g},{hi:g}]):", flush=True)
        for func in FUNCTIONALS:
            tag = " [UNCALIBRATED]" if func in UNCALIBRATED else ""
            cells = []
            for b in bands:
                y = cg[(cg.method == "mst") & (cg.functional == func) & (cg.band == b)
                       & (cg.param >= lo) & (cg.param <= hi)]
                if y.empty:
                    continue
                cells.append(f"{b[:5]}:p={y.cluster_p.median():.3f}"
                             f"/m={y.margin_med.median():+.3f}"
                             f"({int((y.cluster_p < ALPHA).sum())}/{len(y)})")
            print(f"    {func:15s} " + "  ".join(cells) + tag, flush=True)
    else:
        print("\n  --> NO admissible fraction: no f makes both signal bands positive "
              "while the null band stays null. The substrate is NOT fixable on this "
              "grid and the verdict is knob-dependent.", flush=True)
    print(f"\n[a2-analysis] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
