#!/usr/bin/env python3
"""W0-A / A2 — the sparsification stability surface (NOT a best f).

PRE-REGISTERED RULE R2 (frozen in w0a_00_preregistration.md before any number;
reproduced so the script is self-contained):

  1. Parameter-free first. A parameter-free filter is adopted iff its
     (density, verdict) point lies INSIDE the invariance plateau.
  2. Plateau, frozen definition. For band b, the maximal contiguous set of f
     over which the PER-SCALE cohort gate verdict (cleared / not at alpha=0.05,
     per scale, no best-scale collapse) is constant, requiring width >= one
     octave (f_max / f_min >= 2) AND requiring the verdict to be reproduced by
     at least two different density MECHANISMS at matched density.
  3. Fallback: knob-integrated readout -- median over the plateau, plateau
     width stated, no single f ever reported as the pipeline setting.
  4. If no plateau of width >= one octave exists for a band, that band's
     verdict is reported as KNOB-DEPENDENT. It is not rescued by picking an f.

5-POINT CRITICAL PREAMBLE
1. Claim. There is a region of the density knob over which the per-band,
   per-scale cross-phase verdict does not change, so the substrate can be fixed
   without choosing the answer.
2. Null. H0: the verdict is a function of the knob -- no invariant region
   exists, and any stated verdict is an artifact of the chosen f.
3. Strongest plausible alternative the null should control for. That a
   "plateau" is manufactured by the statistic rather than by the graph:
   adjacent f share almost all their edges, so the surface is smooth by
   construction and a flat region carries no information.
4. Does the null control for it, by mechanism? NOT by itself -- and this is the
   design consequence. The sweep deliberately holds the SURROGATE DRAWS FIXED
   across configs within a cell (one matched-strength shuffle of the dense FC is
   re-sparsified by every config), so every config-to-config difference is
   attributable to the backbone and not to surrogate noise; but that also makes
   neighbouring f MORE correlated, not less. The plateau claim is therefore made
   structural, not statistical: it is only asserted where the verdict survives a
   change of the density-selection MECHANISM (mst-union, which forces
   connectivity, vs the plain global threshold, which does not, vs disparity,
   which is local-statistical) at MATCHED density, and it is reported together
   with the structural covariates (density, mean degree, clustering,
   algebraic connectivity, spectral gap, retained weight fraction) that show the
   graph itself changed materially across the plateau. What it CANNOT reject:
   every mechanism in the panel is a weight-ranked filter of the same dense
   |ImCoh| matrix, so a bias common to all weight-ranked filters is invisible
   here; and matched-strength is injected at the FC-matrix stage, so nothing in
   this script tests the coherency estimator, the band split, or session
   nonstationarity (that is lane W0-B).
5. Falsification / limits. Falsified for a band if no contiguous octave-wide
   invariant region exists, or if mst-union and the plain threshold disagree at
   matched density. Limits: n=10 gives a coarse Wilcoxon p-grid; the s-grid is
   the incumbent 16-point logspace(1, 180) and inherits its range choice; the
   plain-threshold arm fragments below the percolation point, and a fragmented
   graph's cophenetic tree is dominated by component structure -- n_components
   is recorded so those cells can be read correctly rather than silently
   averaged in.

Pipeline is identical to the incumbent harness (13_matched_strength_mst020.py):
matched-strength shuffle of the DENSE FC -> sparsify -> combinatorial Laplacian
-> heat kernel -> UPGMA cophenetic at scale s -> cross-phase functional, swept
over s. ``T_probe`` reproduces the incumbent ``rho_sym`` to 1e-16
(w0a_verify_functionals.py), so this is the incumbent gate plus three axes.

FIVE phases, FOUR functionals. The paradigm is transitive inference: task_learn
presents ordered premises, task_test requires judging novel non-adjacent pairs.
A substrate selected only on the standard trace could be wrong for the encoding
and inference-specific components -- they are different functionals of the SAME
cophenetic distances and there is no guarantee their stability plateaus
coincide. So the sweep carries the functional axis rather than assuming it:

  T_probe        standard trace          rho(D_test - D_A,     D_post - D_B)
  T_encode       encoding echo           rho(D_learn - D_A,    D_post - D_B)
  T_probespec    inference-specific      rho(D_test - D_learn, D_post - D_B)
  T_probespec_pe inference | encoding    partial, conditioning on the encoding
                                         reorganisation

all symmetrised over the arbitrary A/B arm assignment and cross-baseline. The
functional axis is nearly free: all Spearmans come from one rank correlation
matrix per scale, and the marginal cost over the four-phase version is the one
extra cophenetic tree (task_learn).

Env: W0A_TRANSFORM = abs | sq   (which imcoh magnitude transform to run on)
     SA_WORKERS, W0A_R (surrogates), W0A_BANDS (comma list),
     W0A_CONFIGS (comma list of config labels; default = the full grid)

Outputs (data/paper_final/w0a_substrate/a2_stability/{transform}/):
  per_patient_scale.csv  config, method, param, functional, patient, band, s,
                         obs_rho, surr_p50/p95, p
  cohort_gate.csv        config, method, param, functional, band, s, gate_p,
                         n_above, obs_med, surr_med, margin_med
  structure.csv          config, patient, band, phase, density, clustering,
                         lambda_2, spectral_gap, weight_fraction, n_components
  config.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import backbone_structure, select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_over_scales, laplacian_eig,
)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.workflow.fc import load_fc_matrix

COHORT = list(PATIENTS_4PHASE)
ALL_BANDS = list(BRAIN_BANDS)
# FIVE phases: transitive inference separates learning the structure
# (task_learn) from applying it (task_test). The phase set is data, not a
# hardcoded tuple -- see lrg_eegfc.utils.fc.heat_multiscale.CROSS_PHASE_ROLES.
PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
FUNCTIONALS = ("T_probe", "T_encode", "T_probespec", "T_probespec_pe")
SGRID = np.logspace(0.0, np.log10(180.0), 16)          # incumbent grid, unchanged
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260825

TRANSFORM = os.environ.get("W0A_TRANSFORM", "abs")
R = int(os.environ.get("W0A_R", 200))
HALVES = ROOT / "data" / "paper_final" / "w0a_substrate" / "halves"
OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a2_stability" / TRANSFORM

# --- the config grid: density knob x mechanism ------------------------------ #
MST_FRACS = (0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.14, 0.20, 0.28, 0.40, 0.55, 0.75, 1.00)
THRESH_FRACS = (0.05, 0.10, 0.20, 0.40)                 # mechanism contrast (may fragment)
DISP_ALPHAS = (0.01, 0.05, 0.10, 0.20, 0.50)

CONFIGS: list[tuple[str, str, float]] = []
CONFIGS += [(f"mst@{f:g}", "mst", f) for f in MST_FRACS]
CONFIGS += [(f"thresh@{f:g}", "thresh", f) for f in THRESH_FRACS]
CONFIGS += [(f"disparity@{a:g}", "disparity", a) for a in DISP_ALPHAS]
CONFIGS += [("tmfg", "tmfg", np.nan), ("perc", "perc", np.nan)]

_WANT = os.environ.get("W0A_CONFIGS", "")
if _WANT:                                       # reduced grid (transform-robustness arm)
    _keep = {c.strip() for c in _WANT.split(",") if c.strip()}
    unknown = _keep - {c[0] for c in CONFIGS}
    if unknown:
        raise ValueError(f"W0A_CONFIGS names no such config: {sorted(unknown)}")
    CONFIGS = [c for c in CONFIGS if c[0] in _keep]


def make_backbone(W, method, param):
    if method == "disparity":
        return select_backbone(W, "disparity", disparity_alpha=param)
    if method in ("mst", "thresh"):
        return select_backbone(W, method, frac=param)
    return select_backbone(W, method)


def load_phase(pat, phase, band):
    """Per-phase FC under the selected transform. A/B are the rest_pre halves."""
    if phase in ("A", "B"):
        W = np.load(HALVES / pat / f"{band}_rest_pre_{phase}_imcoh_{TRANSFORM}.npy")
    else:
        W = load_fc_matrix(pat, phase, band, fc_method=f"imcoh_{TRANSFORM}")
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def per_cell(job):
    """One (patient, band): observed + R shared surrogates across EVERY config."""
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception as exc:                                       # pragma: no cover
        print(f"  [{pat}/{band}] load failed: {exc}", flush=True)
        return None, None
    N = Ws["A"].shape[0]

    # ---- observed + structural covariates, per config ---------------------- #
    obs = {}
    struct_rows = []
    for label, method, param in CONFIGS:
        Bs = {ph: make_backbone(Ws[ph], method, param) for ph in PHASES}
        obs[label] = cross_phase_functionals_over_scales(
            {ph: laplacian_eig(Bs[ph]) for ph in PHASES}, SGRID)
        for ph in PHASES:
            st = backbone_structure(Bs[ph])
            st.update(config=label, method=method, param=param, patient=pat,
                      band=band, phase=ph, N=int(N),
                      weight_fraction=float(Bs[ph].sum() / Ws[ph].sum())
                      if Ws[ph].sum() > 0 else np.nan)
            struct_rows.append(st)

    # ---- surrogates: ONE shuffle per r, re-sparsified by EVERY config ------- #
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = {label: {k: np.full((R, SGRID.size), np.nan) for k in FUNCTIONALS}
            for label, _, _ in CONFIGS}
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        for label, method, param in CONFIGS:
            eig_s = {ph: laplacian_eig(make_backbone(Wsh[ph], method, param))
                     for ph in PHASES}
            so = cross_phase_functionals_over_scales(eig_s, SGRID)
            for k in FUNCTIONALS:
                surr[label][k][r] = so[k]

    rows = []
    for label, method, param in CONFIGS:
        for k in FUNCTIONALS:
            o_all, s_all = obs[label][k], surr[label][k]
            for j, s in enumerate(SGRID):
                col = s_all[:, j]
                col = col[np.isfinite(col)]
                o = o_all[j]
                p = float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan
                rows.append(dict(config=label, method=method, param=param,
                                 functional=k, patient=pat, band=band,
                                 s=float(s), N=int(N),
                                 obs_rho=float(o) if np.isfinite(o) else np.nan,
                                 surr_p50=float(np.nanpercentile(col, 50)) if col.size else np.nan,
                                 surr_p95=float(np.nanpercentile(col, 95)) if col.size else np.nan,
                                 p=p, n_surr=int(col.size)))
    return rows, struct_rows


def cohort_gate(df, bands):
    """Per (config, functional, band, scale) one-sided Wilcoxon on obs - surr_p50."""
    out = []
    for label in df.config.unique():
        d0 = df[df.config == label]
        method = d0.method.iloc[0]
        param = d0.param.iloc[0]
        for func in FUNCTIONALS:
            dF = d0[d0.functional == func]
            for band in bands:
                d1 = dF[dF.band == band]
                for s in np.sort(d1.s.unique()):
                    x = d1[np.isclose(d1.s, s)].dropna(subset=["obs_rho", "surr_p50"])
                    if len(x) < 5:
                        continue
                    d = x.obs_rho.values - x.surr_p50.values
                    try:
                        p = (float(wilcoxon(d, alternative="greater")[1])
                             if np.any(d != 0) else np.nan)
                    except Exception:
                        p = np.nan
                    out.append(dict(config=label, method=method, param=param,
                                    functional=func, band=band, s=float(s), gate_p=p,
                                    n_above=int((x.p < 0.05).sum()), n_pat=len(x),
                                    obs_med=float(x.obs_rho.median()),
                                    surr_med=float(x.surr_p50.median()),
                                    margin_med=float(np.median(d))))
    return pd.DataFrame(out)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="first K patients (timing runs)")
    ap.add_argument("--time-only", action="store_true",
                    help="time one cell at R=5 and extrapolate; do not run the grid")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    bands = os.environ.get("W0A_BANDS", "").split(",") if os.environ.get("W0A_BANDS") else ALL_BANDS
    pats = COHORT[:a.limit] if a.limit else COHORT
    ncpu = int(os.environ.get("SA_WORKERS", 12))

    # warm numba (matched-strength shuffle + TMFG) outside the timed region
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    _ = select_backbone(np.abs(np.random.default_rng(0).standard_normal((12, 12))), "tmfg")

    if a.time_only:
        # Time at TWO surrogate counts and fit cost = fixed + slope*R, so the
        # one-off observed + structural-covariate overhead is not scaled with R.
        global R
        R_full = R
        ts = []
        for r_probe in (4, 16):
            R = r_probe
            t0 = time.time()
            per_cell((0, pats[0], bands[0]))
            ts.append(time.time() - t0)
        slope = (ts[1] - ts[0]) / (16 - 4)
        fixed = ts[0] - 4 * slope
        per_cell_full = fixed + slope * R_full
        n_cells = len(pats) * len(bands)
        est = per_cell_full * n_cells / ncpu
        print(f"[a2/{TRANSFORM}] {len(CONFIGS)} configs; cell cost = {fixed:.1f}s fixed "
              f"(observed + structure) + {slope:.2f}s/surrogate  "
              f"[R=4: {ts[0]:.1f}s, R=16: {ts[1]:.1f}s]", flush=True)
        print(f"[a2/{TRANSFORM}] FULL GRID ESTIMATE at R={R_full}, {n_cells} cells, "
              f"{ncpu} workers: {per_cell_full:.0f}s/cell -> {est/60:.1f} min "
              f"({est/3600:.2f} h)", flush=True)
        return

    print(f"[a2/{TRANSFORM}] {len(pats)*len(bands)} cells x {len(CONFIGS)} configs x "
          f"{len(SGRID)} scales, R={R}, {ncpu} workers -> {OUT}", flush=True)
    print(f"[a2/{TRANSFORM}] configs: {', '.join(c[0] for c in CONFIGS)}", flush=True)
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    t0 = time.time()
    rows, srows = [], []
    with Pool(ncpu) as pool:
        for i, (rl, sl) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
                srows.extend(sl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_patient_scale.csv", index=False)
    pd.DataFrame(srows).to_csv(OUT / "structure.csv", index=False)
    gate = cohort_gate(df, bands)
    gate.to_csv(OUT / "cohort_gate.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="a2_sparsification_stability", transform=TRANSFORM,
        cohort=pats, bands=bands, R=R, swap_factor=SWAP_FACTOR,
        s_grid=[float(x) for x in SGRID], base_seed=BASE_SEED,
        configs=[dict(label=l, method=m, param=(None if not np.isfinite(p) else float(p)))
                 for l, m, p in CONFIGS],
        null="matched-strength on dense FC -> sparsify(config) -> LRG -> rho_sym, per scale; "
             "surrogate draws SHARED across configs within a cell (paired across the knob)",
        rule="w0a_00_preregistration.md R2 (frozen before any number)",
    ), indent=2))
    print(f"\n[a2/{TRANSFORM}] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}", flush=True)

    # ---- terse per-band scale-count summary (full surface is in the CSV) ---- #
    st = pd.DataFrame(srows)
    for func in FUNCTIONALS:
        print(f"\n=== [{func}] n_scales cleared (gate_p < 0.05) of {len(SGRID)}, "
              f"per config x band ===", flush=True)
        print("  " + "config".ljust(16) + "".join(b[:5].rjust(8) for b in bands)
              + "   density", flush=True)
        gF = gate[gate.functional == func]
        for label, _, _ in CONFIGS:
            g = gF[gF.config == label]
            if g.empty:
                continue
            cells = "".join(f"{int((g[g.band == b].gate_p < 0.05).sum()):8d}" for b in bands)
            dens = st[st.config == label].density.median()
            print(f"  {label.ljust(16)}{cells}   {dens:.4f}", flush=True)


if __name__ == "__main__":
    main()
