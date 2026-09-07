#!/usr/bin/env python3
"""Amendment B: does a *contrast between functionals* carry scale structure?

Part A of this lane (``w0s_01``..``w0s_06``) tested scale-locality of one
functional at a time and returned a negative. A contrast is a different object
-- two quantities can each be flat along tau while their difference is not,
because the flat parts cancel -- so the Part-A verdict does not cover it. This
script tests the contrast, under the criteria pre-registered in scope report
``2026-08-31_scale-local-trace-readouts.md`` section 12, written before any
contrast statistic was computed.

Inputs are Lane E's five-phase cells, read-only. Nothing here recomputes FC,
LRG, surrogates or the sham; the contrast is formed from the stored ``obs`` and
``surr`` blocks, whose realization index is shared across the functional axis,
so the pairing survives the difference.

Objects tested (all reported, including the failures)
-----------------------------------------------------
``T_test``, ``T_learn``, ``T_infspec``, ``T_infspec_pe``
    the raw functionals, so that "the contrast decorrelates more than its
    parts" is a number rather than an assertion.
``C_learn_minus_test`` (paired, primary)
    difference formed at raw level, before any surrogate subtraction, so the
    locked ``patient_margin`` contract stays exact rather than approximated.
``C_unpaired`` (sensitivity)
    ``margin[T_learn] - margin[T_test]``, which subtracts two independently
    drawn surrogate medians and discards the pairing.

Nulls
-----
(a) held-out matched-strength realizations of the *same* object -- the object at
    its own noise level, which is the only reference that can tell extra
    information from extra noise;
(b) Lane E's ordered sham -- a fake five-phase arc carved from one rest
    recording with block order preserved and no task -- which is the only
    reference that can see within-recording drift, since matched-strength
    surrogates are redrawn per phase and carry no recording-order information.

Nothing is tested against zero anywhere.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.cohort_gate import (
    axis_cluster_gate,
    axis_reversal_gate,
    effective_tests,
    gate_grid,
)
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from w0s_00_artifacts import CellSet                          # noqa: E402
from w0s_03_scale_local_verdict import _shape_stats, _upper_p  # noqa: E402

ENC = ROOT / "data" / "paper_final" / "lane_e_encinf"
OUT = Path(os.environ.get(
    "W0S_OUT_CONTRAST",
    ROOT / "data" / "paper_final" / "lane_s_scale" / "contrast"))

BANDS = ("delta", "theta", "alpha", "beta", "low_gamma", "high_gamma")
SHAM_SOURCES = ("rest_pre", "rest_post")
#: (name, minuend, subtrahend) for every derived contrast carried through.
CONTRASTS = (("C_learn_minus_test", "T_learn", "T_test"),)
N_PERM = int(os.environ.get("W0S_NPERM", "10000"))
N_DRAWS = int(os.environ.get("W0S_NDRAWS", "100"))
EXPECT_N = 10


class FunctionalCells(CellSet):
    """Lane E's five-phase cells as a :class:`CellSet`, with derived contrasts.

    Lane E stores ``obs (nF, nS, nFunc)`` and ``surr (nF, R, nS, nFunc)``; this
    lane's gate machinery expects ``obs (nF, nS, nRead)`` and
    ``surr (R, nF, nS, nRead)``. The transpose happens once, here, so no
    downstream statistic has to know which lane wrote the file.

    Each contrast is appended as an extra column of the readout axis, formed by
    differencing ``obs`` and ``surr`` along the functional axis *before* any
    knob median or surrogate subtraction. Because realization ``r`` indexes the
    same matched-strength draw for every functional within a cell, the
    difference of the surrogates is the surrogate of the difference, and the
    locked ``patient_margin`` contract stays exact rather than approximated.
    """

    def __init__(self, root, patients, bands, *, source=None, arc="ordered"):
        self.root, self.patients, self.bands = Path(root), list(patients), list(bands)
        self._obs, self._surr, self._diag, self._cache = {}, {}, {}, {}
        self.s = self.fracs = None
        self.readouts, self.N = [], {}
        okey = "obs" if source is None else f"obs_{arc}"
        skey = "surr" if source is None else f"surr_{arc}"
        for b in self.bands:
            for p in self.patients:
                f = (self.root / "cells" /
                     (f"{p}__{b}.npz" if source is None
                      else f"{p}__{b}__{source}.npz"))
                if not f.exists():
                    continue
                z = np.load(f, allow_pickle=False)
                funcs = [str(x) for x in z["funcs"]]
                O, S = z[okey], z[skey]
                if source is not None:
                    O, S = O[0], S[0]
                S = np.moveaxis(S, 1, 0)
                names = list(funcs)
                for name, a, bb in CONTRASTS:
                    ia, ib = funcs.index(a), funcs.index(bb)
                    O = np.concatenate([O, (O[..., ia] - O[..., ib])[..., None]], -1)
                    S = np.concatenate([S, (S[..., ia] - S[..., ib])[..., None]], -1)
                    names.append(name)
                if self.s is None:
                    self.s, self.fracs, self.readouts = z["s"], z["fracs"], names
                self._obs[(p, b)] = O
                self._surr[(p, b)] = S
                self.N[(p, b)] = int(z["N"][0]) if "N" in z.files else -1
        if self.s is None:
            raise FileNotFoundError(f"no cells under {self.root/'cells'}")
        self.ix = {k: i for i, k in enumerate(self.readouts)}
        self.n_scales = int(self.s.size)
        self.R = min(v.shape[0] for v in self._surr.values())
        for k in list(self._surr):
            if self._surr[k].shape[0] > self.R:
                self._surr[k] = self._surr[k][: self.R]

    def surr_sd(self, band, readout):
        """Width of each patient's own null ensemble, ``(K, nS)``.

        The descriptive noise level of the object, measured rather than
        assumed: a contrast of two functionals that share most of their input
        vectors may be *quieter* than either input, not noisier, and that is an
        empirical question about these cells.
        """
        _, S, _ = self.cell(band, readout)
        return np.nanstd(S, axis=1, ddof=1)

    def surr_corr(self, band, a, b):
        """Across-realization correlation of two objects' nulls, ``(K, nS)``.

        How much of the two functionals' noise is common. Near 1 means the
        difference cancels shared noise; near 0 means it adds two independent
        errors. This is the number that decides how much of the "a difference
        is noisier than its parts" caveat applies here.
        """
        _, Sa, _ = self.cell(band, a)
        _, Sb, _ = self.cell(band, b)
        K, R, nS = Sa.shape
        out = np.full((K, nS), np.nan)
        for k in range(K):
            for j in range(nS):
                u, v = Sa[k, :, j], Sb[k, :, j]
                m = np.isfinite(u) & np.isfinite(v)
                if m.sum() > 3 and u[m].std() > 0 and v[m].std() > 0:
                    out[k, j] = np.corrcoef(u[m], v[m])[0, 1]
        return out


def _unpaired_margins(cs, band, a, b):
    """``margin[a] - margin[b]`` -- the pairing-discarding variant, for sensitivity."""
    Ma, labs = cs.margins(band, a)
    Mb, _ = cs.margins(band, b)
    return Ma - Mb, labs


def _per_patient_slopes(M, s):
    """Each patient's own ``Spearman(margin, log s)``, ``(K,)``."""
    M = np.asarray(M, float)
    ok = np.isfinite(M).all(axis=0)
    if ok.sum() < 4:
        return np.full(M.shape[0], np.nan)
    ls = np.log(np.asarray(s, float)[ok])
    return np.array([spearmanr(M[k, ok], ls)[0] for k in range(M.shape[0])])


def _null_bundle(cs, band, readout, s, n_draws):
    """Held-out-realization null for every scale-structure statistic at once.

    Under the null the observation is exchangeable with its own surrogates, so
    promoting draw ``i`` to the observed slot and re-referencing to the rest is
    a draw from this object's null *at this object's own noise level*. That
    last clause is the whole point: a noisier object decorrelates its scales
    and inflates every scale-structure statistic for free, and only this
    comparison can tell extra information from extra noise.
    """
    _, S, _ = cs.cell(band, readout)
    n = min(n_draws, S.shape[1])
    keys = ("n_eff_pr", "slope_abs", "slope_signed", "shape_agreement",
            "mass_pos", "mass_neg", "n_cross", "xstar_mad")
    out = {k: np.full(n, np.nan) for k in keys}
    for i in range(n):
        Mn = cs.heldout_margins(band, readout, i)
        out["n_eff_pr"][i] = effective_tests(Mn)["n_eff_pr"]
        sh = _shape_stats(Mn, s)
        out["slope_abs"][i] = sh["slope_abs"]
        out["slope_signed"][i] = sh["slope_signed"]
        out["shape_agreement"][i] = sh["shape_agreement"]
        rv = axis_reversal_gate(Mn, s, n_perm=0, rng=np.random.default_rng(i))
        out["mass_pos"][i] = rv["mass_pos"]
        out["mass_neg"][i] = rv["mass_neg"]
        out["n_cross"][i] = rv["n_cross"]
        out["xstar_mad"][i] = rv["xstar_mad"]
    return out


def _two_sided_p(obs, null):
    """Upper-tail p on ``|obs - median(null)|`` -- for a statistic with a sign."""
    n = np.asarray(null, float)
    n = n[np.isfinite(n)]
    if not n.size or not np.isfinite(obs):
        return np.nan
    c = float(np.median(n))
    return float((1 + int((np.abs(n - c) >= abs(obs - c)).sum())) / (n.size + 1))


def _paired_vs_sham(M_real, M_sham, s):
    """Per-patient scale structure, real minus that patient's own sham.

    The ordered sham supplies only two cohort-level realizations (one per
    source recording), far too few for a cohort p-value. Paired at the patient
    level it supplies ten, and the question becomes properly powered: is
    patient k's real profile more scale-dependent than patient k's own no-task,
    drift-carrying profile? One-sided Wilcoxon on the paired difference -- a
    difference of two matched quantities, so not a test against zero of a
    positively-biased functional.
    """
    out = dict(n_pair=0, d_slope_abs=np.nan, p_slope_abs=np.nan,
               d_slope_signed=np.nan, p_slope_signed=np.nan,
               real_slope_abs=np.nan, sham_slope_abs=np.nan)
    ra, sa = _per_patient_slopes(M_real, s), _per_patient_slopes(M_sham, s)
    m = np.isfinite(ra) & np.isfinite(sa)
    if m.sum() < 5:
        return out
    out["n_pair"] = int(m.sum())
    out["real_slope_abs"] = float(np.mean(np.abs(ra[m])))
    out["sham_slope_abs"] = float(np.mean(np.abs(sa[m])))
    for key, d in (("slope_abs", np.abs(ra[m]) - np.abs(sa[m])),
                   ("slope_signed", ra[m] - sa[m])):
        out[f"d_{key}"] = float(np.mean(d))
        dd = d[d != 0]
        if dd.size >= 5:
            out[f"p_{key}"] = float(wilcoxon(dd, alternative="greater")[1])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cs = FunctionalCells(ENC / "grid", PATIENTS_4PHASE, BANDS)
    s = np.asarray(cs.s, float)
    objects = list(cs.readouts)
    print(f"[w0s-contrast] {len(objects)} objects x {len(cs.bands)} bands x "
          f"{cs.n_scales} scales | R={cs.R} | "
          f"s in [{s.min():.3g}, {s.max():.3g}]", flush=True)

    sham = {}
    for src in SHAM_SOURCES:
        try:
            sham[src] = FunctionalCells(ENC / "sham", PATIENTS_4PHASE, BANDS,
                                        source=src, arc="ordered")
            bb = sorted({b for _, b in sham[src]._obs})
            print(f"  ordered sham [{src}]: bands {bb} R={sham[src].R}", flush=True)
        except FileNotFoundError:
            print(f"  ordered sham [{src}]: absent", flush=True)

    # -- 1. noise-level descriptives -------------------------------------- #
    nrows = []
    for b in cs.bands:
        if not cs.have(b):
            continue
        sd = {m: float(np.nanmedian(cs.surr_sd(b, m))) for m in objects}
        rho_lt = float(np.nanmedian(cs.surr_corr(b, "T_learn", "T_test")))
        nrows.append(dict(band=b, rho_null_learn_test=rho_lt,
                          **{f"sd_{m}": sd[m] for m in objects},
                          ratio_C_over_test=sd["C_learn_minus_test"] / sd["T_test"],
                          ratio_C_over_learn=sd["C_learn_minus_test"] / sd["T_learn"]))
    pd.DataFrame(nrows).to_csv(OUT / "contrast_noise_levels.csv", index=False)
    print("  noise levels -> contrast_noise_levels.csv", flush=True)

    # -- 2. per-scale gate grid ------------------------------------------- #
    gate_path = OUT / "contrast_gate_grid.csv"
    if gate_path.exists():
        print("  reusing contrast_gate_grid.csv", flush=True)
    else:
        grids = []
        for m in objects:
            cells, labs0 = [], None
            for b in cs.bands:
                if not cs.have(b):
                    continue
                O, S, labs = cs.cell(b, m)
                labs0 = labs
                for j in range(cs.n_scales):
                    cells.append(({"object": m, "band": b, "s": float(s[j])},
                                  O[:, j], S[:, :, j]))
            df = gate_grid(cells, labels=labs0, full=True, n_boot=2000,
                           expect_n=EXPECT_N, calibrate=True,
                           calibration_draws=N_DRAWS, refuse_uncalibrated=True)
            grids.append(df)
            print(f"  gated {m:22s} cleared(q<.05) "
                  f"{int((df.q < .05).sum()):3d}/{len(df)}", flush=True)
        pd.concat(grids, ignore_index=True).to_csv(gate_path, index=False)

    # -- 3. scale structure, per object x band ---------------------------- #
    rows, prof = [], []
    for m in objects + ["C_unpaired"]:
        for b in cs.bands:
            if not cs.have(b):
                continue
            if m == "C_unpaired":
                M, labs = _unpaired_margins(cs, b, "T_learn", "T_test")
                base = "C_learn_minus_test"
            else:
                M, labs = cs.margins(b, m)
                base = m
            e = effective_tests(M)
            sh = _shape_stats(M, s)
            rv = axis_reversal_gate(M, s, n_perm=N_PERM,
                                    rng=np.random.default_rng(11))
            ac = axis_cluster_gate(M, n_perm=N_PERM, rng=np.random.default_rng(7))
            nb = _null_bundle(cs, b, base, s, N_DRAWS)
            r = dict(
                object=m, band=b, n_pat=int(M.shape[0]), n_scales=cs.n_scales,
                margin_med=float(np.nanmedian(M)),
                surr_sd_med=float(np.nanmedian(cs.surr_sd(b, base))),
                n_eff_pr=e["n_eff_pr"], mean_offdiag=e["mean_offdiag"],
                n_eff_null_med=float(np.nanmedian(nb["n_eff_pr"])),
                n_eff_null_p95=float(np.nanpercentile(nb["n_eff_pr"], 95)),
                p_neff=_upper_p(e["n_eff_pr"], nb["n_eff_pr"]),
                slope_abs=sh["slope_abs"],
                slope_abs_null_med=float(np.nanmedian(nb["slope_abs"])),
                p_slope_abs=_upper_p(sh["slope_abs"], nb["slope_abs"]),
                slope_signed=sh["slope_signed"],
                slope_signed_null_med=float(np.nanmedian(nb["slope_signed"])),
                p_slope_signed=_two_sided_p(sh["slope_signed"], nb["slope_signed"]),
                shape_agreement=sh["shape_agreement"],
                shape_null_med=float(np.nanmedian(nb["shape_agreement"])),
                p_shape=_upper_p(sh["shape_agreement"], nb["shape_agreement"]),
                cluster_p=ac["p"], cluster_mass=ac["mass"],
                cluster_lo=(None if ac["cluster"] is None
                            else float(s[ac["cluster"][0]])),
                cluster_hi=(None if ac["cluster"] is None
                            else float(s[ac["cluster"][1]])),
                rev_p_pos=rv["p_pos"], rev_p_neg=rv["p_neg"],
                rev_ordered=rv["ordered"], reversal=rv["reversal"],
                n_cross=rv["n_cross"],
                n_cross_null_med=float(np.nanmedian(nb["n_cross"])),
                p_ncross=_upper_p(rv["n_cross"], nb["n_cross"]),
                xstar_med_log10=rv["xstar_med"], xstar_mad_log10=rv["xstar_mad"],
                xstar_mad_null_med=float(np.nanmedian(nb["xstar_mad"])),
                p_xstar_conc=_upper_p(-rv["xstar_mad"], -nb["xstar_mad"]),
            )
            for src, scs in sham.items():
                if m == "C_unpaired" or not scs.have(b):
                    continue
                Ms, _ = scs.margins(b, base)
                ps = _paired_vs_sham(M, Ms, s)
                rvs = axis_reversal_gate(Ms, s, n_perm=0,
                                         rng=np.random.default_rng(3))
                r.update({f"sham_{src}_{k}": v for k, v in ps.items()})
                r[f"sham_{src}_slope_abs_cohort"] = _shape_stats(Ms, s)["slope_abs"]
                r[f"sham_{src}_n_eff"] = effective_tests(Ms)["n_eff_pr"]
                r[f"sham_{src}_n_cross"] = rvs["n_cross"]
                r[f"sham_{src}_xstar_mad"] = rvs["xstar_mad"]
            rows.append(r)
            zprof = ac["z"]
            for j in range(cs.n_scales):
                prof.append(dict(object=m, band=b, s=float(s[j]),
                                 margin_med=float(np.nanmedian(M[:, j])),
                                 margin_mean=float(np.nanmean(M[:, j])),
                                 z=float(zprof[j]),
                                 frac_pos=float(np.mean(M[:, j] > 0))))
        print(f"  scale structure {m:22s} done", flush=True)

    ss = pd.DataFrame(rows)
    for col, out in (("p_neff", "q_neff"), ("p_slope_abs", "q_slope_abs"),
                     ("p_slope_signed", "q_slope_signed"), ("p_shape", "q_shape"),
                     ("cluster_p", "cluster_q")):
        ss[out] = np.nan
        for m in ss.object.unique():
            sel = (ss.object == m) & ss[col].notna()
            if sel.any():
                ss.loc[sel, out] = np.asarray(
                    bh_fdr(ss.loc[sel, col].to_numpy()), float)
    ss.to_csv(OUT / "contrast_scale_structure.csv", index=False)
    pd.DataFrame(prof).to_csv(OUT / "contrast_per_scale_profiles.csv", index=False)

    # -- 4. pre-registered criteria B1-B3 --------------------------------- #
    crit = []
    for m in ss.object.unique():
        g = ss[ss.object == m]
        b1 = bool((g.q_neff < 0.05).any())
        b2 = bool(g.reversal.any())
        shamcols = [c for c in g.columns if c.endswith("_p_slope_abs")]
        b3 = bool(g[shamcols].lt(0.05).any().any()) if shamcols else False
        crit.append(dict(
            object=m,
            n_eff_range=f"{g.n_eff_pr.min():.2f}-{g.n_eff_pr.max():.2f}",
            n_eff_null_range=(f"{g.n_eff_null_med.min():.2f}-"
                              f"{g.n_eff_null_med.max():.2f}"),
            B1_more_independent_info=b1, B1_min_q=float(g.q_neff.min()),
            B2_reversal=b2,
            B2_bands_both_clusters=int(
                ((g.rev_p_pos < .05) & (g.rev_p_neg < .05)).sum()),
            B3_beats_sham=b3,
            B3_min_p=(float(g[shamcols].min().min()) if shamcols else np.nan),
            PASSES_ALL=bool(b1 and b2 and b3)))
    pd.DataFrame(crit).to_csv(OUT / "contrast_criteria.csv", index=False)
    print(pd.DataFrame(crit).to_string(index=False), flush=True)
    print(f"[w0s-contrast] wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
