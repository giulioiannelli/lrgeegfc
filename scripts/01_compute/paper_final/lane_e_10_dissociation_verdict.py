#!/usr/bin/env python3
r"""Lane E / E5-E8 -- is encoding-versus-inference real, once drift is accounted for?

Everything here is a **difference against a reference that inherits the
construction**. Nothing is tested against zero.

Two references, used for different questions:

  matched-strength   that patient's own strength-preserving surrogate, drawn on
                     the dense FC before sparsification. Controls the graph's
                     strength distribution. It CANNOT see temporal structure,
                     because it shuffles the finished connectivity matrix.
  ordered sham       that patient's own five-phase arc carved out of a resting
                     recording in true temporal order -- same electrodes, same
                     estimator, same duration profile, same backbone, same
                     scale grid, matched-strength applied on top, and no task
                     inside it. This is the reference that reaches session drift
                     and phase adjacency.

Statistics (all per patient, per scale, knob-integrated over the four plateau
fractions before any surrogate comparison):

  margin_ms      obs - median_r surrogate                        [E7 uses this]
  Delta          margin_ms(real) - margin_ms(own ordered sham)   [E5, E6 use this]

and the objects they are applied to: the four canonical functionals, their
drift-controlled counterparts (``d = D_B - D_A`` partialled out), and the
contrast ``C = T_learn - T_test``.

One algebraic fact worth stating because it decides how ``C`` should be read.
The encode/probe role swap maps ``T_learn -> T_test`` exactly, so

    C = T_learn - T_test = T_learn - T_learn(swap)

i.e. the contrast IS the role-swap difference of the encoding functional. It is
therefore **antisymmetric under exchanging the two task blocks**, hence exactly
zero-centred under the null that the two blocks are exchangeable -- no surrogate
required for that part. What exchangeability fails on is duration (``task_test``
is longer in all 10 patients) and position in the session, which is precisely
what the ordered sham measures. So ``C`` referenced to the ordered sham is a
statistic that is calibrated by construction against role, and calibrated
empirically against order and duration.

Gate: ``axis_cluster_gate`` (sign-flip cluster mass over whole per-patient
profiles), collapsing the 16-scale axis to ONE test per band. Family: 6 bands,
BH. ``Delta T_test`` / ``Delta T_learn`` one-sided (the claim is directional);
``Delta C`` two-sided (either sign is a dissociation). Pre-registration:
``lane_e_07_preregistration_addendum.md``.

Outputs -> data/paper_final/lane_e_encinf/verdict2/
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_from_corr_stack as funcs_of,
)
from lrg_eegfc.utils.metrics.cohort_gate import axis_cluster_gate, effective_tests
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, boot_ci_mean
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("LANE_E_OUT",
                           ROOT / "data" / "paper_final" / "lane_e_encinf"))
GRID, SHAM = BASE / "rankcorr_grid" / "cells", BASE / "rankcorr_sham" / "cells"
OUT = BASE / "verdict2"
N_PERM = int(os.environ.get("LANE_E_NPERM", "10000"))
ARMS = ("ordered", "shuffled", "eqdur")

#: statistic -> how to build it from the derived functional dict
PLAIN = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
DRIFT = ("T_test_d", "T_learn_d", "T_infspec_d", "T_infspec_ped")
CONTRASTS = {"C": ("T_learn", "T_test"), "C_d": ("T_learn_d", "T_test_d")}
STATS = PLAIN + DRIFT + tuple(CONTRASTS)


def derive(R: np.ndarray) -> dict:
    """``(..., 8, 8)`` rank-correlation stack -> ``{statistic: (...)}``."""
    f = funcs_of(R)
    out = {k: f[k] for k in PLAIN + DRIFT}
    for name, (a, b) in CONTRASTS.items():
        out[name] = f[a] - f[b]
    return out


def _knob_margin(Robs: np.ndarray, Rsurr: np.ndarray) -> dict:
    """Knob-integrated matched-strength margin per statistic, ``(nS,)``.

    The statistic is averaged over the plateau fractions FIRST (that is what
    knob integration means -- the reported quantity is the fraction-integrated
    one), then referenced to the median of the fraction-integrated surrogate.
    Because one shuffled dense matrix feeds all four fractions, the surrogate
    average is taken over the identical fraction set.
    """
    o = derive(Robs)                       # each (nF, nS)
    s = derive(Rsurr)                      # each (nF, R, nS)
    return {k: np.nanmean(o[k], axis=0)
            - np.nanmedian(np.nanmean(s[k], axis=0), axis=0) for k in o}


def load_real(pats, bands) -> dict:
    """``{(pat, band): {stat: (nS,) margin}}`` plus the raw knob-integrated obs."""
    out, raw = {}, {}
    for p in pats:
        for b in bands:
            f = GRID / f"{p}__{b}.npz"
            if not f.exists():
                continue
            z = np.load(f)
            out[(p, b)] = _knob_margin(z["Robs"].astype(float),
                                       z["Rsurr"].astype(float))
            raw[(p, b)] = {k: np.nanmean(v, axis=0)
                           for k, v in derive(z["Robs"].astype(float)).items()}
    return out, raw


def load_sham(pats, bands, sources) -> dict:
    """``{(pat, band, arm): {stat: (nS,) margin}}``, averaged over reps + sources."""
    acc = {}
    for p in pats:
        for b in bands:
            for src in sources:
                f = SHAM / f"{p}__{b}__{src}.npz"
                if not f.exists():
                    continue
                z = np.load(f)
                for arm in ARMS:
                    ko, ks = f"Robs_{arm}", f"Rsurr_{arm}"
                    if ko not in z.files:
                        continue
                    Ro, Rs = z[ko].astype(float), z[ks].astype(float)
                    for rep in range(Ro.shape[0]):
                        m = _knob_margin(Ro[rep], Rs[rep])
                        acc.setdefault((p, b, arm), []).append(m)
    return {k: {s: np.nanmean([m[s] for m in v], axis=0) for s in v[0]}
            for k, v in acc.items()}


# --------------------------------------------------------------------------- #
def gate(M: np.ndarray, sided: str, rng) -> dict:
    """Axis-cluster gate on a ``(K patients, nS scales)`` margin matrix."""
    if sided == "greater":
        g = axis_cluster_gate(M, n_perm=N_PERM, rng=rng)
        return dict(p=g["p"], mass=g["mass"], cluster=g["cluster"], side="+")
    gp = axis_cluster_gate(M, n_perm=N_PERM, rng=rng)
    gm = axis_cluster_gate(-M, n_perm=N_PERM, rng=rng)
    best = gp if (np.nan_to_num(gp["p"], nan=1.0)
                  <= np.nan_to_num(gm["p"], nan=1.0)) else gm
    p2 = min(2.0 * float(np.nan_to_num(best["p"], nan=1.0)), 1.0)
    return dict(p=p2, mass=best["mass"], cluster=best["cluster"],
                side="+" if best is gp else "-")


def run_family(rows, label, sided, rng, expect_n=10):
    """Gate every band of one statistic, BH within the family, LOO the verdict."""
    out = []
    for band, M in rows:
        if M is None or M.shape[0] < expect_n:
            out.append(dict(family=label, band=band, n_pat=0 if M is None else M.shape[0],
                            p=np.nan, note="incomplete cohort"))
            continue
        g = gate(M, sided, rng)
        prof = np.nanmedian(M, axis=0)
        mu, lo, hi = boot_ci_mean(np.nanmean(M, axis=1), 10000,
                                  np.random.default_rng(7))
        ne = effective_tests(M)
        loo = [gate(np.delete(M, i, axis=0), sided, np.random.default_rng(100 + i))["p"]
               for i in range(M.shape[0])]
        out.append(dict(
            family=label, band=band, n_pat=int(M.shape[0]), p=float(g["p"]),
            side=g["side"], cluster_mass=float(g["mass"]),
            cluster=str(g["cluster"]),
            median_over_scales=float(np.nanmedian(prof)),
            max_abs_scale=float(np.nanmax(np.abs(prof))),
            mean_margin=float(mu), mean_ci_lo=float(lo), mean_ci_hi=float(hi),
            n_eff_pr=float(ne["n_eff_pr"]),
            loo_p_max=float(np.nanmax(loo)), loo_p_min=float(np.nanmin(loo)),
            loo_n_lost=int(np.sum(np.array(loo) >= 0.05))))
    df = pd.DataFrame(out)
    ok = df["p"].notna()
    df.loc[ok, "q"] = bh_fdr(df.loc[ok, "p"].to_numpy())
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="")
    ap.add_argument("--sources", default="rest_pre,rest_post")
    a = ap.parse_args()
    bands = a.bands.split(",") if a.bands else list(BRAIN_BANDS_NAMES)
    pats = list(PATIENTS_4PHASE)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260903)

    real, raw = load_real(pats, bands)
    sham = load_sham(pats, bands, a.sources.split(","))
    have = sorted({b for (_, b) in real})
    bands = [b for b in bands if b in have]
    print(f"[verdict] real cells {len(real)} | sham cells {len(sham)} | bands {bands}",
          flush=True)
    s0 = np.load(next(GRID.glob('*.npz')))["s"]

    # ---- long-form per-patient table (also the hand-off to the scale lane) ---
    rows = []
    for p in pats:
        for b in bands:
            if (p, b) not in real:
                continue
            for j, sv in enumerate(s0):
                r = dict(patient=p, band=b, scale_index=j, s=float(sv))
                for st in STATS:
                    r[f"real_{st}"] = float(real[(p, b)][st][j])
                    r[f"obs_{st}"] = float(raw[(p, b)][st][j])
                    for arm in ARMS:
                        k = (p, b, arm)
                        r[f"sham{arm[0].upper()}_{st}"] = (
                            float(sham[k][st][j]) if k in sham else np.nan)
                    k = (p, b, "ordered")
                    r[f"delta_{st}"] = (float(real[(p, b)][st][j] - sham[k][st][j])
                                        if k in sham else np.nan)
                rows.append(r)
    per_patient = pd.DataFrame(rows)
    per_patient.to_csv(OUT / "per_patient_long.csv", index=False)

    def M_of(col, band, pats_):
        sub = per_patient[per_patient.band == band]
        piv = sub.pivot_table(index="patient", columns="scale_index", values=col)
        piv = piv.reindex([p for p in pats_ if p in piv.index])
        return piv.to_numpy() if len(piv) else None

    # ---- E5: is the persistence headline a drift artifact? ------------------
    e5 = pd.concat([
        run_family([(b, M_of(f"delta_{st}", b, pats)) for b in bands],
                   f"E5 delta_{st}", "greater", rng)
        for st in ("T_test", "T_learn")])
    # ---- E6: do encoding and inference dissociate? --------------------------
    e6 = run_family([(b, M_of("delta_C", b, pats)) for b in bands],
                    "E6 delta_C", "two", rng)
    # ---- E7: does drift control preserve anything? V1 sham, V2 real ---------
    v1 = pd.concat([
        run_family([(b, M_of(f"shamO_{st}", b, pats)) for b in bands],
                   f"E7-V1 sham_ordered {st}", "greater", rng)
        for st in DRIFT + PLAIN])
    v2 = pd.concat([
        run_family([(b, M_of(f"real_{st}", b, pats)) for b in bands],
                   f"E7-V2 real {st}", "greater", rng)
        for st in DRIFT + PLAIN])
    # ---- E6b: contrast under drift control, and against the other arms ------
    e6b = pd.concat([
        run_family([(b, M_of("delta_C_d", b, pats)) for b in bands],
                   "E6b delta_C_d", "two", rng),
        run_family([(b, M_of("real_C", b, pats)) for b in bands],
                   "E6b real_C (matched-strength only)", "two", rng),
        run_family([(b, M_of("shamO_C", b, pats)) for b in bands],
                   "E6b sham_ordered C", "two", rng)])

    for name, df in (("e5_persistence_vs_drift", e5), ("e6_contrast", e6),
                     ("e7_v1_sham_calibration", v1), ("e7_v2_real_retention", v2),
                     ("e6b_contrast_variants", e6b)):
        df.to_csv(OUT / f"{name}.csv", index=False)

    # ---- descriptive: how much of the real margin the ordered sham makes ----
    desc = []
    for b in bands:
        for st in STATS:
            r = np.nanmedian(np.nanmedian(M_of(f"real_{st}", b, pats), axis=0))
            s_ = np.nanmedian(np.nanmedian(M_of(f"shamO_{st}", b, pats), axis=0))
            sh = np.nanmedian(np.nanmedian(M_of(f"shamS_{st}", b, pats), axis=0))
            eq = np.nanmedian(np.nanmedian(M_of(f"shamE_{st}", b, pats), axis=0))
            desc.append(dict(band=b, stat=st, real=r, sham_ordered=s_,
                             sham_shuffled=sh, sham_eqdur=eq, delta=r - s_,
                             pct_sham_reproduces=(100 * s_ / r) if r else np.nan))
    pd.DataFrame(desc).to_csv(OUT / "descriptive_margins.csv", index=False)

    # ---- hand-off CSV for the scale lane, keyed as agreed -------------------
    cols = ["patient", "band", "scale_index", "s"]
    ho = per_patient[cols + ["obs_C", "real_C", "shamO_C", "shamS_C", "shamE_C",
                             "delta_C", "obs_T_learn", "obs_T_test",
                             "real_T_learn", "real_T_test"]].copy()
    ho.rename(columns={"obs_C": "C_real_raw", "real_C": "margin_ms",
                       "shamO_C": "C_sham_ordered", "shamS_C": "C_sham_shuffled",
                       "shamE_C": "C_sham_eqdur"}, inplace=True)
    ho.to_csv(BASE / "verdict" / "contrast_paired.csv", index=False)

    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="lane_e_dissociation_verdict", bands=bands, cohort=pats,
        stats=list(STATS), arms=list(ARMS), n_perm=N_PERM,
        gate="axis_cluster_gate (sign-flip cluster mass, 16 scales -> 1 test)",
        family="6 bands within one statistic, BH",
        reference="matched-strength for E7; real minus own ordered sham for E5/E6",
    ), indent=2))
    print(f"[verdict] -> {OUT}", flush=True)

    for nm, df in (("E5 persistence vs drift", e5), ("E6 contrast", e6)):
        print(f"\n=== {nm} ===")
        print(df[["family", "band", "p", "q", "side", "median_over_scales",
                  "loo_n_lost"]].to_string(index=False))


if __name__ == "__main__":
    main()
