#!/usr/bin/env python3
"""audit_132 — N4 Direction B1: is the SUSTAINED reinstatement CONCENTRATED in OFC?

The whole-brain finding (audit_129 A1 / audit_125 g_TP): rest_post DWELLS in the task
ρ^coph configuration more than rest_pre (10/10, p=0.001) — a sustained state, not
transient flashes. B1 asks the localization question, in the EXACT locked methodology
used for the static β trace (audit_110/112): is that dwelling anatomically concentrated
in the OFC system (the β-trace carrier, N1.3), above each patient's own whole-brain
average system?

5-point critical preamble
  (1) Claim: the sustained post>pre dwelling-in-task-config is concentrated on OFC edges.
  (2) Null: per-patient DEMEANED concentration — Level_S minus the patient's mean Level
      over the 9 a-priori systems; a system clears only if it beats the patient's own
      whole-brain average (removes the global N1 shift; this is audit_110's logic).
  (3) Strongest alternative: OFC contacts are high-strength, and high strength shapes
      cophenetic structure → an apparent OFC concentration that is just hubness.
  (4) Does the null control it? The post−pre contrast already cancels any STATIONARY
      OFC property (OFC strength is present in pre AND post windows equally), and the
      demeaning cancels the global shift. So a false OFC hit needs OFC strength to make
      OFC windows look more task-like in POST SPECIFICALLY — which pure strength cannot
      trivially do. BUT per the mandatory rule, if OFC concentrates here we STILL run a
      matched-strength referee before any claim (this script flags it; does not claim).
  (5) Falsify: OFC demeaned concentration ≤ 0 / not top system / fails BH across systems.

Per patient, band (alpha, beta), system S (endpoint-incidence pairs, audit_112 convention):
  refs (canonical full-phase imcoh_abs):  c_test, c_pre  restricted to keep_S
  per window w (20 s cophenetic cache):   g_S(w) = ρ(c_w|S, c_test|S) − ρ(c_w|S, c_pre|S)
  Level_S = mean_post g_S − mean_pre g_S      (post dwells MORE in task config within S)
  conc_S  = Level_S − mean_over_systems(Level_S)        (demeaned concentration)
Cohort one-sided Wilcoxon (conc_OFC > 0), LOO-max, BH across the 9 systems, per band.
Secondary: within-OFC burstiness (AC1 of g_post vs time-shuffle) — expected null.

Output: data/audit/replay_states/localized_ofc_{per_patient,cohort}.csv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
CACHE = OUT / "cache"
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["alpha", "beta"]
SYSTEMS = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal", "PFC",
           "parietal", "sensorimotor", "occipital"]
MIN_PAIRS = 30


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def load_windows(pat, band, phase):
    f = CACHE / f"{pat}_{band}_{phase}_coph.npz"
    return np.load(f)["coph"].astype(np.float64) if f.exists() else None


def ac1(x):
    x = np.asarray(x, float); x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 else 0.0


def perm_p(x, fn, n, rng):
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def g_within(windows, ct_keep, cp_keep, keep):
    """Per-window task-likeness contrast restricted to keep pairs; ranks once."""
    rt = rankdata(ct_keep); rp = rankdata(cp_keep)
    out = np.empty(windows.shape[0])
    for i, c in enumerate(windows):
        ck = rankdata(c[keep])
        st = np.corrcoef(ck, rt)[0, 1]
        sp = np.corrcoef(ck, rp)[0, 1]
        out[i] = st - sp
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    rows = []
    for pat in COHORT:
        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)
        ndrop = epi_keep_mask(rdf, pat)            # True = keep (non-epi)
        for band in BANDS:
            Wt = load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs")
            Wp = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs")
            win_pre = load_windows(pat, band, "rest_pre")
            win_post = load_windows(pat, band, "rest_post")
            if Wt is None or Wp is None or win_pre is None or win_post is None:
                print(f"  [skip] {pat} {band}: missing inputs"); continue
            N = Wt.shape[0]
            if uv.shape[0] != N or win_pre.shape[1] != N * (N - 1) // 2:
                print(f"  [skip] {pat} {band}: shape mismatch N={N} "
                      f"regions={uv.shape[0]} pairs={win_pre.shape[1]}"); continue
            ct = coph(Wt); cp = coph(Wp)
            iu_i, iu_j = np.triu_indices(N, 1)
            si, sj = uv[iu_i], uv[iu_j]
            for epi_excl in (False, True):
                lvl = {}            # system -> Level_S (finite only)
                ofc_burst = None
                for S in SYSTEMS:
                    keep = (si == S) | (sj == S)
                    if epi_excl:
                        keep = keep & ndrop[iu_i] & ndrop[iu_j]
                    if int(keep.sum()) < MIN_PAIRS:
                        continue
                    g_post = g_within(win_post, ct[keep], cp[keep], keep)
                    g_pre = g_within(win_pre, ct[keep], cp[keep], keep)
                    # NaN-robust: a degenerate (constant) cophenetic sub-vector in a
                    # tiny system yields NaN windows; nanmean them and only keep the
                    # system if its Level is finite (so it cannot poison the demean).
                    lp, lq = np.nanmean(g_post), np.nanmean(g_pre)
                    if np.isfinite(lp) and np.isfinite(lq):
                        lvl[S] = float(lp - lq)
                    if S == "OFC":
                        gp = g_post[np.isfinite(g_post)]
                        if gp.size >= 10:
                            ac_post, p_post = perm_p(gp, ac1, args.n_perm, rng)
                            ofc_burst = (ac_post, p_post)
                if len(lvl) < 3 or "OFC" not in lvl:
                    continue
                mean_lvl = float(np.nanmean(list(lvl.values())))
                for S, v in lvl.items():
                    rows.append(dict(patient=pat, band=band,
                                     epi="exclude" if epi_excl else "include",
                                     system=S, level=v, conc=v - mean_lvl,
                                     ofc_ac1=ofc_burst[0] if (S == "OFC" and ofc_burst) else np.nan,
                                     ofc_burst_p=ofc_burst[1] if (S == "OFC" and ofc_burst) else np.nan))
            print(f"  {pat} {band}: done ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "localized_ofc_per_patient.csv", index=False)

    # cohort: per (band, epi, system) Wilcoxon conc>0, LOO-max, BH across systems
    coh = []
    for band in BANDS:
        for epi in ("include", "exclude"):
            sub = df[(df.band == band) & (df.epi == epi)]
            sysrows = []
            for S in SYSTEMS:
                d = sub[sub.system == S]["conc"].to_numpy()
                d = d[np.isfinite(d)]
                if len(d) < 3:
                    continue
                try:
                    p = float(wilcoxon(d, alternative="greater").pvalue)
                    loo = float(np.nanmax([wilcoxon(np.delete(d, i), alternative="greater").pvalue
                                           for i in range(len(d))]))
                except ValueError:
                    p = loo = np.nan
                r = dict(band=band, epi=epi, system=S, K=len(d),
                         median_conc=float(np.median(d)), n_pos=int((d > 0).sum()),
                         wilcoxon_p=p, loo_max_p=loo, bh_q=np.nan)
                sysrows.append(r)
            if sysrows:
                qs = bh_fdr([r["wilcoxon_p"] for r in sysrows])
                for r, q in zip(sysrows, qs):
                    r["bh_q"] = q
            coh += sysrows
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "localized_ofc_cohort.csv", index=False)

    print("\n=== B1 — sustained-reinstatement CONCENTRATION (demeaned, per band) ===")
    for band in BANDS:
        for epi in ("include", "exclude"):
            blk = cohdf[(cohdf.band == band) & (cohdf.epi == epi)].sort_values("wilcoxon_p")
            if blk.empty:
                continue
            print(f"\n  [{band} | epi-{epi}]  (median demeaned concentration, post−pre)")
            for _, r in blk.iterrows():
                star = "*" if r.bh_q < 0.05 else ("·" if r.wilcoxon_p < 0.05 else " ")
                print(f"    {r.system:16s} conc={r.median_conc:+.4f} "
                      f"{r.n_pos:2d}/{r.K} W_p={r.wilcoxon_p:.3f} "
                      f"LOO={r.loo_max_p:.3f} q={r.bh_q:.3f}{star}")
    # OFC burstiness secondary
    ob = df[(df.system == "OFC") & np.isfinite(df.ofc_burst_p)]
    if len(ob):
        print(f"\n  OFC within-system burstiness: "
              f"{int((ob.ofc_burst_p < 0.05).sum())}/{len(ob)} patient×band cells p<0.05 "
              f"(time-shuffle); median AC1={ob.ofc_ac1.median():+.3f}")
    print(f"\n[audit_132] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
