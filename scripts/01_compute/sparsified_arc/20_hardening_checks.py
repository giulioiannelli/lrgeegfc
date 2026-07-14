#!/usr/bin/env python3
"""Two hardening checks for the 'only-multiscale cognitive feature' headline.

CHECK 1 (controls-ladder robustness): is beta ENCODING (T_learn) really cophenetic-only?
  Reproduce the cohort gate from per_cell (Wilcoxon(obs-surr_p50, greater), n=10), then
  stress it: leave-one-patient-out (does raw_fc ever cross .05? does cophenetic stay <.05?),
  n_above, and the null-separation magnitude. Reads controls_ladder{_TAG}/per_cell.csv so
  it can run on the R=200 cache now and the R=1000 rerun later (SA_TAG=beta_R1000).

CHECK 2 (scale-span flatness): is the beta trace / encoding 'scale-invariant' (constant
  effect across scales) or merely 'significant at every scale'? Per-scale effect
  (obs - surr_p50) over 16 scales: Friedman across scales (H0 flat), coefficient of
  variation of the per-scale median, and a per-patient log-s trend (Wilcoxon of Spearman
  rho vs 0). Verdict: scale-invariant (Friedman ns + low CV + no trend) vs spans-all-scales.

Output: data/sparsified_arc/hardening/{ladder_robustness,scale_flatness}.csv
"""
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, friedmanchisquare, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

TAG = os.environ.get("SA_TAG", "")               # "" -> R=200 all-band; "beta_R1000" -> rerun
LAD = ROOT / f"data/sparsified_arc/controls_ladder{('_'+TAG) if TAG else ''}/per_cell.csv"
MS = ROOT / "data/sparsified_arc/ms_mst020/per_patient_scale.csv"
ENC = ROOT / "data/sparsified_arc/enc_inf_arc_mst020"
OUT = ROOT / "data/sparsified_arc/hardening"


def gate(diff):
    diff = np.asarray(diff, float); diff = diff[np.isfinite(diff)]
    if diff.size < 5 or not np.any(diff != 0):
        return np.nan
    return float(wilcoxon(diff, alternative="greater")[1])


def check1_ladder():
    d = pd.read_csv(LAD)
    band = "beta"
    descs = ["coph_meso", "coph_taumin", "geodesic", "clustering", "raw_fc",
             "resistance", "strength"]
    func = "T_learn"   # encoding
    rows = []
    for desc in descs:
        x = d[(d.band == band) & (d.descriptor == desc)]
        if x.empty:
            continue
        obs = x[f"{func}_obs"].values
        surr = x[f"{func}_surr_p50"].values
        diff = obs - surr
        g = gate(diff)
        n_above = int(np.sum(diff > 0))
        # LOO
        loo = [gate(np.delete(diff, i)) for i in range(len(diff))]
        loo = [p for p in loo if np.isfinite(p)]
        loo_arr = np.array(loo) if loo else np.array([np.nan])
        rows.append(dict(descriptor=desc, functional=func, band=band, n=len(diff),
                         gate_p=g, n_above=f"{n_above}/{len(diff)}",
                         med_obs=float(np.median(obs)), med_surr=float(np.median(surr)),
                         med_effect=float(np.median(diff)),
                         loo_min=float(np.nanmin(loo_arr)), loo_max=float(np.nanmax(loo_arr)),
                         loo_med=float(np.nanmedian(loo_arr)),
                         loo_n_sig=f"{int(np.sum(loo_arr<0.05))}/{len(loo)}"))
    return pd.DataFrame(rows)


def _flatness(pat_scale_effect, s_vals, label):
    """pat_scale_effect: (n_pat, n_scale) effect matrix; s_vals: scales."""
    E = np.asarray(pat_scale_effect, float)
    keep = np.all(np.isfinite(E), axis=1)
    E = E[keep]; n_pat, n_s = E.shape
    per_scale_med = np.median(E, axis=0)
    per_scale_gate = np.array([gate(E[:, j]) for j in range(n_s)])
    n_fire = int(np.sum(per_scale_gate < 0.05))
    # Friedman across scales (H0: same distribution every scale)
    try:
        fr_stat, fr_p = friedmanchisquare(*[E[:, j] for j in range(n_s)])
    except Exception:
        fr_stat, fr_p = np.nan, np.nan
    cv = float(np.std(per_scale_med) / np.mean(per_scale_med)) if np.mean(per_scale_med) != 0 else np.nan
    # per-patient log-s trend
    rhos = [spearmanr(np.log(s_vals), E[i])[0] for i in range(n_pat)]
    rhos = [r for r in rhos if np.isfinite(r)]
    trend_p = float(wilcoxon(rhos, alternative="two-sided")[1]) if len(rhos) >= 5 and np.any(np.array(rhos) != 0) else np.nan
    verdict = ("scale-invariant" if (np.isfinite(fr_p) and fr_p > 0.05 and cv < 0.35
                                     and (not np.isfinite(trend_p) or trend_p > 0.05))
               else "spans-all-scales (magnitude varies)")
    return dict(target=label, n_pat=n_pat, n_scales=n_s, n_fire=f"{n_fire}/{n_s}",
                friedman_p=float(fr_p), cv_median=cv, trend_rho_med=float(np.median(rhos)),
                trend_p=trend_p, med_effect=float(np.median(per_scale_med)), verdict=verdict), per_scale_med, per_scale_gate, s_vals


def check2_flatness():
    out_rows, detail = [], []
    # (a) TRACE beta/alpha/delta from ms_mst020 (beta scale-invariant vs alpha single-scale)
    ms = pd.read_csv(MS)
    for band in ["beta", "alpha", "delta"]:
        b = ms[ms.band == band].copy()
        s_vals = np.sort(b.s.unique())
        pats = sorted(b.patient.unique())
        E = np.full((len(pats), len(s_vals)), np.nan)
        for i, p in enumerate(pats):
            for j, s in enumerate(s_vals):
                r = b[(b.patient == p) & np.isclose(b.s, s)]
                if not r.empty:
                    E[i, j] = r.obs_rho.iloc[0] - r.surr_p50.iloc[0]
        rec, med, gp, sv = _flatness(E, s_vals, f"{band}_trace")
        out_rows.append(rec)
        for j in range(len(sv)):
            detail.append(dict(target=f"{band}_trace", s=float(sv[j]), med_effect=float(med[j]), gate_p=float(gp[j])))
    # (b) ENCODING beta from enc_inf npz
    files = sorted(glob.glob(str(ENC / "beta" / "*.npz")))
    s_grid = np.load(files[0])["s"]
    E2 = np.full((len(files), len(s_grid)), np.nan)
    for i, f in enumerate(files):
        dd = np.load(f)
        E2[i] = dd["T_learn__obs"] - dd["T_learn__surr_p50"]
    rec2, med2, gp2, sv2 = _flatness(E2, s_grid, "beta_encoding")
    out_rows.append(rec2)
    for j in range(len(sv2)):
        detail.append(dict(target="beta_encoding", s=float(sv2[j]), med_effect=float(med2[j]), gate_p=float(gp2[j])))
    return pd.DataFrame(out_rows), pd.DataFrame(detail)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"=== CHECK 1 — controls-ladder robustness (beta ENCODING) [tag={TAG or 'R200'}] ===")
    c1 = check1_ladder(); c1.to_csv(OUT / f"ladder_robustness{('_'+TAG) if TAG else ''}.csv", index=False)
    print(c1.to_string(index=False))
    print("\n  READ: raw_fc gate_p + LOO range vs .05 ; cophenetic loo_n_sig should be full.\n")

    print("=== CHECK 2 — scale-span flatness (beta trace + encoding) ===")
    c2, det = check2_flatness()
    c2.to_csv(OUT / "scale_flatness.csv", index=False); det.to_csv(OUT / "scale_flatness_detail.csv", index=False)
    print(c2.to_string(index=False))
    print("\n  per-scale median effect (beta_trace / beta_encoding):")
    for tgt in ["beta_trace", "beta_encoding"]:
        dt = det[det.target == tgt].sort_values("s")
        cells = " ".join(f"{m:.2f}{'*' if p<0.05 else ' '}" for m, p in zip(dt.med_effect, dt.gate_p))
        print(f"    {tgt:14s}: {cells}")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
