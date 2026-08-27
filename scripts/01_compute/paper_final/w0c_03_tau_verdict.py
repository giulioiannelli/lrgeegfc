#!/usr/bin/env python3
r"""W0-C part 2: what the diffusion scale actually buys, decomposed into three claims.

The three are tested and reported separately because they have been conflated,
and because they can come out differently -- which they do.

(a) DETECTION -- does any scale carry cross-phase trace that the dense raw FC
    does not already carry? Gate ``coph|raw`` (the cophenetic geometry with its
    dependence on the raw edges removed) at every scale, under two
    residualizers of very different strength, and gate the reverse ``raw|coph``
    for asymmetry. A negative under the WEAKER residualizer is the informative
    one: if even removing only the monotone part of raw kills the incremental
    trace, no residualization artifact is available as an excuse.

(b) CHARACTERIZATION -- is the shape of the trace-vs-scale profile
    band-discriminative ACROSS PATIENTS? Per-patient margin profiles, z-scored
    over scale so only shape survives, then leave-one-patient-out
    nearest-centroid classification of band identity against a
    within-patient label-permutation null. Plus paired per-patient contrasts of
    profile centroid and profile flatness. A cohort Friedman that fails to
    reject is reported for continuity and is *not* treated as evidence for
    invariance; the equivalence question is answered by stating the smallest
    flatness difference the design could have detected.

(c) SELECTION -- the hierarchy rejects what raw accepts. Verified from the gate
    grid, then attacked mechanistically: the edge-locality hypothesis says the
    rejected bands' trace is carried by a few strong pairs that hierarchical
    coarse-graining absorbs. Probed three ways -- concentration of the per-pair
    contributions (Gini, top-share), ablation of the top contributors, and how
    much of the raw per-pair reorganization the hierarchy can express at each
    scale. Tested as PAIRED per-patient contrasts, because the previous attempt
    tested a cross-band correlation at n = 6 bands and could not have detected
    anything.

(d) The scale axis gets a unit. Effective resolved components ``N_eff(s)``,
    effective communication neighbourhood ``m(s)``, and the physical reach
    ``ell(s)`` in mm, so no scale is described in words without a number.

Reads the cells from ``w0c_01_gate_and_tau_grid``.
Outputs (data/paper_final/w0c_gate_tau/):
  detection_grid.csv, characterization.csv, selection.csv, scale_units.csv
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, spearmanr, wilcoxon

from lrg_eegfc.utils.metrics.cohort_gate import DESCRIPTIVE_ALPHA, gate_grid
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
BASE = Path(os.environ.get("W0C_BASE", ROOT / "data" / "paper_final" / "w0c_gate_tau"))
GRID = BASE / "grid"

import importlib.util as _ilu                                          # noqa: E402
_spec = _ilu.spec_from_file_location(
    "w0c_gate_verdict", Path(__file__).with_name("w0c_02_gate_verdict.py"))
_gv = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_gv)
load_cells, cells_for = _gv.load_cells, _gv.cells_for


def _paired(x, y):
    """Paired one-sided Wilcoxon ``x > y`` and ``x < y`` on finite pairs."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.allclose(x[m], y[m]):
        return np.nan, np.nan, int(m.sum())
    g = float(wilcoxon(x[m], y[m], alternative="greater").pvalue)
    l = float(wilcoxon(x[m], y[m], alternative="less").pvalue)
    return g, l, int(m.sum())


# --------------------------------------------------------------------------- #
def detection(D_all):
    print("\n" + "=" * 78, flush=True)
    print("(a) DETECTION -- does the hierarchy see anything raw does not?", flush=True)
    print("=" * 78, flush=True)
    rows = []
    for mi, mname in ((1, "coph_res_lin"), (2, "coph_res_np"), (3, "raw_res_lin")):
        D = D_all[mi]
        g = gate_grid(cells_for(D), labels=D["pats"], rng=np.random.default_rng(0))
        rows.append(g)
        print(f"\n  {mname}  (whole-grid BH over {len(g)} cells; * = q < "
              f"{DESCRIPTIVE_ALPHA})", flush=True)
        print("    band        cleared   q_min     s@q_min   p_min    margin_med@q_min",
              flush=True)
        for b in D["bands"]:
            x = g[g.band == b]
            if x.q.notna().sum() == 0:
                print(f"    {b:11s}  all cells degenerate", flush=True)
                continue
            r = x.loc[x.q.idxmin()]
            print(f"    {b:11s} {int((x.q < DESCRIPTIVE_ALPHA).sum()):2d}/{len(x):2d}    "
                  f"{'*' if r.q < DESCRIPTIVE_ALPHA else ' '}{float(r.q):6.4f}  "
                  f"{float(r.s):8.2f}  {float(r.p):7.4f}  {float(r.margin_med):+.4f}",
                  flush=True)
    det = pd.concat(rows, ignore_index=True)
    det.to_csv(BASE / "detection_grid.csv", index=False)
    return det


# --------------------------------------------------------------------------- #
def _profiles(D):
    """``(nB, K, nS)`` per-patient margin profiles for the cophenetic trace."""
    return D["obs"] - np.nanmedian(D["surr"], axis=3)


def _shape(M):
    """z-score each profile over scale: amplitude out, shape in."""
    mu = np.nanmean(M, axis=-1, keepdims=True)
    sd = np.nanstd(M, axis=-1, keepdims=True)
    return np.divide(M - mu, sd, out=np.full_like(M, np.nan), where=sd > 0)


def _unit(X):
    """Row-wise mean-centred, unit-norm copy (so a dot product is a correlation)."""
    Y = X - X.mean(axis=-1, keepdims=True)
    n = np.linalg.norm(Y, axis=-1, keepdims=True)
    return np.divide(Y, n, out=np.zeros_like(Y), where=n > 0)


def _loo_nearest_centroid(S):
    """LOO nearest-centroid accuracy of band identity from profile shape.

    ``S`` is ``(nB, K, nS)``. For each held-out patient the band centroids are
    rebuilt from the other patients, and each of that patient's band profiles is
    assigned to the centroid it correlates with most. Vectorised: the leave-one-
    out centroid is the full sum minus the held-out row, so no inner loop.
    """
    nB, K, _ = S.shape
    ok = np.isfinite(S).all(axis=-1)                                 # (nB, K)
    Sf = np.where(ok[..., None], S, 0.0)
    tot = Sf.sum(axis=1, keepdims=True)                              # (nB, 1, nS)
    cnt = ok.sum(axis=1, keepdims=True)[..., None]                   # (nB, 1, 1)
    cen = (tot - Sf) / np.maximum(cnt - ok[..., None], 1)            # (nB, K, nS)
    U = _unit(Sf)                                                    # (nB, K, nS)
    C = _unit(cen)
    # corr[b, k, c] = profile of band b in patient k vs centroid of band c in k
    corr = np.einsum("bks,cks->bkc", U, C)
    pred = np.argmax(np.where(np.isfinite(corr), corr, -np.inf), axis=2)  # (nB, K)
    conf = np.zeros((nB, nB), int)
    hits = total = 0
    for b in range(nB):
        for k in range(K):
            if not ok[b, k]:
                continue
            conf[b, pred[b, k]] += 1
            hits += int(pred[b, k] == b)
            total += 1
    return hits / total if total else np.nan, conf, total


def characterization(D, out):
    print("\n" + "=" * 78, flush=True)
    print("(b) CHARACTERIZATION -- is the scale profile band-discriminative?", flush=True)
    print("=" * 78, flush=True)
    bands, pats, s = D["bands"], D["pats"], D["s"]
    M = _profiles(D)                                                 # (nB, K, nS)
    S = _shape(M)

    acc, conf, tot = _loo_nearest_centroid(S)
    rng = np.random.default_rng(0)
    n_perm = 2000
    null = np.full(n_perm, np.nan)
    for i in range(n_perm):
        Sp = np.empty_like(S)
        for k in range(len(pats)):
            Sp[:, k, :] = S[rng.permutation(len(bands)), k, :]
        null[i] = _loo_nearest_centroid(Sp)[0]
    p_perm = (1 + int(np.sum(null >= acc))) / (n_perm + 1)
    print(f"\n  LOO nearest-centroid band classification from profile SHAPE:", flush=True)
    print(f"    accuracy {acc:.3f} on {tot} profiles (chance {1/len(bands):.3f}); "
          f"within-patient label-permutation p = {p_perm:.4f} "
          f"(null mean {null.mean():.3f}, 95th pct {np.percentile(null,95):.3f})",
          flush=True)
    print("    confusion (rows = true band, cols = predicted):", flush=True)
    print("      " + "".join(f"{b[:5]:>7s}" for b in bands), flush=True)
    for i, b in enumerate(bands):
        print(f"      {b:11s}" + "".join(f"{conf[i,j]:7d}" for j in range(len(bands))),
              flush=True)

    # per-patient profile summaries
    ls = np.log(s)
    rows = []
    for ib, b in enumerate(bands):
        for ip, p in enumerate(pats):
            m = M[ib, ip]
            w = np.clip(m, 0.0, None)
            cen = float(np.sum(w * ls) / np.sum(w)) if np.sum(w) > 0 else np.nan
            fin = m[np.isfinite(m)]
            cv = float(np.std(fin) / np.mean(fin)) if fin.size and np.mean(fin) > 0 else np.nan
            rng_rel = (float((np.max(fin) - np.min(fin)) / np.max(fin))
                       if fin.size and np.max(fin) > 0 else np.nan)
            has = np.isfinite(m).any()
            rows.append(dict(band=b, patient=p, centroid_log_s=cen,
                             s_centroid=float(np.exp(cen)) if np.isfinite(cen) else np.nan,
                             cv=cv, rel_range=rng_rel,
                             margin_max=float(np.nanmax(m)) if has else np.nan,
                             s_at_max=float(s[int(np.nanargmax(m))]) if has else np.nan,
                             n_pos=int(np.nansum(m > 0))))
    prof = pd.DataFrame(rows)
    prof.to_csv(out, index=False)

    print("\n  per-patient profile summaries (cohort median [IQR]):", flush=True)
    print("    band        s_centroid          CV of margin        s at max margin",
          flush=True)
    for b in bands:
        x = prof[prof.band == b]
        print(f"    {b:11s} {x.s_centroid.median():6.2f} "
              f"[{x.s_centroid.quantile(.25):5.2f},{x.s_centroid.quantile(.75):6.2f}]  "
              f"{x.cv.median():6.2f} [{x.cv.quantile(.25):5.2f},{x.cv.quantile(.75):6.2f}] "
              f"  {x.s_at_max.median():7.2f}", flush=True)

    # A flatness statistic defined for every patient. CV and relative range both
    # require a positive mean / max margin and silently drop patients; the rank
    # correlation between the margin and log s does not, and it is signed, so it
    # separates "flat" from "tuned" without discarding anyone.
    slope = {}
    for ib, b in enumerate(bands):
        v = []
        for ip in range(len(pats)):
            m = M[ib, ip]
            ok = np.isfinite(m)
            v.append(spearmanr(ls[ok], m[ok])[0] if ok.sum() >= 5 else np.nan)
        slope[b] = np.array(v, float)
    print("\n  scale-dependence of the margin profile, defined for every patient:",
          flush=True)
    print("    rho(margin, log s) per patient -- 0 = flat, |rho| large = tuned",
          flush=True)
    print("    band        median rho   median |rho|   n", flush=True)
    for b in bands:
        v = slope[b]
        print(f"    {b:11s} {np.nanmedian(v):+8.3f}     {np.nanmedian(np.abs(v)):8.3f}   "
              f"{int(np.isfinite(v).sum())}", flush=True)
    print("\n    paired |rho| contrasts vs beta (is beta flatter?):", flush=True)
    for b in bands:
        if b == "beta":
            continue
        g, l, n = _paired(np.abs(slope[b]), np.abs(slope["beta"]))
        print(f"      {b:11s} |rho| {np.nanmedian(np.abs(slope[b])):.3f} vs beta "
              f"{np.nanmedian(np.abs(slope['beta'])):.3f}  n={n}  "
              f"p(band more scale-tuned than beta)={g:.4f}", flush=True)

    print("\n  paired per-patient contrasts (alpha vs beta):", flush=True)
    a = prof[prof.band == "alpha"].set_index("patient")
    bt = prof[prof.band == "beta"].set_index("patient")
    common = [p for p in pats if p in a.index and p in bt.index]
    for col, lab in (("centroid_log_s", "profile centroid (log s)"),
                     ("cv", "profile CV (scale-dependence)"),
                     ("rel_range", "relative range")):
        x = a.loc[common, col].to_numpy(float)
        y = bt.loc[common, col].to_numpy(float)
        g, l, n = _paired(x, y)
        print(f"    {lab:32s} alpha {np.nanmedian(x):+7.3f} vs beta "
              f"{np.nanmedian(y):+7.3f}  n={n}  p(a>b)={g:.4f}  p(a<b)={l:.4f}",
              flush=True)

    print("\n  cohort Friedman across scales, per band (NOT evidence for invariance):",
          flush=True)
    for ib, b in enumerate(bands):
        X = M[ib]
        ok = np.isfinite(X).all(axis=1)
        if ok.sum() >= 3:
            st, pv = friedmanchisquare(*[X[ok, j] for j in range(X.shape[1])])
            print(f"    {b:11s} chi2={st:8.2f}  p={pv:.4f}  (n={int(ok.sum())} patients, "
                  f"{X.shape[1]} scales)", flush=True)

    # how large a scale-dependence could this design have detected?
    print("\n  power of the invariance claim (what a null Friedman could NOT rule out):",
          flush=True)
    ib_b = bands.index("beta")
    Xb = M[ib_b][np.isfinite(M[ib_b]).all(axis=1)]
    sd_within = float(np.median(np.std(Xb, axis=1)))
    sd_between = float(np.median(np.std(Xb, axis=0)))
    rng2 = np.random.default_rng(3)
    for eff in (0.25, 0.5, 0.75, 1.0, 1.5):
        hits = 0
        n_sim = 400
        for _ in range(n_sim):
            base = rng2.normal(0, sd_between, size=(Xb.shape[0], 1))
            trend = eff * sd_within * np.linspace(-0.5, 0.5, Xb.shape[1])[None, :]
            noise = rng2.normal(0, sd_within, size=Xb.shape)
            Y = base + trend + noise
            hits += int(friedmanchisquare(*[Y[:, j] for j in range(Y.shape[1])]).pvalue < 0.05)
        print(f"    monotone scale trend of {eff:.2f} x within-patient SD "
              f"-> Friedman rejects {hits/n_sim:5.1%} of the time", flush=True)
    return prof


# --------------------------------------------------------------------------- #
def selection(D, desc, gate_coph, gate_raw):
    print("\n" + "=" * 78, flush=True)
    print("(c) SELECTION -- what the hierarchy rejects, and why", flush=True)
    print("=" * 78, flush=True)
    bands, pats = D["bands"], D["pats"]

    # --- the selection itself, tested rather than inferred ----------------- #
    # "raw clears and the hierarchy does not" compares two verdicts, and a
    # difference between a significant and a non-significant result is not
    # itself significant. The claim is tested directly here: within each
    # patient, is the raw margin larger than the cophenetic margin, and does
    # that gap differ between a band the hierarchy keeps and one it rejects?
    m_raw = D["obs_raw"] - np.nanmedian(D["surr_raw"], axis=2)        # (nB, K)
    m_coph = D["obs"] - np.nanmedian(D["surr"], axis=3)               # (nB, K, nS)
    print("\n  DOES THE HIERARCHY REALLY DISCARD ANYTHING? paired within patient.",
          flush=True)
    print("  gap(s) = raw margin - cophenetic margin, per patient. Positive means", flush=True)
    print("  the raw representation registers more of the trace than the hierarchy.",
          flush=True)
    print("    band        gap at s=1.04   p(gap>0)   gap at s=6.40   p(gap>0)   "
          "median over all s", flush=True)
    sel_rows = []
    idx_sel = [int(np.argmin(np.abs(D["s"] - v))) for v in (1.0, 6.4)]
    for ib, b in enumerate(bands):
        gaps = m_raw[ib][:, None] - m_coph[ib]                        # (K, nS)
        cells = []
        for j in idx_sel:
            gj = gaps[:, j]
            gj = gj[np.isfinite(gj)]
            pj = (float(wilcoxon(gj, alternative="greater").pvalue)
                  if gj.size >= 3 and not np.allclose(gj, 0) else np.nan)
            cells.append((float(np.median(gj)), pj))
        med_all = float(np.nanmedian(gaps))
        sel_rows.append(dict(band=b, gap_s1=cells[0][0], p_s1=cells[0][1],
                             gap_s6=cells[1][0], p_s6=cells[1][1], gap_med=med_all))
        print(f"    {b:11s} {cells[0][0]:+8.4f}      {cells[0][1]:7.4f}   "
              f"{cells[1][0]:+8.4f}      {cells[1][1]:7.4f}   {med_all:+8.4f}",
              flush=True)
    print("\n  and the interaction -- is the gap LARGER for a rejected band than for beta?",
          flush=True)
    ibb = bands.index("beta")
    for ib, b in enumerate(bands):
        if b == "beta":
            continue
        for j, sv in zip(idx_sel, (1.04, 6.40)):
            d = ((m_raw[ib] - m_coph[ib][:, j]) - (m_raw[ibb] - m_coph[ibb][:, j]))
            d = d[np.isfinite(d)]
            pj = (float(wilcoxon(d, alternative="greater").pvalue)
                  if d.size >= 3 and not np.allclose(d, 0) else np.nan)
            print(f"    {b:11s} vs beta at s={sv:5.2f}: median gap difference "
                  f"{np.median(d):+.4f}  p(larger than beta)={pj:.4f}  n={d.size}",
                  flush=True)
    pd.DataFrame(sel_rows).to_csv(BASE / "selection_gap.csv", index=False)

    print("\n  raw vs hierarchy, side by side (the two verdicts, for reference):",
          flush=True)
    print("    band        raw q     hierarchy: cleared scales   verdict", flush=True)
    for b in bands:
        rq = float(gate_raw[gate_raw.band == b].q.iloc[0])
        x = gate_coph[gate_coph.band == b]
        n = int((x.q < DESCRIPTIVE_ALPHA).sum())
        raw_ok = rq < DESCRIPTIVE_ALPHA
        v = ("kept" if raw_ok and n > 0 else
             "REJECTED by hierarchy" if raw_ok and n == 0 else
             "hierarchy-only" if not raw_ok and n > 0 else "null in both")
        print(f"    {b:11s} {'*' if raw_ok else ' '}{rq:6.4f}   {n:2d}/{len(x):2d}"
              f"                     {v}", flush=True)

    # ---- edge-locality ---------------------------------------------------- #
    d = desc.copy()
    abl_cols = [c for c in d.columns if c.startswith("abl_") and c != "abl_full"]
    for c in abl_cols:
        d["ret_" + c[4:]] = np.where(d.abl_full > 0.02, d[c] / d.abl_full, np.nan)
    print("\n  edge-locality of the RAW trace (cohort median):", flush=True)
    print("    band        Gini    top1%   top5%   retained after deleting top "
          "0.1%/1%/5% of pairs", flush=True)
    for b in bands:
        x = d[d.band == b]
        print(f"    {b:11s} {x.gini_raw.median():.3f}  {x.top1pct_raw.median():.3f}  "
              f"{x.top5pct_raw.median():.3f}   {x['ret_0.001'].median():+.3f} / "
              f"{x['ret_0.01'].median():+.3f} / {x['ret_0.05'].median():+.3f}", flush=True)

    print("\n  paired per-patient contrasts vs beta (the edge-locality test).", flush=True)
    print("    The hypothesis predicts a REJECTED band is MORE concentrated than beta", flush=True)
    print("    (larger top-1% share) and LESS robust to deleting its top pairs", flush=True)
    print("    (smaller retained fraction).", flush=True)
    piv = {c: d.pivot_table(index="patient", columns="band", values=c)
           for c in ("gini_raw", "top1pct_raw", "ret_0.01", "ret_0.05")}
    print("    band vs beta   top1% median (band/beta)  p(band>beta)  "
          "retained median (band/beta)  p(band<beta)", flush=True)
    for other in bands:
        if other == "beta":
            continue
        t1, t2 = piv["top1pct_raw"], piv["ret_0.01"]
        if other not in t1.columns:
            continue
        x1, y1 = t1[other].to_numpy(float), t1["beta"].to_numpy(float)
        x2, y2 = t2[other].to_numpy(float), t2["beta"].to_numpy(float)
        g1, _, n1 = _paired(x1, y1)
        _, l2, n2 = _paired(x2, y2)
        print(f"    {other:11s}    {np.nanmedian(x1):.3f} / {np.nanmedian(y1):.3f}"
              f"           {g1:8.4f}      {np.nanmedian(x2):+.3f} / "
              f"{np.nanmedian(y2):+.3f}          {l2:8.4f}  (n={n1}/{n2})", flush=True)

    # ---- hierarchy absorption --------------------------------------------- #
    print("\n  how much of the RAW per-pair task reorganization the hierarchy can "
          "express (Spearman^2):", flush=True)
    s = D["s"]
    idx = [int(np.argmin(np.abs(s - v))) for v in (0.05, 1.0, 5.6, 30.0, 180.0)]
    print("    band        " + "".join(f"s={s[i]:<8.2f}" for i in idx), flush=True)
    rows = []
    for ib, b in enumerate(bands):
        v = np.nanmedian(D["absorb_r2"][ib], axis=0)
        print(f"    {b:11s} " + "".join(f"{v[i]:<10.3f}" for i in idx), flush=True)
        for js in range(s.size):
            rows.append(dict(band=b, s=float(s[js]),
                             absorb_r2_med=float(np.nanmedian(D["absorb_r2"][ib, :, js])),
                             n_eff_med=float(np.nanmedian(D["n_eff"][ib, :, js])),
                             m_comm_med=float(np.nanmedian(D["m_comm"][ib, :, js]))))
    print("\n  paired per-patient absorption contrasts at each landmark scale "
          "(beta vs low_gamma):", flush=True)
    ibb, ibg = bands.index("beta"), bands.index("low_gamma")
    for i in idx:
        x = D["absorb_r2"][ibb, :, i]
        y = D["absorb_r2"][ibg, :, i]
        g, l, n = _paired(x, y)
        print(f"    s={s[i]:7.2f}  beta {np.nanmedian(x):.3f} vs low_gamma "
              f"{np.nanmedian(y):.3f}  n={n}  p(beta>low_gamma)={g:.4f}", flush=True)
    pd.DataFrame(rows).to_csv(BASE / "selection.csv", index=False)
    d.to_csv(BASE / "edge_locality.csv", index=False)
    return d


# --------------------------------------------------------------------------- #
def scale_units(D):
    print("\n" + "=" * 78, flush=True)
    print("(d) THE SCALE AXIS, IN UNITS", flush=True)
    print("=" * 78, flush=True)
    s = D["s"]
    N = np.nanmedian([np.nanmedian(D["n_eff"][:, ip, 0]) for ip in range(len(D["pats"]))])
    print(f"\n  cohort median over all (patient, band); backbone node count ~ {N:.0f}",
          flush=True)
    print("       s     N_eff (resolved components)   m (communication neighbourhood)",
          flush=True)
    rows = []
    for js, sv in enumerate(s):
        ne = float(np.nanmedian(D["n_eff"][:, :, js]))
        mc = float(np.nanmedian(D["m_comm"][:, :, js]))
        rows.append(dict(s=float(sv), n_eff_med=ne, m_comm_med=mc))
        if js % 3 == 0 or js == s.size - 1:
            print(f"    {sv:8.2f}          {ne:8.1f}                     {mc:8.1f}",
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(BASE / "scale_units.csv", index=False)
    return df


def main():
    D_all = {mi: load_cells(mi, nm) for mi, nm in
             ((0, "coph"), (1, "coph_res_lin"), (2, "coph_res_np"), (3, "raw_res_lin"))}
    D = D_all[0]
    gate_coph = pd.read_csv(BASE / "gate_grid.csv")
    gate_raw = pd.read_csv(BASE / "gate_raw.csv")
    desc = pd.read_csv(GRID / "descriptive.csv")

    detection(D_all)
    characterization(D, BASE / "characterization.csv")
    selection(D, desc, gate_coph, gate_raw)
    scale_units(D)
    print(f"\n[w0c-tau] -> {BASE}", flush=True)


if __name__ == "__main__":
    main()
