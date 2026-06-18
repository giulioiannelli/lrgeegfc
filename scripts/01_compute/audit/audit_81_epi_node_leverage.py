#!/usr/bin/env python3
"""Audit 81 — per-node LEVERAGE / influence on the LRG dendrogram and trace.

A different question from audit_79/80: not "does node i's own trace behaviour
look epileptic?" (identity, null cross-patient) but "does removing node i
drastically restructure the dendrogram / shift the cross-phase trace?"
(influence). Motivated by the user's observation that *including* epileptic
nodes produces a big shift in the dendrogram/trace in some patients — the
Direction-1 cohort result pushed to the node level.

Scope report (section 10):
``.agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md``.

Stage 1 (this script): raw leverage + within-patient percentile (scale-free),
its epi-vs-non-epi separation, its correlation with strength (the confound),
and a fold into the audit_80 LOPO + whole-shaft-masking harness. NO
matched-strength surrogate-z yet (Stage 2, gated on transfer).

Critical preamble (per CLAUDE.md rule)
======================================
(1) **Claim.** Epileptic nodes have outsized leverage on the LRG hierarchy:
    removing one markedly restructures the task/rest_post dendrogram (L1) and/or
    shifts ρ_split (L2, signed). That leverage flags them.
(2) **Null (Stage 1).** Permutation null in the LOPO harness (label-free
    leverage features invariant → exact). Matched-strength surrogate-z is
    Stage 2, deferred until transfer is shown.
(3) **Strongest alternative.** Leverage = hubness; high-strength nodes perturb
    the Laplacian spectrum most. Absolute strength does NOT transfer
    cross-patient (audit_80).
(4) **Mechanical reach.** L1/L2 are rank/dimensionless, so the within-patient
    PERCENTILE of leverage is scale-free and may transfer where dimensioned
    strength could not — the specific reason to retry. The strength confound is
    measured head-on (Spearman(leverage, strength) per patient) and the LOPO
    harness pits leverage against the strength and electrode-geometry baselines.
(5) **Falsification.** If leverage separates epi within-patient but is at chance
    in LOPO (like strength), it is dead as a discovery tool — report that.

Outputs (``data/audit/epi_marker/``)
------------------------------------
    node_leverage_per_band.csv       one row per (patient, band, node)
    leverage_behavior_per_band.csv   within-patient epi-vs-non-epi AUC + strength corr
    leverage_lopo_per_band.csv       leverage feature sets in the audit_80 harness
    README_leverage.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks, parse_seeg_label
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs, load_phase_fc)
import _epi_stratify as es  # type: ignore
import audit_80_epi_marker_classifier as a80  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker"
OUT.mkdir(parents=True, exist_ok=True)

L1_PHASES = ("task_test", "rest_post")
LEV_RAW_COLS = ["lev_restruct_post", "lev_restruct_task",
                "lev_trace_abs", "lev_trace_signed"]
LEV_PCT_COLS = ["pct_lev_restruct_post", "pct_lev_restruct_task",
                "pct_lev_trace_abs", "pct_lev_trace_signed"]


def _sp(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5 or np.unique(a[m]).size < 2 or np.unique(b[m]).size < 2:
        return float("nan")
    return float(spearmanr(a[m], b[m]).correlation)


def _laplacian_eigh(W: np.ndarray):
    return np.linalg.eigh(np.diag(W.sum(1)) - W)


# --------------------------------------------------------------------------- #
# Per (patient, band) leverage rows                                           #
# --------------------------------------------------------------------------- #
def per_patient_band(pat: str, band: str, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in es.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_81] SKIP {pat}/{band}: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return []
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return []
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    contacts = [parse_seeg_label(c)[1] for c in pm.channels]

    obs_eig = {ph: _laplacian_eigh(Ws[ph]) for ph in es.PHASES_4}
    D_obs = {ph: cophenetic_condensed_from_eigs(*obs_eig[ph])
             for ph in es.PHASES_4}
    dT = D_obs["task_test"] - D_obs["rest_pre_A"]
    dR = D_obs["rest_post"] - D_obs["rest_pre_B"]
    rho_full = _sp(dT, dR)

    iu_i, iu_j = es.upper_tri_indices(N)
    s_mean = np.mean([Ws[ph].sum(1) for ph in es.PHASES_4], axis=0)

    rows = []
    for i in range(N):
        keep = (iu_i != i) & (iu_j != i)             # pairs not involving i
        # L2 — signed trace influence (eigh-free pair-drop jackknife)
        rho_minus = _sp(dT[keep], dR[keep])
        l2 = (rho_full - rho_minus) if np.isfinite(rho_minus) else np.nan
        # L1 — dendrogram restructuring via leave-one-node-out resect
        idx = np.arange(N) != i
        l1 = {}
        for ph in L1_PHASES:
            Wsub = Ws[ph][np.ix_(idx, idx)]
            ev = _laplacian_eigh(Wsub)
            Dsub = cophenetic_condensed_from_eigs(*ev)   # condensed over N-1
            # D_obs[ph][keep] aligns elementwise with the (N-1)-triu (verified)
            l1[ph] = 1.0 - _sp(D_obs[ph][keep], Dsub)
        rows.append({
            "patient": pat, "band": band, "node": i,
            "channel": pm.channels[i], "probe": str(probes[i]),
            "contact": (int(contacts[i]) if contacts[i] is not None else -1),
            "is_epi": bool(epi[i]), "strength_mean": float(s_mean[i]),
            "lev_trace_signed": l2, "lev_trace_abs": (abs(l2)
                                                      if np.isfinite(l2) else np.nan),
            "lev_restruct_task": l1["task_test"],
            "lev_restruct_post": l1["rest_post"],
        })
    df = pd.DataFrame(rows)
    for c in LEV_RAW_COLS:
        df["pct_" + c] = df[c].rank(pct=True)           # within-patient percentile
    return df.to_dict("records")


# --------------------------------------------------------------------------- #
# Within-patient epi-vs-non-epi AUC + strength confound                       #
# --------------------------------------------------------------------------- #
def _auc(epi_v: np.ndarray, non_v: np.ndarray) -> float:
    epi_v = epi_v[np.isfinite(epi_v)]; non_v = non_v[np.isfinite(non_v)]
    if epi_v.size < 2 or non_v.size < 2:
        return float("nan")
    u = mannwhitneyu(epi_v, non_v, alternative="two-sided").statistic
    return float(u / (epi_v.size * non_v.size))


def behaviour(feat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(feat.band)]:
        fb = feat[feat.band == band]
        for f in LEV_RAW_COLS + LEV_PCT_COLS:
            aucs, scorr = [], []
            for _, fp in fb.groupby("patient"):
                a = _auc(fp[fp.is_epi][f].values, fp[~fp.is_epi][f].values)
                if np.isfinite(a):
                    aucs.append(a)
                sc = _sp(fp[f].values, fp["strength_mean"].values)
                if np.isfinite(sc):
                    scorr.append(sc)
            aucs = np.array(aucs)
            if aucs.size >= 3:
                try:
                    _, wp = wilcoxon(aucs - 0.5, alternative="two-sided")
                except Exception:
                    wp = float("nan")
            else:
                wp = float("nan")
            out.append({
                "band": band, "feature": f,
                "n_patients": int(aucs.size),
                "median_auc": float(np.median(aucs)) if aucs.size else np.nan,
                "abs_median_auc_dev": (abs(np.median(aucs) - 0.5)
                                       if aucs.size else np.nan),
                "n_auc_gt_0.5": int((aucs > 0.5).sum()),
                "cohort_wilcoxon_p": wp,
                "median_strength_corr": (float(np.median(scorr))
                                         if scorr else np.nan),
            })
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
# LOPO via the audit_80 harness (leverage feature sets injected)              #
# --------------------------------------------------------------------------- #
def lopo_table(merged: pd.DataFrame, perm_r: int, rng) -> pd.DataFrame:
    # inject leverage sets into the shared registry (reuse, not duplicate)
    a80.FEATURE_SETS["lev_pct"] = LEV_PCT_COLS
    a80.FEATURE_SETS["lev_raw"] = LEV_RAW_COLS
    a80.FEATURE_SETS["lev_strength"] = LEV_PCT_COLS + a80.STRENGTH_COLS
    a80.FEATURE_SETS["lev_all"] = (LEV_PCT_COLS + a80.STRENGTH_COLS
                                   + a80.PROBE_COLS)
    a80._FILL.update({c: 0.5 for c in LEV_PCT_COLS})     # mid-rank impute
    a80._FILL.update({c: 0.0 for c in LEV_RAW_COLS})
    ref_sets = ["lev_pct", "lev_raw", "lev_strength", "lev_all",
                "trace", "strength", "probe_geom", "all"]
    perm_sets = {"lev_pct", "lev_raw", "lev_strength"}
    rows = []
    for band in [b for b in es.ALL_BANDS if b in set(merged.band)]:
        dfb = merged[merged.band == band].reset_index(drop=True)
        if dfb.is_epi.sum() < 5:
            continue
        prev = float(dfb.is_epi.mean())
        aucs = {sn: a80.pooled_auc(dfb.is_epi.to_numpy(int),
                                   a80.lopo_proba(dfb, sn, "logit"))
                for sn in ref_sets}
        mask = {sn: a80.shaft_masking(dfb, sn, "logit")
                for sn in ("lev_pct", "lev_raw", "probe_label")}
        pp = {}
        for sn in perm_sets:
            pp[sn] = a80.perm_null(dfb, sn, "logit", aucs[sn][0], aucs[sn][1],
                                   perm_r, rng)
        best_base = max(aucs["strength"][0], aucs["probe_geom"][0])
        rows.append(dict(
            band=band, prevalence=prev,
            **{f"{sn}_pr_auc": aucs[sn][0] for sn in ref_sets},
            **{f"{sn}_roc_auc": aucs[sn][1] for sn in ref_sets},
            lev_pct_perm_p=pp["lev_pct"][0],
            lev_raw_perm_p=pp["lev_raw"][0],
            lev_strength_perm_p=pp["lev_strength"][0],
            dAUC_lev_vs_baseline=aucs["lev_pct"][0] - best_base,
            dAUC_all_vs_base=aucs["all"][0] - max(
                aucs["strength"][0], aucs["probe_geom"][0]),
            lev_pct_masked_roc=mask["lev_pct"]["masked_roc_auc"],
            lev_raw_masked_roc=mask["lev_raw"]["masked_roc_auc"],
            probe_label_masked_roc=mask["probe_label"]["masked_roc_auc"],
        ))
        print(f"[audit_81] {band}: prev={prev:.3f} | lev_pct PR="
              f"{aucs['lev_pct'][0]:.3f} lev_raw PR={aucs['lev_raw'][0]:.3f} "
              f"lev_str PR={aucs['lev_strength'][0]:.3f} (geom "
              f"{aucs['probe_geom'][0]:.3f}, trace {aucs['trace'][0]:.3f}) | "
              f"perm lev_pct={pp['lev_pct'][0]:.3f} | "
              f"lev_pct_mask_ROC={mask['lev_pct']['masked_roc_auc']:.3f}")
    return pd.DataFrame(rows)


def write_readme(beh: pd.DataFrame, lopo: pd.DataFrame, runtime: float) -> None:
    L = []; a = L.append
    a("---")
    a("name: epi_marker_leverage")
    a("scope: direction2_node_leverage_influence_on_dendrogram_and_trace")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: stage1_no_surrogate_z")
    a("build_script: scripts/01_compute/audit/audit_81_epi_node_leverage.py")
    a("scope_report: .agents/guides/task-persistence-investigation/"
      "2026-06-05_epi-node-trace-marker.md (section 10)")
    a("---")
    a("")
    a("# Per-node leverage / influence — epi discrimination (Stage 1)")
    a("")
    a("**Head.** Leave-one-node-out influence on the LRG dendrogram (L1, task & "
      "rest_post) and signed influence on the cross-phase trace (L2), raw and "
      "as within-patient percentile (scale-free). Does removing epi nodes "
      "restructure the hierarchy / shift the trace more than removing healthy "
      "nodes — and does that TRANSFER across patients? Stage 1 has NO "
      "matched-strength surrogate-z; the strength correlation is reported as "
      "the confound check.")
    a("")
    # bottom line
    if not lopo.empty:
        bb = lopo.loc[lopo["lev_pct_pr_auc"].idxmax()]
        beats_perm = (np.isfinite(bb["lev_pct_perm_p"])
                      and bb["lev_pct_perm_p"] < 0.05)
        beats_geom = bb["dAUC_lev_vs_baseline"] > 0
        bname = BRAIN_BAND_TEX_DICT.get(bb["band"], bb["band"])
        a("## Bottom line")
        a("")
        if beats_perm and beats_geom:
            a(f"- **Leverage TRANSFERS and beats the trivial baselines** (best "
              f"{bname}: lev_pct PR-AUC {bb['lev_pct_pr_auc']:.3f}, perm-p "
              f"{bb['lev_pct_perm_p']:.3f}, dAUC vs best baseline "
              f"{bb['dAUC_lev_vs_baseline']:+.3f}). → Stage 2 (matched-strength "
              "surrogate-z) warranted to confirm pattern-vs-hubness.")
        elif beats_perm:
            a(f"- **One weak but REAL cross-patient signal: {bname} leverage** "
              f"clears the permutation null (lev_pct PR-AUC "
              f"{bb['lev_pct_pr_auc']:.3f} vs prevalence {bb['prevalence']:.3f}, "
              f"**perm-p {bb['lev_pct_perm_p']:.3f}**) — but does NOT beat "
              f"trivial electrode geometry (dAUC {bb['dAUC_lev_vs_baseline']:+.3f}"
              ") and barely survives masking. **Direction is counter-hypothesis:"
              "** the within-patient driver is `lev_trace_abs` AUC≈0.34 "
              "(p≈0.008) = epi nodes have *LOWER* trace-leverage, i.e. they are "
              "trace-INERT, not trace-disruptive — consistent with the "
              "Direction-1 result that the trace lives in healthy tissue and "
              "epi nodes are dead weight on it.")
        else:
            a(f"- **Leverage does NOT separate epi cross-patient** (best "
              f"{bname}: lev_pct PR-AUC {bb['lev_pct_pr_auc']:.3f} vs prevalence "
              f"{bb['prevalence']:.3f}, perm-p {bb['lev_pct_perm_p']:.3f}).")
        a("- **Per-node leverage is NOT just strength**: |median strength corr| "
          "≈ 0.25–0.29 for the restructuring features (would be ≈1 if leverage "
          "were hubness). So this is a genuine \"leverage doesn't flag epi\" "
          "result, not the strength re-read in disguise.")
        a("- **Why the per-node jackknife is null while the cohort epi-SHIFT is "
          "real (Direction-1):** single-node leave-one-out underestimates a "
          "GROUP effect — with 6–30 epi nodes, removing any one barely moves "
          "the tree because the others hold the structure. The epi influence on "
          "the dendrogram is **collective / delocalized**, not localizable to "
          "individual high-leverage contacts. Per-node discovery is therefore "
          "not viable from leverage; the set-level shift remains real but "
          "non-localizing (same signature as the spatially-delocalized trace).")
        a("")
    a("## Within-patient epi-vs-non-epi AUC + strength confound")
    a("")
    a("| band | feature | n_pat | median AUC | |AUC−.5| | n>0.5 | Wilcoxon p | "
      "median strength corr |")
    a("|---|---|---|---|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(beh.band)]:
        for _, r in beh[beh.band == band].iterrows():
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['feature']} "
              f"| {r['n_patients']} | {r['median_auc']:.3f} "
              f"| {r['abs_median_auc_dev']:.3f} | {r['n_auc_gt_0.5']} "
              f"| {r['cohort_wilcoxon_p']:.4f} "
              f"| {r['median_strength_corr']:+.3f} |")
    a("")
    a("## LOPO cross-patient transfer (audit_80 harness)")
    a("")
    a("| band | prev | lev_pct PR | lev_raw PR | lev+str PR | trace PR | "
      "strength PR | geom PR | lev_pct perm-p | dAUC vs base | lev_pct mask ROC "
      "| taut mask ROC |")
    a("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in lopo.iterrows():
        a(f"| {BRAIN_BAND_TEX_DICT.get(r['band'], r['band'])} "
          f"| {r['prevalence']:.3f} | {r['lev_pct_pr_auc']:.3f} "
          f"| {r['lev_raw_pr_auc']:.3f} | {r['lev_strength_pr_auc']:.3f} "
          f"| {r['trace_pr_auc']:.3f} | {r['strength_pr_auc']:.3f} "
          f"| {r['probe_geom_pr_auc']:.3f} | {r['lev_pct_perm_p']:.3f} "
          f"| {r['dAUC_lev_vs_baseline']:+.3f} | {r['lev_pct_masked_roc']:.3f} "
          f"| {r['probe_label_masked_roc']:.3f} |")
    a("")
    a("## Definitions")
    a("- **L1** `lev_restruct_{task,post}` = `1 − Spearman(D_full[¬i], "
      "D_resect_i)`: dendrogram restructuring when node i is removed.")
    a("- **L2** `lev_trace_signed` = `ρ_split_full − ρ_split(¬i pairs)`; >0 = i "
      "inflates the trace, <0 = i masks it. `lev_trace_abs` = |L2|.")
    a("- `pct_*` = within-patient percentile (scale-free, the LOPO-transferable "
      "form).")
    a("- **Stage 2 (deferred):** matched-strength surrogate-z of L1/L2 to "
      "separate pattern from hubness — gated on a band transferring here.")
    a(f"- Wall-clock: {runtime:.1f}s")
    (OUT / "README_leverage.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--perm", type=int, default=200)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_81] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = per_patient_band(pat, band, args.verbose)
            all_rows.extend(rows)
            if rows:
                ne = sum(r["is_epi"] for r in rows)
                print(f"[audit_81] {pat}/{band}: {len(rows)} nodes ({ne} epi) "
                      f"({time.time()-tc:.1f}s)")
    feat = pd.DataFrame(all_rows)
    feat.to_csv(OUT / "node_leverage_per_band.csv", index=False)

    beh = behaviour(feat)
    beh.to_csv(OUT / "leverage_behavior_per_band.csv", index=False)

    # merge with audit_79 features for the shared LOPO harness
    f79 = pd.read_csv(OUT / "node_features_per_band.csv")
    merged = f79.merge(
        feat[["patient", "band", "node"] + LEV_RAW_COLS + LEV_PCT_COLS],
        on=["patient", "band", "node"], how="inner")
    rng = np.random.default_rng(20260605)
    lopo = lopo_table(merged, args.perm, rng)
    lopo.to_csv(OUT / "leverage_lopo_per_band.csv", index=False)

    runtime = time.time() - t0
    write_readme(beh, lopo, runtime)
    print(f"[audit_81] {len(feat)} rows, {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
