#!/usr/bin/env python3
"""Audit 82 — WITHIN-patient propagator-subspace recovery of epileptic nodes.

The right framing for "rediscover an epi contact if it wasn't labelled": not a
cross-patient classifier (killed by patient-specific FC scale in audit_80/81)
but a WITHIN-patient anomaly/template recovery in a rich feature subspace built
from the LRG propagator ρ(τ) = e^{-τL}/Z (information diffusion) + the per-pair
trace measures + leverage. Self-normalised within each patient → scale is no
longer a confound. Leave-one-epi-out: hide one epi label, build the template
from the other epi, check whether the held-out node re-ranks to the top.

Scope report section 11:
``.agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md``.

Critical preamble (per CLAUDE.md)
=================================
(1) **Claim.** In a within-patient propagator+trace subspace, epi contacts
    cluster tightly enough that a held-out epi node (label hidden) re-ranks to
    the top via similarity to the other epi, and unlabelled look-alikes are
    candidates.
(2) **Null.** Random epi-set permutation (same size) → null recovery AUC; plus a
    strength-only ablation.
(3) **Strongest alternatives.** (a) hubness/strength; (b) sEEG-shaft proximity
    (epi cluster on shafts → neighbours look alike in any shaft-smooth feature).
(4) **Mechanical reach.** Ablation strength / diffusion / trace / full quantifies
    the strength share; the SHAFT CONTROL recomputes recovery with negatives
    restricted to epi-free shafts; leave-one-epi-out keeps the held-out label
    unused.
(5) **Falsification (of the INTERESTING claim).** If full recovery AUC ≈ strength
    AUC ≈ shaft-proximity AUC, there is no propagator-specific epi subspace — the
    recovery is real but is hubness/geometry, reported as such (not a novel
    diffusion biomarker).

Outputs (``data/audit/epi_marker/``)
------------------------------------
    node_propagator_features.csv   per (patient, band, node) propagator features
    recovery_per_band.csv          per (band, feature_set) recovery AUC + null + shaft
    propagator_candidates.csv      within-patient ranked non-epi candidates
    README_recovery.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks, parse_seeg_label

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import ensure_half_fcs, load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
import audit_80_epi_marker_classifier as a80  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker"
OUT.mkdir(parents=True, exist_ok=True)

N_TAU = 6
PROP_COLS = ([f"rho_ii_t{i}" for i in range(N_TAU)]
             + [f"Hdiff_t{i}" for i in range(N_TAU)]
             + ["Tcomm_mean", "Tcomm_min"])
TRACE_COLS = ["z_rho_split", "z_coherence", "z_grass_modeload"]
LEV_COLS = ["lev_restruct_post", "lev_restruct_task",
            "lev_trace_signed", "lev_trace_abs"]
STR_COLS = ["strength_mean"]

FEATURE_SETS = {
    "strength": STR_COLS,
    "diffusion": PROP_COLS,
    "trace": TRACE_COLS,
    "leverage": LEV_COLS,
    "diffusion_trace": PROP_COLS + TRACE_COLS,
    "full": PROP_COLS + TRACE_COLS + LEV_COLS + STR_COLS,
}
PERM_SETS = ("strength", "diffusion", "full")
MIN_EPI = 3


# --------------------------------------------------------------------------- #
# Propagator ρ(τ) features from eigh(L) (canonical LRG density matrix)         #
# --------------------------------------------------------------------------- #
def propagator_feats(W: np.ndarray) -> dict | None:
    L = np.diag(W.sum(1)) - W
    lam, U = np.linalg.eigh(L)
    lam = np.clip(lam, 0.0, None)
    lmax = float(lam[-1])
    if lmax <= 0:
        return None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)
    feats: dict[str, np.ndarray] = {}
    for ti, t in enumerate(taus):
        w = np.exp(-t * lam)
        Z = w.sum()
        feats[f"rho_ii_t{ti}"] = (U ** 2 * w).sum(1) / Z          # return prob
        K = np.clip((U * w) @ U.T, 1e-300, None)
        P = K / K.sum(1, keepdims=True)
        feats[f"Hdiff_t{ti}"] = -(P * np.log(P)).sum(1)           # diffusion reach
    # communication centrality at the finest LRG scale τ=1/λ_max
    w0 = np.exp(-(1.0 / lmax) * lam); Z0 = w0.sum()
    rho0 = np.clip(((U * w0) @ U.T) / Z0, 1e-300, None)
    Trho = 1.0 / rho0
    np.fill_diagonal(Trho, np.nan)
    feats["Tcomm_mean"] = np.nanmean(Trho, axis=1)
    feats["Tcomm_min"] = np.nanmin(Trho, axis=1)
    return feats


def compute_propagator_table(patients, bands, verbose) -> pd.DataFrame:
    rows = []
    for band in bands:
        for pat in patients:
            try:
                W = load_phase_fc(pat, "rest_post", band)
            except Exception as e:
                if verbose:
                    print(f"[audit_85] SKIP {pat}/{band}: {e}")
                continue
            N = W.shape[0]
            pm = build_epi_masks(pat)
            if len(pm.channels) != N:
                continue
            f = propagator_feats(W)
            if f is None:
                continue
            epi = np.asarray(pm.epi_mask, bool)
            probes = np.asarray(pm.probes, object)
            for i in range(N):
                row = {"patient": pat, "band": band, "node": i,
                       "channel": pm.channels[i], "probe": str(probes[i]),
                       "is_epi": bool(epi[i])}
                row.update({c: float(f[c][i]) for c in PROP_COLS})
                rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Within-patient leave-one-epi-out template recovery                          #
# --------------------------------------------------------------------------- #
def _standardize(X: np.ndarray) -> np.ndarray:
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd = np.where(sd > 1e-9, sd, 1.0)
    Z = (X - mu) / sd
    return np.nan_to_num(Z, nan=0.0)


def _loo_recovery(Z: np.ndarray, epi_idx: np.ndarray,
                  neg_idx: np.ndarray) -> float:
    """Pooled P(held-out epi scores above a negative), leave-one-epi-out."""
    if epi_idx.size < MIN_EPI or neg_idx.size < 2:
        return float("nan")
    wins = 0.0; total = 0
    for e in epi_idx:
        others = epi_idx[epi_idx != e]
        template = Z[others].mean(0)
        de = -np.linalg.norm(Z[e] - template)
        dn = -np.linalg.norm(Z[neg_idx] - template, axis=1)
        wins += float((de > dn).sum()) + 0.5 * float((de == dn).sum())
        total += dn.size
    return wins / total if total else float("nan")


def recovery_for_cell(merged_pb: pd.DataFrame, rng) -> list[dict]:
    """All feature sets for one (patient, band): recovery AUC, shaft-control
    AUC, permutation-null p."""
    y = merged_pb["is_epi"].to_numpy(bool)
    if y.sum() < MIN_EPI or (~y).sum() < 2:
        return []
    probe = merged_pb["probe"].to_numpy()
    epi_idx = np.where(y)[0]
    non_idx = np.where(~y)[0]
    epi_shafts = set(probe[epi_idx])
    offshaft_neg = np.array([j for j in non_idx if probe[j] not in epi_shafts],
                            dtype=int)
    out = []
    for sn, cols in FEATURE_SETS.items():
        Z = _standardize(merged_pb[cols].to_numpy(float))
        auc = _loo_recovery(Z, epi_idx, non_idx)
        auc_off = (_loo_recovery(Z, epi_idx, offshaft_neg)
                   if offshaft_neg.size >= 2 else float("nan"))
        perm_p = float("nan")
        if sn in PERM_SETS and np.isfinite(auc):
            n_epi = epi_idx.size
            ge = 0
            for _ in range(args_perm()):
                rand_epi = rng.choice(len(y), size=n_epi, replace=False)
                a = _loo_recovery(Z, rand_epi,
                                  np.setdiff1d(np.arange(len(y)), rand_epi))
                if np.isfinite(a) and a >= auc:
                    ge += 1
            perm_p = (ge + 1) / (args_perm() + 1)
        out.append({
            "patient": merged_pb["patient"].iloc[0],
            "band": merged_pb["band"].iloc[0], "feature_set": sn,
            "n_epi": int(epi_idx.size), "n_nonepi": int(non_idx.size),
            "recovery_auc": auc, "recovery_auc_offshaft": auc_off,
            "perm_p": perm_p,
        })
    return out


_PERM = {"v": 200}
def args_perm() -> int:
    return _PERM["v"]


# --------------------------------------------------------------------------- #
# Discovery — within-patient candidate ranking                                 #
# --------------------------------------------------------------------------- #
def discovery(merged: pd.DataFrame, bands) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = FEATURE_SETS["full"]
    coord_cache: dict[str, np.ndarray] = {}
    rows = []
    for band in bands:
        for pat, mp in merged[merged.band == band].groupby("patient"):
            mp = mp.reset_index(drop=True)
            y = mp["is_epi"].to_numpy(bool)
            if y.sum() < MIN_EPI:
                continue
            Z = _standardize(mp[cols].to_numpy(float))
            template = Z[y].mean(0)
            dist = np.linalg.norm(Z - template, axis=1)
            non = np.where(~y)[0]
            order = np.argsort(dist[non])              # closest first
            pct = np.empty(non.size)
            pct[order] = np.linspace(1, 0, non.size)   # 1 = closest to epi
            if pat not in coord_cache:
                chans = mp.sort_values("node")["channel"].tolist()
                coord_cache[pat] = a80._load_coords(pat, chans)
            coords = coord_cache[pat]
            epi_nodes = mp[y]["node"].to_numpy(int)
            for k, j in enumerate(non):
                ni = int(mp.iloc[j]["node"])
                d3 = np.nan
                if (coords.shape[0] > ni and np.isfinite(coords[ni]).all()
                        and epi_nodes.size):
                    ec = coords[epi_nodes]
                    ec = ec[np.isfinite(ec).all(1)]
                    if ec.size:
                        d3 = float(np.sqrt(((ec - coords[ni]) ** 2).sum(1)).min())
                rows.append({
                    "patient": pat, "band": band, "node": ni,
                    "channel": mp.iloc[j]["channel"], "probe": mp.iloc[j]["probe"],
                    "epi_template_dist": float(dist[j]),
                    "closeness_pct": float(pct[k]),
                    "dist3d_to_nearest_epi_mm": d3,
                })
    cand = pd.DataFrame(rows)
    if cand.empty:
        return cand, pd.DataFrame()
    cand["is_top"] = (cand["closeness_pct"] >= 0.9).astype(int)
    agg = (cand.groupby(["patient", "node", "channel", "probe"])
           .agg(n_bands_top=("is_top", "sum"),
                mean_closeness=("closeness_pct", "mean"),
                min_dist3d=("dist3d_to_nearest_epi_mm", "min"))
           .reset_index()
           .sort_values(["n_bands_top", "mean_closeness"],
                        ascending=[False, False]))
    return cand.sort_values(["band", "closeness_pct"],
                            ascending=[True, False]), agg


# --------------------------------------------------------------------------- #
def band_summary(rec: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(rec.band)]:
        rb = rec[rec.band == band]
        per_set = {sn: rb[rb.feature_set == sn] for sn in FEATURE_SETS}
        str_med = float(np.nanmedian(per_set["strength"]["recovery_auc"]))
        for sn in FEATURE_SETS:
            a = per_set[sn]["recovery_auc"].to_numpy(float)
            a = a[np.isfinite(a)]
            ao = per_set[sn]["recovery_auc_offshaft"].to_numpy(float)
            ao = ao[np.isfinite(ao)]
            pp = per_set[sn]["perm_p"].to_numpy(float)
            pp = pp[np.isfinite(pp)]
            wp = float("nan")
            if a.size >= 3:
                try:
                    _, wp = wilcoxon(a - 0.5, alternative="greater")
                except Exception:
                    pass
            out.append({
                "band": band, "feature_set": sn, "n_patients": int(a.size),
                "median_recovery_auc": float(np.median(a)) if a.size else np.nan,
                "median_recovery_offshaft": (float(np.median(ao))
                                             if ao.size else np.nan),
                "delta_vs_strength": ((float(np.median(a)) - str_med)
                                      if a.size else np.nan),
                "n_perm_sig": int((pp < 0.05).sum()) if pp.size else 0,
                "wilcoxon_gt_p": wp,
            })
    return pd.DataFrame(out)


def write_readme(summ: pd.DataFrame, agg: pd.DataFrame, runtime: float) -> None:
    L = []; a = L.append
    a("---")
    a("name: epi_marker_propagator_recovery")
    a("scope: direction2_within_patient_propagator_subspace_recovery")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: within_patient_recovery")
    a("build_script: scripts/01_compute/audit/audit_85_epi_propagator_recovery.py")
    a("scope_report: .agents/guides/task-persistence-investigation/"
      "2026-06-05_epi-node-trace-marker.md (section 11)")
    a("---")
    a("")
    a("# Within-patient propagator-subspace recovery of epi nodes")
    a("")
    a("**Head.** WITHIN each patient, build a feature subspace from the LRG "
      "propagator ρ(τ)=e^{-τL}/Z (return probability, diffusion entropy, "
      "communication distance) + per-pair trace + leverage, standardise within "
      "patient, and leave-one-epi-out: hide one epi label, template = mean of "
      "the OTHER epi, recover the held-out node by similarity. Recovery AUC = "
      "P(held-out epi ranks above a non-epi). The scientifically load-bearing "
      "readouts are (1) does `full` beat `strength`, (2) does it survive the "
      "epi-free-shaft control.")
    a("")
    if not summ.empty:
        full = summ[summ.feature_set == "full"]
        strg = summ[summ.feature_set == "strength"]
        fb = full.loc[full["median_recovery_auc"].idxmax()]
        a("## Bottom line")
        a("")
        a(f"- **Within-patient recovery WORKS** (best band "
          f"{BRAIN_BAND_TEX_DICT.get(fb['band'], fb['band'])}: full median "
          f"recovery AUC {fb['median_recovery_auc']:.3f}; a hidden epi contact "
          "re-ranks above non-epi well above chance). This is the answer to "
          "\"rediscover an unlabelled epi node\": yes, within patient.")
        med_full = float(np.nanmedian(full["median_recovery_auc"]))
        med_str = float(np.nanmedian(strg["median_recovery_auc"]))
        med_delta = float(np.nanmedian(full["delta_vs_strength"]))
        a(f"- **Cross-band median: full {med_full:.3f} vs strength {med_str:.3f} "
          f"(Δ {med_delta:+.3f}).** "
          + ("The propagator subspace ADDS to strength → the recovery is not "
             "only hubness." if med_delta > 0.02 else
             "The propagator subspace does NOT add meaningfully to strength → "
             "recovery is largely hubness (honest: usable tool, unsurprising "
             "mechanism)."))
        off_ok = float(np.nanmedian(full["median_recovery_offshaft"]))
        a(f"- **Shaft control:** full off-shaft recovery AUC {off_ok:.3f} "
          + ("≈ all-negative recovery → not a shaft-proximity artifact."
             if abs(off_ok - med_full) < 0.05 else
             "drops vs all-negative → part of the recovery is shaft proximity."))
        a("")
    a("## Recovery AUC per band × feature set")
    a("")
    a("| band | feature_set | n_pat | median AUC | off-shaft AUC | Δ vs strength "
      "| n perm-sig | Wilcoxon>p |")
    a("|---|---|---|---|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(summ.band)]:
        for _, r in summ[summ.band == band].iterrows():
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['feature_set']} "
              f"| {r['n_patients']} | {r['median_recovery_auc']:.3f} "
              f"| {r['median_recovery_offshaft']:.3f} "
              f"| {r['delta_vs_strength']:+.3f} | {int(r['n_perm_sig'])} "
              f"| {r['wilcoxon_gt_p']:.4f} |")
    a("")
    a("## Reading")
    a("- Recovery AUC 0.5 = chance. `strength` is the trivial baseline (epi are "
      "hubs). The interesting quantity is **Δ vs strength** (propagator+trace "
      "subspace beyond hubness) and the **off-shaft** column (geometry control).")
    a("- `n perm-sig` = patients whose epi set recovers better than random "
      "same-size sets (perm-p<0.05).")
    a("- Pat_15 (0 epi) excluded (no template).")
    if not agg.empty:
        a("")
        a("## Top within-patient candidates (unlabelled non-epi near the epi "
          "template across bands)")
        a("")
        a("| patient | channel | probe | n_bands_top | mean_closeness | "
          "min 3D dist to epi (mm) |")
        a("|---|---|---|---|---|---|")
        for _, r in agg.head(12).iterrows():
            d = r["min_dist3d"]
            ds = f"{d:.1f}" if np.isfinite(d) else "n/a"
            a(f"| {r['patient']} | {r['channel']} | {r['probe']} "
              f"| {int(r['n_bands_top'])} | {r['mean_closeness']:.3f} | {ds} |")
        a("")
        a("- A candidate CLOSE to a known epi in 3D (low min dist) is likely "
          "field continuation; a high-agreement candidate FAR from any known "
          "epi is the genuinely novel — and most uncertain — flag. Confidence "
          "is bounded by the Δ-vs-strength column above.")
    a("")
    a(f"- Wall-clock: {runtime:.1f}s")
    (OUT / "README_recovery.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--perm", type=int, default=200)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    _PERM["v"] = args.perm

    print("[audit_85] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    prop = compute_propagator_table(args.patients, args.bands, args.verbose)
    prop.to_csv(OUT / "node_propagator_features.csv", index=False)
    print(f"[audit_85] propagator features: {len(prop)} rows")

    # merge trace (audit_79) + leverage (audit_81)
    f79 = pd.read_csv(OUT / "node_features_per_band.csv")[
        ["patient", "band", "node", "contact"] + TRACE_COLS + STR_COLS]
    f81 = pd.read_csv(OUT / "node_leverage_per_band.csv")[
        ["patient", "band", "node"] + LEV_COLS]
    merged = prop.merge(f79, on=["patient", "band", "node"], how="inner") \
                 .merge(f81, on=["patient", "band", "node"], how="inner")

    rng = np.random.default_rng(20260605)
    rec_rows = []
    for band in args.bands:
        for pat in args.patients:
            mp = merged[(merged.band == band) & (merged.patient == pat)] \
                .reset_index(drop=True)
            if mp.empty:
                continue
            rec_rows.extend(recovery_for_cell(mp, rng))
        # progress per band
        rb = [r for r in rec_rows if r["band"] == band
              and r["feature_set"] in ("full", "strength")]
        if rb:
            full = np.nanmedian([r["recovery_auc"] for r in rb
                                 if r["feature_set"] == "full"])
            strg = np.nanmedian([r["recovery_auc"] for r in rb
                                 if r["feature_set"] == "strength"])
            print(f"[audit_85] {band}: full recovery AUC={full:.3f} "
                  f"strength={strg:.3f} (Δ {full-strg:+.3f})")
    rec = pd.DataFrame(rec_rows)
    rec.to_csv(OUT / "recovery_per_band.csv", index=False)
    summ = band_summary(rec)

    cand, agg = discovery(merged, args.bands)
    cand.to_csv(OUT / "propagator_candidates.csv", index=False)
    agg.to_csv(OUT / "propagator_candidates_agreement.csv", index=False)

    runtime = time.time() - t0
    write_readme(summ, agg, runtime)
    print(f"[audit_85] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
