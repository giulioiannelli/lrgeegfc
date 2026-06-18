#!/usr/bin/env python3
"""Audit 88 — signed (Kunegis) Laplacian eigenmodes as an epi-node marker.

Uses the SIGN of ImCoh (the lead/lag information that imcoh_abs throws away). The
signed band-averaged ImCoh ``A`` is exactly antisymmetric (``A = -A^T``; verified
``max|A + A^T| = 0``). To analyse it with a REAL, symmetric, positive-semidefinite
signed Laplacian we keep the mixed signs by mirroring the upper triangle into a
SYMMETRIC signed adjacency ``B`` (``B_ij = B_ji = A_ij`` for ``i<j``); ~half of
``B``'s entries stay negative. Then

    L_s = D_bar - B,    D_bar_ii = sum_j |B_ij|        (Kunegis signed Laplacian)

is real-symmetric and PSD, because
    x^T L_s x = sum_{i<j} |B_ij| (x_i - sign(B_ij) x_j)^2 >= 0.
``eigh`` therefore returns real eigenvalues and REAL eigenvectors whose squared
entries localise each node in the signed (frustration) modes. Question: do
epileptic contacts localise in the low / high signed-Laplacian modes BEYOND node
strength (hubness)?

NB on gauge. The sign of ``B`` rides on node ordering (which member of a pair is
``i<j``). Any *physically* meaningful orientation (put ``+`` on the leader) makes
every entry positive and collapses ``B`` to ``|A|`` -- the plain magnitude
Laplacian we already analyse -- so mixed signs require an arbitrary convention. We
therefore add a CHANNEL-RELABEL robustness check: recompute recovery under random
node relabelings; a real signal survives relabeling, labeling noise does not.
(This replaces the magnetic Laplacian ``L_H = D_bar - iA``, removed at user
request.)

Critical preamble
=================
(1) Claim: epi contacts localise in the signed-Laplacian eigensubspace beyond
    node strength.
(2) Null: chance recovery = 0.5 (hide-and-seek AUC, hidden epi vs healthy); per-
    band cohort Wilcoxon vs 0.5; AND the paired test signed-subspace vs hubness.
(3) Strongest alternative: it is hubness (the |ImCoh| degree); or it is gauge
    noise (the sign convention is arbitrary).
(4) Reach: hubness is the head-to-head baseline (all-contacts recovery, fair
    Fisher scorer, audit_87 design); the relabel check addresses the gauge
    alternative. Depth / within-shaft is NOT used -- retracted: depth is not a
    deployable marker and only "recovers" hidden epi via on-shaft contiguity.
(5) Falsification: if signed-subspace recovery <= hubness, or collapses under
    relabeling, the ImCoh sign carries no epi marker beyond hubness -- report it.

Outputs (``data/audit/epi_signed_laplacian/``)
    node_signed_laplacian_features.csv      per (patient, band, node)
    signed_laplacian_behavior_per_band.csv  epi-vs-non AUC + recovery vs hubness + gauge
    README.md
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
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.workflow.fc import load_fc_matrix

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_87_epi_masked_recovery import _recovery_auc, _standardize  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_signed_laplacian"
OUT.mkdir(parents=True, exist_ok=True)

ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08",
            "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
K_LOW = 5
MIN_EPI = 6
R_GAUGE = 8

SL_FEATURES = ["signed_strength", "abs_signed_strength",
               "sl_lowmode_loc", "sl_himode_loc"]
SUBSPACE = ["sl_lowmode_loc", "sl_himode_loc"]          # the eigenmode subspace
FEATURE_SETS = {
    "signed_subspace": SUBSPACE,
    "signed_all": SL_FEATURES,
    "hubness": ["strength_abs"],
}


def _mirror_upper(A: np.ndarray) -> np.ndarray:
    """Symmetric signed adjacency: mirror the upper triangle WITHOUT sign flip."""
    iu = np.triu_indices_from(A, 1)
    B = np.zeros_like(A)
    B[iu] = A[iu]
    return B + B.T


def signed_laplacian_feats(A: np.ndarray):
    """Per-node features from the real symmetric PSD signed Laplacian of A."""
    if not np.allclose(A, -A.T, atol=1e-6):
        A = 0.5 * (A - A.T)                              # enforce antisymmetry
    B = _mirror_upper(A)
    absdeg = np.abs(B).sum(1)
    if not np.isfinite(absdeg).all() or absdeg.max() <= 0:
        return None, None
    Ls = np.diag(absdeg) - B                             # Kunegis signed Laplacian
    w, U = np.linalg.eigh(Ls)                            # symmetric -> real w, real U
    N = A.shape[0]
    K = min(K_LOW, N - 1)
    sig = B.sum(1)
    feats = {
        "signed_strength": sig,                          # net signed degree (gauge-dep.)
        "abs_signed_strength": np.abs(sig),
        "sl_lowmode_loc": (U[:, :K] ** 2).sum(1),        # low signed modes
        "sl_himode_loc": (U[:, N - K:] ** 2).sum(1),     # high (frustration) modes
        "strength_abs": absdeg,                          # reference = hubness (|ImCoh| deg)
    }
    return feats, w


def feats_in_frame(A: np.ndarray, perm: np.ndarray):
    """SL features recomputed in a relabeled frame, mapped back to node order."""
    Ap = A[np.ix_(perm, perm)]
    f, _ = signed_laplacian_feats(Ap)
    if f is None:
        return None
    inv = np.argsort(perm)
    return {k: v[inv] for k, v in f.items()}


def _auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan")
    return float(mannwhitneyu(a, b, alternative="two-sided").statistic
                 / (a.size * b.size))


def _med_p(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if v.size < 3:
        return (float(np.median(v)) if v.size else np.nan), np.nan
    try:
        _, p = wilcoxon(v - 0.5, alternative="greater")
    except Exception:
        p = np.nan
    return float(np.median(v)), p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=ALL_BANDS)
    ap.add_argument("--patients", nargs="+", default=PATIENTS)
    args = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(20260606)
    gauge_rng = np.random.default_rng(20260607)

    rows, psd_min = [], []
    for band in args.bands:
        for pat in args.patients:
            A = load_fc_matrix(pat, "rest_post", band, "imcoh")
            if A is None:
                continue
            A = np.asarray(A, float)
            N = A.shape[0]
            pm = build_epi_masks(pat)
            if len(pm.channels) != N:
                continue
            f, w = signed_laplacian_feats(A)
            if f is None:
                continue
            psd_min.append(float(w[0]))
            epi = np.asarray(pm.epi_mask, bool)
            probes = np.asarray(pm.probes, object)
            for i in range(N):
                rows.append({"patient": pat, "band": band, "node": i,
                             "is_epi": bool(epi[i]), "probe": str(probes[i]),
                             **{k: float(f[k][i]) for k in f}})
    feat = pd.DataFrame(rows)
    feat.to_csv(OUT / "node_signed_laplacian_features.csv", index=False)
    if psd_min:
        print(f"[audit_88] signed Laplacian min eigenvalue across cells: "
              f"{min(psd_min):.3e} (>= 0 confirms PSD)")

    summ = []
    for band in args.bands:
        fb = feat[feat.band == band]
        row = {"band": band}
        # per-feature epi-vs-non AUC (oriented) + strength confound
        for col in SL_FEATURES + ["strength_abs"]:
            devs, corrs = [], []
            for _, fp in fb.groupby("patient"):
                if fp.is_epi.sum() < 2 or (~fp.is_epi).sum() < 2:
                    continue
                a = _auc(fp[fp.is_epi][col].values, fp[~fp.is_epi][col].values)
                if np.isfinite(a):
                    devs.append(max(a, 1 - a))
                sc = spearmanr(fp[col].values, fp["strength_abs"].values).correlation
                if np.isfinite(sc):
                    corrs.append(sc)
            devs = np.array(devs)
            wp = np.nan
            if devs.size >= 3:
                try:
                    _, wp = wilcoxon(devs - 0.5, alternative="greater")
                except Exception:
                    pass
            row[f"{col}_auc_oriented"] = float(np.median(devs)) if devs.size else np.nan
            row[f"{col}_auc_p"] = wp
            row[f"{col}_strength_corr"] = (float(np.median(np.abs(corrs)))
                                           if corrs else np.nan)

        # all-contacts hide-and-seek recovery per feature set + paired vs hubness + gauge
        rec = {k: [] for k in FEATURE_SETS}
        paired = {"signed_subspace": [], "signed_all": []}
        rec_gauge = {"signed_subspace": [], "signed_all": []}
        rec_gauge_std = {"signed_subspace": [], "signed_all": []}
        for pat, fp in fb.groupby("patient"):
            fp = fp.reset_index(drop=True)
            y = fp.is_epi.to_numpy(bool)
            if y.sum() < MIN_EPI:
                continue
            epi_idx = np.where(y)[0]
            healthy = np.where(~y)[0]
            here = {}
            for sn, cols in FEATURE_SETS.items():
                Z = _standardize(fp[cols].to_numpy(float))
                r = _recovery_auc(Z, epi_idx, healthy, rng)
                rec[sn].append(r); here[sn] = r
            for sn in ("signed_subspace", "signed_all"):
                if np.isfinite(here[sn]) and np.isfinite(here["hubness"]):
                    paired[sn].append(here[sn] - here["hubness"])
            # gauge: recompute SL feats under random relabelings of this patient's A
            A = np.asarray(load_fc_matrix(pat, "rest_post", band, "imcoh"), float)
            gl = {"signed_subspace": [], "signed_all": []}
            for _ in range(R_GAUGE):
                perm = gauge_rng.permutation(A.shape[0])
                ff = feats_in_frame(A, perm)
                if ff is None:
                    continue
                fdf = pd.DataFrame(ff)
                for sn in ("signed_subspace", "signed_all"):
                    Z = _standardize(fdf[FEATURE_SETS[sn]].to_numpy(float))
                    gl[sn].append(_recovery_auc(Z, epi_idx, healthy, gauge_rng))
            for sn in ("signed_subspace", "signed_all"):
                arr = np.array([x for x in gl[sn] if np.isfinite(x)])
                if arr.size:
                    rec_gauge[sn].append(float(np.median(arr)))
                    rec_gauge_std[sn].append(float(np.std(arr)))

        for sn in FEATURE_SETS:
            m, p = _med_p(rec[sn])
            row[f"rec_{sn}_median"] = m
            row[f"rec_{sn}_p_gt_chance"] = p
        for sn in ("signed_subspace", "signed_all"):
            pv = np.array([x for x in paired[sn] if np.isfinite(x)])
            pp = np.nan
            if pv.size >= 3:
                try:
                    _, pp = wilcoxon(pv, alternative="greater")
                except Exception:
                    pass
            row[f"{sn}_minus_hubness_median"] = float(np.median(pv)) if pv.size else np.nan
            row[f"{sn}_beats_hubness_p"] = pp
            row[f"rec_{sn}_gauge_median"] = (float(np.median(rec_gauge[sn]))
                                             if rec_gauge[sn] else np.nan)
            row[f"rec_{sn}_gauge_withinpat_std"] = (float(np.median(rec_gauge_std[sn]))
                                                    if rec_gauge_std[sn] else np.nan)
        summ.append(row)
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "signed_laplacian_behavior_per_band.csv", index=False)

    # README
    L = ["---", "name: epi_signed_laplacian",
         "scope: direction2_signed_imcoh_kunegis_laplacian",
         "era: COHORT_N10 / IMCOH (signed, antisymmetric)",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_88_signed_laplacian.py",
         "scope_report: .agents/guides/task-persistence-investigation/2026-06-05_signed-magnetic-laplacian-epi.md",
         "---", "",
         "# Signed (Kunegis) Laplacian epi features (within-patient, vs hubness)",
         "",
         "**Head.** Real symmetric PSD signed Laplacian `L_s = D_bar - B` of the "
         "SIGNED ImCoh (`B` = upper-triangle-mirrored antisymmetric `A`, mixed "
         "signs kept). Eigenmode-localisation features `sl_lowmode_loc` / "
         f"`sl_himode_loc` (participation in the {K_LOW} lowest / highest signed "
         "modes), plus net signed degree. Head-to-head with node strength "
         "(hubness) on the all-contacts hide-and-seek recovery, with a channel-"
         "relabel gauge robustness check. Depth / within-shaft NOT used.", "",
         "## Per-feature epi-vs-non-epi separation (within-patient, oriented AUC)",
         "",
         "| band | feature | oriented AUC | p>chance | |strength corr| |",
         "|---|---|---|---|---|"]
    for band in [b for b in ALL_BANDS if b in set(summ.band)]:
        r = summ[summ.band == band].iloc[0]
        for col in SL_FEATURES + ["strength_abs"]:
            L.append(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {col} "
                     f"| {r[f'{col}_auc_oriented']:.3f} | {r[f'{col}_auc_p']:.4f} "
                     f"| {r[f'{col}_strength_corr']:.3f} |")
    L += ["", "## All-contacts recovery — signed subspace vs hubness (+ gauge)", "",
          "| band | signed_subspace | p>chance | hubness | Δ(sub−hub) | beats-hub p "
          "| gauge median | gauge within-pat std |",
          "|---|---|---|---|---|---|---|---|"]
    for band in [b for b in ALL_BANDS if b in set(summ.band)]:
        r = summ[summ.band == band].iloc[0]
        L.append(
            f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
            f"| {r['rec_signed_subspace_median']:.3f} "
            f"| {r['rec_signed_subspace_p_gt_chance']:.4f} "
            f"| {r['rec_hubness_median']:.3f} "
            f"| {r['signed_subspace_minus_hubness_median']:+.3f} "
            f"| {r['signed_subspace_beats_hubness_p']:.4f} "
            f"| {r['rec_signed_subspace_gauge_median']:.3f} "
            f"| {r['rec_signed_subspace_gauge_withinpat_std']:.3f} |")
    L += ["", "## Reading",
          "- `oriented AUC` 0.5 = no separation. `|strength corr|` near 1 ⇒ the "
          "feature is the |ImCoh| degree (hubness) in disguise.",
          "- `signed_subspace` recovery = can the signed-Laplacian eigenmode "
          "subspace re-find hidden epi among ALL healthy contacts. The marker is "
          "real beyond hubness iff `signed_subspace > chance` AND `Δ(sub−hub) > 0` "
          "(beats-hub p < 0.05).",
          "- `gauge median` = same recovery under random channel relabelings; "
          "`gauge within-pat std` small ⇒ signal is gauge-robust, large ⇒ the "
          "mixed-sign convention was carrying noise. If the canonical recovery far "
          "exceeds the gauge median, the result is an artifact of the index-order "
          "sign convention, not the ImCoh sign.",
          f"- splits/patient {200}; relabelings/patient {R_GAUGE}; min epi {MIN_EPI}; "
          f"wall-clock {time.time()-t0:.1f}s"]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    print("[audit_88] all-contacts recovery (median over patients):")
    for band in [b for b in ALL_BANDS if b in set(summ.band)]:
        r = summ[summ.band == band].iloc[0]
        print(f"  {band:10s}: signed_subspace={r['rec_signed_subspace_median']:.3f} "
              f"(p{r['rec_signed_subspace_p_gt_chance']:.3f}) vs hubness "
              f"{r['rec_hubness_median']:.3f} | Δ={r['signed_subspace_minus_hubness_median']:+.3f} "
              f"(beats-hub p{r['signed_subspace_beats_hubness_p']:.3f}) | "
              f"gauge med={r['rec_signed_subspace_gauge_median']:.3f} "
              f"std={r['rec_signed_subspace_gauge_withinpat_std']:.3f}")
    print(f"[audit_88] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
