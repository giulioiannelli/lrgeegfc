#!/usr/bin/env python3
"""Audit 90 — C5: spatial-matched random-subset null for the epi diffusion community.

This is the ONE control left before the audit_89 finding (epileptic contacts form
a strength-independent diffusion community, ρ_EE elevated) is airtight. audit_89
already closed strength (C3, strength-matched random subsets on the real graph).
The remaining alternative: epi contacts sit in ONE compact anatomical region, and
spatial clustering — not an epi-specific diffusion property — inflates ρ_EE.

C5 draws random node subsets matched to the epi set on BOTH:
  - node strength (strength-quintile profile preserved EXACTLY, audit_89 C3 idiom),
  - spatial spread (radius of gyration in implant (x,y,z) matched to the epi set).
The matching is constructive: start from a strength-quintile-matched subset, then
do WITHIN-quintile swaps that drive its radius of gyration toward the epi set's —
so the null fixes strength exactly and matches spatial compactness, and any
remaining ρ_EE elevation is neither hubness nor spatial clustering.

Also emits per-patient epi GEOMETRY diagnostics (V2): n_epi, radius of gyration,
n shafts, hemisphere span — to explain the audit_89 non-responders (Pat_10, Pat_07).

Critical preamble
=================
(1) Claim: epi–epi diffusion (ρ_EE) is elevated beyond what an equal-size group of
    equal-strength contacts with the SAME spatial compactness would show — i.e. the
    diffusion community is not an artifact of epi sitting in one region.
(2) Null: strength-quintile-preserving, radius-of-gyration-matched random subsets on
    the REAL graph (R_SUB draws); z and upper-tail p of the epi off-diagonal mean ρ.
(3) Strongest alternative: epi are spatially clustered → proximity mediates ρ_EE.
(4) Mechanical reach: matching radius of gyration forces null groups to be as
    spatially compact as epi; preserving the strength-quintile profile fixes the
    strength channel. Survival ⇒ not compactness, not hubness. CANNOT control a
    NON-compact but still privileged spatial pattern (cross-shaft C4 covers the
    on-shaft part). ImCoh already suppresses zero-lag volume conduction, so
    "nearby ⇒ high FC" is not assumed.
(5) Falsification: if z_spat_EE ≈ 0 (obs within the spatial null spread) in every
    band, the co-diffusion is spatial clustering, not an epi diffusion community —
    and the audit_89 positives reduce to anatomy.

Outputs (``data/audit/epi_propagator_blocks/``)
    spatial_null_per_patient.csv   per (patient, band, tau): z_spat_EE + geometry
    spatial_null_cohort.csv        per (band, tau)
    epi_geometry_per_patient.csv   per (patient, band): V2 outlier diagnostics
    README_spatial_null.md
"""
from __future__ import annotations

import argparse
import sys
import time
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
from audit_89_epi_propagator_blocks import subset_offdiag_means  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_propagator_blocks"
OUT.mkdir(parents=True, exist_ok=True)

R_SUB = 400           # spatial+strength-matched random subsets
N_SBINS = 5           # strength quantile bins (same as audit_89 C3)
N_SWAPS = 600         # within-quintile swaps to match radius of gyration
RG_TOL = 0.15         # report match quality vs this fractional tolerance


def radius_of_gyration(coords: np.ndarray, idx: np.ndarray) -> float:
    """RMS distance of a node subset to its own centroid (spatial spread)."""
    c = coords[idx]
    return float(np.sqrt(((c - c.mean(0)) ** 2).sum(1).mean()))


def spatial_strength_matched_subsets(strength, coords, epi_mask, n_draws, rng,
                                     n_swaps=N_SWAPS):
    """Random subsets matched to the epi set on strength-quintile profile (exact,
    by construction) AND spatial radius of gyration (driven by within-quintile
    swaps). Only coord-valid nodes are eligible. Returns (subsets, target_Rg,
    achieved_Rg_list)."""
    n = strength.size
    ranks = np.argsort(np.argsort(strength))
    bins = np.minimum((ranks * N_SBINS) // n, N_SBINS - 1)
    coord_ok = np.isfinite(coords).all(1)
    epi_idx = np.where(epi_mask)[0]
    epi_ok = epi_idx[coord_ok[epi_idx]]
    if epi_ok.size < MIN_EPI:
        return [], np.nan, []
    target = radius_of_gyration(coords, epi_ok)
    pools = {b: np.where((bins == b) & coord_ok)[0] for b in range(N_SBINS)}
    epi_count = Counter(bins[epi_ok].tolist())
    n_sub = int(sum(epi_count.values()))
    subsets, achieved = [], []
    for _ in range(n_draws):
        chosen, binof, ok = set(), {}, True
        for b, k in epi_count.items():
            pool = pools[b]
            if pool.size < k:
                ok = False
                break
            pick = rng.choice(pool, size=k, replace=False)
            for p in pick.tolist():
                chosen.add(p)
                binof[p] = b
        if not ok or len(chosen) != n_sub:
            continue
        cur = radius_of_gyration(coords, np.fromiter(chosen, int))
        members = list(chosen)
        for _ in range(n_swaps):
            m = int(rng.choice(members))
            b = binof[m]
            c = int(rng.choice(pools[b]))
            if c in chosen:
                continue
            trial = chosen - {m} | {c}
            r_new = radius_of_gyration(coords, np.fromiter(trial, int))
            if abs(r_new - target) < abs(cur - target):
                chosen, cur = trial, r_new
                del binof[m]
                binof[c] = b
                members = list(chosen)
        subsets.append(np.fromiter(chosen, int))
        achieved.append(cur)
    return subsets, target, achieved


def epi_geometry(pat, epi_mask, probe, coords, regions):
    """V2 diagnostic: epi-set geometry that may explain (non-)response."""
    epi_idx = np.where(epi_mask)[0]
    coord_ok = np.isfinite(coords).all(1)
    epi_ok = epi_idx[coord_ok[epi_idx]]
    g = {
        "patient": pat,
        "n_epi": int(epi_mask.sum()),
        "n_epi_coord": int(epi_ok.size),
        "n_epi_shafts": int(len(set(probe[epi_idx]))),
        "rg_epi": radius_of_gyration(coords, epi_ok) if epi_ok.size >= 2 else np.nan,
        "rg_all": radius_of_gyration(coords, np.where(coord_ok)[0]),
    }
    hemis = regions["hemisphere"].to_numpy()[epi_idx]
    hcount = Counter(h for h in hemis if h in ("L", "R"))
    g["epi_hemi"] = "/".join(f"{k}{v}" for k, v in sorted(hcount.items())) or "?"
    sysn = regions["system"].to_numpy()[epi_idx]
    g["epi_n_systems"] = int(len({s for s in sysn if s and s != "unknown"}))
    g["rg_epi_frac"] = g["rg_epi"] / g["rg_all"] if g["rg_all"] else np.nan
    return g


def per_patient_band(pat, band, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_90] SKIP {pat}/{band}: {e}")
        return None, None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None, None
    epi = np.asarray(pm.epi_mask, bool)
    probe = np.asarray(pm.probes, object)
    if epi.sum() < MIN_EPI or (~epi).sum() < 2:
        return None, None
    reg = load_channel_regions(pat)
    if len(reg) != N:
        if verbose:
            print(f"[audit_90] SKIP {pat}/{band}: regions N {len(reg)} != FC {N}")
        return None, None
    coords = reg[["x", "y", "z"]].to_numpy(float)

    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    lmax = float(lam_o[-1])
    if lmax <= 0:
        return None, None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)

    strength = W.sum(1)
    rng = np.random.default_rng(20260608 + sum(ord(c) for c in pat + band))
    subsets, target, achieved = spatial_strength_matched_subsets(
        strength, coords, epi, R_SUB, rng)

    coord_ok = np.isfinite(coords).all(1)
    epi_ok = np.where(epi)[0][coord_ok[np.where(epi)[0]]]
    obs_epi = subset_offdiag_means(lam_o, U_o, epi_ok, taus)
    if subsets:
        null_means = np.array([subset_offdiag_means(lam_o, U_o, s, taus)
                               for s in subsets])
    else:
        null_means = np.full((1, N_TAU), np.nan)

    geo = epi_geometry(pat, epi, probe, coords, reg)
    geo["band"] = band
    geo["rg_null_med"] = float(np.median(achieved)) if achieved else np.nan
    geo["rg_match_frac"] = (abs(geo["rg_null_med"] - target) / target
                            if achieved and target else np.nan)
    geo["n_spat_subsets"] = len(subsets)

    rows = []
    for ti in range(N_TAU):
        nm = null_means[:, ti]
        nm = nm[np.isfinite(nm)]
        oe = obs_epi[ti]
        row = {"patient": pat, "band": band, "tau_idx": ti,
               "tau": float(taus[ti]), "n_epi": int(epi_ok.size),
               "obs_EE": oe, "rg_epi": target,
               "rg_null_med": geo["rg_null_med"],
               "n_spat_subsets": len(subsets)}
        if nm.size >= 10 and np.isfinite(oe) and nm.std() > 0:
            row["z_spat_EE"] = (oe - nm.mean()) / nm.std()
            row["p_spat_EE"] = float((nm >= oe).mean())
        else:
            row["z_spat_EE"] = np.nan
            row["p_spat_EE"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows), geo


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for ti in sorted(df[df.band == band].tau_idx.unique()):
            d = df[(df.band == band) & (df.tau_idx == ti)]
            z = d["z_spat_EE"].to_numpy()
            z = z[np.isfinite(z)]
            wp = np.nan
            if z.size >= 3:
                try:
                    _, wp = wilcoxon(z, alternative="greater")
                except Exception:
                    pass
            out.append({
                "band": band, "tau_idx": int(ti), "n_patients": len(d),
                "med_z_spat_EE": float(np.median(z)) if z.size else np.nan,
                "wilcoxon_p_spat_EE": wp,
                "n_sig_spat_EE": int((d["p_spat_EE"].to_numpy() < 0.05).sum()),
            })
    return pd.DataFrame(out)


def write_readme(coh, geo_df, runtime):
    L = ["---", "name: epi_propagator_spatial_null",
         "scope: C5_spatial_matched_random_subset_null_epi_diffusion_community",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_90_epi_propagator_spatial_null.py",
         "---", "",
         "# C5 — spatial-matched null: is the epi diffusion community spatial "
         "clustering?", "",
         "**Head.** The audit_89 epi diffusion community (ρ_EE elevated beyond "
         "strength) is here tested against random groups matched to the epi set on "
         "BOTH node strength (quintile profile preserved) AND spatial compactness "
         "(radius of gyration in implant (x,y,z)). `z_spat_EE > 0` (small p) ⇒ epi "
         "co-diffuse MORE than an equal-size, equal-strength, equally-compact random "
         "group — so the community is not an artifact of epi sitting in one region. "
         f"Multiscale over {N_TAU} τ. Cohort n={int(coh.n_patients.max()) if not coh.empty else 0} "
         "(Pat_15 has 0 epi).", "",
         "## Cohort C5 (τ chosen to max the spatial-null z)", "",
         "| band | τ_idx | med z_spat_EE | Wilcoxon p | n sig/N |",
         "|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        cb = coh[coh.band == band].dropna(subset=["med_z_spat_EE"])
        if cb.empty:
            continue
        r = cb.loc[cb.med_z_spat_EE.idxmax()]
        L.append(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {int(r.tau_idx)} "
                 f"| {r.med_z_spat_EE:+.2f} | {r.wilcoxon_p_spat_EE:.4f} "
                 f"| {int(r.n_sig_spat_EE)}/{int(r.n_patients)} |")
    L += ["", "## Reading",
          "- This is the JOINT strength+spatial control. Combined with audit_89 C3 "
          "(strength alone) and C4 (cross-shaft), a band that survives here is "
          "neither hubness nor anatomical clustering.",
          "- Numbers per τ in `spatial_null_cohort.csv`; per patient in "
          "`spatial_null_per_patient.csv` (outliers read directly).",
          "- Spatial-match quality (achieved vs target radius of gyration) in "
          "`epi_geometry_per_patient.csv` (`rg_match_frac`; target tol "
          f"{RG_TOL:.2f}).", "",
          "## V2 — epi-set geometry (outlier diagnostics)", "",
          "| patient | n_epi | n_shafts | epi hemi | rg_epi/rg_all | rg match |",
          "|---|---|---|---|---|---|"]
    if geo_df is not None and not geo_df.empty:
        g1 = geo_df.drop_duplicates("patient")
        for _, r in g1.iterrows():
            L.append(f"| {r.patient} | {int(r.n_epi)} | {int(r.n_epi_shafts)} "
                     f"| {r.epi_hemi} | {r.rg_epi_frac:.2f} "
                     f"| {r.rg_match_frac:.2f} |")
    L += ["", f"- Wall-clock: {runtime:.1f}s"]
    (OUT / "README_spatial_null.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    frames, geos = [], []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            res, geo = per_patient_band(pat, band, args.verbose)
            if res is not None:
                frames.append(res)
                geos.append(geo)
                print(f"[audit_90] {pat}/{band}: n_epi={int(res.n_epi.iloc[0])} "
                      f"rg_match={geo['rg_match_frac']:.2f} "
                      f"({time.time()-tc:.1f}s)")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "spatial_null_per_patient.csv", index=False)
    geo_df = pd.DataFrame(geos)
    geo_df.to_csv(OUT / "epi_geometry_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "spatial_null_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, geo_df, runtime)

    print("\n[audit_90] C5 spatial-matched null (τ max z):")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        cb = coh[coh.band == band].dropna(subset=["med_z_spat_EE"])
        if cb.empty:
            continue
        r = cb.loc[cb.med_z_spat_EE.idxmax()]
        print(f"  {band:10s}: z_spat_EE {r.med_z_spat_EE:+.2f} "
              f"(p{r.wilcoxon_p_spat_EE:.3f}, {int(r.n_sig_spat_EE)}/{int(r.n_patients)})")
    print(f"[audit_90] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
