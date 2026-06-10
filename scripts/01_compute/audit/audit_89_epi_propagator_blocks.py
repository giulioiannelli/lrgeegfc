#!/usr/bin/env python3
"""Audit 89 — propagator BLOCK structure: do epileptic nodes form a preferential
diffusion community, beyond node strength?

Scope pivot (user, 2026-06-08): stop scoring nodes by strength/hubness (trivial,
known). Use the LRG propagator ρ(τ) = e^{-τL}/Z as a *relational* lens — partition
node PAIRS by epileptic membership and ask whether diffusion preferentially stays
WITHIN the epileptic set:

    ρ_EE  = mean ρ_ij over epi–epi pairs
    ρ_EN  = mean ρ_ij over epi–nonepi pairs   (cross / interface)
    ρ_NN  = mean ρ_ij over nonepi–nonepi pairs

The whole question is whether ρ_EE is elevated BEYOND what each node's strength
already forces — so every block mean is compared to the matched-strength surrogate
ensemble (R=200, 4-cycle ±δ, cached rest_post eigs, seed 20260511), which fixes
each node's strength exactly while rewiring the topology. Multiscale: a τ-grid
from 1/λ_max to 10/λ_max. Per band.

Critical preamble
=================
(1) Claim: epileptic nodes form a strength-INDEPENDENT preferential diffusion
    block (ρ_EE elevated vs surrogate), band-specifically.
(2) Null: matched-strength surrogate — recompute ρ_EE / ρ_EN / ρ_NN on each
    surrogate (same epi node-set, strengths fixed) → z and upper-tail p.
(3) Strongest alternatives: (a) epi are hubs → ρ high everywhere incl. EE;
    (b) epi cluster on the same sEEG shafts → on-shaft proximity inflates ρ_EE;
    (c) global connectivity scale.
(4) Mechanical reach: matched-strength fixes each node's strength exactly, so
    z_EE>0 means epi–epi diffusion exceeds the strength expectation (kills a,c);
    the CROSS-PROBE-only EE block (epi pairs on different shafts) kills (b).
(5) Falsification: if z_EE ≈ 0 (obs within surrogate spread) in every band, the
    propagator block structure is just strength + geometry — no diffusion marker.

Outputs (``data/audit/epi_propagator_blocks/``)
    propagator_blocks_per_patient.csv   per (patient, band, tau)
    propagator_blocks_cohort.csv        per (band, tau)
    README.md
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
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_propagator_blocks"
OUT.mkdir(parents=True, exist_ok=True)

R_SUB = 400          # strength-matched random subsets for the specificity null
N_SBINS = 5          # strength quantile bins for matching


def subset_offdiag_means(lam, U, idx, taus):
    """Per-τ off-diagonal mean ρ within an arbitrary node subset `idx`."""
    Us = U[idx]
    k = idx.size
    iu = np.triu_indices(k, 1)
    out = np.empty(len(taus))
    for ti, t in enumerate(taus):
        w = np.exp(-t * lam); Z = w.sum()
        rho = (Us * w) @ Us.T / Z
        out[ti] = rho[iu].mean()
    return out


def strength_matched_subsets(strength, epi_mask, n_draws, rng):
    """Random node subsets matched to the epi set in size AND strength profile
    (stratified by strength quantile). Sampled from ALL nodes (mild epi overlap
    is conservative — it raises the null). Returns a list of index arrays."""
    n = strength.size
    ranks = np.argsort(np.argsort(strength))
    bins = np.minimum((ranks * N_SBINS) // n, N_SBINS - 1)
    epi_count = Counter(bins[epi_mask].tolist())
    pools = {b: np.where(bins == b)[0] for b in range(N_SBINS)}
    n_epi = int(epi_mask.sum())
    subsets = []
    for _ in range(n_draws):
        chosen, used = [], set()
        for b, k in epi_count.items():
            pool = pools[b]
            take = min(k, pool.size)
            if take:
                pick = rng.choice(pool, size=take, replace=False)
                chosen.extend(pick.tolist()); used.update(pick.tolist())
        short = n_epi - len(chosen)
        if short > 0:
            rest = np.array([x for x in range(n) if x not in used])
            if rest.size >= short:
                chosen.extend(rng.choice(rest, size=short, replace=False).tolist())
        if len(chosen) == n_epi:
            subsets.append(np.array(chosen))
    return subsets


def block_means(lam, U, epi_mask, probe, taus):
    """Per-τ (ρ_EE, ρ_EN, ρ_NN, ρ_EE_xprobe) off-diagonal block means.

    EE / EE_xprobe use the small epi×epi sub-propagator directly; EN / NN use
    O(N) eigenvector projections (e^T ρ e etc.) so surrogate loops stay cheap.
    """
    lam = np.clip(np.asarray(lam, float), 0.0, None)
    e = epi_mask.astype(float)
    n_e = int(epi_mask.sum()); n_n = int((~epi_mask).sum())
    if n_e < 2 or n_n < 2:
        return None
    Ue = U[epi_mask]                              # (n_e, N)
    p = Ue.sum(0)                                 # U^T e
    q = U.sum(0)                                  # U^T 1
    ep = probe[epi_mask]
    xpb = ep[:, None] != ep[None, :]
    iu = np.triu_indices(n_e, 1)
    xmask = xpb[iu]
    rows = []
    for t in taus:
        w = np.exp(-t * lam); Z = w.sum()
        rho_ee = (Ue * w) @ Ue.T / Z             # (n_e, n_e)
        ee_off = rho_ee[iu]
        mean_EE = float(ee_off.mean()) if ee_off.size else np.nan
        mean_EEx = float(ee_off[xmask].mean()) if xmask.any() else np.nan
        D_epi = float(np.trace(rho_ee))
        S_EE = float((w * p * p).sum() / Z)
        S_Eall = float((w * p * q).sum() / Z)
        S_all = float((w * q * q).sum() / Z)
        S_NN = S_all - 2 * S_Eall + S_EE
        D_non = 1.0 - D_epi
        mean_EN = (S_Eall - S_EE) / (n_e * n_n)
        mean_NN = (S_NN - D_non) / (n_n * (n_n - 1)) if n_n > 1 else np.nan
        rows.append((mean_EE, mean_EN, mean_NN, mean_EEx))
    return rows


def per_patient_band(pat, band, n_surr, swap, verbose):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_89] SKIP {pat}/{band}: {e}")
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probe = np.asarray(pm.probes, object)
    if epi.sum() < MIN_EPI or (~epi).sum() < 2:
        return None
    lam_o, U_o = np.linalg.eigh(np.diag(W.sum(1)) - W)
    lmax = float(lam_o[-1])
    if lmax <= 0:
        return None
    taus = np.geomspace(1.0 / lmax, 10.0 / lmax, N_TAU)
    obs = block_means(lam_o, U_o, epi, probe, taus)
    if obs is None:
        return None
    n_ee = int(epi.sum() * (epi.sum() - 1) // 2)
    ep = probe[epi]
    n_eex = int((ep[:, None] != ep[None, :])[np.triu_indices(epi.sum(), 1)].sum())

    # strength-matched random-subset specificity null (on the REAL graph):
    # is the epi group more diffusion-coherent than an arbitrary equal-size,
    # equal-strength group? This is the control matched-strength rewiring cannot
    # give (rewiring destroys ALL community structure, so any coherent subset
    # beats it; this asks whether epi beat a comparable random group).
    strength = W.sum(1)
    epi_idx = np.where(epi)[0]
    rng_sub = np.random.default_rng(20260608 + sum(ord(c) for c in pat + band))
    subsets = strength_matched_subsets(strength, epi, R_SUB, rng_sub)
    obs_epi_mean = subset_offdiag_means(lam_o, U_o, epi_idx, taus)
    if subsets:
        null_means = np.array([subset_offdiag_means(lam_o, U_o, s, taus)
                               for s in subsets])
    else:
        null_means = np.full((1, N_TAU), np.nan)

    evals, evecs = load_or_compute_eigs_at_path(
        es.surr_eig_path("full", pat, band, "rest_post", n_surr, swap),
        W, n_surr, swap, es.cell_rng(pat, band, "rest_post", "full"),
        verbose=verbose)
    R = evals.shape[0]
    # surrogate block means: (R, N_TAU, 4)
    surr = np.full((R, N_TAU, 4), np.nan)
    for r in range(R):
        if not np.isfinite(evals[r]).all():
            continue
        bm = block_means(evals[r], evecs[r], epi, probe, taus)
        if bm is not None:
            surr[r] = np.asarray(bm)

    rows = []
    labels = ["EE", "EN", "NN", "EEx"]
    for ti in range(N_TAU):
        row = {"patient": pat, "band": band, "tau_idx": ti,
               "tau": float(taus[ti]), "n_epi": int(epi.sum()),
               "n_ee_pairs": n_ee, "n_ee_xprobe_pairs": n_eex}
        for bi, lab in enumerate(labels):
            o = obs[ti][bi]
            s = surr[:, ti, bi]; s = s[np.isfinite(s)]
            row[f"obs_{lab}"] = o
            if s.size >= 5 and np.isfinite(o) and s.std() > 0:
                row[f"z_{lab}"] = (o - s.mean()) / s.std()
                row[f"p_{lab}"] = float((s >= o).mean())   # upper tail: epi diffuse MORE
            else:
                row[f"z_{lab}"] = np.nan; row[f"p_{lab}"] = np.nan
        # contrast EE - NN (paired across surrogates)
        oc = obs[ti][0] - obs[ti][2]
        sc = surr[:, ti, 0] - surr[:, ti, 2]; sc = sc[np.isfinite(sc)]
        if sc.size >= 5 and np.isfinite(oc) and sc.std() > 0:
            row["z_EE_minus_NN"] = (oc - sc.mean()) / sc.std()
            row["p_EE_minus_NN"] = float((sc >= oc).mean())
        else:
            row["z_EE_minus_NN"] = np.nan; row["p_EE_minus_NN"] = np.nan
        row["obs_EE_minus_NN"] = oc
        # strength-matched random-subset null (specificity, on real graph)
        nm = null_means[:, ti]; nm = nm[np.isfinite(nm)]
        oe = obs_epi_mean[ti]
        if nm.size >= 10 and np.isfinite(oe) and nm.std() > 0:
            row["z_rand_EE"] = (oe - nm.mean()) / nm.std()
            row["p_rand_EE"] = float((nm >= oe).mean())
        else:
            row["z_rand_EE"] = np.nan; row["p_rand_EE"] = np.nan
        row["n_rand_subsets"] = int(nm.size)
        rows.append(row)
    return pd.DataFrame(rows)


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for ti in sorted(df[df.band == band].tau_idx.unique()):
            d = df[(df.band == band) & (df.tau_idx == ti)]
            row = {"band": band, "tau_idx": int(ti), "n_patients": len(d)}
            for lab in ["EE", "EN", "NN", "EEx", "EE_minus_NN", "rand_EE"]:
                z = d[f"z_{lab}"].to_numpy(); z = z[np.isfinite(z)]
                row[f"med_z_{lab}"] = float(np.median(z)) if z.size else np.nan
                wp = np.nan
                if z.size >= 3:
                    try:
                        _, wp = wilcoxon(z, alternative="greater")
                    except Exception:
                        pass
                row[f"wilcoxon_p_{lab}"] = wp
                pcol = d[f"p_{lab}"].to_numpy()
                row[f"n_sig_{lab}"] = int(np.sum(pcol < 0.05))
            # raw descriptive block means
            for lab in ["EE", "EN", "NN", "EEx"]:
                o = d[f"obs_{lab}"].to_numpy(); o = o[np.isfinite(o)]
                row[f"med_obs_{lab}"] = float(np.median(o)) if o.size else np.nan
            out.append(row)
    return pd.DataFrame(out)


def write_readme(coh, runtime):
    L = ["---", "name: epi_propagator_blocks",
         "scope: direction3_propagator_block_diffusion_epi_vs_nonepi",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_89_epi_propagator_blocks.py",
         "---", "",
         "# Propagator block structure — do epileptic nodes diffuse preferentially "
         "to each other, beyond strength?", "",
         "**Head.** LRG propagator ρ(τ)=e^{-τL}/Z partitioned by epileptic "
         "membership: ρ_EE (epi–epi), ρ_EN (interface), ρ_NN (healthy–healthy). "
         "Each block mean is z-scored against the matched-strength surrogate "
         "(strengths fixed, topology rewired), so z_EE>0 = epi–epi diffusion "
         "exceeds the strength expectation. `EEx` = epi–epi restricted to pairs on "
         "DIFFERENT shafts (kills the on-shaft-proximity confound). Multiscale over "
         f"{N_TAU} τ from 1/λmax to 10/λmax.", "",
         "## Epi-specific diffusion enrichment (cohort, τ chosen to max the "
         "random-subset z)", "",
         "Two honest epi-specific measures + the inflated raw for context. "
         "`z_rand_EE` = epi vs strength-matched random groups on the REAL graph "
         "(the decisive specificity control). `z(EE−NN)` = epi block vs healthy "
         "block, both vs matched-strength rewiring. `raw z_EE` is inflated by "
         "general community structure (rises with frequency) — context only.", "",
         "| band | τ_idx | **z_rand_EE** | p | n sig/N | **z(EE−NN)** | p | "
         "z_EEx | raw z_EE (inflated) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        cb = coh[coh.band == band]
        cb_valid = cb.dropna(subset=["med_z_rand_EE"])
        if cb_valid.empty:
            continue
        r = cb_valid.loc[cb_valid.med_z_rand_EE.idxmax()]
        L.append(
            f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {int(r.tau_idx)} "
            f"| {r.med_z_rand_EE:+.2f} | {r.wilcoxon_p_rand_EE:.4f} "
            f"| {int(r.n_sig_rand_EE)}/{int(r.n_patients)} "
            f"| {r.med_z_EE_minus_NN:+.2f} | {r.wilcoxon_p_EE_minus_NN:.4f} "
            f"| {r.med_z_EEx:+.2f} | {r.med_z_EE:+.2f} |")
    L += ["", "## Reading",
          "- **`z_rand_EE` is the decisive measure**: epi within-group diffusion "
          "vs strength-matched RANDOM groups on the real graph. >0 (small p) ⇒ "
          "epileptic contacts diffuse to each other more than an arbitrary equal-"
          "size, equal-strength set of contacts — i.e. a genuine, strength-"
          "independent, epi-SPECIFIC diffusion community. Hubness cannot produce "
          "this (strength is matched in the null).",
          "- `z(EE−NN)` is the second epi-specific measure: epi block enriched "
          "MORE than the healthy block, both vs matched-strength rewiring.",
          "- `raw z_EE` (vs matched-strength rewiring, absolute) is INFLATED: "
          "rewiring destroys all community structure, so any coherent subset beats "
          "it; raw z_EE rises with frequency (sparser graph), so it is context "
          "only, NOT the epi claim. `z_EEx` = same on cross-shaft epi pairs.",
          "- Per-τ detail in `propagator_blocks_cohort.csv`; per-patient in "
          "`propagator_blocks_per_patient.csv`.",
          f"- Wall-clock: {runtime:.1f}s"]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    frames = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            res = per_patient_band(pat, band, args.n_surrogates,
                                   args.swap_factor, args.verbose)
            if res is not None:
                frames.append(res)
                print(f"[audit_89] {pat}/{band}: n_epi={int(res.n_epi.iloc[0])} "
                      f"({time.time()-tc:.1f}s)")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "propagator_blocks_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "propagator_blocks_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, runtime)

    print("\n[audit_89] epi-SPECIFIC diffusion enrichment (τ max z_rand_EE):")
    print("  band       : z_rand_EE (vs random group) | z(EE-NN) | raw z_EE")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        cb = coh[coh.band == band].dropna(subset=["med_z_rand_EE"])
        if cb.empty:
            continue
        r = cb.loc[cb.med_z_rand_EE.idxmax()]
        print(f"  {band:10s}: {r.med_z_rand_EE:+.2f} (p{r.wilcoxon_p_rand_EE:.3f}, "
              f"{int(r.n_sig_rand_EE)}/{int(r.n_patients)}) | "
              f"{r.med_z_EE_minus_NN:+.2f} (p{r.wilcoxon_p_EE_minus_NN:.3f}) | "
              f"{r.med_z_EE:+.2f}")
    print(f"[audit_89] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
