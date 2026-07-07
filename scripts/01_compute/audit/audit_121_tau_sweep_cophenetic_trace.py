"""audit_121 — τ-sensitivity of the cophenetic cohort trace.

Scope: `.agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md`

Every cophenetic task-trace verdict was computed at the single finest diffusion
time τ = 1/λ_max. This sweeps τ = α/λ_max (per phase, so α=1 reproduces the
locked cache bit-exactly) over α ∈ [0.5, 40] — finer-than-current → genuine LRG
mesoscale → past-Fiedler collapse — and asks whether the cohort trace
(`T_d`, `ρ_split`) and its band ranking are τ-flat (robust) or τ-structured.

Observed-statistic only; no matched-strength null here (that is the follow-up
for any positive claim at α≠1). D(τ) reconstructs exactly from cached eigenpairs
(no expm, no FC reload). The α=1 reconstruction is asserted equal to the cached
`ultrametric_matrix` per (patient, band, phase).

Outputs:
  data/audit/tau_sweep_trace/per_patient.csv
  data/audit/tau_sweep_trace/cohort.csv
"""
from __future__ import annotations

import warnings
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.notebook import move_to_rootf
from lrgsglib.core import compute_normalized_linkage, extract_ultrametric_matrix

move_to_rootf(pathname="lrgeegfc")

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HALVES_CACHE = Path("data/cache/imcoh_lrg_halves")
OUT_DIR = Path("data/audit/tau_sweep_trace")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Shared α-grid: α=1 present (anchor), 0.5 = finer than current, 40 = past
# Fiedler (collapse) for every cell. τ_phase = α / λmax_phase.
ALPHA = np.unique(np.concatenate([
    [0.5, 0.7, 0.85, 1.0],
    np.geomspace(1.0, 40.0, 26),
]))
M = len(ALPHA)
VAR_FLOOR = 1e-10          # below this, D_coph is ~constant → ρ^coph unstable
_EMPTY_G: dict[int, nx.Graph] = {}


def _giant(n: int) -> nx.Graph:
    if n not in _EMPTY_G:
        _EMPTY_G[n] = nx.empty_graph(n)
    return _EMPTY_G[n]


def load_cell(patient: str, band: str, phase: str, halves: bool = False):
    root = HALVES_CACHE if halves else IMCOH_LRG_CACHE
    sub = f"{phase}" if not halves else phase  # phase already e.g. rest_pre_A
    p = root / patient / f"{band}_{sub}_lrg_imcoh-abs.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    ev = np.asarray(z["eigenvalues"], float)
    U = np.asarray(z["eigenvectors"], float)
    n = int(z["n_nodes"])
    ultra = np.asarray(z["ultrametric_matrix"], float)  # condensed, canonical
    return dict(ev=ev, U=U, n=n, ultra=ultra, lam_max=float(np.max(ev)),
                lam_gap=float(np.min(ev[ev > 1e-10])))


def D_raw(ev, U, tau):
    diag = np.exp(-tau * ev)
    Z = float(diag.sum())
    rho = (U * diag) @ U.T / Z
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    fin = np.isfinite(D)
    if not fin.all():
        D = np.where(fin, D, np.nanmax(D[fin]) if fin.any() else 1e12)
    return D


def D_coph(D, n):
    Zlink, _, _ = compute_normalized_linkage(squareform(D), _giant(n))
    return extract_ultrametric_matrix(Zlink, n)  # square


def triu(M):
    return M[np.triu_indices(M.shape[0], k=1)]


def rho_coph(a, b):
    """Spearman with zero-variance guard (collapsed geometry → nan)."""
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = spearmanr(a, b).statistic
    return float(r)


def process(patient: str, band: str):
    full = {ph: load_cell(patient, band, ph)
            for ph in ("rest_pre", "task_test", "rest_post")}
    if any(v is None for v in full.values()):
        return None
    halves = {h: load_cell(patient, band, h, halves=True)
              for h in ("rest_pre_A", "rest_pre_B")}
    have_halves = all(v is not None for v in halves.values())

    # Anchor check (α=1, each phase own 1/λmax) — bit-exact vs cache.
    for ph, c in full.items():
        u1_raw = D_raw(c["ev"], c["U"], 1.0 / c["lam_max"])
        u1 = squareform(D_coph(u1_raw, c["n"]))
        if not np.allclose(u1, c["ultra"], atol=1e-8):
            raise AssertionError(
                f"anchor mismatch {patient}/{band}/{ph}: "
                f"max|Δ|={np.max(np.abs(u1 - c['ultra'])):.2e}")

    lam_max_g = float(np.exp(np.mean(np.log(
        [full[p]["lam_max"] for p in full]))))
    lam_gap_g = float(np.exp(np.mean(np.log(
        [full[p]["lam_gap"] for p in full]))))
    alpha_F = lam_max_g / lam_gap_g

    # canonical cophenetic vectors (α=1) for self-similarity reference
    u_ref = {ph: c["ultra"] for ph, c in full.items()}

    rows = []
    for ai, alpha in enumerate(ALPHA):
        # per-phase τ = α / λmax_phase ; cophenetic + raw triu vectors
        dco, draw = {}, {}
        for ph, c in full.items():
            Dr = D_raw(c["ev"], c["U"], alpha / c["lam_max"])
            draw[ph] = triu(Dr)
            dco[ph] = squareform(D_coph(Dr, c["n"]))
        hv = {}
        if have_halves:
            for h, c in halves.items():
                Dr = D_raw(c["ev"], c["U"], alpha / c["lam_max"])
                hv[(h, "raw")] = triu(Dr)
                hv[(h, "coph")] = squareform(D_coph(Dr, c["n"]))

        # cophenetic legs + T_d (locked: T_d = ρ(task,post) − ρ(pre,task))
        r_pt_c = rho_coph(dco["rest_pre"], dco["task_test"])
        r_tp_c = rho_coph(dco["task_test"], dco["rest_post"])
        Td_c = r_tp_c - r_pt_c
        r_pt_r = rho_coph(draw["rest_pre"], draw["task_test"])
        r_tp_r = rho_coph(draw["task_test"], draw["rest_post"])
        Td_r = r_tp_r - r_pt_r

        # split-baseline ρ_split = Spearman(task−preA, post−preB)
        # placebo ρ_indep = Spearman(preA−preB, task−post): manifestly trace-free
        #   (baseline-split noise vs a task/post difference, NO shared phase term).
        #   ≈0 at all τ unless coarsening manufactures spurious correlation between
        #   independent difference vectors → the low-rank-shared-mode artifact test.
        rs_c = rs_r = ri_c = ri_r = np.nan
        if have_halves:
            rs_c = rho_coph(dco["task_test"] - hv[("rest_pre_A", "coph")],
                            dco["rest_post"] - hv[("rest_pre_B", "coph")])
            rs_r = rho_coph(draw["task_test"] - hv[("rest_pre_A", "raw")],
                            draw["rest_post"] - hv[("rest_pre_B", "raw")])
            ri_c = rho_coph(hv[("rest_pre_A", "coph")] - hv[("rest_pre_B", "coph")],
                            dco["task_test"] - dco["rest_post"])
            ri_r = rho_coph(hv[("rest_pre_A", "raw")] - hv[("rest_pre_B", "raw")],
                            draw["task_test"] - draw["rest_post"])

        # geometry diagnostics (dco[ph] and cached u_ref[ph] are both condensed,
        # same scipy ordering → compare directly; == 1.0 at α=1 by construction)
        gself = {ph: rho_coph(dco[ph], u_ref[ph]) for ph in full}
        var_pre = float(np.nanvar(dco["rest_pre"]))

        rows.append(dict(
            patient=patient, band=band, alpha_idx=ai, alpha=float(alpha),
            alpha_fiedler=alpha_F,
            T_d_coph=Td_c, T_d_raw=Td_r,
            rho_split_coph=rs_c, rho_split_raw=rs_r,
            rho_indep_coph=ri_c, rho_indep_raw=ri_r,
            rho_pre_task_coph=r_pt_c, rho_task_post_coph=r_tp_c,
            G_self_pre=gself["rest_pre"], G_self_task=gself["task_test"],
            G_self_post=gself["rest_post"],
            var_pre_coph=var_pre,
            collapsed=int(var_pre < VAR_FLOOR),
        ))
    return rows


def main():
    print(f"τ-sweep: {len(PATIENTS)} patients × {len(BANDS)} bands × {M} α-points "
          f"(α∈[{ALPHA[0]:.2f}, {ALPHA[-1]:.1f}], α=1 anchored)")
    all_rows = []
    for band in BANDS:
        for pat in PATIENTS:
            r = process(pat, band)
            if r:
                all_rows.extend(r)
            else:
                print(f"  SKIP {pat} {band} (missing cache)")
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_DIR / "per_patient.csv", index=False)
    print(f"  anchor checks passed for all loaded cells (α=1 == cache)")

    # cohort aggregation per (band, alpha_idx)
    def w1(x):
        x = x[np.isfinite(x)]
        if x.size < 3 or np.allclose(x, 0):
            return np.nan
        try:
            return float(wilcoxon(x, alternative="greater").pvalue)
        except Exception:
            return np.nan

    crows = []
    for (band, ai), g in df.groupby(["band", "alpha_idx"]):
        tdc = g["T_d_coph"].values
        tdr = g["T_d_raw"].values
        rsc = g["rho_split_coph"].values
        ric = g["rho_indep_coph"].values
        crows.append(dict(
            band=band, alpha_idx=int(ai), alpha=float(g["alpha"].iloc[0]),
            median_alpha_fiedler=float(np.median(g["alpha_fiedler"])),
            median_T_d_coph=float(np.nanmedian(tdc)),
            n_trace_coph=int(np.nansum(tdc > 0)),
            n_finite_coph=int(np.isfinite(tdc).sum()),
            wilcox_p_coph=w1(tdc),
            median_T_d_raw=float(np.nanmedian(tdr)),
            n_trace_raw=int(np.nansum(tdr > 0)),
            wilcox_p_raw=w1(tdr),
            median_rho_split_coph=float(np.nanmedian(rsc)),
            n_pos_split=int(np.nansum(rsc > 0)),
            wilcox_p_split=w1(rsc),
            median_rho_indep_coph=float(np.nanmedian(ric)),
            median_abs_rho_indep_coph=float(np.nanmedian(np.abs(ric))),
            median_G_self=float(np.nanmedian(g["G_self_task"].values)),
            frac_collapsed=float(np.mean(g["collapsed"].values)),
        ))
    cohort = pd.DataFrame(crows)
    cohort.to_csv(OUT_DIR / "cohort.csv", index=False)

    # print: anchor row (α=1) + a coarse row per band
    a1 = int(np.argmin(np.abs(ALPHA - 1.0)))
    print("\n=== α = 1 anchor (== locked canonical verdict) ===")
    print(" band        med T_d_coph  n_trace  wilcox_p   med ρ_split  n_pos")
    for band in BANDS:
        r = cohort[(cohort.band == band) & (cohort.alpha_idx == a1)]
        if r.empty:
            continue
        r = r.iloc[0]
        print(f" {band:<11} {r.median_T_d_coph:+.4f}     {r.n_trace_coph}/10"
              f"     {r.wilcox_p_coph:.4f}    {r.median_rho_split_coph:+.4f}"
              f"   {r.n_pos_split}/10")

    print("\n=== β τ-response: TRACE (T_d, ρ_split) vs PLACEBO (ρ_indep) vs α ===")
    print("   ρ_indep is trace-free → must stay ≈0; if it tracks ρ_split at coarse")
    print("   τ, the ρ_split coarse-τ growth is a shared-low-rank artifact.\n")
    bsub = cohort[cohort.band == "beta"].sort_values("alpha_idx")
    print(" alpha   T_d_coph  ρ_split   ρ_indep(placebo)  G_self  frac_coll  (αF≈%.1f)"
          % bsub["median_alpha_fiedler"].iloc[0])
    for _, r in bsub.iterrows():
        mark = " <α=1" if r.alpha_idx == a1 else ("  <Fiedler"
               if abs(r.alpha - r.median_alpha_fiedler) < 0.15 * r.median_alpha_fiedler else "")
        print(f" {r.alpha:6.2f}  {r.median_T_d_coph:+.4f}   {r.median_rho_split_coph:+.4f}"
              f"    {r.median_rho_indep_coph:+.4f}        {r.median_G_self:+.3f}"
              f"   {r.frac_collapsed:.2f}{mark}")

    print(f"\nWrote {OUT_DIR}/per_patient.csv  ({len(df)} rows)")
    print(f"Wrote {OUT_DIR}/cohort.csv")


if __name__ == "__main__":
    main()
