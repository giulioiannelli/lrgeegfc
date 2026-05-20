"""τ-sweep of ρ_split on raw `D(τ) = 1/ρ̂(τ)` for β.

The single-τ raw-D result at τ = 1/λ_max gave "trace everywhere"
(δ, α, β, γ_l all within-baseline significant), making ρ_split nearly
redundant with raw FC. This script tests whether SOME OTHER τ —
larger than 1/λ_max, probing coarser LRG scales — gives band-selective
behaviour (β-only or β-strongest).

For each (patient, band) we sweep τ over a log-spaced grid from
1/λ_max (finest LRG scale, current paper convention) to 1/λ_min^+
(coarsest LRG scale, where λ_min^+ is the smallest non-zero eigenvalue
of L, the Fiedler-related scale). At each τ:

    ρ̂(τ) = V diag(exp(-τλ)) V^T / Σ exp(-τλ)
    D(τ) = 1/ρ̂(τ),   max-symmetrized, zero diag
    Δ_task = D_taskT(τ) - D_rsPre_A(τ)
    Δ_rest = D_rsPost(τ) - D_rsPre_B(τ)
    ρ_split(τ) = Spearman(Δ_task, Δ_rest)

Output: `data/preprint/rho_split_raw_D/tau_sweep_per_patient.csv`
        `data/preprint/rho_split_raw_D/tau_sweep_cohort.csv`
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.notebook import move_to_rootf

move_to_rootf(pathname="lrgeegfc")


PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
N_TAU = 25
HALVES_CACHE = Path("data/cache/imcoh_lrg_halves")
OUT_DIR = Path("data/preprint/rho_split_raw_D")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_eig(p: Path):
    z = np.load(p, allow_pickle=True)
    return np.asarray(z["eigenvalues"], dtype=float), np.asarray(z["eigenvectors"], dtype=float)


def D_at_tau(eigvals: np.ndarray, eigvecs: np.ndarray, tau: float) -> np.ndarray:
    diag_exp = np.exp(-tau * eigvals)
    Z = float(np.sum(diag_exp))
    rho = (eigvecs * diag_exp) @ eigvecs.T / Z
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    finite = np.isfinite(D)
    if not finite.all():
        cap = np.nanmax(D[finite]) if finite.any() else 1e6
        D = np.where(finite, D, cap)
    return D


def triu_vec(M): return M[np.triu_indices(M.shape[0], k=1)]


def per_cell(patient: str, band: str, n_tau: int = N_TAU) -> list[dict] | None:
    paths = {
        "preA": HALVES_CACHE / patient / f"{band}_rest_pre_A_lrg_imcoh-abs.npz",
        "preB": HALVES_CACHE / patient / f"{band}_rest_pre_B_lrg_imcoh-abs.npz",
        "task": IMCOH_LRG_CACHE / patient / f"{band}_task_test_lrg_imcoh-abs.npz",
        "post": IMCOH_LRG_CACHE / patient / f"{band}_rest_post_lrg_imcoh-abs.npz",
    }
    for p in paths.values():
        if not p.exists():
            return None
    eigs = {k: load_eig(p) for k, p in paths.items()}

    # τ grid: from 1/λ_max to 1/λ_min+, log-spaced. Use the rsPre_A spectrum
    # to define the patient×band grid (use the geometric-mean of phases for
    # robustness)
    lam_maxes = []
    lam_mins_pos = []
    for k in ("preA", "preB", "task", "post"):
        ev = eigs[k][0]
        lam_maxes.append(float(np.max(ev)))
        pos = ev[ev > 1e-10]
        lam_mins_pos.append(float(np.min(pos)) if pos.size else float("nan"))
    lam_max_geom = float(np.exp(np.mean(np.log(lam_maxes))))
    lam_min_geom = float(np.exp(np.mean(np.log([x for x in lam_mins_pos if np.isfinite(x)]))))

    tau_min = 1.0 / lam_max_geom
    tau_max = 1.0 / lam_min_geom
    if tau_max <= tau_min:
        return None
    taus = np.geomspace(tau_min, tau_max, n_tau)

    rows = []
    for tau in taus:
        D = {k: D_at_tau(eigs[k][0], eigs[k][1], tau) for k in eigs}
        dt = triu_vec(D["task"]) - triu_vec(D["preA"])
        dr = triu_vec(D["post"]) - triu_vec(D["preB"])
        rho, _ = spearmanr(dt, dr)
        rows.append({
            "patient": patient, "band": band,
            "tau": float(tau),
            "tau_over_tau_min": float(tau / tau_min),
            "rho_split": float(rho),
        })
    return rows


print("Running τ-sweep ...")
all_rows = []
for band in BANDS:
    for pat in PATIENTS:
        cell = per_cell(pat, band)
        if cell:
            all_rows.extend(cell)
        else:
            print(f"  SKIP {pat} {band}")

df = pd.DataFrame(all_rows)
df.to_csv(OUT_DIR / "tau_sweep_per_patient.csv", index=False)

# Cohort: for each (band, tau_index), aggregate ρ_split across patients
df["tau_idx"] = df.groupby(["patient", "band"]).cumcount()
cohort_rows = []
for (band, tau_idx), grp in df.groupby(["band", "tau_idx"]):
    rhos = grp["rho_split"].values
    rhos = rhos[np.isfinite(rhos)]
    if rhos.size < 2:
        continue
    try:
        wp = float(wilcoxon(rhos, alternative="greater").pvalue)
    except Exception:
        wp = float("nan")
    cohort_rows.append({
        "band": band, "tau_idx": int(tau_idx),
        "tau_median":  float(np.median(grp["tau"].values)),
        "tau_over_tau_min_median": float(np.median(grp["tau_over_tau_min"].values)),
        "median_rho_split": float(np.median(rhos)),
        "n_pos": int((rhos > 0).sum()),
        "wilcoxon_p_one_sided": wp,
    })

cohort = pd.DataFrame(cohort_rows)
cohort.to_csv(OUT_DIR / "tau_sweep_cohort.csv", index=False)

# Print: pivot to (band × tau_idx) summary of cohort median + n_pos + p
print("\n=== Cohort ρ_split(τ) sweep — per band ===")
for band in BANDS:
    sub = cohort[cohort["band"] == band].sort_values("tau_idx")
    if sub.empty:
        continue
    print(f"\n--- {band} ---")
    print(f"  tau_idx  tau/tau_min     median_rho   n_pos    p")
    for _, r in sub.iterrows():
        sig = "*" if r["wilcoxon_p_one_sided"] < 0.05 else " "
        print(f"   {int(r['tau_idx']):>3}    {r['tau_over_tau_min_median']:>8.3f}    "
              f"{r['median_rho_split']:+.4f}    {int(r['n_pos'])}/10    "
              f"{r['wilcoxon_p_one_sided']:.4f}{sig}")

print(f"\nWrote {OUT_DIR}/tau_sweep_per_patient.csv")
print(f"Wrote {OUT_DIR}/tau_sweep_cohort.csv")
