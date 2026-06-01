"""audit_73 — comparative methodology audit of cross-phase pair-trace measures.

Four per-pair cross-phase measures
    {ρ_Spearman, ρ_Pearson, s_TR (asymmetric slope), R² = Pearson²}
on three LRG primitives
    {ρ̂(τ_max), D(τ_max) = 1/ρ̂, D_coph = cophenet(UPGMA(D(τ_max)))}
under the split-baseline + matched-strength null, over the full n=10 cohort
× 6 bands.

This is a COMPARATIVE METHODOLOGY AUDIT, not a manuscript-replacement run.
The current §5.3 / §3.2 headline is `ρ_split^coph` = ρ_Spearman on D_coph
(audit_63). Here we add the magnitude-aware symmetric (ρ_Pearson), the
asymmetric through-origin regression slope `s_TR = ⟨Δ_task, Δ_rest⟩/‖Δ_task‖²`
("fraction of the task-induced per-pair shift recovered in rsPost"), and its
companion R², and we run all four on all three primitives so the user can see
which (measure, primitive) cells reproduce the β-band selectivity that
`ρ_split^coph` carries — on more principled, units-bearing footing.

NAMING — the regression slope is `s_TR` (math) / `slope` (code, CSV column);
NEVER `β`, which is reserved for the 13–30 Hz band. See
`.agents/guides/04_rules/never-always-list.md` and
`feedback_no_beta_for_regression_slope`.

Scope report (read first):
  `.agents/guides/task-persistence-investigation/2026-05-30_asymmetric-pair-trace-regression.md`

Split-baseline convention (identical to audit_63 / preprint_05):
  Δ_task = X(task_test) − X(rest_pre_A)
  Δ_rest = X(rest_post) − X(rest_pre_B)
with rest_pre_A ⊥ rest_pre_B (independent rsPre halves) so the two
delta-vectors do not share a baseline-noise term.

Null: per-patient observed scalar vs R=200 matched-strength surrogate scalars
(4-cycle ±δ strength-preserving rewiring, seed 20260511, swap_factor 20). The
surrogate eigendecomposition ensemble is already on disk for all 10 × 6 × 4
cells (`data/cache/matched_strength_surrogate_lrg/`); this script loads it and
recomputes the four measures on the three primitives per surrogate. Per
`feedback_matched_strength_mandatory`, the matched-strength null is re-run for
`s_TR` here — ρ_split's clearance does NOT transfer to a magnitude-aware
measure.

Output (LONG format — chosen over the 48-wide schema for trivial downstream
pivoting; resolves scope §10 open-question #1):
  data/audit/pair_trace_measures_comparison/per_patient_per_band.csv
      one row per (patient, band, primitive, measure) = 10×6×3×4 = 720 rows
  data/audit/pair_trace_measures_comparison/cohort_summary.csv
      one row per (band, primitive, measure) = 6×3×4 = 72 rows, with
      cohort-paired one-sided Wilcoxon (obs > surrogate-median) + LOO
  data/audit/pair_trace_measures_comparison/README.md  (auto lab notes)

CLI:
  python scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py \
      [--patients Pat_05] [--bands beta] [--max-surr 200] \
      [--out-suffix _dryrun] [-v]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr, wilcoxon

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import CACHE_ROOT, IMCOH_LRG_CACHE
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.utils.metrics.hypothesis import (
    loo_sensitivity,
    regression_slope_through_origin,
)
from lrg_eegfc.utils.metrics.tree import cophenet_matrix
from lrg_eegfc.utils.surrogate.matched_strength import load_or_compute_surrogate_eigs
from lrg_eegfc.workflow.fc import load_fc_matrix

move_to_rootf(pathname="lrgeegfc")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BANDS = ["alpha", "beta", "low_gamma", "delta", "theta", "high_gamma"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
PRIMITIVES = ("rhohat", "D", "Dcoph")
MEASURES = ("spearman", "pearson", "slope", "rsq")

N_SURROGATES = 200
SWAP_FACTOR = 20
SEED = 20260511

HALVES_LRG_CACHE = Path("data/cache/imcoh_lrg_halves")
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
OUT_DIR = Path("data/audit/pair_trace_measures_comparison")


# ---------------------------------------------------------------------------
# Primitive reconstruction from a Laplacian eigendecomposition
# ---------------------------------------------------------------------------
def primitives_from_eig(
    eigvals: np.ndarray, eigvecs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(ρ̂(τ_max), D(τ_max))`` full N×N matrices from one eig pair.

    Mirrors `lrg_ultrametric_condensed` (audit_63) / `D_raw_from_eig`
    (preprint_05) so the D_coph path reproduces `ρ_split^coph` exactly.
    """
    lam_max = float(np.max(eigvals))
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    Z = float(np.sum(diag_exp))
    rho_hat = (eigvecs * diag_exp) @ eigvecs.T / Z
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho_hat
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    finite = np.isfinite(D)
    if not finite.all():
        cap = np.nanmax(D[finite]) if finite.any() else 1e6
        D = np.where(finite, D, cap)
    return rho_hat, D


def primitive_pair_vectors(
    eigs_by_phase: dict[str, tuple[np.ndarray, np.ndarray]],
    iu: tuple[np.ndarray, np.ndarray],
) -> dict[str, dict[str, np.ndarray]]:
    """Build per-pair upper-triangle vectors for all 3 primitives × 4 phases.

    Returns ``{primitive: {phase: vec[m]}}``. The D_coph vector uses the
    condensed cophenet (same row-major upper-triangle ordering as
    ``np.triu_indices``), so it aligns pair-for-pair with the ρ̂ / D vectors.
    """
    out = {p: {} for p in PRIMITIVES}
    for phase, (ev, vc) in eigs_by_phase.items():
        rho_hat, D = primitives_from_eig(ev, vc)
        out["rhohat"][phase] = rho_hat[iu]
        out["D"][phase] = D[iu]
        Zc = linkage(squareform(D, checks=False), method="average")
        out["Dcoph"][phase] = cophenet_matrix(Zc, condensed=True)
    return out


def four_measures(dT: np.ndarray, dR: np.ndarray) -> dict[str, float]:
    """ρ_Spearman, ρ_Pearson, s_TR slope, R² on a pair of shift vectors."""
    sp, _ = spearmanr(dT, dR)
    pe, _ = pearsonr(dT, dR)
    slope, rsq = regression_slope_through_origin(dT, dR)
    return {
        "spearman": float(sp),
        "pearson": float(pe),
        "slope": float(slope),
        "rsq": float(rsq),
    }


def scalars_from_pair_vectors(
    pv: dict[str, dict[str, np.ndarray]]
) -> dict[str, dict[str, float]]:
    """Apply the split-baseline deltas + four measures to all 3 primitives."""
    out = {}
    for prim in PRIMITIVES:
        dT = pv[prim]["task_test"] - pv[prim]["rest_pre_A"]
        dR = pv[prim]["rest_post"] - pv[prim]["rest_pre_B"]
        out[prim] = four_measures(dT, dR)
    return out


# ---------------------------------------------------------------------------
# Observed + surrogate loaders
# ---------------------------------------------------------------------------
def _load_eig(p: Path) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(p, allow_pickle=True)
    return (
        np.asarray(z["eigenvalues"], dtype=float),
        np.asarray(z["eigenvectors"], dtype=float),
    )


def observed_eigs(pat: str, band: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {
        "rest_pre_A": _load_eig(
            HALVES_LRG_CACHE / pat / f"{band}_rest_pre_A_lrg_imcoh-abs.npz"
        ),
        "rest_pre_B": _load_eig(
            HALVES_LRG_CACHE / pat / f"{band}_rest_pre_B_lrg_imcoh-abs.npz"
        ),
        "task_test": _load_eig(
            IMCOH_LRG_CACHE / pat / f"{band}_task_test_lrg_imcoh-abs.npz"
        ),
        "rest_post": _load_eig(
            IMCOH_LRG_CACHE / pat / f"{band}_rest_post_lrg_imcoh-abs.npz"
        ),
    }


def load_W_for_surrogate(pat: str, phase: str, band: str) -> np.ndarray:
    """FC matrix used only for the surrogate-ensemble shape check (cache hit)."""
    if phase == "rest_pre_A":
        W = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_A_imcoh_abs.npy")
    elif phase == "rest_pre_B":
        W = np.load(HALVES_FC_CACHE / pat / f"{band}_rest_pre_B_imcoh_abs.npy")
    elif phase == "task_test":
        W = load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs")
    elif phase == "rest_post":
        W = load_fc_matrix(pat, "rest_post", band, fc_method="imcoh_abs")
    else:
        raise ValueError(phase)
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


# ---------------------------------------------------------------------------
# Per-cell pipeline
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, max_surr: int, verbose: bool) -> list[dict] | None:
    t0 = time.time()
    try:
        obs_e = observed_eigs(pat, band)
    except FileNotFoundError as e:
        if verbose:
            print(f"[{pat} {band}] SKIP obs: {e}")
        return None

    N = obs_e["rest_pre_A"][1].shape[0]
    if any(e[1].shape[0] != N for e in obs_e.values()):
        if verbose:
            print(f"[{pat} {band}] SKIP: phase node-count mismatch")
        return None
    iu = np.triu_indices(N, k=1)
    n_pairs = iu[0].size

    obs_scalars = scalars_from_pair_vectors(primitive_pair_vectors(obs_e, iu))

    # --- surrogate ensemble (cache hit expected) ---
    rng = np.random.default_rng(SEED)
    surr_e = {}
    for phase in PHASES:
        W = load_W_for_surrogate(pat, phase, band)
        evals, evecs = load_or_compute_surrogate_eigs(
            pat, band, phase, W,
            n_surr=N_SURROGATES, swap_factor=SWAP_FACTOR, seed=SEED, rng=rng,
            fc_method="imcoh_abs", verbose=False,
        )
        surr_e[phase] = (evals, evecs)

    R = min(max_surr, N_SURROGATES)
    surr_scalars = {p: {m: np.full(R, np.nan) for m in MEASURES} for p in PRIMITIVES}
    for r in range(R):
        eigs_r = {ph: (surr_e[ph][0][r], surr_e[ph][1][r]) for ph in PHASES}
        if any(
            not (np.all(np.isfinite(ev)) and np.all(np.isfinite(vc)))
            for ev, vc in eigs_r.values()
        ):
            continue
        try:
            sc = scalars_from_pair_vectors(primitive_pair_vectors(eigs_r, iu))
        except Exception:
            continue
        for prim in PRIMITIVES:
            for m in MEASURES:
                surr_scalars[prim][m][r] = sc[prim][m]

    # --- emit one long-format row per (primitive, measure) ---
    rows = []
    for prim in PRIMITIVES:
        for m in MEASURES:
            obs = obs_scalars[prim][m]
            surr = surr_scalars[prim][m]
            surr = surr[np.isfinite(surr)]
            if surr.size == 0:
                continue
            surr_mean = float(np.mean(surr))
            surr_std = float(np.std(surr, ddof=1)) if surr.size > 1 else np.nan
            z = (obs - surr_mean) / surr_std if surr_std and surr_std > 0 else np.nan
            rows.append({
                "patient": pat, "band": band,
                "primitive": prim, "measure": m,
                "N_nodes": int(N), "n_pairs": int(n_pairs),
                "n_surr_finite": int(surr.size),
                "obs": float(obs),
                "surr_mean": surr_mean,
                "surr_p50": float(np.median(surr)),
                "surr_std": surr_std,
                "surr_p5": float(np.quantile(surr, 0.05)),
                "surr_p95": float(np.quantile(surr, 0.95)),
                "obs_z": float(z),
                "obs_p_one_sided": float(np.mean(surr >= obs)),
            })
    if verbose:
        d = obs_scalars["Dcoph"]
        print(
            f"[{pat} {band}] Dcoph: spearman={d['spearman']:+.4f} "
            f"pearson={d['pearson']:+.4f} s_TR={d['slope']:+.4f} "
            f"R²={d['rsq']:.4f}  ({time.time() - t0:.1f}s)"
        )
    return rows


# ---------------------------------------------------------------------------
# Cohort aggregation
# ---------------------------------------------------------------------------
def _wilcoxon_greater(diff: np.ndarray) -> float:
    d = np.asarray(diff, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < 3 or np.allclose(d, 0.0):
        return np.nan
    try:
        return float(wilcoxon(d, alternative="greater").pvalue)
    except ValueError:
        return np.nan


def cohort_summary(df_pp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in df_pp["band"].unique():
        for prim in PRIMITIVES:
            for m in MEASURES:
                sub = df_pp[
                    (df_pp.band == band)
                    & (df_pp.primitive == prim)
                    & (df_pp.measure == m)
                ].copy()
                if sub.empty:
                    continue
                obs = sub["obs"].to_numpy(float)
                surr = sub["surr_p50"].to_numpy(float)
                diff = obs - surr
                labels = sub["patient"].tolist()
                p_w = _wilcoxon_greater(diff)
                med_obs = float(np.median(obs))
                med_surr = float(np.median(surr))
                ratio = med_obs / med_surr if med_surr != 0 else np.inf
                n_above = int((sub["obs_p_one_sided"] < 0.05).sum())
                loo = loo_sensitivity(diff, test_fn=_wilcoxon_greater, labels=labels)
                rows.append({
                    "band": band, "primitive": prim, "measure": m,
                    "n": len(sub),
                    "obs_median": med_obs,
                    "surr_median": med_surr,
                    "ratio_obs_over_surr": ratio,
                    "n_above_p05": n_above,
                    "wilcoxon_p_greater": p_w,
                    "loo_max_p": loo["worst"],
                    "loo_argmax_patient": loo["worst_patient"],
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--patients", default=",".join(PATIENTS_4PHASE))
    ap.add_argument("--bands", default=",".join(BANDS))
    ap.add_argument("--max-surr", type=int, default=N_SURROGATES,
                    help="cap surrogates (dry-run speed); default 200")
    ap.add_argument("--out-suffix", default="",
                    help="suffix on output CSV names (dry-run isolation)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    patients = [p.strip() for p in args.patients.split(",") if p.strip()]
    bands = [b.strip() for b in args.bands.split(",") if b.strip()]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for band in bands:
        print(f"\n=== {band} ===")
        for pat in patients:
            r = per_cell(pat, band, args.max_surr, args.verbose)
            if r:
                all_rows.extend(r)

    if not all_rows:
        print("No rows produced — check caches.")
        return

    df_pp = pd.DataFrame(all_rows)
    pp_path = OUT_DIR / f"per_patient_per_band{args.out_suffix}.csv"
    df_pp.to_csv(pp_path, index=False)
    print(f"\nWrote {pp_path}  ({len(df_pp)} rows)")

    df_co = cohort_summary(df_pp)
    co_path = OUT_DIR / f"cohort_summary{args.out_suffix}.csv"
    df_co.to_csv(co_path, index=False)
    print(f"Wrote {co_path}  ({len(df_co)} rows)")

    # headline: D_coph rows across bands/measures
    print("\n=== D_coph cohort (matched-strength Wilcoxon p_greater) ===")
    head = df_co[df_co.primitive == "Dcoph"].pivot(
        index="band", columns="measure", values="wilcoxon_p_greater"
    )
    print(head.to_string(float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
