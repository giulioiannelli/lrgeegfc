#!/usr/bin/env python3
"""Audit 40 -- tau-sweep CTM rho_split on per-phase intrinsic [tau_min, tau*].

For each (patient, band, u) with u in linspace(0, 1, 8), and each phase
phi in {TT, RPost, RPreA, RPreB, RPostA, RPostB}:

    tau_phi(u) = tau_min(phi) * (tau*(phi) / tau_min(phi)) ** u
    tau_min(phi) = 1 / lambda_max(phi)
    tau*(phi)    = argmax_tau C(tau; eigvals(phi))   in [tau_min, 1/lambda_2]

Each phase is evaluated at its own intrinsic time scale at the same
normalized fraction u, so the comparison axis u is comparable across
patients and bands.

CTM-style scalars per (patient, band, u):
    rho_split(p, b, u) = Spearman( triu(D^TT(tau_TT(u)) - D^RPreA(tau_RPreA(u))),
                                    triu(D^RPost(tau_RPost(u)) - D^RPreB(tau_RPreB(u))) )
    rho_drift(p, b, u) = Spearman( triu(D^RPreB - D^RPreA), triu(D^RPostB - D^RPostA) )
                          (matches official h2e_split_half / CTM Run C definition:
                           within-session noise consistency between RPre and RPost halves)

Bands (first pass): alpha, beta, low_gamma (the controlled CTM bands).

Outputs
-------
data/audit/ctm_tau_sweep/tau_sweep_rho.csv
data/audit/ctm_tau_sweep/tau_star_per_pb.csv
data/outputs/figures/section_5_lrg_trace/ctm_tau_sweep/cohort.pdf
data/outputs/figures/section_5_lrg_trace/ctm_tau_sweep/per_patient.pdf
data/outputs/figures/section_5_lrg_trace/ctm_tau_sweep/c_tau_sanity.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.diagnostics import (
    _entropy_at_tau, _compute_propagator, _ultrametric_from_propagator,
)

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
OUT = ROOT / "data" / "audit" / "ctm_tau_sweep"
FIG = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "ctm_tau_sweep"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["alpha", "beta", "low_gamma"]
N_U = 8
BAND_TEX = {"alpha": r"$\alpha$", "beta": r"$\beta$", "low_gamma": r"$\gamma_l$"}

PHASE_TAGS = ["TT", "RPost", "RPreA", "RPreB", "RPostA", "RPostB"]


def find_tau_min_star(eigenvalues, n_grid=600, pad_low=0.5, pad_high=3.0):
    """Return tau_min = 1/lambda_max and tau* = argmax C(tau) over a wide grid.

    The standard [tau_min, 1/lambda_2] window is too narrow for our
    weight-heterogeneous fully-connected (continuous spectrum) case --
    C(tau) often peaks at tau > 1/lambda_2. Take the global argmax over a
    grid extended pad_high decades past 1/lambda_2, with sanity check that
    the peak is interior (not at the grid boundary).

    Returns dict with tau_min, tau_star, tau_grid, S, C, lam2, lam_max,
    plus a 'boundary_hit' flag set True if argmax landed at the upper
    grid edge (suggests we should pad further).
    """
    eig = np.maximum(np.asarray(eigenvalues, dtype=float), 0.0)
    lam2 = float(eig[1]) if eig[1] > 1e-10 else 1e-6
    lam_max = float(eig[-1])
    tau_min = 1.0 / lam_max
    tau_ref = 1.0 / lam2
    log_tau = np.linspace(
        np.log10(tau_min) - pad_low,
        np.log10(tau_ref) + pad_high,
        n_grid,
    )
    tau_grid = 10 ** log_tau
    S = np.array([_entropy_at_tau(eig, t) for t in tau_grid])
    C = -np.gradient(S, log_tau)
    # Restrict argmax to tau >= tau_min so we don't pick the low-tau noise
    valid = tau_grid >= tau_min
    idx_in = np.argmax(C[valid])
    tau_star = float(tau_grid[valid][idx_in])
    boundary_hit = bool(idx_in == valid.sum() - 1)
    return dict(tau_min=tau_min, tau_star=tau_star,
                tau_grid=tau_grid, S=S, C=C, lam2=lam2, lam_max=lam_max,
                boundary_hit=boundary_hit)


def D_at_tau(eigenvalues, eigenvectors, tau):
    K = _compute_propagator(eigenvalues, eigenvectors, tau)
    return _ultrametric_from_propagator(K)


def load_phases(pat, band):
    """Load all 6 LRG phases needed for the τ-sweep CTM."""
    out = {}
    out["TT"] = load_lrg_result(pat, "task_test", band, "imcoh_abs")
    out["RPost"] = load_lrg_result(pat, "rest_post", band, "imcoh_abs")
    out["RPreA"] = load_lrg_result(pat, "rest_pre_A", band, "imcoh_abs",
                                   cache_root=HALVES_CACHE)
    out["RPreB"] = load_lrg_result(pat, "rest_pre_B", band, "imcoh_abs",
                                   cache_root=HALVES_CACHE)
    out["RPostA"] = load_lrg_result(pat, "rest_post_A", band, "imcoh_abs",
                                    cache_root=HALVES_CACHE)
    out["RPostB"] = load_lrg_result(pat, "rest_post_B", band, "imcoh_abs",
                                    cache_root=HALVES_CACHE)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    rows = []
    tau_rows = []
    c_tau_records = []  # for sanity figure
    u_grid = np.linspace(0.0, 1.0, N_U)

    for pat in PATIENTS:
        print(f"=== {pat} ===")
        for band in BANDS:
            phases = load_phases(pat, band)
            if any(p is None or p.eigenvalues is None or p.eigenvectors is None
                   for p in phases.values()):
                missing = [k for k, v in phases.items()
                           if v is None or v.eigenvalues is None]
                print(f"  [{band}] skip — missing {missing}")
                continue

            # Per-phase intrinsic [tau_min, tau*]
            tau_info = {tag: find_tau_min_star(phases[tag].eigenvalues)
                        for tag in PHASE_TAGS}
            for tag, info in tau_info.items():
                tau_rows.append(dict(patient=pat, band=band, phase=tag,
                                     tau_min=info["tau_min"],
                                     tau_star=info["tau_star"],
                                     ratio_log10=np.log10(
                                         info["tau_star"] / info["tau_min"]),
                                     lambda_2=info["lam2"],
                                     lambda_max=info["lam_max"],
                                     boundary_hit=int(info["boundary_hit"])))
                c_tau_records.append(dict(
                    patient=pat, band=band, phase=tag,
                    tau=info["tau_grid"], S=info["S"], C=info["C"],
                    tau_min=info["tau_min"], tau_star=info["tau_star"]))

            N = phases["TT"].eigenvalues.shape[0]
            iu = np.triu_indices(N, k=1)

            for ui, u in enumerate(u_grid):
                tau_at_u = {
                    tag: tau_info[tag]["tau_min"]
                    * (tau_info[tag]["tau_star"] / tau_info[tag]["tau_min"]) ** u
                    for tag in PHASE_TAGS
                }
                D = {tag: D_at_tau(phases[tag].eigenvalues,
                                   phases[tag].eigenvectors,
                                   tau_at_u[tag])
                     for tag in PHASE_TAGS}

                d_task = (D["TT"] - D["RPreA"])[iu]
                d_rest = (D["RPost"] - D["RPreB"])[iu]
                # Official ρ_null_drift: within-session half-difference, RPre vs RPost
                d_drift_a = (D["RPreB"] - D["RPreA"])[iu]
                d_drift_b = (D["RPostB"] - D["RPostA"])[iu]

                rs, _ = spearmanr(d_task, d_rest)
                rd, _ = spearmanr(d_drift_a, d_drift_b)

                row = dict(patient=pat, band=band, u_idx=int(ui), u=float(u),
                           rho_split=float(rs) if not np.isnan(rs) else 0.0,
                           rho_drift=float(rd) if not np.isnan(rd) else 0.0,
                           gap=float(rs - rd) if not (
                               np.isnan(rs) or np.isnan(rd)) else 0.0)
                for tag in PHASE_TAGS:
                    row[f"tau_{tag}"] = tau_at_u[tag]
                rows.append(row)
            print(f"  [{band}] done")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tau_sweep_rho.csv", index=False)
    pd.DataFrame(tau_rows).to_csv(OUT / "tau_star_per_pb.csv", index=False)

    # Cohort plot
    fig, axes = plt.subplots(1, len(BANDS), figsize=(13.5, 4), sharey=True)
    for ax, band in zip(axes, BANDS):
        sub = df[df["band"] == band]
        for pat in PATIENTS:
            psub = sub[sub["patient"] == pat].sort_values("u")
            if not psub.empty:
                ax.plot(psub["u"], psub["rho_split"], color="#9ec6e5",
                        lw=0.7, alpha=0.55)
                ax.plot(psub["u"], psub["rho_drift"], color="#cccccc",
                        lw=0.5, alpha=0.45, ls="--")
        agg = (sub.groupby("u_idx")
                  .agg(med_split=("rho_split", "median"),
                       med_drift=("rho_drift", "median"),
                       q25_split=("rho_split", lambda s: float(np.nanpercentile(s, 25))),
                       q75_split=("rho_split", lambda s: float(np.nanpercentile(s, 75))),
                       u=("u", "first"))
                  .reset_index())
        ax.fill_between(agg["u"], agg["q25_split"], agg["q75_split"],
                        color="#1a4f73", alpha=0.18, lw=0)
        ax.plot(agg["u"], agg["med_split"], color="#1a4f73", lw=2.2,
                label=r"$\rho_{\mathrm{split}}$ (cohort median)")
        ax.plot(agg["u"], agg["med_drift"], color="#444", lw=1.4, ls="--",
                label=r"$\rho_{\mathrm{drift}}$ (cohort median)")
        ax.axhline(0, color="0.5", lw=0.6, ls=":")
        ax.set_xlabel(r"$u = (\log\tau - \log\tau_{\min})/(\log\tau^* - \log\tau_{\min})$")
        ax.set_title(BAND_TEX[band])
        if band == "alpha":
            ax.set_ylabel(r"$\rho$ (Spearman, triu-pair-shift agreement)")
            ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "cohort.pdf")
    plt.close(fig)

    # Per-patient small-multiples
    fig, axes = plt.subplots(len(PATIENTS), len(BANDS),
                              figsize=(10, 17), sharex=True, sharey=True)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            psub = df[(df["patient"] == pat) & (df["band"] == band)].sort_values("u")
            if not psub.empty:
                ax.plot(psub["u"], psub["rho_split"], color="#1a4f73", lw=1.2,
                        marker="o", ms=3)
                ax.plot(psub["u"], psub["rho_drift"], color="#888", lw=0.9,
                        ls="--")
                ax.axhline(0, color="0.5", lw=0.5, ls=":")
            if i == 0:
                ax.set_title(BAND_TEX[band])
            if j == 0:
                ax.set_ylabel(pat, fontsize=8, rotation=0, ha="right", va="center")
            if i == len(PATIENTS) - 1:
                ax.set_xlabel(r"$u$")
    fig.tight_layout()
    fig.savefig(FIG / "per_patient.pdf")
    plt.close(fig)

    # C(τ) sanity figure: 4 representative cells
    sample_cells = [("Pat_02", "beta"), ("Pat_06", "alpha"),
                    ("Pat_05", "low_gamma"), ("Pat_15", "beta")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    show_phases = ["TT", "RPost", "RPreA", "RPostA"]
    color_map = {"TT": "#1a4f73", "RPost": "#d62728",
                 "RPreA": "#888", "RPostA": "#2ca02c"}
    for ax, (pat_s, band_s) in zip(axes.flat, sample_cells):
        any_drawn = False
        for rec in c_tau_records:
            if (rec["patient"] == pat_s and rec["band"] == band_s
                    and rec["phase"] in show_phases):
                ax.semilogx(rec["tau"], rec["C"], lw=1.0,
                            color=color_map[rec["phase"]],
                            label=rec["phase"], alpha=0.9)
                ax.axvline(rec["tau_star"], color=color_map[rec["phase"]],
                           lw=0.5, ls="--", alpha=0.6)
                any_drawn = True
        if any_drawn:
            ax.axvspan(
                10 ** np.log10(min(r["tau_min"] for r in c_tau_records
                                   if r["patient"] == pat_s and r["band"] == band_s
                                   and r["phase"] in show_phases)),
                10 ** np.log10(max(r["tau_star"] for r in c_tau_records
                                   if r["patient"] == pat_s and r["band"] == band_s
                                   and r["phase"] in show_phases)),
                color="#fff7e6", alpha=0.4, zorder=-1)
        ax.set_title(f"{pat_s}  {BAND_TEX[band_s]}")
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(r"$C(\tau) = -dS/d\log\tau$")
        ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(FIG / "c_tau_sanity.pdf")
    plt.close(fig)

    # Per-band cohort summary print
    print("\n=== per-band cohort summary ===")
    for band in BANDS:
        sub = df[df["band"] == band]
        agg = (sub.groupby("u_idx")
                  .agg(med_split=("rho_split", "median"),
                       med_drift=("rho_drift", "median"),
                       n_pos=("rho_split", lambda s: int((s > 0).sum())),
                       u=("u", "first"))
                  .reset_index())
        print(f"\n[{band}]")
        print(agg.to_string(index=False))

    print(f"\nWrote outputs to {OUT}\nFigures to {FIG}")


if __name__ == "__main__":
    main()
