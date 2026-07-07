#!/usr/bin/env python
"""audit_142 — Cross-phase LRG thermodynamic observables (FIRST-PASS SIGNAL PROBE).

HEAD
----
Does learning a relational (transitive-inference) task reorganize the brain's
INTRINSIC multiscale architecture — shifting its characteristic integration
scale (the C(tau) specific-heat peak) and/or changing its effective
dimensionality (the von Neumann entropy S(tau)) — and does that change PERSIST
into post-task rest?  This mines the LRG framework's NATIVE thermodynamic
observables cross-phase (Villegas et al., Nat. Phys. 2023), which the team's
cophenetic/Grassmann results have never touched.  STATIC (phase-averaged)
measure, so it sidesteps the timescale problem.

5-POINT CRITICAL PREAMBLE (per house rule)
------------------------------------------
1. CLAIM (to PROBE, not confirm): task learning shifts C(tau) peak location
   and/or S(tau) effective-dimensionality, and the shift persists rest_pre ->
   rest_post, esp. in alpha/beta.
2. NULL: none run here. This is the compute+plot+evaluate first pass. Per
   house rule we do NOT pre-register an acceptance threshold and do NOT run the
   matched-strength null yet (later, only if observed signal exists).
3. STRONGEST ALTERNATIVE the eventual null must control: electrode count N and
   raw connectivity STRENGTH differ across phases; both move S(tau)/C(tau)
   trivially. We sidestep N by NEVER comparing absolute peak locations across
   patients (each patient is its own control; only WITHIN-patient cross-phase
   deltas). Strength is NOT controlled here -> any positive signal is
   UNVERIFIED until matched-strength.
4. WHAT THIS PASS CANNOT DO: it cannot distinguish a genuine architectural
   reorganization from a strength-driven one, and cannot reject the null
   (no null run). It only answers: is there ANY consistent cohort cross-phase
   direction worth taking to the null stage?
5. FALSIFY: no consistent sign across the cohort (n_consistent ~ 5/10) on the
   persistence contrast for the primary observables/bands => dead end.

OBSERVABLES (per patient, per phase, per band)
----------------------------------------------
Recomputed from `eigenvalues` via the canonical
`workflow.diagnostics.compute_entropy_curve` on the per-patient resolution
window [1/lambda_max, 1/lambda_2] (denser/cleaner than the legacy cached
fixed-grid `entropy_*`; the cached `entropy_C` is on a fixed [1e-3,1e5] grid
that does not align with the resolution window).

  log10_tau_star      : log10(tau) at the dominant C(tau) peak  [characteristic scale]
  log10_tau_star_norm : log10(tau* / tau_min) = log10(tau* * lambda_max)
                        [dimensionless within-patient characteristic scale]
  C_max               : max of C(tau) in window  [sharpness of scale separation]
  N_peaks             : # prominent C(tau) peaks
  S_mean              : mean of S(tau)/ln(N) over the resolution window
                        [effective-dimensionality proxy; HIGH = more dims/less integrated]
  integration         : 1 - S_mean  [HIGH = more integrated]
  S_at_tau_min        : S(tau)/ln(N) at the fine canonical scale tau_min
                        (recorded for completeness; saturates ~1, weak discriminator)

CROSS-PHASE CONTRASTS (sign convention: positive = "moves like the task / persists")
------------------------------------------------------------------------------------
  encode      : X(task_learn) - X(rest_pre)
  inference   : X(task_test)  - X(task_learn)
  persist     : X(rest_post)  - X(rest_pre)        [KEY]
  arc         : 1 if |X(rest_post)-X(task_test)| < |X(rest_pre)-X(task_test)|
                else 0  -> rest_post moved TOWARD the task value [KEY, count /10]

Cohort: per band x observable x contrast -> median delta, n_consistent /10,
Wilcoxon (one-sided 'greater' from lrg_eegfc + two-sided from scipy). Lead with
direction + consistency; p secondary (house rule).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.diagnostics import (
    compute_entropy_curve,
    compute_susceptibility_diagnostics,
)
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z
from lrg_eegfc.visuals.styles import use_lrg_style

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
BANDS = ["alpha", "beta", "delta", "theta", "low_gamma", "high_gamma"]
FC = "imcoh_abs"

OBS = ["log10_tau_star", "log10_tau_star_norm", "C_max",
       "N_peaks", "S_mean", "integration", "S_at_tau_min"]


# ---------------------------------------------------------------------------
def per_phase_scalars(eigenvalues: np.ndarray) -> dict:
    """All per-phase observables from the eigenvalue spectrum."""
    ev = np.maximum(np.asarray(eigenvalues, dtype=float), 0.0)
    N = len(ev)
    ec = compute_entropy_curve(ev)
    sd = compute_susceptibility_diagnostics(ev)

    Snorm = ec["S"] / np.log(N)
    in_win = (ec["tau"] >= ec["tau_min"]) & (ec["tau"] <= ec["tau_max"])
    S_mean = float(Snorm[in_win].mean())
    i_tmin = int(np.argmin(np.abs(ec["tau"] - ec["tau_min"])))
    S_at_tau_min = float(Snorm[i_tmin])

    tau_star = sd["tau_star1"]
    log10_tau_star = float(np.log10(tau_star)) if tau_star > 0 else np.nan
    log10_tau_star_norm = (
        float(np.log10(tau_star * ec["lambda_max"]))
        if tau_star > 0 else np.nan
    )
    return dict(
        N=N,
        log10_tau_star=log10_tau_star,
        log10_tau_star_norm=log10_tau_star_norm,
        C_max=float(sd["C_max"]),
        N_peaks=int(sd["N_peaks"]),
        S_mean=S_mean,
        integration=1.0 - S_mean,
        S_at_tau_min=S_at_tau_min,
    )


def collect_curves_and_scalars():
    """Returns (rows_df, curves) where curves[(pat,band,phase)] holds the
    common-grid S/C for figures."""
    rows = []
    curves = {}
    for pat in PATIENTS:
        for band in BANDS:
            for ph in PHASES:
                try:
                    r = load_lrg_result(pat, ph, band, FC)
                except Exception as e:  # noqa: BLE001
                    print(f"  MISS {pat} {band} {ph}: {e}")
                    continue
                ev = np.maximum(np.asarray(r.eigenvalues, float), 0.0)
                sc = per_phase_scalars(ev)
                sc.update(patient=pat, band=band, phase=ph)
                rows.append(sc)
                ec = compute_entropy_curve(ev)
                curves[(pat, band, ph)] = dict(
                    log10_tau=ec["log10_tau"],
                    Snorm=ec["S"] / np.log(len(ev)),
                    C=ec["C"],
                    tau_min=ec["tau_min"], tau_max=ec["tau_max"],
                )
    return pd.DataFrame(rows), curves


# ---------------------------------------------------------------------------
def cross_phase_contrasts(df: pd.DataFrame) -> pd.DataFrame:
    """Per (patient, band) -> the 4 contrasts for each observable."""
    out = []
    for (pat, band), g in df.groupby(["patient", "band"]):
        gp = g.set_index("phase")
        if not all(ph in gp.index for ph in PHASES):
            continue
        for obs in OBS:
            x_pre = gp.loc["rest_pre", obs]
            x_lrn = gp.loc["task_learn", obs]
            x_tst = gp.loc["task_test", obs]
            x_pst = gp.loc["rest_post", obs]
            arc = int(abs(x_pst - x_tst) < abs(x_pre - x_tst))
            out.append(dict(
                patient=pat, band=band, observable=obs,
                encode=x_lrn - x_pre,
                inference=x_tst - x_lrn,
                persist=x_pst - x_pre,
                arc=arc,
            ))
    return pd.DataFrame(out)


def cohort_table(contrasts: pd.DataFrame) -> pd.DataFrame:
    """Per band x observable x contrast -> median, n_consistent/10, Wilcoxon."""
    out = []
    for band in BANDS:
        for obs in OBS:
            sub = contrasts[(contrasts.band == band)
                            & (contrasts.observable == obs)]
            n = len(sub)
            if n == 0:
                continue
            # arc is a binary count, treat separately
            for con in ["encode", "inference", "persist"]:
                x = sub[con].to_numpy(dtype=float)
                x = x[np.isfinite(x)]
                med = float(np.median(x)) if len(x) else np.nan
                n_pos = int((x > 0).sum())
                n_neg = int((x < 0).sum())
                # consistency = patients agreeing with the SIGN OF THE MEDIAN
                if med > 0:
                    n_cons = n_pos
                elif med < 0:
                    n_cons = n_neg
                else:
                    n_cons = max(n_pos, n_neg)
                # one-sided 'greater' (positive direction = library helper)
                z_g, p_g = wilcoxon_z(x)
                # two-sided
                xnz = x[x != 0.0]
                if len(xnz) >= 3:
                    try:
                        p_two = float(stats.wilcoxon(
                            xnz, alternative="two-sided",
                            method="approx").pvalue)
                    except Exception:  # noqa: BLE001
                        p_two = np.nan
                else:
                    p_two = np.nan
                out.append(dict(
                    band=band, observable=obs, contrast=con,
                    n=len(x), median=med,
                    n_pos=n_pos, n_neg=n_neg, n_consistent=n_cons,
                    wilcox_z_greater=z_g, p_one_sided_greater=p_g,
                    p_two_sided=p_two,
                ))
            # arc: binary fraction moving toward task
            a = sub["arc"].to_numpy(dtype=float)
            a = a[np.isfinite(a)]
            n_toward = int(a.sum())
            # binomial test vs 0.5
            try:
                p_bin = float(stats.binomtest(
                    n_toward, len(a), 0.5, alternative="greater").pvalue)
            except Exception:  # noqa: BLE001
                p_bin = np.nan
            out.append(dict(
                band=band, observable=obs, contrast="arc_toward_task",
                n=len(a), median=float(np.mean(a)),
                n_pos=n_toward, n_neg=len(a) - n_toward,
                n_consistent=n_toward,
                wilcox_z_greater=np.nan, p_one_sided_greater=np.nan,
                p_two_sided=p_bin,
            ))
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
def make_figures(curves: dict, outdir):
    """Per-band cohort-mean S(tau) and C(tau), one line per phase, on a
    within-patient z-aligned common grid."""
    use_lrg_style()
    phase_colors = {
        "rest_pre": "#4477AA", "task_learn": "#EE6677",
        "task_test": "#CCBB44", "rest_post": "#228833",
    }
    # Common grid in NORMALIZED log10(tau): align each patient on
    # log10(tau) - log10(tau_min)  (so 0 = fine canonical scale tau_min).
    grid = np.linspace(-0.5, 1.5, 200)  # tau_min..~tau_max span is ~0.6 dec

    for band in BANDS:
        fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
        for ph in PHASES:
            S_stack, C_stack = [], []
            for pat in PATIENTS:
                key = (pat, band, ph)
                if key not in curves:
                    continue
                c = curves[key]
                x = c["log10_tau"] - np.log10(c["tau_min"])
                S_i = np.interp(grid, x, c["Snorm"], left=np.nan, right=np.nan)
                C_i = np.interp(grid, x, c["C"], left=np.nan, right=np.nan)
                S_stack.append(S_i)
                C_stack.append(C_i)
            if not S_stack:
                continue
            S_m = np.nanmean(np.vstack(S_stack), axis=0)
            C_m = np.nanmean(np.vstack(C_stack), axis=0)
            axes[0].plot(grid, S_m, color=phase_colors[ph], label=ph, lw=1.6)
            axes[1].plot(grid, C_m, color=phase_colors[ph], label=ph, lw=1.6)
        axes[0].set_xlabel(r"$\log_{10}(\tau/\tau_{\min})$")
        axes[0].set_ylabel(r"$S(\tau)/\ln N$  (cohort mean)")
        axes[1].set_xlabel(r"$\log_{10}(\tau/\tau_{\min})$")
        axes[1].set_ylabel(r"$C(\tau)=-dS/d\log\tau$  (cohort mean)")
        for ax in axes:
            ax.axvline(0.0, color="0.6", lw=0.8, ls=":")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center",
                   bbox_to_anchor=(0.5, -0.04), ncol=4, frameon=False)
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        fp = outdir / f"entropy_specheat_{band}.pdf"
        fig.savefig(fp, bbox_inches="tight")
        plt.close(fig)
        print(f"  fig -> {fp}")


# ---------------------------------------------------------------------------
def main():
    from pathlib import Path
    repo_root = Path(__file__).resolve().parents[3]
    outdir = repo_root / "data" / "audit" / "entropy_crossphase"
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    print("[1/4] collecting curves + scalars ...")
    df, curves = collect_curves_and_scalars()
    df.to_csv(outdir / "per_phase_scalars.csv", index=False)

    print("[2/4] cross-phase contrasts ...")
    contrasts = cross_phase_contrasts(df)
    contrasts.to_csv(outdir / "per_patient_contrasts.csv", index=False)

    print("[3/4] cohort table ...")
    cohort = cohort_table(contrasts)
    cohort = cohort.sort_values(
        ["contrast", "band", "observable"]).reset_index(drop=True)
    cohort.to_csv(outdir / "cohort_summary.csv", index=False)

    print("[4/4] figures ...")
    make_figures(curves, figdir)

    # ---- terminal summary: lead with persistence + arc in alpha/beta ----
    print("\n==== KEY CELLS (persist + arc, alpha/beta) ====")
    key = cohort[
        (cohort.contrast.isin(["persist", "arc_toward_task"]))
        & (cohort.band.isin(["alpha", "beta"]))
    ].copy()
    with pd.option_context("display.max_rows", None,
                           "display.width", 200):
        print(key[["band", "observable", "contrast", "n",
                   "median", "n_consistent",
                   "p_one_sided_greater", "p_two_sided"]].to_string(index=False))

    print(f"\nWrote: {outdir}/cohort_summary.csv")
    return cohort


if __name__ == "__main__":
    main()
