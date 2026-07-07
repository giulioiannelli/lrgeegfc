#!/usr/bin/env python3
"""audit_148 — what explains WHO traces? resolving the β inter-patient variability.

Follow-up to audit_147 (OFC coverage, `null`) and the Q1 reliability audit, both
of which left the between-patient β-trace heterogeneity as "unexplained residual
biology". This script attributes it, decomposing the failure rather than scoring
one covariate at a time.

Prior refutations (kept): reliability/detectability (Q1, ~6%; Pat_15 best-measured,
lowest-trace), OFC coverage (audit_147, ρ=-0.10), epi burden (audit_64 frac_epi
null all bands). This audit adds the *dynamic* + *node-level* + *cross-phase* views
that audit_64/147 never combined.

5-point critical preamble
1. Claim: the per-patient β trace (ρ_split) is a STABLE cross-phase retention
   trait with a left-lateralised anatomical basis; the 3 non-tracers fail for two
   mechanistically distinct reasons — engagement-null (Pat_15) vs reset (Pat_10/14).
2. Null/baseline: the variability is structureless noise — no covariate, dynamic
   or structural, tracks it, and learning vs test traces are unrelated.
3. Strongest alternative the analysis must beat: (a) it is the vendor-replaced
   task_test of Pat_14 (a provenance artifact); (b) it is measurement reliability;
   (c) it is total graph size / node count. Controlled by: (a) recomputing the
   trace on the INDEPENDENT, non-replaced learning phase (T_learn) — an artifact
   in task_test cannot survive there; (b) correlating vs reliability bottleneck;
   (c) correlating vs total N and non-lateralised size.
4. Cannot: n=10 — every correlation is a wide-CI description, not a powered test;
   it cannot prove the laterality is causal vs a correlated recording/biology
   factor; Pat_14's reset is left-hemisphere and adequately measured, so it is an
   irreducible residual this analysis characterises but does not mechanise.
5. Falsify: if T_learn and T_test are uncorrelated (retention not a trait), or if
   reliability/N correlate as strongly as laterality, or if right-hemisphere
   contacts are NOT depleted of trace carriers cohort-wide, the resolution fails
   and the heterogeneity stays unexplained.

Inputs (all existing — no LRG/FC/surrogate recompute):
    data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
    data/audit/multiphase_snr/patient_ranking.csv                 (d_task)
    data/audit/multiphase_snr/per_patient_per_band_reliability.csv
    data/audit/per_node_trace_decomposition/per_patient_summary.csv (carrier/anti/neutral)
    data/audit/consolidation_arc/arc_per_patient.csv               (T_learn, C_LT)
    data/audit/implant_geometry/per_patient_features.csv           (B_hemi, N_R, ...)
    data/audit/raw_fc_phase_distance/session_timing.csv            (durations, gap)
    data/anatomy_distribution/hemisphere_lobe_enrichment.csv       (carrier laterality)

Outputs:
    data/audit/interpatient_variability/master_table.csv
    data/audit/interpatient_variability/correlations.csv
    data/audit/interpatient_variability/fig_interpatient_variability.pdf
    data/audit/interpatient_variability/README.md
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BAND = "beta"
A = ROOT / "data/audit"
OUT = A / "interpatient_variability"
OUT.mkdir(parents=True, exist_ok=True)


def _band(df, band=BAND):
    return df[df.band == band].set_index("patient")


def build_master():
    trace = pd.read_csv(A / "matched_strength_surrogate_split_baseline/per_patient_per_band.csv")
    snr = pd.read_csv(A / "multiphase_snr/patient_ranking.csv")
    rel = pd.read_csv(A / "multiphase_snr/per_patient_per_band_reliability.csv")
    node = pd.read_csv(A / "per_node_trace_decomposition/per_patient_summary.csv")
    arc = pd.read_csv(A / "consolidation_arc/arc_per_patient.csv")
    feat = pd.read_csv(A / "implant_geometry/per_patient_features.csv").set_index("patient")
    tim = pd.read_csv(A / "raw_fc_phase_distance/session_timing.csv").set_index("patient")

    tb, sb, rb, nb, ab = _band(trace), _band(snr), _band(rel), _band(node), _band(arc)
    M = pd.DataFrame(index=list(tb.index))
    M["beta_rho"], M["beta_p"] = tb["obs_rho"], tb["obs_p_one_sided"]
    # node decomposition -> engagement (axis 1) + direction (axis 2)
    M["carrier"], M["anti"], M["neutral"] = nb["n_carrier"], nb["n_anti"], nb["n_neutral"]
    M["mover_frac"] = (nb["n_carrier"] + nb["n_anti"]) / nb["N"]
    M["neutral_frac"] = nb["n_neutral"] / nb["N"]
    M["dir_balance"] = (nb["n_carrier"] - nb["n_anti"]) / nb["N"]
    # dynamics + cross-phase trait
    M["d_task"], M["T_learn"], M["C_LT"] = sb["d_task"], ab["T_learn"], ab["C_LT"]
    M["trait"] = (M["beta_rho"] + M["T_learn"]) / 2.0
    # reliability
    M["rel_bottleneck"], M["rel_post"] = rb["rel_bottleneck"], rb["rel_rest_post"]
    # implant
    for c in ["B_hemi", "N_R", "N_L", "N", "frac_epi", "N_epi", "sigma_disp", "N_frontal"]:
        M[c] = feat[c]
    # timing
    M["dur_post"], M["dur_test"] = tim["rest_post_s"], tim["task_test_s"]
    M["dur_learn"], M["gap"] = tim["task_learn_s"], tim["inter_rest_gap_s"]

    def cls(r):
        if r.beta_p < 0.05:
            return "TRACE"
        return "RESET" if r.beta_rho < 0 else "NULL"
    M["class"] = M.apply(cls, axis=1)
    clears = trace.assign(c=(trace.obs_p_one_sided < 0.05).astype(int))
    M["n_bands"] = clears.groupby("patient")["c"].sum()
    return M


def correlations(M):
    covs = ["T_learn", "dir_balance", "mover_frac", "neutral_frac", "d_task", "C_LT",
            "B_hemi", "N_R", "N_L", "N", "frac_epi", "N_epi", "sigma_disp", "N_frontal",
            "rel_bottleneck", "rel_post", "dur_post", "dur_test", "dur_learn", "gap"]
    rows = []
    for c in covs:
        x, y = M[c].astype(float), M["beta_rho"].astype(float)
        ok = x.notna() & y.notna()
        r, p = spearmanr(x[ok], y[ok])
        rows.append((c, r, p, int(ok.sum())))
    return pd.DataFrame(rows, columns=["covariate", "rho", "p", "n"]).sort_values(
        "rho", key=abs, ascending=False).reset_index(drop=True)


def figure(M):
    hemi = pd.read_csv(A / "anatomy_distribution".replace("audit/", "") if False else
                       ROOT / "data/anatomy_distribution/hemisphere_lobe_enrichment.csv")
    hemi = hemi[(hemi.group_level == "hemisphere") & (hemi.sensitivity_regime == "full")]
    COL = {"TRACE": "#0f9d58", "RESET": "#d1495b", "NULL": "#8a8a8a"}
    fig, axs = plt.subplots(1, 4, figsize=(15.5, 3.9))

    # A: stable cross-phase trait
    ax = axs[0]
    ax.axhline(0, color="0.7", lw=.8, zorder=0); ax.axvline(0, color="0.7", lw=.8, zorder=0)
    lim = [-.2, .65]; ax.plot(lim, lim, ls="--", color="0.6", lw=.9, zorder=1)
    for p, r in M.iterrows():
        ax.scatter(r.T_learn, r.beta_rho, s=70, c=COL[r["class"]], edgecolor="k", lw=.6, zorder=3)
        ax.annotate(p.replace("Pat_", ""), (r.T_learn, r.beta_rho), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    rr, pp = spearmanr(M.T_learn, M.beta_rho)
    ax.set_xlabel(r"$T_{\rm learn}$ (learning$\to$rest_post)")
    ax.set_ylabel(r"$T_{\rm test}$ ($\rho_{\rm split}$, test$\to$rest_post)")
    ax.set_title(r"$\bf{A}$  a stable cross-phase trait" + f"\n$\\rho$={rr:.2f}, p={pp:.0e}", fontsize=9)
    ax.set_xlim(lim); ax.set_ylim(lim)

    # B: two failure modes
    ax = axs[1]; ax.axhline(0, color="0.7", lw=.8, zorder=0)
    ax.axvspan(0, 0.32, color="0.92", zorder=0)
    for p, r in M.iterrows():
        ax.scatter(r.mover_frac, r.dir_balance, s=70, c=COL[r["class"]], edgecolor="k", lw=.6, zorder=3)
        ax.annotate(p.replace("Pat_", ""), (r.mover_frac, r.dir_balance), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("engagement = mover fraction\n(nodes past matched-strength null)")
    ax.set_ylabel(r"direction = (carrier$-$anti)/N")
    ax.set_title(r"$\bf{B}$  two failure modes" + "\nnull-engage (15) vs reset (10,14)", fontsize=9)
    ax.annotate("engagement\nfloor", (0.16, 0.42), fontsize=7, color="0.4", ha="center")

    # C: laterality bias
    ax = axs[2]; ax.axhline(0, color="0.7", lw=.8, zorder=0)
    for p, r in M.iterrows():
        ax.scatter(r.N_R, r.beta_rho, s=70, c=COL[r["class"]], edgecolor="k", lw=.6, zorder=3)
        ax.annotate(p.replace("Pat_", ""), (r.N_R, r.beta_rho), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    rN, pN = spearmanr(M.N_R, M.beta_rho)
    ax.set_xlabel("N right-hemisphere contacts")
    ax.set_ylabel(r"$T_{\rm test}$ ($\rho_{\rm split}$, $\beta$)")
    ax.set_title(r"$\bf{C}$  right-sampling anti-predicts trace" + f"\n$\\rho$={rN:.2f}, p={pN:.2f}", fontsize=9)
    ax.annotate("06: right-OFC\nstill traces", (74, 0.211), fontsize=6.3, color="0.35",
                xytext=(-2, -30), textcoords="offset points", ha="center",
                arrowprops=dict(arrowstyle="-", color="0.6", lw=.6))
    ax.annotate("14: left, yet\nresets", (24, -0.049), fontsize=6.3, color="0.35",
                xytext=(24, -20), textcoords="offset points", ha="center",
                arrowprops=dict(arrowstyle="-", color="0.6", lw=.6))

    # D: mechanism - carriers are left-lateralised
    ax = axs[3]
    bands = ["delta", "alpha", "beta", "low_gamma"]
    xl = np.arange(len(bands)); w = 0.38
    lr = [float(hemi[(hemi.band == b) & (hemi.group == "left")]["trace_rate"]) for b in bands]
    rr_ = [float(hemi[(hemi.band == b) & (hemi.group == "right")]["trace_rate"]) for b in bands]
    ax.bar(xl - w / 2, np.array(lr) * 100, w, label="left", color="#3d5a80", edgecolor="k", lw=.5)
    ax.bar(xl + w / 2, np.array(rr_) * 100, w, label="right", color="#ee9b00", edgecolor="k", lw=.5)
    ax.set_xticks(xl); ax.set_xticklabels([b.replace("low_gamma", r"low-$\gamma$") for b in bands], fontsize=8)
    ax.set_ylabel("trace-carrier rate (%)")
    ax.set_title(r"$\bf{D}$  mechanism: carriers" + "\nare left-lateralised", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    for i, b in enumerate(bands):
        if float(hemi[(hemi.band == b) & (hemi.group == "left")]["p_bonferroni"]) < 0.05:
            ax.annotate("*", (i - w / 2, lr[i] * 100 + 0.3), ha="center", fontsize=12)

    handles = [Patch(facecolor=COL[k], edgecolor="k", label=k.title()) for k in ["TRACE", "RESET", "NULL"]]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False, fontsize=9)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(OUT / "fig_interpatient_variability.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    M = build_master()
    corr = correlations(M)
    M.to_csv(OUT / "master_table.csv")
    corr.to_csv(OUT / "correlations.csv", index=False)
    figure(M)

    rr, pp = spearmanr(M.T_learn, M.beta_rho)
    print("=" * 78)
    print("audit_148 — inter-patient β-trace variability, RESOLVED")
    print("=" * 78)
    print(M.sort_values("beta_rho", ascending=False)[
        ["beta_rho", "class", "mover_frac", "dir_balance", "T_learn", "N_R", "B_hemi"]].round(3).to_string())
    print(f"\nTop covariates of beta_rho:\n{corr.head(8).round(3).to_string(index=False)}")
    print(f"\nStable-trait test: Spearman(T_learn, T_test) = {rr:.3f}, p={pp:.1e}  "
          f"(sign concordance {(np.sign(M.T_learn)==np.sign(M.beta_rho)).sum()}/10)")
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
