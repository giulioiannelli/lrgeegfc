#!/usr/bin/env python3
"""Audit 91 — eigenmode localization: does the task-trace live in localized
modes (Q1), and do localized modes trap on epileptic nodes (Q2)?

Scope report (with full 5-point preamble):
``.agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md``.

All nulls are the degree/strength-preserving matched-strength surrogate whose
Laplacian eigendecomposition is ALREADY cached (audit_63/66 family); the
observed eigenvectors are one ``eigh(L = D − W)`` per (patient, band, phase),
identical in construction to the surrogate eigendecomposition, so observed and
null PR / mode-mass are directly comparable.

Critical preamble (per CLAUDE.md rule — full version in the scope report)
========================================================================
(1) Claim. (S) Observed FC-Laplacian modes are more localized than a
    strength-matched graph. (Q1) The α per-pair cophenetic trace lives in
    LOCALIZED, NON-LEADING modes while the leading (extended) subspace is
    unchanged → explains the α Grassmann-null / cophenet-positive dissociation.
    (Q2) Localized modes concentrate mass on epileptic nodes beyond strength
    (Anderson-trap analogue).
(2) Null. Matched-strength surrogate (R=200, cached eigvecs). Same epi indices
    on surrogate eigvecs for Q2 (preserves the epi set's strength profile).
(3) Strongest alternative. High-strength nodes accumulate mode mass with no
    localization/trapping mechanism; Q2 also confounded by sEEG-shaft spatial
    autocorrelation.
(4) Mechanical reach. The surrogate preserves per-node strength EXACTLY, so an
    observed PR below null / m^E above null cannot be a degree re-read. It does
    NOT control electrode-shaft geometry (Q2) → a negative Q2 is conclusive, a
    positive Q2 is only a candidate. Q1 is answered by observed-only geometry
    (leading-vs-bulk PR, PR-vs-displacement rank corr); the null only certifies
    the modes are localized at all.
(5) Falsification. Q1 dies if leading modes are not more extended than bulk OR
    displacement does not concentrate on localized modes. Q2 dies if m^E does
    not beat the strength null cohort-wide → report clean negative.

Outputs
-------
``data/audit/eigenmode_localization/``
    localization_per_cell.csv   one row per (patient, band, phase)
    q1_displacement_per_band.csv per (band) Q1 leading/bulk + PR-displacement
    q2_epi_mass_per_band.csv     per (band) Q2 epi mode-mass vs strength null
    README.md
    figures/fig_01_pr_obs_vs_null.pdf
    figures/fig_02_q1_pr_vs_displacement.pdf
    figures/fig_03_q2_epi_mass.pdf
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.metrics.hypothesis import surrogate_p_value
from lrg_eegfc.utils.metrics.spectral import (
    node_set_mode_mass,
    participation_number,
    subspace_displacement,
)
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
)
import _epi_stratify as es  # type: ignore


OUT = ROOT / "data" / "audit" / "eigenmode_localization"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Q1 leading-subspace cut: K slowest nontrivial modes (the block the
# leading-subspace Grassmann probe lives on). Reported across a few K.
K_LEAD = (5, 10, 20)
K_LEAD_PRIMARY = 10


# ---------------------------------------------------------------------------
# Eigendecomposition helpers (observed == surrogate construction)
# ---------------------------------------------------------------------------
def obs_eigvecs(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(eigvals[N], eigvecs[N, N])`` of L = D − W, ascending. Col 0 = zero."""
    L = np.diag(W.sum(axis=1)) - W
    return np.linalg.eigh(L)


def _nontrivial(evals: np.ndarray, evecs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Drop the trivial constant zero-mode (column 0)."""
    return evals[1:], evecs[:, 1:]


# ---------------------------------------------------------------------------
# (S + Q1) per-cell localization + leading/bulk + epi mass
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, phase: str, n_surr: int, swap_factor: int,
             verbose: bool) -> dict | None:
    try:
        W = load_phase_fc(pat, phase, band)
    except Exception as e:
        if verbose:
            print(f"[audit_91] SKIP {pat}/{band}/{phase}: {e}")
        return None
    N = W.shape[0]
    ev, V = obs_eigvecs(W)
    lam, Vn = _nontrivial(ev, V)
    lam_hat = lam / ev[-1]
    pr = participation_number(Vn)                       # (N-1,)

    # observed localization summaries
    med_pr = float(np.nanmedian(pr))
    min_pr = float(np.nanmin(pr))
    # leading-vs-bulk (Q1 structural)
    lead_bulk = {}
    for K in K_LEAD:
        if Vn.shape[1] > K + 3:
            lead_bulk[f"pr_lead_K{K}"] = float(np.nanmedian(pr[:K]))
            lead_bulk[f"pr_bulk_K{K}"] = float(np.nanmedian(pr[K:]))
        else:
            lead_bulk[f"pr_lead_K{K}"] = np.nan
            lead_bulk[f"pr_bulk_K{K}"] = np.nan

    # epi mode-mass (Q2). NOTE: the mean over ALL modes of m^E_k is exactly
    # |E|/N for ANY orthonormal basis (rows are unit-norm: Σ_k v_k(i)²=1), so a
    # spectrum-averaged epi mass is basis-invariant and cannot detect a trap.
    # An Anderson-style trap is a property of the LOCALIZED modes: we therefore
    # measure (i) the epi mass of the most-localized decile of modes (where a
    # trap must show) and (ii) the single largest epi-mass mode — neither is
    # basis-invariant, so observed and strength-null genuinely differ.
    pm = build_epi_masks(pat)
    epi_ok = len(pm.epi_mask) == N
    epi = np.asarray(pm.epi_mask, dtype=bool) if epi_ok else None
    n_epi = int(epi.sum()) if epi_ok else 0
    if epi_ok and n_epi > 0:
        mE_modes = node_set_mode_mass(Vn, epi)          # (N-1,)
        loc_sel = pr <= np.nanpercentile(pr, 10)        # most-localized decile
        mE_loc_obs = float(np.nanmean(mE_modes[loc_sel])) if loc_sel.any() else np.nan
        mE_max_obs = float(np.nanmax(mE_modes))
        frac_epi = n_epi / N
        # SAME-SIZE RANDOM-NODE-SET control on the OBSERVED eigvecs: disambiguates
        # a genuine epi trap from the localization-extreme artifact (observed
        # modes are localized, so max-of-any-mass beats the extended null). If
        # epi m^E_max is unexceptional among random sets of size |E|, the
        # strength-null m^E_max excess is NOT epi-specific.
        rng_ctrl = np.random.default_rng(int(pat.split("_")[-1]) * 1000 + N)
        rand_max = np.array([
            float(node_set_mode_mass(Vn, rng_ctrl.choice(N, n_epi, replace=False)).max())
            for _ in range(200)])
        p_epi_max_rand = surrogate_p_value(mE_max_obs, rand_max, tail="upper")
    else:
        mE_loc_obs = np.nan
        mE_max_obs = np.nan
        frac_epi = np.nan
        p_epi_max_rand = np.nan

    # surrogate null (cached eigvecs)
    evals_s, evecs_s = load_or_compute_eigs_at_path(
        es.surr_eig_path("full", pat, band, phase, n_surr, swap_factor),
        W, n_surr, swap_factor,
        es.cell_rng(pat, band, phase, "full"), verbose=verbose)
    med_pr_surr = []
    mE_loc_surr = []
    mE_max_surr = []
    for r in range(evals_s.shape[0]):
        if not np.isfinite(evals_s[r]).all():
            continue
        _, Vsn = _nontrivial(evals_s[r], evecs_s[r])
        prs = participation_number(Vsn)
        med_pr_surr.append(np.nanmedian(prs))
        if epi_ok and n_epi > 0:
            mEs = node_set_mode_mass(Vsn, epi)
            sel = prs <= np.nanpercentile(prs, 10)
            mE_loc_surr.append(float(np.nanmean(mEs[sel])) if sel.any() else np.nan)
            mE_max_surr.append(float(np.nanmax(mEs)))
    med_pr_surr = np.asarray(med_pr_surr)
    mE_loc_surr = np.asarray(mE_loc_surr)
    mE_max_surr = np.asarray(mE_max_surr)

    # localization p: observed PR BELOW null (lower tail = more localized)
    p_loc = surrogate_p_value(med_pr, med_pr_surr, tail="lower")
    # epi-trap p: observed epi mass of localized modes ABOVE null (upper tail)
    p_epi_loc = (surrogate_p_value(mE_loc_obs, mE_loc_surr, tail="upper")
                 if epi_ok and n_epi > 0 else np.nan)
    p_epi_max = (surrogate_p_value(mE_max_obs, mE_max_surr, tail="upper")
                 if epi_ok and n_epi > 0 else np.nan)

    row = {
        "patient": pat, "band": band, "phase": phase, "N": N,
        "med_pr_obs": med_pr, "med_pr_surr": float(np.median(med_pr_surr))
        if med_pr_surr.size else np.nan,
        "pr_ratio_obs_over_surr": (med_pr / float(np.median(med_pr_surr)))
        if med_pr_surr.size and np.median(med_pr_surr) > 0 else np.nan,
        "min_pr_obs": min_pr, "p_localized": p_loc,
        "n_epi": n_epi, "frac_epi": frac_epi,
        "mE_loc_obs": mE_loc_obs,
        "mE_loc_surr": float(np.nanmedian(mE_loc_surr)) if mE_loc_surr.size else np.nan,
        "p_epi_trap_loc": p_epi_loc,
        "mE_max_obs": mE_max_obs,
        "mE_max_surr": float(np.nanmedian(mE_max_surr)) if mE_max_surr.size else np.nan,
        "p_epi_trap_max": p_epi_max,
        "p_epi_max_vs_random_nodeset": p_epi_max_rand,
        **lead_bulk,
    }
    return row


# ---------------------------------------------------------------------------
# (Q1) cross-phase PR-vs-displacement (observed-only geometry)
# ---------------------------------------------------------------------------
def q1_displacement(pat: str, band: str, verbose: bool) -> dict | None:
    """Per-patient Spearman(PR_rest_post, displacement(rest_pre→rest_post)).

    ρ > 0 ⇒ EXTENDED modes (high PR) displace more (trace carried by extended
    modes); ρ < 0 ⇒ LOCALIZED modes carry the cross-phase change (the Q1
    hypothesis). Uses the full rest_pre / rest_post FCs (not the split
    halves) — this is the cross-phase trace direction.
    """
    try:
        Wp = load_phase_fc(pat, "rest_pre", band)
        Wq = load_phase_fc(pat, "rest_post", band)
    except Exception as e:
        if verbose:
            print(f"[audit_91] Q1 SKIP {pat}/{band}: {e}")
        return None
    if Wp.shape != Wq.shape:
        return None
    _, Vp = obs_eigvecs(Wp)
    evq, Vq = obs_eigvecs(Wq)
    _, Vpn = _nontrivial(np.zeros(Vp.shape[0]), Vp)     # drop col 0
    lamq, Vqn = _nontrivial(evq, Vq)
    disp = subspace_displacement(Vpn, Vqn)              # per rest_post mode
    pr = participation_number(Vqn)
    rho, _ = spearmanr(pr, disp)
    return {
        "patient": pat, "band": band,
        "spearman_pr_disp": float(rho),
        "med_disp_localized": float(np.nanmedian(disp[pr <= np.nanmedian(pr)])),
        "med_disp_extended": float(np.nanmedian(disp[pr > np.nanmedian(pr)])),
    }


# ---------------------------------------------------------------------------
# Cohort aggregation
# ---------------------------------------------------------------------------
def cohort_q1(cell: pd.DataFrame, disp: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(cell.band)]:
        cb = cell[(cell.band == band) & (cell.phase == "rest_post")]
        db = disp[disp.band == band]
        row = {"band": band}
        # leading vs bulk (count patients where leading more localized)
        for K in K_LEAD:
            lead = cb[f"pr_lead_K{K}"].values
            bulk = cb[f"pr_bulk_K{K}"].values
            m = np.isfinite(lead) & np.isfinite(bulk)
            row[f"n_lead_more_localized_K{K}"] = int((lead[m] < bulk[m]).sum())
            row[f"n_defined_K{K}"] = int(m.sum())
        # PR-displacement
        rhos = db["spearman_pr_disp"].values
        rhos = rhos[np.isfinite(rhos)]
        row["median_spearman_pr_disp"] = float(np.median(rhos)) if rhos.size else np.nan
        row["n_rho_negative"] = int((rhos < 0).sum())   # localized carry trace
        row["n_rho_positive"] = int((rhos > 0).sum())   # extended carry trace
        row["n_patients"] = int(rhos.size)
        if rhos.size >= 3:
            try:
                # two-sided: is the median displacement-localization link != 0
                _, wp = wilcoxon(rhos, alternative="two-sided")
            except Exception:
                wp = np.nan
            row["wilcoxon_p_rho"] = float(wp)
        else:
            row["wilcoxon_p_rho"] = np.nan
        out.append(row)
    return pd.DataFrame(out)


def cohort_q2(cell: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(cell.band)]:
        cb = cell[(cell.band == band) & (cell.phase == "rest_post")]
        cb = cb[cb.n_epi > 0]
        row = {"band": band}
        for tag, ocol, scol, pcol in (
            ("loc", "mE_loc_obs", "mE_loc_surr", "p_epi_trap_loc"),
            ("max", "mE_max_obs", "mE_max_surr", "p_epi_trap_max"),
        ):
            p = cb[pcol].values
            p = p[np.isfinite(p)]
            mE_o = cb[ocol].values
            mE_s = cb[scol].values
            m = np.isfinite(mE_o) & np.isfinite(mE_s)
            if m.sum() >= 3:
                try:
                    _, wp = wilcoxon(mE_o[m] - mE_s[m], alternative="greater")
                except Exception:
                    wp = np.nan
            else:
                wp = np.nan
            row[f"n_patients_epi"] = int(m.sum())
            row[f"n_trap_{tag}_p_lt_0.05"] = int((p < 0.05).sum())
            row[f"median_mE_{tag}_obs"] = float(np.median(mE_o[m])) if m.any() else np.nan
            row[f"median_mE_{tag}_surr"] = float(np.median(mE_s[m])) if m.any() else np.nan
            row[f"cohort_wilcoxon_p_{tag}"] = float(wp)
        # random-node-set control (observed eigvecs): is epi m^E_max special vs
        # random sets of size |E|? If not, the strength-null m^E_max excess is
        # the localization-extreme artifact, not an epi trap.
        pr_rand = cb["p_epi_max_vs_random_nodeset"].values
        pr_rand = pr_rand[np.isfinite(pr_rand)]
        row["median_p_max_vs_random"] = float(np.median(pr_rand)) if pr_rand.size else np.nan
        row["n_max_special_vs_random_p_lt_0.05"] = int((pr_rand < 0.05).sum())
        out.append(row)
    return pd.DataFrame(out)


def cohort_localization(cell: pd.DataFrame) -> pd.DataFrame:
    """Structural claim: observed PR below strength null, per band (rest_post)."""
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(cell.band)]:
        cb = cell[(cell.band == band) & (cell.phase == "rest_post")]
        ratio = cb["pr_ratio_obs_over_surr"].values
        ratio = ratio[np.isfinite(ratio)]
        p = cb["p_localized"].values
        p = p[np.isfinite(p)]
        out.append({
            "band": band, "n_patients": int(ratio.size),
            "median_pr_ratio": float(np.median(ratio)) if ratio.size else np.nan,
            "n_localized_p_lt_0.05": int((p < 0.05).sum()),
        })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def make_figures(cell: pd.DataFrame, disp: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()

    bands = [b for b in es.ALL_BANDS if b in set(cell.band)]
    rp = cell[cell.phase == "rest_post"]

    # Fig 1: observed vs null median PR per band (paired per patient)
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    x = np.arange(len(bands))
    for j, band in enumerate(bands):
        cb = rp[rp.band == band]
        o = cb["med_pr_obs"].values
        s = cb["med_pr_surr"].values
        ax.scatter(np.full_like(o, x[j] - 0.13), o, s=14, color="C3",
                   alpha=0.8, zorder=3, label="observed" if j == 0 else None)
        ax.scatter(np.full_like(s, x[j] + 0.13), s, s=14, color="C0",
                   alpha=0.8, zorder=3, label="strength null" if j == 0 else None)
        for oi, si in zip(o, s):
            ax.plot([x[j] - 0.13, x[j] + 0.13], [oi, si], color="0.7",
                    lw=0.5, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                       rotation=30, ha="right")
    ax.set_ylabel(r"median participation number $\overline{PR}$")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_01_pr_obs_vs_null.pdf", transparent=True)
    plt.close(fig)

    # Fig 2: Q1 — per-patient Spearman(PR, displacement) per band
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    for j, band in enumerate(bands):
        db = disp[disp.band == band]
        r = db["spearman_pr_disp"].values
        r = r[np.isfinite(r)]
        ax.scatter(np.full_like(r, x[j]) + np.random.uniform(-0.08, 0.08, r.size),
                   r, s=16, color="C2", alpha=0.8, zorder=3)
        ax.plot([x[j] - 0.2, x[j] + 0.2], [np.median(r)] * 2, color="k",
                lw=1.4, zorder=4)
    ax.axhline(0.0, color="0.4", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                       rotation=30, ha="right")
    ax.set_ylabel(r"$\rho_{\mathrm{Spearman}}(PR,\ \mathrm{displacement})$")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_02_q1_pr_vs_displacement.pdf", transparent=True)
    plt.close(fig)

    # Fig 3: Q2 — epi mass of the MOST-LOCALIZED modes, observed vs null
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    rpe = rp[rp.n_epi > 0]
    for j, band in enumerate(bands):
        cb = rpe[rpe.band == band]
        o = cb["mE_loc_obs"].values
        s = cb["mE_loc_surr"].values
        ax.scatter(np.full_like(o, x[j] - 0.13), o, s=14, color="C3",
                   alpha=0.8, zorder=3, label="observed" if j == 0 else None)
        ax.scatter(np.full_like(s, x[j] + 0.13), s, s=14, color="C0",
                   alpha=0.8, zorder=3, label="strength null" if j == 0 else None)
        for oi, si in zip(o, s):
            ax.plot([x[j] - 0.13, x[j] + 0.13], [oi, si], color="0.7",
                    lw=0.5, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in bands],
                       rotation=30, ha="right")
    ax.set_ylabel(r"epi mass of localized modes $m^{E}_{\mathrm{loc}}$")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_03_q2_epi_mass.pdf", transparent=True)
    plt.close(fig)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(loc: pd.DataFrame, q1: pd.DataFrame, q2: pd.DataFrame,
                 runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: eigenmode_localization")
    a("scope: q1_trace_in_localized_modes + q2_epi_eigenmode_trap")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: complete_both_negative")
    a("build_script: scripts/01_compute/audit/audit_91_eigenmode_localization.py")
    a("scope_report: .agents/guides/task-persistence-investigation/2026-05-08_epi-eigenmode-localization.md")
    a("---")
    a("")
    a("# Eigenmode localization — Q1 (trace) and Q2 (epi trap)")
    a("")
    a("**Head.** Observed FC-Laplacian eigenmodes are strongly localized "
      "relative to a degree/strength-matched null (3-6x, 10/10 both alpha and "
      "beta) — but this localization explains NEITHER the alpha "
      "Grassmann/cophenet dissociation (Q1) NOR an epileptic-node trap (Q2). "
      "Q1 is doubly falsified: the leading (slowest) modes are the MOST "
      "localized (not extended), and cross-phase mode displacement INCREASES "
      "with participation number (extended modes carry the trace). Q2 is a "
      "clean negative: epi mode-mass does not beat the strength null in any "
      "band. All nulls are the cached matched-strength surrogate (R=200).")
    a("")
    a("## Structural localization (observed PR vs strength null, rest_post)")
    a("")
    a("| band | n | median PR ratio (obs/null) | n localized p<.05 |")
    a("|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(loc.band)]:
        r = loc[loc.band == band].iloc[0]
        a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {int(r['n_patients'])} "
          f"| {r['median_pr_ratio']:.3f} | {int(r['n_localized_p_lt_0.05'])} |")
    a("")
    a("Ratio < 1 = observed more localized than degree predicts. Lower-tail "
      "surrogate p (PR below null).")
    a("")
    a("## Q1 — does the trace live in localized non-leading modes? (rest_post)")
    a("")
    a("| band | lead<bulk K10 | median rho(PR,disp) | rho>0 (extended carry) "
      "| rho<0 (localized carry) | Wilcoxon p |")
    a("|---|---|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(q1.band)]:
        r = q1[q1.band == band].iloc[0]
        a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
          f"| {int(r['n_lead_more_localized_K10'])}/{int(r['n_defined_K10'])} "
          f"| {r['median_spearman_pr_disp']:+.3f} "
          f"| {int(r['n_rho_positive'])}/{int(r['n_patients'])} "
          f"| {int(r['n_rho_negative'])}/{int(r['n_patients'])} "
          f"| {r['wilcoxon_p_rho']:.4f} |")
    a("")
    a("- `lead<bulk`: patients where the 10 slowest nontrivial modes are MORE "
      "localized (lower PR) than the bulk. The Q1 hypothesis needed leading = "
      "EXTENDED (lead>bulk); the data give lead<bulk.")
    a("- `rho(PR,disp)`: per-patient Spearman of participation number vs "
      "cross-phase displacement of rest_post modes. rho>0 ⇒ EXTENDED modes "
      "displace more (trace carried by extended modes). The Q1 hypothesis "
      "needed rho<0 (localized modes carry the change).")
    a("- **Verdict: Q1 FALSIFIED.** The alpha dissociation is not a "
      "localized-vs-extended mode split.")
    a("")
    a("## Q2 — do localized modes trap on epileptic nodes? (rest_post)")
    a("")
    a("Two non-basis-invariant statistics (the mean over ALL modes of epi mass "
      "is identically |E|/N for any orthonormal basis, so it is vacuous and is "
      "NOT used): epi mass of the most-localized mode decile `m^E_loc`, and the "
      "single largest epi-mass mode `m^E_max`. Null = SAME epi indices on the "
      "strength-matched surrogate eigvecs (epi strength profile preserved). "
      "Upper-tail p (trap = more mass than degree).")
    a("")
    a("| band | n (epi) | n trap_loc p<.05 | m^E_loc obs/null | Wilcoxon p_loc "
      "| n trap_max p<.05 | m^E_max obs/null | Wilcoxon p_max |")
    a("|---|---|---|---|---|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(q2.band)]:
        r = q2[q2.band == band].iloc[0]
        a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
          f"| {int(r['n_patients_epi'])} "
          f"| {int(r['n_trap_loc_p_lt_0.05'])} "
          f"| {r['median_mE_loc_obs']:.3f}/{r['median_mE_loc_surr']:.3f} "
          f"| {r['cohort_wilcoxon_p_loc']:.4f} "
          f"| {int(r['n_trap_max_p_lt_0.05'])} "
          f"| {r['median_mE_max_obs']:.3f}/{r['median_mE_max_surr']:.3f} "
          f"| {r['cohort_wilcoxon_p_max']:.4f} |")
    a("")
    a("**Disambiguation of the `m^E_max` firing.** `m^E_max` exceeds the "
      "STRENGTH null because the observed modes are localized (the structural "
      "finding): the max of ANY per-mode mass beats the extended-null max. The "
      "random-node-set control (same |E|, OBSERVED eigvecs) shows the epi set "
      "is NOT special:")
    a("")
    a("| band | median epi m^E_max upper-p vs random sets | n special p<.05 |")
    a("|---|---|---|")
    for band in [b for b in es.ALL_BANDS if b in set(q2.band)]:
        r = q2[q2.band == band].iloc[0]
        a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} "
          f"| {r['median_p_max_vs_random']:.3f} "
          f"| {int(r['n_max_special_vs_random_p_lt_0.05'])} |")
    a("")
    a("- **Verdict: Q2 CLEAN NEGATIVE.** The Anderson-trap statistic "
      "(`m^E_loc`, epi mass of localized modes) does not beat the strength "
      "null; the `m^E_max` strength-null excess is the localization-extreme "
      "artifact (epi behaves like a random set of its size on the observed "
      "modes). Consistent with the established spatial-delocalization result "
      "(per-patient localization null 0-1/10, anatomy retracted).")
    a("")
    a("## Provenance")
    a("- Observed eigvecs: eigh(L=D-W) per (patient,band,phase), identical to "
      "the surrogate construction.")
    a("- Null: matched-strength R=200 cached eigvecs "
      "(`load_or_compute_eigs_at_path`, seed 20260511, full-graph cache).")
    a("- Library: `participation_number`, `node_set_mode_mass`, "
      "`subspace_displacement` (utils.metrics.spectral); `surrogate_p_value` "
      "(utils.metrics.hypothesis).")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `localization_per_cell.csv` — one row per (patient, band, phase)")
    a("- `q1_displacement_per_band.csv` — per-patient PR-displacement Spearman")
    a("- `q2_epi_mass_per_band.csv` — Q2 cohort table")
    a("- `figures/fig_01_pr_obs_vs_null.pdf` — observed vs null PR per band")
    a("- `figures/fig_02_q1_pr_vs_displacement.pdf` — Q1 PR-displacement")
    a("- `figures/fig_03_q2_epi_mass.pdf` — Q2 epi mode-mass")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    # per-cell localization/Q2 need the cached matched-strength ensemble, which
    # exists only for the split-baseline phase vocabulary (no plain rest_pre).
    ap.add_argument("--phases", nargs="+", default=list(es.PHASES_4))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_91] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    cell_rows: list[dict] = []
    disp_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            for phase in args.phases:
                r = per_cell(pat, band, phase, args.n_surrogates,
                             args.swap_factor, args.verbose)
                if r is not None:
                    cell_rows.append(r)
            d = q1_displacement(pat, band, args.verbose)
            if d is not None:
                disp_rows.append(d)
            print(f"[audit_91] {pat}/{band} ({time.time()-tc:.1f}s)")
    runtime = time.time() - t0

    cell = pd.DataFrame(cell_rows)
    disp = pd.DataFrame(disp_rows)
    cell.to_csv(OUT / "localization_per_cell.csv", index=False)
    disp.to_csv(OUT / "q1_displacement_per_band.csv", index=False)

    loc = cohort_localization(cell)
    q1 = cohort_q1(cell, disp)
    q2 = cohort_q2(cell)
    q1.to_csv(OUT / "q1_cohort_per_band.csv", index=False)
    q2.to_csv(OUT / "q2_epi_mass_per_band.csv", index=False)

    if not args.no_figures:
        make_figures(cell, disp)
    write_readme(loc, q1, q2, runtime)

    print("\n=== STRUCTURAL (PR obs/null ratio, rest_post) ===")
    print(loc.to_string(index=False))
    print("\n=== Q1 (trace in localized modes?) ===")
    print(q1[["band", "n_lead_more_localized_K10", "n_defined_K10",
              "median_spearman_pr_disp", "n_rho_positive",
              "n_rho_negative"]].to_string(index=False))
    print("\n=== Q2 (epi trap?) loc=Anderson-stat vs strength null; "
          "max_vs_random=artifact control ===")
    print(q2[["band", "n_patients_epi", "n_trap_loc_p_lt_0.05",
              "cohort_wilcoxon_p_loc", "n_trap_max_p_lt_0.05",
              "median_p_max_vs_random",
              "n_max_special_vs_random_p_lt_0.05"]].to_string(index=False))
    print(f"\n[audit_91] {len(cell)} cells, {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
