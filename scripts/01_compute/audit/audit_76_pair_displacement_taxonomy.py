#!/usr/bin/env python3
"""Audit 76 — per-pair cross-phase displacement portrait (decompose rho_split).

Decomposes the validated rho_split^coph (audit_63, matched-strength
resilient) into its per-pair anatomy. rho_split is just the rank
correlation of a scatter whose every dot is ONE pair (i,j):

    x = Delta_task(i,j) = D_coph_task(i,j) - D_coph_preA(i,j)
    y = Delta_rest(i,j) = D_coph_post(i,j) - D_coph_preB(i,j)

FINDING (beta, 2026-06-03): a per-pair HARD partition (anchor/reset/
trace/anti by magnitude threshold) does NOT describe the cophenetic
trace. Diagnostic shows rho_split lives in the BULK of small movers
(rho among bottom-50% |dT| ~= rho_full; rho among top-10% movers is
often lower or negative). The trace is a broad, low-amplitude,
sign-coherent drift, not a distinct subset of reorganizing pairs. So
the characterization is done at the (patient, band) CELL level via two
threshold-free quantities, and the dots are shown as a cloud colored by
sign-alignment (no magnitude bins):

    reorg_magnitude  m = median|Delta_task| / median|Delta_base|
                         (how much task moved the tree vs baseline noise)
    coherence        c = fraction of pairs with sign(dT)==sign(dR)
    rho_split        = the validated rank correlation (tilt of the cloud)

Cell readings that emerge cleanly at beta:
    TRACE   : m high, c high, rho>0        (e.g. Pat_03/06 m~10 c~0.99)
    ANCHOR  : m ~< 1 (task barely moves tree), rho~0   (Pat_15 m=0.47)
    RANDOM  : m ~> 1 but c~0.5, rho~0       (Pat_10/14 - incoherent;
                                             trace hides in Grassmann subspace)

NOT a new statistic: the cohort trace claim remains rho_split (already
matched-strength validated). This is the descriptive decomposition of
WHERE that validated signal lives. Pipeline is bit-exact with audit_63
(its load_phase_fc + lrg_ultrametric_condensed imported directly);
rho_split reproduces the locked per_patient_per_band.csv (integrity print).

Usage: python audit_76_pair_displacement_taxonomy.py [--band beta]
Outputs data/audit/pair_displacement_taxonomy/:
    figures/{band}_pair_clouds.pdf      - 10-panel per-pair cloud, sign-colored
    per_patient_{band}.csv              - per-patient m, c, rho, bulk/top split
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()

# Import audit_63's validated FC->cophenet pipeline verbatim (no fork).
_spec = importlib.util.spec_from_file_location(
    "audit_63", ROOT / "scripts" / "01_compute" / "audit"
    / "audit_63_split_baseline_surrogate.py")
a63 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a63)

OUT = ROOT / "data" / "audit" / "pair_displacement_taxonomy"
FIG = OUT / "figures"
VALIDATED_CSV = (ROOT / "data" / "audit"
                 / "matched_strength_surrogate_split_baseline"
                 / "per_patient_per_band.csv")

C_ALIGN = "#c1121f"    # sign(dT)==sign(dR): moved same way task & rest
C_CROSS = "#1f77b4"    # opposite sign: task and rest disagree
NOISE_PCT = 95.0       # |Delta_base| percentile, used only as a visual ruler


def phase_cophenet(pat: str, phase: str, band: str) -> np.ndarray:
    W = a63.load_phase_fc(pat, phase, band)
    return a63.lrg_ultrametric_condensed(W)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="beta")
    args = ap.parse_args()
    band = args.band

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    a63.ensure_half_fcs  # noqa: B018  (cache already populated; load_phase_fc reads it)

    ref = pd.read_csv(VALIDATED_CSV)
    ref = ref[ref.band == band].set_index("patient")["obs_rho"].to_dict()

    use_lrg_style()
    fig, axes = plt.subplots(2, 5, figsize=(15, 6.2), sharex=True, sharey=True)
    axes = axes.ravel()

    rows = []
    print(f"\n================ {band}  per-pair displacement portrait ================")
    print(f"{'patient':<8}{'rho':>7}{'(valid)':>9}{'  m_task':>8}{'  m_rest':>8}"
          f"{'  coher':>7}{'  rho_bulk':>10}{'  rho_top':>9}   reading")
    for ax, pat in zip(axes, a63.COHORT):
        D_preA = phase_cophenet(pat, "rest_pre_A", band)
        D_preB = phase_cophenet(pat, "rest_pre_B", band)
        D_task = phase_cophenet(pat, "task_test", band)
        D_post = phase_cophenet(pat, "rest_post", band)

        dT = D_task - D_preA
        dR = D_post - D_preB
        dBase = D_preA - D_preB
        rho = float(spearmanr(dT, dR).statistic)

        # threshold-free cell descriptors: task & rest displacement in
        # baseline-noise units, plus where the rank signal lives (bulk vs top).
        base = max(np.median(np.abs(dBase)), 1e-12)
        m_task = float(np.median(np.abs(dT)) / base)
        m_rest = float(np.median(np.abs(dR)) / base)
        c = float(np.mean(np.sign(dT) == np.sign(dR)))
        order = np.argsort(np.abs(dT))
        bulk = order[:len(order) // 2]
        top = order[int(0.9 * len(order)):]
        rho_bulk = float(spearmanr(dT[bulk], dR[bulk]).statistic)
        rho_top = float(spearmanr(dT[top], dR[top]).statistic)

        # post-hoc cell reading (descriptive, not a gate). Two axes: did it move
        # (m_task, m_rest vs baseline noise) and did it persist (rho, validated).
        moved_task = m_task >= 1.0
        moved_rest = m_rest >= 1.0
        if not moved_task and not moved_rest:
            reading = "anchor"
        elif moved_task and not moved_rest:
            reading = "reset"
        elif not moved_task and moved_rest:
            reading = "rest-drift"
        else:
            reading = "trace" if rho >= 0.15 else "random"
        if reading in ("anchor", "reset", "rest-drift") and rho >= 0.15:
            reading = "weak-trace"      # faint but rank-coherent drift

        rv = ref.get(pat, np.nan)
        flag = "" if abs(rho - rv) < 1e-3 else f"  <-MISMATCH {rho-rv:+.2e}"
        print(f"{pat:<8}{rho:>7.3f}{rv:>9.3f}{m_task:>8.2f}{m_rest:>8.2f}{c:>7.2f}"
              f"{rho_bulk:>10.3f}{rho_top:>9.3f}   {reading}{flag}")
        rows.append(dict(patient=pat, band=band, rho_split=rho, rho_validated=rv,
                         m_task=m_task, m_rest=m_rest, coherence=c,
                         rho_bulk=rho_bulk, rho_top=rho_top, reading=reading))

        # --- threshold-free cloud, symlog in baseline-noise units, sign-colored ---
        t = float(np.percentile(np.abs(dBase), NOISE_PCT)) or float(np.std(dBase)) or 1.0
        x, y = dT / t, dR / t
        align = np.sign(dT) == np.sign(dR)
        ax.scatter(x[~align], y[~align], s=2.0, c=C_CROSS, alpha=0.25, linewidths=0)
        ax.scatter(x[align], y[align], s=2.0, c=C_ALIGN, alpha=0.25, linewidths=0)
        lim = 60.0
        ax.plot([-lim, lim], [-lim, lim], color="0.4", lw=0.6, ls="--", zorder=0)
        ax.axhline(0, color="0.7", lw=0.5, zorder=0)
        ax.axvline(0, color="0.7", lw=0.5, zorder=0)
        ax.axhspan(-1, 1, color="0.85", alpha=0.4, zorder=0)
        ax.axvspan(-1, 1, color="0.85", alpha=0.4, zorder=0)
        ax.set_xscale("symlog", linthresh=1.0)
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.text(0.04, 0.90, pat.replace("Pat_", "P"), transform=ax.transAxes,
                fontsize=9, fontweight="bold")

    for ax in axes[5:]:
        ax.set_xlabel(r"$\Delta_{\mathrm{task}}$  (baseline-noise units)")
    for ax in (axes[0], axes[5]):
        ax.set_ylabel(r"$\Delta_{\mathrm{rest}}$")

    handles = [Patch(facecolor=C_ALIGN, label="aligned (same sign task & rest)"),
               Patch(facecolor=C_CROSS, label="crossed (opposite sign)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out_pdf = FIG / f"{band}_pair_clouds.pdf"
    fig.savefig(out_pdf, transparent=True, bbox_inches="tight")
    plt.close(fig)

    df = pd.DataFrame(rows)
    csv = OUT / f"per_patient_{band}.csv"
    df.to_csv(csv, index=False)

    counts = df.reading.value_counts()
    print(f"\ncohort ({band}): "
          + ", ".join(f"{n} {k}" for k, n in counts.items()))
    print(f"\nWrote:\n  {out_pdf}\n  {csv}")


if __name__ == "__main__":
    main()
