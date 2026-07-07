#!/usr/bin/env python3
r"""audit_156 — does the LRG cophenet resolve structure raw FC CANNOT? (double dissociation)

CRITICAL PREAMBLE (per feedback_critical_null_preamble):
  (1) CLAIM: the LRG cophenetic trace measures persistence of the COMMUNITY
      HIERARCHY, a phenomenon that raw pairwise FC (a local, non-multiscale
      representation) cannot resolve; conversely raw FC responds to hierarchy-
      neutral edge drift that the cophenet correctly filters. => the two are
      NON-EQUIVALENT representations, not "cophenet = smoothed raw".
  (2) NULL / straw-man to kill: "cophenet is just raw FC after a monotone
      transform, so it measures the same thing". If true, raw-trace and
      cophenet-trace would TRACK each other across every regime.
  (3) STRONGEST ALTERNATIVE the test must beat: a single construction can be
      cherry-picked. So we SWEEP the nuisance (edge noise / drift amplitude) and
      require the dissociation to be a ROBUST REGIME, not one lucky point; and we
      show BOTH directions (H: cophenet-only, E: raw-only).
  (4) WHAT WOULD FALSIFY the claim: if raw-trace and cophenet-trace moved
      together in BOTH regimes (no gap), cophenet would add nothing over raw.
  (5) LIMITATIONS: SBM blocks are idealised; real FC is not exactly block-
      structured. This proves the phenomenon is POSSIBLE and non-resolvable by
      raw; the REAL data placement (data/audit/raw_vs_multiscale: beta traces on
      BOTH probes and cophenet SHARPENS it 6->7/10; theta traces on raw +0.117
      7/10 but cophenet KILLS it -0.040 2/10) shows the regimes are REALISED --
      beta ~ Regime H, theta ~ Regime E. It does NOT claim the real beta trace is
      raw-invisible (it is not); it claims the DISCRIMINATION (beta vs theta) and
      the community-level persistence are non-raw-resolvable.

Split-baseline trace (matches rho_sym): rho(X) = Spearman(X_task - X_preA,
X_post - X_preB) over condensed pairs, X in {raw W, cophenet(W)}.

Writes: data/reports/rho_sym_band_map/fig_raw_vs_cophenet_resolving_power.pdf
        data/audit/raw_vs_cophenet_resolving_power/sweeps.csv
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.metrics.node_localization import canonical_cophenet
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()
OUT = ROOT / "data/reports/rho_sym_band_map/fig_raw_vs_cophenet_resolving_power.pdf"
CSV = ROOT / "data/audit/raw_vs_cophenet_resolving_power/sweeps.csv"

N, K, P_IN, P_OUT, REAL = 60, 4, 0.6, 0.15, 80
C_RAW, C_COPH = "#9467bd", "#1a9a4f"


def _sbm(assign, sigma, rng):
    same = assign[:, None] == assign[None, :]
    m = np.where(same, P_IN, P_OUT)
    n = rng.normal(0, sigma, (N, N)); n = (n + n.T) / 2
    W = np.clip(m + n, 0, None); np.fill_diagonal(W, 0.0)
    return (W + W.T) / 2


def _trace(pa, tk, pb, po):
    return spearmanr(tk - pa, po - pb).correlation


def _raw(W): return squareform(W, checks=False)
def _coph(W): return canonical_cophenet(W)


def regime_H(sigmas, rng):
    """Task REORGANIZES blocks and rest_post KEEPS them; edges redrawn each phase."""
    P_pre = rng.integers(0, K, N); P_task = P_pre.copy()
    mv = rng.choice(N, int(0.35 * N), replace=False); P_task[mv] = rng.integers(0, K, len(mv))
    raw_m, raw_s, coph_m, coph_s = [], [], [], []
    for sg in sigmas:
        R, C = [], []
        for _ in range(REAL):
            pa, pb = _sbm(P_pre, sg, rng), _sbm(P_pre, sg, rng)
            tk, po = _sbm(P_task, sg, rng), _sbm(P_task, sg, rng)
            R.append(_trace(_raw(pa), _raw(tk), _raw(pb), _raw(po)))
            C.append(_trace(_coph(pa), _coph(tk), _coph(pb), _coph(po)))
        raw_m.append(np.mean(R)); raw_s.append(np.std(R))
        coph_m.append(np.mean(C)); coph_s.append(np.std(C))
    return map(np.array, (raw_m, raw_s, coph_m, coph_s))


def regime_E(amps, rng):
    """Blocks STATIC; a persistent WITHIN-BLOCK zero-mean edge rewiring (hierarchy-neutral)."""
    P = rng.integers(0, K, N); same = P[:, None] == P[None, :]
    raw_m, raw_s, coph_m, coph_s = [], [], [], []
    for amp in amps:
        R, C = [], []
        for _ in range(REAL):
            Dp = rng.normal(0, amp, (N, N)); Dp = (Dp + Dp.T) / 2
            Dp[~same] = 0.0; np.fill_diagonal(Dp, 0.0)
            pa, pb = _sbm(P, 0.15, rng), _sbm(P, 0.15, rng)
            tk = np.clip(_sbm(P, 0.15, rng) + Dp, 0, None)
            po = np.clip(_sbm(P, 0.15, rng) + Dp, 0, None)
            R.append(_trace(_raw(pa), _raw(tk), _raw(pb), _raw(po)))
            C.append(_trace(_coph(pa), _coph(tk), _coph(pb), _coph(po)))
        raw_m.append(np.mean(R)); raw_s.append(np.std(R))
        coph_m.append(np.mean(C)); coph_s.append(np.std(C))
    return map(np.array, (raw_m, raw_s, coph_m, coph_s))


def _panel(ax, x, rm, rs, cm, cs, xlabel, tag, shade_where, note):
    se = 1.0 / np.sqrt(REAL)
    ax.axhline(0, color="0.7", lw=0.8, ls="--", zorder=1)
    ax.plot(x, rm, "-o", color=C_RAW, lw=2, ms=4, label="raw FC trace", zorder=3)
    ax.fill_between(x, rm - rs * se, rm + rs * se, color=C_RAW, alpha=0.18, zorder=2)
    ax.plot(x, cm, "-o", color=C_COPH, lw=2, ms=4, label="LRG cophenet trace", zorder=3)
    ax.fill_between(x, cm - cs * se, cm + cs * se, color=C_COPH, alpha=0.18, zorder=2)
    m = shade_where
    ax.fill_between(x[m], np.minimum(rm, cm)[m], np.maximum(rm, cm)[m],
                    color="0.75", alpha=0.35, zorder=1)
    ax.set_xlabel(xlabel); ax.set_ylabel(r"cross-phase trace  (Spearman $\rho$)")
    ax.text(0.03, 0.96, tag, transform=ax.transAxes, fontsize=12, va="top", fontweight="bold")
    ax.text(0.5, -0.30, note, transform=ax.transAxes, fontsize=7.6, ha="center", va="top",
            color="0.25")
    ax.legend(loc="upper right", frameon=False, fontsize=8)


def main():
    rng = np.random.default_rng(20260706)
    sigmas = np.array([0.05, 0.10, 0.15, 0.22, 0.30, 0.40, 0.50, 0.60])
    amps = np.array([0.08, 0.15, 0.22, 0.30, 0.40, 0.50, 0.65, 0.80])
    hrm, hrs, hcm, hcs = regime_H(sigmas, rng)
    erm, ers, ecm, ecs = regime_E(amps, rng)

    print("Regime H (hierarchy persists, edges noisy):")
    for s, r, c in zip(sigmas, hrm, hcm):
        print(f"  sigma={s:.2f}  raw={r:+.3f}  coph={c:+.3f}  gap(coph-raw)={c-r:+.3f}")
    print("Regime E (hierarchy-neutral within-block drift):")
    for a, r, c in zip(amps, erm, ecm):
        print(f"  amp={a:.2f}  raw={r:+.3f}  coph={c:+.3f}  gap(raw-coph)={r-c:+.3f}")

    CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.concat([
        pd.DataFrame(dict(regime="H", x=sigmas, raw=hrm, raw_sd=hrs, coph=hcm, coph_sd=hcs)),
        pd.DataFrame(dict(regime="E", x=amps, raw=erm, raw_sd=ers, coph=ecm, coph_sd=ecs)),
    ]).to_csv(CSV, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    _panel(axes[0], sigmas, hrm, hrs, hcm, hcs, r"edge noise  $\sigma$", r"$\mathbf{A}$",
           hcm - hrm > 0.08,
           "H: task reorganizes the community hierarchy, rest_post keeps it, edges\n"
           "redrawn each phase.  raw collapses with noise; cophenet resolves the\n"
           "persistent hierarchy — the $\\beta$-like regime (raw+cophenet trace).")
    _panel(axes[1], amps, erm, ers, ecm, ecs, r"within-block drift amplitude", r"$\mathbf{B}$",
           erm - ecm > 0.06,
           "E: blocks static, a hierarchy-neutral within-block edge pattern persists.\n"
           "raw sees the drift; cophenet is blind until it becomes structural — the\n"
           "$\\theta$-like regime (raw traces +0.117, cophenet kills it to −0.040).")
    fig.subplots_adjust(bottom=0.30, wspace=0.24)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\nwrote {OUT}\nwrote {CSV}")


if __name__ == "__main__":
    main()
