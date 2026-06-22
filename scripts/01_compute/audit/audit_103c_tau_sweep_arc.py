#!/usr/bin/env python3
"""audit_103c — τ-sweep of the consolidation arc (observed-statistic diagnostic).

Scope + 5-point preamble:
`.agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md`.

Question: every N2 number (encoding echo `T_learn`, inference-specific
`T_infspec·e`, the α/β dissociation) was computed at the SINGLE finest diffusion
time τ = 1/λ_max. Is the decomposition scale-robust, and is inference-specific
consolidation fine-scale (τ≈1/λ_max) or mesoscale (τ≫1/λ_max)?

This is an OBSERVED-ONLY sweep (no matched-strength surrogate — that is mandatory
only for a NEW positive claim at some τ≠1/λ_max, deferred to a follow-up). Pure
reuse of audit_103 (loaders + functionals); the only new primitive is the
cophenetic at a parameterized τ. Each phase Laplacian is eigendecomposed ONCE and
the heat kernel e^{−τL} is reformed per τ from the cached spectrum → the sweep is
cheap. Correctness anchor: at τ_mult=1 the cophenetic is byte-identical to
audit_63's `lrg_ultrametric_condensed`, so `T_*` reproduce audit_103 exactly.

τ is a per-phase multiple of that phase's own 1/λ_max (matching the canonical,
which uses each phase's own λ_max). τ_mult ∈ [1, 10] (the SOZ-marker work showed
the informative slow-diffusion regime ≈ 4–10/λ_max in this dataset).

Output: data/audit/consolidation_arc/arc_tau_sweep.csv (one row per pat×band×τ).
"""
from __future__ import annotations

import argparse
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial.distance import squareform
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_103_cophenetic_consolidation_arc import (  # type: ignore
    ARC_PHASES, BANDS, COHORT, OUT, _obs_functionals,
    ensure_half_fcs, load_phase_fc,
)

TAU_MULTS = np.logspace(0.0, 1.0, 13)   # 1 .. 10 × (1/λ_max), log-spaced
ARC_CSV = OUT / "arc_per_patient.csv"   # canonical (τ_mult=1) for the correctness anchor


def _eig_L(W: np.ndarray):
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return np.linalg.eigh(L)             # (evals ascending, evecs)


def _cophenet_at_taumult(evals: np.ndarray, evecs: np.ndarray, tau_mult: float):
    """Cophenetic condensed vector at τ = tau_mult/λ_max. At tau_mult=1 this is
    byte-identical to audit_63.lrg_ultrametric_condensed (correctness anchor)."""
    lam_max = evals[-1]
    tau = tau_mult / lam_max
    diag_exp = np.exp(-tau * evals)
    rho = (evecs * diag_exp) @ evecs.T
    rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        Trho = 1.0 / rho
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)
    finite = np.isfinite(Trho)
    frac_finite = float(finite.mean())
    if not finite.all():
        cap = np.nanmax(Trho[finite]) if finite.any() else 1e6
        Trho = np.where(finite, Trho, cap)
    cond = squareform(Trho, checks=False)
    Z = linkage(cond, method="average")
    return cophenet(Z), frac_finite


def main() -> None:
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, BANDS)
        for band in BANDS:
            try:
                Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
                eigs = {ph: _eig_L(Ws[ph]) for ph in ARC_PHASES}   # once per phase
            except Exception as exc:  # noqa: BLE001
                print(f"[audit_103c] FAIL load {pat} {band}: {type(exc).__name__}: {exc}")
                continue
            for tm in TAU_MULTS:
                D, fr_min = {}, 1.0
                for ph in ARC_PHASES:
                    d, fr = _cophenet_at_taumult(*eigs[ph], tm)
                    D[ph] = d
                    fr_min = min(fr_min, fr)
                fn = _obs_functionals(D)
                rows.append({"patient": pat, "band": band, "tau_mult": float(tm),
                             **fn, "finite_frac": fr_min})
        print(f"[audit_103c] {pat} done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "arc_tau_sweep.csv", index=False)
    print(f"\n[audit_103c] wrote {OUT/'arc_tau_sweep.csv'} ({len(df)} rows)\n")

    # --- Correctness anchor: τ_mult=1 reproduces audit_103 arc_per_patient.csv ---
    if ARC_CSV.exists():
        ref = pd.read_csv(ARC_CSV).set_index(["patient", "band"])
        cur = df[np.isclose(df["tau_mult"], 1.0)].set_index(["patient", "band"])
        diffs = []
        for k in ("T_learn", "T_infspec_pe"):
            for idx, r in cur.iterrows():
                if idx in ref.index and k in ref.columns:
                    diffs.append(abs(r[k] - float(ref.loc[idx, k])))
        md = max(diffs) if diffs else np.nan
        status = "PASS" if md < 1e-9 else ("CLOSE" if md < 1e-4 else "FAIL")
        print(f"[ANCHOR] τ_mult=1 vs audit_103: max|Δ|={md:.2e} → {status}")

    # --- Cohort τ-response: median + sign + observed-only one-sided Wilcoxon vs 0 ---
    def _wilcox_gt0(v):
        v = v[np.isfinite(v)]
        if v.size < 3 or np.allclose(v, 0):
            return np.nan
        try:
            return float(wilcoxon(v, alternative="greater").pvalue)
        except ValueError:
            return np.nan

    print("\n=== τ-RESPONSE (cohort median; n_pos/10; observed-only Wilcoxon>0 p) ===")
    print("    (NOT the matched-strength verdict — a shape diagnostic; positive=trace)\n")
    for band in ("alpha", "beta", "low_gamma", "theta", "delta", "high_gamma"):
        sub = df[df["band"] == band]
        if sub.empty:
            continue
        print(f"--- {band} ---")
        print(f"  {'τ/λmax':>7} | {'T_learn med':>12} {'p>0':>7} | "
              f"{'T_infspec·e med':>15} {'p>0':>7} | {'min finite':>10}")
        for tm in TAU_MULTS:
            s = sub[np.isclose(sub["tau_mult"], tm)]
            tl, ti = s["T_learn"].values, s["T_infspec_pe"].values
            print(f"  {tm:>7.2f} | {np.median(tl):>+12.3f} {_wilcox_gt0(tl):>7.3f} | "
                  f"{np.median(ti):>+15.3f} {_wilcox_gt0(ti):>7.3f} | "
                  f"{s['finite_frac'].min():>10.3f}")
        print()
    print("[read] FLAT across τ → N2 verdicts scale-robust. RISING T_infspec·e with "
          "τ → inference is a MESOSCALE phenomenon (slower diffusion integrates "
          "multi-step relational paths). Watch whether α stays below β (dissociation "
          "holds) or converges. A surrogate null is required before calling any "
          "τ≠1 peak a result.")


def make_figure() -> None:
    """τ-response curves (cohort median ± IQR) for T_learn and T_infspec·e, read
    from the cached sweep CSV. β/α emphasized; other bands greyed. PDF, vector."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()

    df = pd.read_csv(OUT / "arc_tau_sweep.csv")
    taus = np.sort(df["tau_mult"].unique())
    # emphasis palette (saturated, no near-white): beta, alpha highlighted.
    EMPH = {"beta": "#1f77b4", "alpha": "#d62728"}
    GREY = "0.70"
    order = ["beta", "alpha", "low_gamma", "theta", "delta", "high_gamma"]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharex=True)
    for ax, key, tag, ttl in (
        (axes[0], "T_learn", "a", r"encoding echo  $T_{\mathrm{learn}}$"),
        (axes[1], "T_infspec_pe", "b",
         r"inference-specific  $T_{\mathrm{inf{\cdot}e}}$"),
    ):
        for band in order:
            sub = df[df["band"] == band]
            med = np.array([np.median(sub[np.isclose(sub.tau_mult, t)][key]) for t in taus])
            c = EMPH.get(band, GREY)
            lw = 2.2 if band in EMPH else 1.0
            z = 3 if band in EMPH else 1
            ax.plot(taus, med, color=c, lw=lw, zorder=z,
                    label=band.replace("_", "-") if band in EMPH or ax is axes[1] else None)
            if band in EMPH:
                q1 = np.array([np.percentile(sub[np.isclose(sub.tau_mult, t)][key], 25) for t in taus])
                q3 = np.array([np.percentile(sub[np.isclose(sub.tau_mult, t)][key], 75) for t in taus])
                ax.fill_between(taus, q1, q3, color=c, alpha=0.15, zorder=z - 1, lw=0)
        ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
        ax.axvline(1.0, color="0.4", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\tau \,/\, (1/\lambda_{\max})$")
        ax.set_ylabel("cohort median concordance")
        ax.text(0.04, 0.93, tag, transform=ax.transAxes, fontweight="bold", va="top")
        ax.set_title(ttl)

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=len(labels), frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    outdir = OUT / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "arc_tau_sweep.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"[audit_103c] wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", action="store_true",
                    help="make the τ-response figure from the cached sweep CSV (no recompute)")
    args = ap.parse_args()
    if args.figure:
        make_figure()
    else:
        main()
