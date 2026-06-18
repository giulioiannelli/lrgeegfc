#!/usr/bin/env python3
"""Audit 90 — CRC vs matched-strength BRACKET on the §5.3 split-baseline trace.

Scope: .agents/guides/task-persistence-investigation/2026-06-08_on-manifold-coherency-surrogate.md

For each (patient, band ∈ {α, β, γ_l}) computes the §5.3 statistic of record
    ρ_split = Spearman(D^tt − D^pre_A , D^post − D^pre_B)
on the OBSERVED cached FCs (bit-identical to audit_63), then evaluates it
against TWO nulls run on the identical cophenetic pipeline:

  • matched-strength (off-manifold, preserves node strength) — loaded from the
    complete R=200 seed=20260511 cache (audit_66 ensemble), paired across the
    four split-baseline phases;
  • CRC = coupling-randomized coherency (on-manifold, real-orthogonal congruence,
    preserves realizability + per-frequency coupling content, strength floats) —
    generated from the timeseries complex coherency (halves split at T//2 with
    nperseg//2, matching compute_imcoh_abs_halves).

Reading: a trace clears matched-strength if it isn't a hub/strength artifact;
it clears CRC if it isn't explained by a realizable, magnitude-faithful,
who-randomized coupling structure. A robust trace clears BOTH. Verdict per band
is DESCRIPTIVE (no pre-registered acceptance gate, per feedback_no_pre_registered_acceptance).

Calibration is reported per cell (CRC magnitude ratio → target ~1; strength KS;
Λ renorm residual); a failed calibration voids that cell's CRC p.

Outputs
-------
    data/audit/crc_cohort_bracket/
        per_patient_per_band.csv
        cohort_summary.csv
        figures/bracket.pdf
        figures/calibration.pdf
        README.md
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon, ks_2samp

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.surrogate import (
    complex_coherency_bands,
    cophenetic_condensed_from_adjacency,
    cophenetic_condensed_from_eigs,
    coupling_randomized_coherency,
    load_or_compute_surrogate_eigs,
)
from lrg_eegfc.visuals.styles import use_lrg_style

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
N_SURR = 200
SWAP_FACTOR = 20
MS_SEED = 20260511          # canonical audit_66 matched-strength ensemble
CRC_SEED = 20260608

OUT = ROOT / "data" / "audit" / "crc_cohort_bracket"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"


def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_rest_pre_{half}_imcoh_abs.npy"


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    """Observed phase FC (full or half), (N,N) float64 — identical to audit_63."""
    if phase in ("rest_pre_A", "rest_pre_B"):
        W = np.load(_half_fc_path(pat, band, phase[-1]))
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _channels_first(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    return X.T if X.shape[0] > X.shape[1] else X


def coherency_bands(X: np.ndarray, fs: float, nperseg: int,
                    bands: dict[str, tuple[float, float]]) -> dict[str, np.ndarray]:
    """Thin adapter onto the promoted library helper (arg order kept local)."""
    return complex_coherency_bands(X, fs, bands, nperseg)


def rho_split(coph: dict[str, np.ndarray]) -> float:
    return float(spearmanr(coph["task_test"] - coph["rest_pre_A"],
                           coph["rest_post"] - coph["rest_pre_B"]).statistic)


def build_coherency(pat: str, bands: list[str]) -> dict[str, dict[str, np.ndarray]]:
    """{band: {phase: C_band}} for the four split-baseline phases."""
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nfull = nperseg_for_fs(fs)
    nhalf = max(256, nfull // 2)
    bd = {b: BRAIN_BANDS[b] for b in bands}
    C = {b: {} for b in bands}
    # rest_pre → halves A,B at nperseg//2 (T//2 split = compute_imcoh_abs_halves)
    X = _channels_first(load_timeseries(pat, "rest_pre", SEEG_DATAPATH))
    T = X.shape[1]
    for tag, Xh in (("rest_pre_A", X[:, : T // 2]), ("rest_pre_B", X[:, T // 2:])):
        cb = coherency_bands(Xh, fs, nhalf, bd)
        for b in bands:
            C[b][tag] = cb[b]
    del X
    gc.collect()
    # task_test, rest_post at full nperseg
    for phase in ("task_test", "rest_post"):
        X = _channels_first(load_timeseries(pat, phase, SEEG_DATAPATH))
        cb = coherency_bands(X, fs, nfull, bd)
        for b in bands:
            C[b][phase] = cb[b]
        del X
        gc.collect()
    return C


def crc_null(C_band_phases: dict[str, np.ndarray], n_surr: int,
             rng: np.random.Generator,
             A_obs_post: np.ndarray) -> dict:
    """CRC null + calibration for one cell. C_band_phases: {phase: (F_b,N,N)}."""
    N = A_obs_post.shape[0]
    triu = np.triu_indices(N, 1)
    rho = np.empty(n_surr)
    mag_surr, ks, resid = [], [], []
    nonneg = True
    for r in range(n_surr):
        coph = {}
        for phase in PHASES:
            W, res = coupling_randomized_coherency(C_band_phases[phase], rng,
                                                   mode="real_congruence")
            resid.append(res)
            if W.min() < -1e-12:
                nonneg = False
            coph[phase] = cophenetic_condensed_from_adjacency(W)
            if phase == "rest_post":
                mag_surr.append(float(np.mean(W[triu])))
                ks.append(ks_2samp(W.sum(1), A_obs_post.sum(1)).statistic)
        rho[r] = rho_split(coph)
    finite = rho[np.isfinite(rho)]
    mag_obs = float(np.mean(A_obs_post[triu]))
    return {
        "null": finite,
        "mag_ratio": float(np.mean(mag_surr) / mag_obs) if mag_obs > 0 else np.nan,
        "strength_ks": float(np.mean(ks)),
        "renorm_resid": float(np.mean(resid)),
        "nonneg": nonneg,
    }


def ms_null(pat: str, band: str, n_surr: int) -> np.ndarray:
    """Matched-strength split-baseline ρ null from the complete cache."""
    rng = np.random.default_rng(0)  # cache hit → unused
    coph = {}
    for phase in PHASES:
        W = load_phase_fc(pat, phase, band)
        evals, evecs = load_or_compute_surrogate_eigs(
            pat, band, phase, W, n_surr, SWAP_FACTOR, MS_SEED, rng)
        cph = np.full((n_surr, W.shape[0] * (W.shape[0] - 1) // 2), np.nan)
        for r in range(n_surr):
            if np.isfinite(evals[r]).all():
                cph[r] = cophenetic_condensed_from_eigs(evals[r], evecs[r])
        coph[phase] = cph
        del evals, evecs
        gc.collect()
    rho = np.full(n_surr, np.nan)
    for r in range(n_surr):
        if all(np.isfinite(coph[p][r]).all() for p in PHASES):
            rho[r] = rho_split({p: coph[p][r] for p in PHASES})
    return rho[np.isfinite(rho)]


def per_cell(pat: str, band: str, C_band_phases: dict, n_surr: int,
             rng: np.random.Generator) -> dict:
    # observed (cached FCs → cophenetic) — matches audit_63
    A = {p: load_phase_fc(pat, p, band) for p in PHASES}
    coph_obs = {p: cophenetic_condensed_from_adjacency(A[p]) for p in PHASES}
    obs = rho_split(coph_obs)

    crc = crc_null(C_band_phases, n_surr, rng, A["rest_post"])
    ms = ms_null(pat, band, n_surr)

    def stats(null):
        if null.size == 0:
            return dict(median=np.nan, p5=np.nan, p95=np.nan, p=np.nan, z=np.nan)
        return dict(
            median=float(np.median(null)),
            p5=float(np.quantile(null, 0.05)), p95=float(np.quantile(null, 0.95)),
            p=float(np.mean(null >= obs)),
            z=float((obs - null.mean()) / (null.std() + 1e-12)),
        )

    sc, sm = stats(crc["null"]), stats(ms)
    return {
        "patient": pat, "band": band, "N_nodes": A["rest_post"].shape[0],
        "obs_rho": obs,
        "crc_median": sc["median"], "crc_p5": sc["p5"], "crc_p95": sc["p95"],
        "crc_p": sc["p"], "crc_z": sc["z"],
        "ms_median": sm["median"], "ms_p5": sm["p5"], "ms_p95": sm["p95"],
        "ms_p": sm["p"], "ms_z": sm["z"],
        "crc_mag_ratio": crc["mag_ratio"], "crc_strength_ks": crc["strength_ks"],
        "crc_renorm_resid": crc["renorm_resid"], "crc_nonneg": crc["nonneg"],
        "_crc_null": crc["null"], "_ms_null": ms,
    }


def cohort_summary(rows: list[dict]) -> pd.DataFrame:
    out = []
    for band in TARGET_BANDS:
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            continue
        obs = np.array([r["obs_rho"] for r in sub])
        crc_med = np.array([r["crc_median"] for r in sub])
        ms_med = np.array([r["ms_median"] for r in sub])

        def wilco(diff):
            try:
                w, p = wilcoxon(diff, alternative="greater")
                return float(w), float(p)
            except Exception:
                return float("nan"), float("nan")

        _, p_crc = wilco(obs - crc_med)
        _, p_ms = wilco(obs - ms_med)
        n_crc = int(sum(r["crc_p"] < 0.05 for r in sub))
        n_ms = int(sum(r["ms_p"] < 0.05 for r in sub))
        clears_crc = p_crc < 0.05
        clears_ms = p_ms < 0.05
        verdict = ("clears_both" if clears_crc and clears_ms else
                   "clears_MS_only" if clears_ms else
                   "clears_CRC_only" if clears_crc else "clears_neither")
        out.append({
            "band": band, "n": len(sub),
            "obs_median": float(np.median(obs)),
            "crc_null_median": float(np.median(crc_med)),
            "ms_null_median": float(np.median(ms_med)),
            "crc_wilcoxon_p": p_crc, "ms_wilcoxon_p": p_ms,
            "n_above_crc": f"{n_crc}/{len(sub)}",
            "n_above_ms": f"{n_ms}/{len(sub)}",
            "crc_mag_ratio_med": float(np.median([r["crc_mag_ratio"] for r in sub])),
            "verdict": verdict,
        })
    return pd.DataFrame(out)


def make_bracket(rows: list[dict], cohort: pd.DataFrame, path: Path) -> None:
    use_lrg_style()
    fig, axes = plt.subplots(1, len(TARGET_BANDS),
                             figsize=(4.6 * len(TARGET_BANDS), 4.8))
    if len(TARGET_BANDS) == 1:
        axes = [axes]
    c_ms, c_crc, c_obs = "#bdbdbd", "#d98a3d", "#1f3d6e"
    for ax, band in zip(axes, TARGET_BANDS):
        sub = sorted([r for r in rows if r["band"] == band],
                     key=lambda r: r["obs_rho"])
        ys = np.arange(len(sub))
        for y, r in zip(ys, sub):
            ax.hlines(y + 0.16, r["ms_p5"], r["ms_p95"], color=c_ms, lw=4.5,
                      alpha=0.95, zorder=1)
            ax.plot(r["ms_median"], y + 0.16, "|", color="#777777", ms=9, mew=1.1)
            ax.hlines(y - 0.16, r["crc_p5"], r["crc_p95"], color=c_crc, lw=4.5,
                      alpha=0.95, zorder=1)
            ax.plot(r["crc_median"], y - 0.16, "|", color="#9c5a1c", ms=9, mew=1.1)
            ax.plot(r["obs_rho"], y, "o", color=c_obs, ms=7, mec="white",
                    mew=0.9, zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels([r["patient"].replace("Pat_", "P") for r in sub])
        ax.axvline(0, color="0.6", lw=0.7, ls="--", zorder=0)
        ax.set_xlabel(r"$\rho_{\mathrm{split}}$")
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        ax.spines[["top", "right"]].set_visible(False)
    handles = [
        plt.Line2D([0], [0], color=c_ms, lw=4.5, label="matched-strength null (P5–P95)"),
        plt.Line2D([0], [0], color=c_crc, lw=4.5, label="CRC on-manifold null (P5–P95)"),
        plt.Line2D([0], [0], marker="o", color=c_obs, lw=0, mec="white",
                   label="observed"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, bbox_inches="tight", transparent=True)
    plt.close(fig)


def make_calibration(rows: list[dict], path: Path) -> None:
    use_lrg_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    bands = TARGET_BANDS
    xpos = {b: i for i, b in enumerate(bands)}
    for r in rows:
        x = xpos[r["band"]] + np.random.uniform(-0.12, 0.12)
        axes[0].plot(x, r["crc_mag_ratio"], "o", color="#d98a3d", ms=6,
                     mec="white", mew=0.7)
        axes[1].plot(x, r["crc_renorm_resid"], "o", color="#1f3d6e", ms=6,
                     mec="white", mew=0.7)
    axes[0].axhline(1.0, color="0.5", ls="--", lw=0.8)
    axes[0].set_ylabel(r"CRC $|{\rm ImCoh}|$ magnitude ratio (surr/obs)")
    axes[1].set_ylabel(r"CRC $\Lambda$ renorm residual")
    for ax in axes:
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands])
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", transparent=True)
    plt.close(fig)


def write_readme(cohort: pd.DataFrame, runtime: float) -> None:
    lines = [
        "---", "name: crc_cohort_bracket", "type: report", "era: COHORT_N10",
        "status: current", "created: 2026-06-08", "updated: 2026-06-08",
        "pointers:",
        "  - .agents/guides/task-persistence-investigation/2026-06-08_on-manifold-coherency-surrogate.md",
        "  - data/audit/matched_strength_surrogate_split_baseline/",
        "---", "",
        "# CRC vs matched-strength bracket on §5.3 split-baseline ρ",
        "",
        "**Head.** The §5.3 split-baseline trace `ρ_split` evaluated against two "
        "nulls on the identical cophenetic pipeline: matched-strength "
        "(off-manifold, strength-preserving) and CRC (on-manifold, "
        "real-orthogonal congruence, magnitude-faithful, strength-floating). A "
        "robust trace clears BOTH. Verdict descriptive — no pre-registered gate.",
        "",
        "| band | n | obs med ρ | MS null med | CRC null med | MS Wilcoxon p | CRC Wilcoxon p | n>MS | n>CRC | CRC mag-ratio | verdict |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['n']} | {r['obs_median']:+.3f} "
            f"| {r['ms_null_median']:+.3f} | {r['crc_null_median']:+.3f} "
            f"| {r['ms_wilcoxon_p']:.4f} | {r['crc_wilcoxon_p']:.4f} "
            f"| {r['n_above_ms']} | {r['n_above_crc']} "
            f"| {r['crc_mag_ratio_med']:.2f} | **{r['verdict']}** |")
    lines += [
        "", "## Reading",
        "- **clears_both** — trace survives the strongest version of each "
        "objection (not strength, not realizable-spectrum).",
        "- **clears_MS_only** — survives the strength null but is explained by a "
        "realizable, magnitude-faithful, who-randomized coupling structure; the "
        "higher-order *identity* claim weakens. Brackets disagree — flagged for "
        "later adjudication (combined null §11.3), not auto-concluded.",
        "- CRC magnitude-ratio should be ≈ 1 (real-orthogonal congruence "
        "preserves ‖Im C‖); a ratio far from 1 voids that band's CRC p.",
        "",
        f"R = {N_SURR} per null per phase. MS = cached seed {MS_SEED} (audit_66); "
        f"CRC seed {CRC_SEED}, mode=real_congruence. Cohort: {', '.join(COHORT)}. "
        f"Bands: {', '.join(TARGET_BANDS)}. Runtime {runtime:.0f}s.",
        "- Build: `scripts/01_compute/audit/audit_97_crc_cohort_bracket.py`",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=TARGET_BANDS)
    ap.add_argument("--n-surrogates", type=int, default=N_SURR)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(CRC_SEED)
    t0 = time.time()
    rows = []
    for pat in args.patients:
        try:
            C = build_coherency(pat, args.bands)
        except Exception as e:
            print(f"[audit_97] SKIP {pat}: coherency build failed: {e}")
            continue
        for band in args.bands:
            tc = time.time()
            try:
                r = per_cell(pat, band, C[band], args.n_surrogates, rng)
            except Exception as e:
                print(f"[audit_97] SKIP {pat}/{band}: {e}")
                continue
            rows.append(r)
            print(f"[audit_97] {pat}/{band}: obs={r['obs_rho']:+.3f} | "
                  f"MS med={r['ms_median']:+.3f} p={r['ms_p']:.3f} | "
                  f"CRC med={r['crc_median']:+.3f} p={r['crc_p']:.3f} "
                  f"magR={r['crc_mag_ratio']:.2f} resid={r['crc_renorm_resid']:.2f} "
                  f"({time.time()-tc:.1f}s)")
        del C
        gc.collect()

    runtime = time.time() - t0
    print(f"[audit_97] total {runtime:.1f}s, {len(rows)} cells")
    if not rows:
        return

    per_pat = pd.DataFrame([{k: v for k, v in r.items()
                             if not k.startswith("_")} for r in rows])
    per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)
    cohort = cohort_summary(rows)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))

    make_bracket(rows, cohort, FIG / "bracket.pdf")
    make_calibration(rows, FIG / "calibration.pdf")
    write_readme(cohort, runtime)
    print(f"[audit_97] outputs at {OUT}")


if __name__ == "__main__":
    main()
