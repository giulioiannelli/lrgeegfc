#!/usr/bin/env python3
"""Audit 99 — coordinated cross-phase null (SB-CRC) on the §5.3 trace.

Scope: .agents/guides/task-persistence-investigation/2026-06-11_coordinated-cross-phase-null.md

Both nulls we already run (matched-strength, CRC) randomize EACH PHASE
INDEPENDENTLY, so under the §5.3 split-baseline protocol both collapse to a
ρ_split floor of ≈ 0. Clearing that floor only proves *some* cross-phase
structure exists — it does NOT separate task-specific persistence from a
task-independent STABLE COUPLING BACKBONE shared across all phases. SB-CRC
closes that gap:

    B(f)   = Σ_φ w_φ C_φ(f)          # shared backbone (equal-condition weights)
    Δ_φ(f) = C_φ(f) − B(f)           # per-phase deviation
    C̃_φ(f) = renorm( B(f) + O_φ Δ_φ(f) O_φᵀ )   # rotate ONLY the deviation, INDEPENDENT O_φ
    W̃_φ    = mean_f |Im C̃_φ(f)|

Because the cophenetic map is nonlinear, a shared backbone does NOT cancel in
the Δ_task / Δ_rest differences, so the SB-CRC floor LIFTS OFF ZERO to exactly
the shared-backbone contribution. The observed ρ_split is thereby decomposed:

    obs ≈ [independent-CRC floor ≈ 0] + [SB-CRC floor − 0]  (shared backbone)
                                       + [obs − SB-CRC floor] (task-specific)

Reading: obs clears the SB-CRC floor → task-specific persistence beyond the
stable backbone; obs sits INSIDE the SB-CRC floor → the trace is a
task-INDEPENDENT stable-fingerprint effect (a §5.3 headline reframe). This null
can bite the headline — that is its job.

Common frequency grid (scope §5.1): all four phases computed at the FULL
nperseg (rest_pre split at T//2 but with the full nperseg, not nperseg//2), so
the per-frequency backbone B(f) is well-defined. Consequence: the within-pipeline
observed differs slightly from the canonical audit_63 value (cached nperseg//2
halves); BOTH are reported, and the floor + within-pipeline observed share the
grid so the test is internally exact.

Runs whatever ``--patients`` / ``--bands`` are given (default = full n=10 cohort
× {beta, alpha, low_gamma}, R=200). Per cell it emits within-pipeline obs,
canonical obs (audit_63 cross-ref), the SB-CRC floor, the independent-CRC
reference floor, the obs = backbone + task-specific decomposition, and
calibration (magnitude ratio, renorm + PSD residual, backbone energy fraction).
Cohort verdict per band is DESCRIPTIVE (no pre-registered acceptance gate, per
feedback_no_pre_registered_acceptance): a paired one-sided Wilcoxon asks whether
observed exceeds the per-patient SB-CRC floor across the cohort.

The off-manifold caveat (the deviation rotation B + OΔOᵀ is not strictly PSD) was
quantified in the pilot: projecting each surrogate to the nearest PSD matrix
changes the band-averaged |ImCoh| readout by ρ=0.99 (mean per-pair diff 0.004),
so the floor is robust to it; magnitude stays faithful (ratio ~1.1).

Outputs
-------
    data/audit/coordinated_cross_phase_null/
        per_cell.csv
        cohort_summary.csv          (when > 1 patient)
        figures/decomposition.pdf
        README.md
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

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
    coupling_randomized_coherency,
    deviation_rotated_coherency,
    shared_backbone_deviation,
)
from lrg_eegfc.visuals.styles import use_lrg_style

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["beta", "alpha", "low_gamma"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
# Backbone weighting schemes for B(f) = Σ_φ w_φ C_φ(f):
#   equal_condition — the two rest_pre halves count ONCE between them, so each of
#     the three CONDITIONS (rest, task, post) carries weight 1/3. CONSERVATIVE:
#     B absorbs structure common to task+post (incl. persistent task structure),
#     so this can HIDE a real trace → an upper bound on the backbone contribution.
#   baseline — B = rest_pre only (the pre-task fingerprint, task/post excluded).
#     The PRIMARY test of task-specific persistence: deviations are measured FROM
#     baseline, the backbone cannot peek at the task, so obs > floor cleanly means
#     task and post deviations-from-baseline are aligned (= the trace).
WEIGHT_SCHEMES = {
    "equal_condition": {"rest_pre_A": 1 / 6, "rest_pre_B": 1 / 6,
                        "task_test": 1 / 3, "rest_post": 1 / 3},
    "baseline": {"rest_pre_A": 1 / 2, "rest_pre_B": 1 / 2,
                 "task_test": 0.0, "rest_post": 0.0},
}
N_SURR = 200
SEED = 20260611
N_PSD_DRAWS = 5                      # draws on which to compute the PSD residual

OUT = ROOT / "data" / "audit" / "coordinated_cross_phase_null"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"


# --------------------------------------------------------------------------- #
#  loaders (canonical cached observed = audit_63; full-grid coherency = SB-CRC) #
# --------------------------------------------------------------------------- #
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_rest_pre_{half}_imcoh_abs.npy"


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    """Canonical observed phase FC (cached nperseg//2 halves) — audit_63."""
    if phase in ("rest_pre_A", "rest_pre_B"):
        W = np.load(_half_fc_path(pat, band, phase[-1]))
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    return 0.5 * (np.clip(W, 0.0, 1.0) + np.clip(W, 0.0, 1.0).T)


def _channels_first(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    return X.T if X.shape[0] > X.shape[1] else X


def build_coherency_common_grid(pat: str, bands: list[str]) -> dict[str, dict[str, np.ndarray]]:
    """{band: {phase: C_band}} — all four phases at the FULL nperseg (common grid)."""
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nfull = nperseg_for_fs(fs)
    bd = {b: BRAIN_BANDS[b] for b in bands}
    C = {b: {} for b in bands}
    X = _channels_first(load_timeseries(pat, "rest_pre", SEEG_DATAPATH))
    T = X.shape[1]
    for tag, Xh in (("rest_pre_A", X[:, : T // 2]), ("rest_pre_B", X[:, T // 2:])):
        cb = complex_coherency_bands(Xh, fs, bd, nfull)         # FULL nperseg on the halves
        for b in bands:
            C[b][tag] = cb[b]
    del X
    gc.collect()
    for phase in ("task_test", "rest_post"):
        X = _channels_first(load_timeseries(pat, phase, SEEG_DATAPATH))
        cb = complex_coherency_bands(X, fs, bd, nfull)
        for b in bands:
            C[b][phase] = cb[b]
        del X
        gc.collect()
    return C


def adjacency_from_coherency(C_band: np.ndarray) -> np.ndarray:
    """Observed |ImCoh| adjacency mean_f |Im C(f)| (symmetric, ≥0, zero-diag)."""
    A = np.abs(np.imag(C_band)).mean(axis=0)
    np.fill_diagonal(A, 0.0)
    return A


def rho_split(coph: dict[str, np.ndarray]) -> float:
    return float(spearmanr(coph["task_test"] - coph["rest_pre_A"],
                           coph["rest_post"] - coph["rest_pre_B"]).statistic)


# --------------------------------------------------------------------------- #
#  nulls                                                                        #
# --------------------------------------------------------------------------- #
def sbcrc_null(C_band_phases: dict[str, np.ndarray], n_surr: int,
               rng: np.random.Generator, A_full_post: np.ndarray,
               weights: dict) -> dict:
    """Coordinated SB-CRC floor + calibration. Shared B, independent O_φ per phase."""
    N = A_full_post.shape[0]
    triu = np.triu_indices(N, 1)
    B, Delta = shared_backbone_deviation(C_band_phases, weights)
    rho = np.empty(n_surr)
    mag_surr, resid_all, psd_all = [], [], []
    nonneg = True
    for r in range(n_surr):
        coph = {}
        for phase in PHASES:
            # psd_project=True (default): on-manifold surrogate, magnitude-stable,
            # psd residual (pre-projection min eigenvalue) returned for free.
            W, res, psd = deviation_rotated_coherency(B, Delta[phase], rng)
            resid_all.append(res)
            if np.isfinite(psd):
                psd_all.append(psd)
            if W.min() < -1e-12:
                nonneg = False
            coph[phase] = cophenetic_condensed_from_adjacency(W)
            if phase == "rest_post":
                mag_surr.append(float(np.mean(W[triu])))
        rho[r] = rho_split(coph)
    mag_obs = float(np.mean(A_full_post[triu]))
    bb_energy = float(np.linalg.norm(B)
                      / np.mean([np.linalg.norm(C_band_phases[p]) for p in PHASES]))
    return {
        "null": rho[np.isfinite(rho)],
        "mag_ratio": float(np.mean(mag_surr) / mag_obs) if mag_obs > 0 else np.nan,
        "renorm_resid": float(np.mean(resid_all)),
        "psd_resid": float(np.min(psd_all)) if psd_all else np.nan,
        "nonneg": nonneg,
        "backbone_energy_frac": bb_energy,
    }


def independent_crc_floor(C_band_phases: dict[str, np.ndarray], n_surr: int,
                          rng: np.random.Generator) -> np.ndarray:
    """Independent-CRC floor on the SAME common grid (rotate B+Δ together ⇒ ≈ 0)."""
    rho = np.empty(n_surr)
    for r in range(n_surr):
        coph = {}
        for phase in PHASES:
            W, _ = coupling_randomized_coherency(C_band_phases[phase], rng,
                                                 mode="real_congruence")
            coph[phase] = cophenetic_condensed_from_adjacency(W)
        rho[r] = rho_split(coph)
    return rho[np.isfinite(rho)]


def _stats(null: np.ndarray, obs: float) -> dict:
    if null.size == 0:
        return dict(median=np.nan, p5=np.nan, p95=np.nan, p=np.nan, z=np.nan)
    return dict(
        median=float(np.median(null)),
        p5=float(np.quantile(null, 0.05)), p95=float(np.quantile(null, 0.95)),
        p=float(np.mean(null >= obs)),
        z=float((obs - null.mean()) / (null.std() + 1e-12)),
    )


def per_cell(pat: str, band: str, C_band_phases: dict, n_surr: int,
             rng: np.random.Generator, weights: dict) -> dict:
    # within-pipeline observed (full common grid)
    A_full = {p: adjacency_from_coherency(C_band_phases[p]) for p in PHASES}
    obs_full = rho_split({p: cophenetic_condensed_from_adjacency(A_full[p]) for p in PHASES})
    # canonical observed (cached nperseg//2 halves) = audit_63
    A_canon = {p: load_phase_fc(pat, p, band) for p in PHASES}
    obs_canon = rho_split({p: cophenetic_condensed_from_adjacency(A_canon[p]) for p in PHASES})

    sb = sbcrc_null(C_band_phases, n_surr, rng, A_full["rest_post"], weights)
    crc0 = independent_crc_floor(C_band_phases, n_surr, rng)

    ssb = _stats(sb["null"], obs_full)
    s0 = _stats(crc0, obs_full)
    shared_part = ssb["median"] - s0["median"]          # backbone contribution
    task_part = obs_full - ssb["median"]                # task-specific contribution
    return {
        "patient": pat, "band": band, "N_nodes": A_full["rest_post"].shape[0],
        "obs_full": obs_full, "obs_canonical": obs_canon,
        "indep_crc_floor": s0["median"], "indep_crc_p5": s0["p5"], "indep_crc_p95": s0["p95"],
        "sbcrc_floor": ssb["median"], "sbcrc_p5": ssb["p5"], "sbcrc_p95": ssb["p95"],
        "sbcrc_p": ssb["p"], "sbcrc_z": ssb["z"],
        "shared_backbone_part": shared_part, "task_specific_part": task_part,
        "sbcrc_mag_ratio": sb["mag_ratio"], "sbcrc_renorm_resid": sb["renorm_resid"],
        "sbcrc_psd_resid": sb["psd_resid"], "sbcrc_nonneg": sb["nonneg"],
        "backbone_energy_frac": sb["backbone_energy_frac"],
        "_sb_null": sb["null"], "_crc0_null": crc0,
    }


# --------------------------------------------------------------------------- #
#  cohort summary + figure                                                      #
# --------------------------------------------------------------------------- #
def cohort_summary(rows: list[dict], bands: list[str]) -> pd.DataFrame:
    """Per band: the decomposition + paired one-sided Wilcoxon (obs > SB-CRC floor)."""
    out = []
    for band in bands:
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            continue
        obs = np.array([r["obs_full"] for r in sub])
        sb_floor = np.array([r["sbcrc_floor"] for r in sub])
        crc0 = np.array([r["indep_crc_floor"] for r in sub])
        try:
            _, p_coord = wilcoxon(obs - sb_floor, alternative="greater")
        except ValueError:
            p_coord = float("nan")
        n_above = int(sum(r["sbcrc_p"] < 0.05 for r in sub))
        out.append({
            "band": band, "n": len(sub),
            "obs_median": float(np.median(obs)),
            "indep_crc_floor_median": float(np.median(crc0)),
            "sbcrc_floor_median": float(np.median(sb_floor)),
            "shared_backbone_median": float(np.median([r["shared_backbone_part"] for r in sub])),
            "task_specific_median": float(np.median([r["task_specific_part"] for r in sub])),
            "coord_wilcoxon_p": float(p_coord),
            "n_above_sbcrc": f"{n_above}/{len(sub)}",
            "mag_ratio_median": float(np.median([r["sbcrc_mag_ratio"] for r in sub])),
            "verdict": ("task_specific" if p_coord < 0.05 else "backbone_explained"),
        })
    return pd.DataFrame(out)


def make_decomposition(rows: list[dict], bands: list[str], path: Path) -> None:
    """Per-band strip: each patient's observed dot over the SB-CRC floor bar, with
    the independent-CRC floor (≈0) as the reference tick. Dot right of the bar =
    task-specific beyond the stable backbone."""
    use_lrg_style()
    bands = [b for b in bands if any(r["band"] == b for r in rows)]
    fig, axes = plt.subplots(1, len(bands), figsize=(4.6 * len(bands), 4.8),
                             squeeze=False)
    c0, csb, cobs = "#bdbdbd", "#1f6e4d", "#1f3d6e"
    for ax, band in zip(axes[0], bands):
        sub = sorted([r for r in rows if r["band"] == band],
                     key=lambda r: r["obs_full"])
        ys = np.arange(len(sub))
        for y, r in zip(ys, sub):
            ax.hlines(y + 0.16, r["indep_crc_p5"], r["indep_crc_p95"], color=c0,
                      lw=4.0, alpha=0.95, zorder=1)
            ax.hlines(y - 0.16, r["sbcrc_p5"], r["sbcrc_p95"], color=csb,
                      lw=4.5, alpha=0.95, zorder=1)
            ax.plot(r["sbcrc_floor"], y - 0.16, "|", color="#0d3d28", ms=9, mew=1.1)
            ax.plot(r["obs_full"], y, "o", color=cobs, ms=7, mec="white",
                    mew=0.9, zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels([r["patient"].replace("Pat_", "P") for r in sub])
        ax.axvline(0, color="0.6", lw=0.7, ls="--", zorder=0)
        ax.set_xlabel(r"$\rho_{\mathrm{split}}$")
        ax.set_title(BRAIN_BAND_TEX_DICT[band])
        ax.spines[["top", "right"]].set_visible(False)
    handles = [
        plt.Line2D([0], [0], color=c0, lw=4.0, label="independent-CRC floor (≈0)"),
        plt.Line2D([0], [0], color=csb, lw=4.5, label="SB-CRC floor — shared backbone (P5–P95)"),
        plt.Line2D([0], [0], marker="o", color=cobs, lw=0, mec="white",
                   label="observed"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, bbox_inches="tight", transparent=True)
    plt.close(fig)


def write_readme(cohort: pd.DataFrame, rows: list[dict], runtime: float,
                 backbone: str, weights: dict, out_dir: Path) -> None:
    lines = [
        "---", "name: coordinated_cross_phase_null", "type: report",
        "era: COHORT_N10", "status: current", "created: 2026-06-11",
        "updated: 2026-06-11", "pointers:",
        "  - .agents/guides/task-persistence-investigation/2026-06-11_coordinated-cross-phase-null.md",
        "  - data/audit/crc_cohort_bracket/",
        "---", "",
        f"# SB-CRC coordinated cross-phase null — cohort ({backbone} backbone)",
        "",
        f"**Backbone = `{backbone}`** ("
        + ("PRIMARY: B = rest_pre only, task/post excluded — the backbone cannot "
           "peek at the task, so obs > floor cleanly means task-specific persistence."
           if backbone == "baseline" else
           "CONSERVATIVE upper bound: B absorbs task+post common structure, so it can "
           "hide a real trace.") + ")",
        "",
        "**Head.** Decomposes the §5.3 `ρ_split` into a shared-backbone part and a "
        "task-specific part by holding the cross-phase common coherency `B` fixed "
        "across all surrogate phases and rotating only each phase's deviation "
        "(independent `O_φ`). The matched-strength and CRC nulls randomize each "
        "phase independently, so their floor is ≈0 (the independent-CRC column); "
        "SB-CRC's floor lifts to the shared-backbone contribution. The coordinated "
        "test asks whether observed exceeds the SB-CRC floor — i.e. whether the "
        "trace is task-specific beyond a task-independent stable fingerprint.",
        "",
        "| band | n | obs med | indep-CRC floor | SB-CRC floor | shared part | task part | coord Wilcoxon p | n>SB-CRC | mag-ratio | verdict |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['n']} | {r['obs_median']:+.3f} "
            f"| {r['indep_crc_floor_median']:+.3f} | {r['sbcrc_floor_median']:+.3f} "
            f"| {r['shared_backbone_median']:+.3f} | {r['task_specific_median']:+.3f} "
            f"| {r['coord_wilcoxon_p']:.4f} | {r['n_above_sbcrc']} "
            f"| {r['mag_ratio_median']:.2f} | **{r['verdict']}** |")
    lines += [
        "", "## Reading",
        "- **task_specific** — observed exceeds the SB-CRC floor across the cohort "
        "(coord Wilcoxon p < 0.05): the §5.3 trace is task-specific beyond a "
        "task-independent stable backbone.",
        "- **backbone_explained** — observed sits inside the SB-CRC floor: the "
        "cross-phase ρ_split is substantially a stable-fingerprint effect. FLAG for "
        "adjudication, do NOT auto-reframe the manuscript.",
        "- `shared part` = SB-CRC floor − independent floor (the stable-backbone "
        "contribution); `task part` = observed − SB-CRC floor.",
        "- Off-manifold caveat is cosmetic (pilot: PSD-projected readout ρ=0.99); "
        "magnitude faithful (`mag-ratio` ≈ 1, real-orthogonal congruence preserves "
        "‖Im Δ‖).",
        "",
        "## Per-cell",
        "| patient | band | obs(full) | obs(canon) | SB-CRC floor | shared | task | SB-CRC p | magR | PSD resid |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['patient']} | {r['band']} | {r['obs_full']:+.3f} "
            f"| {r['obs_canonical']:+.3f} | {r['sbcrc_floor']:+.3f} "
            f"| {r['shared_backbone_part']:+.3f} | {r['task_specific_part']:+.3f} "
            f"| {r['sbcrc_p']:.3f} | {r['sbcrc_mag_ratio']:.2f} "
            f"| {r['sbcrc_psd_resid']:.1e} |")
    lines += [
        "",
        f"R = {N_SURR} per null per phase, seed {SEED}, backbone={backbone} weights "
        f"{weights}. Common full-nperseg grid. PSD projection ON. Runtime {runtime:.0f}s.",
        "- Build: `scripts/01_compute/audit/audit_99_coordinated_null.py "
        f"--backbone {backbone}`",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=TARGET_BANDS)
    ap.add_argument("--n-surrogates", type=int, default=N_SURR)
    ap.add_argument("--backbone", choices=list(WEIGHT_SCHEMES), default="baseline",
                    help="shared-backbone weighting (baseline = primary; "
                         "equal_condition = conservative upper bound)")
    args = ap.parse_args()
    weights = WEIGHT_SCHEMES[args.backbone]

    out_dir = OUT / args.backbone
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    rows = []
    for pat in args.patients:
        try:
            C = build_coherency_common_grid(pat, args.bands)
        except Exception as e:
            print(f"[audit_99] SKIP {pat}: coherency build failed: {e}")
            continue
        for band in args.bands:
            tc = time.time()
            try:
                r = per_cell(pat, band, C[band], args.n_surrogates, rng, weights)
            except Exception as e:
                print(f"[audit_99] SKIP {pat}/{band}: {e}")
                continue
            rows.append(r)
            print(f"[audit_99] [{args.backbone}] {pat}/{band}: obs_full={r['obs_full']:+.3f} "
                  f"(canon={r['obs_canonical']:+.3f}) | indep floor="
                  f"{r['indep_crc_floor']:+.3f} | SB-CRC floor={r['sbcrc_floor']:+.3f} "
                  f"p={r['sbcrc_p']:.3f} | shared={r['shared_backbone_part']:+.3f} "
                  f"task={r['task_specific_part']:+.3f} | magR={r['sbcrc_mag_ratio']:.2f} "
                  f"psd={r['sbcrc_psd_resid']:.1e} ({time.time()-tc:.1f}s)")
        del C
        gc.collect()

    runtime = time.time() - t0
    print(f"[audit_99] [{args.backbone}] total {runtime:.1f}s, {len(rows)} cells")
    if not rows:
        return
    per = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")}
                        for r in rows])
    per.to_csv(out_dir / "per_cell.csv", index=False)
    cohort = cohort_summary(rows, args.bands)
    cohort.to_csv(out_dir / "cohort_summary.csv", index=False)
    make_decomposition(rows, args.bands, fig_dir / "decomposition.pdf")
    write_readme(cohort, rows, runtime, args.backbone, weights, out_dir)
    print(cohort.to_string(index=False))
    print(f"[audit_99] outputs at {out_dir}")


if __name__ == "__main__":
    main()
