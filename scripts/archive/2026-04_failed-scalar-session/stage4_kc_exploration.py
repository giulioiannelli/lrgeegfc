#!/usr/bin/env python3
"""Stage 4 — Kendall-Colijn per-band λ-sweep characterization.

Not a pass/fail test; an exploration. For every (patient, band, phase pair)
we compute the KC distance over a fine λ grid, cache the underlying
(m, M) vectors once per tree for speed, then ask three questions per band:

1. **Does KC see the difference between phases at all?**
   Compare the three raw inter-phase distances
   ``d(rest_pre, rest_post)``, ``d(task_test, rest_post)``,
   ``d(rest_pre, task_test)`` across λ. If KC is informative, these
   should differ systematically; if it's noise-dominated, they collapse.

2. **How does each band behave across λ?**
   Per-band curves of rank-biserial for two directional contrasts:
   - H2a: ``s = d(rpre, rpost) − d(tt, rpost)`` (positive = rpost closer to tt)
   - H2b: ``s = d(rpre, tt) − d(tt, rpost)`` (positive = task-post stays close to
     task, while task-pre diverges)
   and of ``n_positive/9`` per band. Watch for bands whose rb-peak sits at
   distinct λ regions — a topology-favoring band vs a heights-favoring band
   vs a band that's noise at all λ is a real distinguishing characterization.

3. **Is the effect robust across patients per-band?**
   Per-patient per-band heatmap of H2a contrast at λ ∈ {0, 0.5, 1}. Lets you
   see which patients drive which bands' trends at each λ.

No pre-registered pass criterion. Outputs are descriptive; narrative in
`data/reports/imcoh_vi/stage4_kc_exploration.md` lays out what each band
does across λ.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import REPORTS_ROOT
from lrg_eegfc.utils.metrics import kc_vectors
from lrg_eegfc.workflow.lrg import load_lrg_result


# --- shared stats helpers ---
_shared_path = ROOT / "scripts" / "01_compute" / "_shared.py"
_spec = importlib.util.spec_from_file_location("_shared_stage4", _shared_path)
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)
wilcoxon_z = _shared.wilcoxon_z
rank_biserial = _shared.rank_biserial
boot_ci_mean = _shared.boot_ci_mean
bh_fdr = _shared.bh_fdr


FC_METHOD = "imcoh_abs"
PATIENTS = [p for p in PATIENTS_4PHASE if p != "Pat_14"]
PHASE_PAIRS = [
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
    ("rest_pre", "task_test"),
]
LAMBDAS = np.round(np.linspace(0.0, 1.0, 21), 3)  # 21 points, 0.05 step
OUT_DIR = REPORTS_ROOT / "imcoh_vi"
FIG_DIR = OUT_DIR / "figures"


def load_vectors() -> dict:
    """Per (patient, phase, band), return the (m, M) KC vectors."""
    cache = {}
    needed_phases = sorted({p for pair in PHASE_PAIRS for p in pair})
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for phase in needed_phases:
                lrg = load_lrg_result(pat, phase, band, FC_METHOD)
                if lrg is None:
                    print(f"  MISSING: {pat} {phase} {band}")
                    continue
                m, M = kc_vectors(lrg.linkage_matrix)
                cache[(pat, phase, band)] = (m, M)
    return cache


def kc_distance_from_cache(cache: dict, key1, key2, lam: float,
                            m_norm: dict, M_norm: dict) -> float:
    """KC L2 distance from cached normalized vectors. m_norm/M_norm hold
    the cohort-normalization factors per (patient, band) pair.
    """
    m1, M1 = cache[key1]
    m2, M2 = cache[key2]
    # Normalization: pool across the two trees being compared.
    # We use per-(patient, band) normalization so that inter-phase
    # comparisons within a (patient, band) are on a common scale.
    m_max, M_max = m_norm[key1[::2]], M_norm[key1[::2]]  # (patient, band)
    m1n, m2n = m1 / m_max, m2 / m_max
    M1n, M2n = M1 / M_max, M2 / M_max
    v1 = (1 - lam) * m1n + lam * M1n
    v2 = (1 - lam) * m2n + lam * M2n
    return float(np.linalg.norm(v1 - v2))


def build_norms(cache: dict) -> tuple[dict, dict]:
    """Per (patient, band), pool m_max and M_max across all 4 phases."""
    m_max: dict[tuple[str, str], float] = {}
    M_max: dict[tuple[str, str], float] = {}
    for (pat, phase, band), (m, M) in cache.items():
        k = (pat, band)
        m_max[k] = max(m_max.get(k, 0), float(m.max()))
        M_max[k] = max(M_max.get(k, 0), float(M.max()))
    # Avoid zeros
    for k in m_max:
        if m_max[k] == 0:
            m_max[k] = 1.0
        if M_max[k] == 0:
            M_max[k] = 1e-12
    return m_max, M_max


def compute_distances(cache: dict, m_norm, M_norm) -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for p1, p2 in PHASE_PAIRS:
                k1, k2 = (pat, p1, band), (pat, p2, band)
                if k1 not in cache or k2 not in cache:
                    continue
                for lam in LAMBDAS:
                    d = kc_distance_from_cache(cache, k1, k2, float(lam),
                                                m_norm, M_norm)
                    rows.append({
                        "patient": pat, "band": band,
                        "phase1": p1, "phase2": p2,
                        "lam": float(lam), "d": d,
                    })
    return pd.DataFrame(rows)


def build_contrasts(dist_df: pd.DataFrame) -> pd.DataFrame:
    """Per (patient, band, lam): H2a = d(pre,post)-d(tt,post); H2b = d(pre,tt)-d(tt,post)."""
    rows = []
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for lam in LAMBDAS:
                sub = dist_df[(dist_df["patient"] == pat) &
                              (dist_df["band"] == band) &
                              (dist_df["lam"].round(3) == round(float(lam), 3))]
                def pick(p1, p2):
                    r = sub[(sub["phase1"] == p1) & (sub["phase2"] == p2)]
                    return float(r["d"].values[0]) if not r.empty else np.nan
                d_pre_post = pick("rest_pre", "rest_post")
                d_tt_post = pick("task_test", "rest_post")
                d_pre_tt = pick("rest_pre", "task_test")
                if not (np.isfinite(d_pre_post) and np.isfinite(d_tt_post)
                        and np.isfinite(d_pre_tt)):
                    continue
                rows.append({
                    "patient": pat, "band": band, "lam": float(lam),
                    "d_pre_post": d_pre_post, "d_tt_post": d_tt_post,
                    "d_pre_tt": d_pre_tt,
                    "H2a": d_pre_post - d_tt_post,
                    "H2b": d_pre_tt - d_tt_post,
                })
    return pd.DataFrame(rows)


def cohort_stats_per_band_lam(contrast_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in BRAIN_BANDS_NAMES:
        for lam in LAMBDAS:
            sub = contrast_df[(contrast_df["band"] == band) &
                              (contrast_df["lam"].round(3) == round(float(lam), 3))]
            for contrast_name in ("H2a", "H2b"):
                s = sub[contrast_name].values.astype(float)
                s = s[np.isfinite(s)]
                n = int(len(s))
                n_pos = int((s > 0).sum())
                z, p = wilcoxon_z(s)
                rb = rank_biserial(s)
                rows.append({
                    "band": band, "lam": float(lam),
                    "contrast": contrast_name, "n": n, "n_positive": n_pos,
                    "median": float(np.median(s)) if n else np.nan,
                    "z": z, "p": p, "rb": rb,
                })
    return pd.DataFrame(rows)


# ----- Figures ------------------------------------------------------------

def plot_lambda_curves(stats_df: pd.DataFrame, out: Path) -> None:
    """Per band: rb and n_positive as a function of λ for H2a and H2b."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[i // 3, i % 3]
        for contrast_name, color, marker in [("H2a", "C0", "o"), ("H2b", "C3", "^")]:
            s = stats_df[(stats_df["band"] == band) &
                         (stats_df["contrast"] == contrast_name)].sort_values("lam")
            ax.plot(s["lam"], s["rb"], marker=marker, color=color,
                    lw=1.5, markersize=4, label=f"{contrast_name}")
            # highlight points with p<0.05 (uncorrected) as filled circles
            sig = s[s["p"] < 0.05]
            if not sig.empty:
                ax.scatter(sig["lam"], sig["rb"], s=80, facecolors="none",
                           edgecolors=color, lw=1.8, zorder=5)
        ax.axhline(0, color="k", lw=0.6, ls="--")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-1.05, 1.05)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=11)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8, loc="best")
        if i >= 3:
            ax.set_xlabel(r"$\lambda$ (0 = topology, 1 = heights)")
        if i % 3 == 0:
            ax.set_ylabel("rank-biserial")
    fig.suptitle(
        f"KC(λ) rank-biserial per band — H2a vs H2b contrasts\n"
        f"|ImCoh|, n={len(PATIENTS)}; open circles = p<0.05 uncorrected",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_distance_distributions(dist_df: pd.DataFrame, out: Path) -> None:
    """Per band: boxplots of the 3 phase-pair distances at λ=0, 0.5, 1.0."""
    lam_picks = [0.0, 0.5, 1.0]
    pair_labels = {
        ("rest_pre", "rest_post"): "pre–post",
        ("task_test", "rest_post"): "tt–post",
        ("rest_pre", "task_test"): "pre–tt",
    }
    fig, axes = plt.subplots(len(BRAIN_BANDS_NAMES), len(lam_picks),
                             figsize=(11, 2.2 * len(BRAIN_BANDS_NAMES)),
                             sharex=False, sharey=False)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        for j, lam in enumerate(lam_picks):
            ax = axes[i, j]
            sub = dist_df[(dist_df["band"] == band) &
                          (dist_df["lam"].round(3) == round(lam, 3))]
            box_data = []
            labels = []
            for p1, p2 in PHASE_PAIRS:
                vals = sub[(sub["phase1"] == p1) & (sub["phase2"] == p2)]["d"].values
                box_data.append(vals)
                labels.append(pair_labels[(p1, p2)])
            bp = ax.boxplot(box_data, labels=labels, showfliers=False,
                            patch_artist=True, widths=0.6)
            colors = ["#cccccc", "#f19b8f", "#7fafd6"]
            for patch, c in zip(bp["boxes"], colors):
                patch.set_facecolor(c)
                patch.set_alpha(0.7)
            # Strip of individual points
            for k, vals in enumerate(box_data):
                ax.scatter(np.full(len(vals), k + 1) + np.random.uniform(-0.08, 0.08, len(vals)),
                           vals, s=10, color="k", alpha=0.5, zorder=3)
            if j == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], rotation=0,
                              ha="right", va="center", fontsize=11)
            if i == 0:
                ax.set_title(f"λ={lam}", fontsize=10)
            ax.grid(alpha=0.3, axis="y")
            ax.tick_params(axis="x", labelsize=8)
            ax.tick_params(axis="y", labelsize=7)
    fig.suptitle(
        f"KC distance distributions per band × phase pair × λ\n"
        f"|ImCoh|, n={len(PATIENTS)} (points = patients)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_per_patient_detail(contrast_df: pd.DataFrame, out: Path) -> None:
    """Per band: H2a per patient across λ (heatmap). Rows = patients."""
    fig, axes = plt.subplots(len(BRAIN_BANDS_NAMES), 1,
                             figsize=(9, 1.6 * len(BRAIN_BANDS_NAMES)),
                             sharex=True)
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[i]
        mat = np.full((len(PATIENTS), len(LAMBDAS)), np.nan)
        for r, pat in enumerate(PATIENTS):
            for c, lam in enumerate(LAMBDAS):
                sub = contrast_df[(contrast_df["band"] == band) &
                                  (contrast_df["patient"] == pat) &
                                  (contrast_df["lam"].round(3) == round(float(lam), 3))]
                if not sub.empty:
                    mat[r, c] = float(sub["H2a"].values[0])
        # Symmetric color range per band
        amax = np.nanmax(np.abs(mat)) or 1e-9
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
                       vmin=-amax, vmax=amax,
                       extent=[LAMBDAS[0] - 0.025, LAMBDAS[-1] + 0.025,
                               len(PATIENTS) - 0.5, -0.5])
        ax.set_yticks(range(len(PATIENTS)))
        ax.set_yticklabels(PATIENTS, fontsize=8)
        ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], rotation=0, ha="right",
                      va="center", fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.7, label="H2a (d_pre_post − d_tt_post)")
    axes[-1].set_xlabel(r"$\lambda$ (0 = topology, 1 = heights)")
    fig.suptitle(f"Per-patient H2a contrast across λ, per band\n"
                 f"|ImCoh|, n={len(PATIENTS)}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


# ----- Narrative ---------------------------------------------------------

def write_report(stats_df: pd.DataFrame, out: Path) -> None:
    lines = [
        "# Stage 4 — KC per-band λ-sweep characterization",
        "",
        f"Date: 2026-04-24. FC: `{FC_METHOD}`. Patients: n={len(PATIENTS)} "
        f"(Pat_14 excluded). 6 bands × {len(LAMBDAS)} λ values × 3 phase pairs.",
        "",
        "## What this asks",
        "",
        "Not a pre-registered pass/fail test. An exploration of whether the "
        "Kendall-Colijn distance sees structure in our dendrograms, and whether "
        "the bands distinguish themselves along the λ axis (pure topology at "
        "λ=0, pure heights at λ=1, a blend in between).",
        "",
        "Contrasts reported:",
        "- **H2a** = d(rest_pre, rest_post) − d(task_test, rest_post); "
        "positive ⇒ rest_post closer to task_test than to rest_pre.",
        "- **H2b** = d(rest_pre, task_test) − d(task_test, rest_post); "
        "positive ⇒ rest_post stays near task_test while rest_pre diverges.",
        "",
        "Per-patient normalization: (m, M) vectors pooled per (patient, band) "
        "across the 4 phases before blending. Keeps the intra-(patient, band) "
        "phase-pair distances on a comparable scale.",
        "",
        "## Per-band λ-extrema summary",
        "",
        "For each band, find the λ where rb is maximized for H2a and H2b; "
        "report rb, n_positive/9, and raw Wilcoxon p at that λ. Also report "
        "rb at λ=0 (topology) and λ=1 (heights) for reference.",
        "",
    ]
    # Build per-band table
    summary_rows = []
    for band in BRAIN_BANDS_NAMES:
        row = {"band": BRAIN_BAND_TEX_DICT[band]}
        for contrast in ("H2a", "H2b"):
            sub = stats_df[(stats_df["band"] == band) &
                           (stats_df["contrast"] == contrast)].set_index("lam")
            if sub.empty:
                continue
            best_lam = float(sub["rb"].astype(float).idxmax())
            best = sub.loc[best_lam]
            row[f"{contrast}_best_lam"] = best_lam
            row[f"{contrast}_best_rb"] = float(best["rb"])
            row[f"{contrast}_best_npos"] = int(best["n_positive"])
            row[f"{contrast}_best_p"] = float(best["p"])
            # End-points
            if 0.0 in sub.index:
                row[f"{contrast}_rb_lam0"] = float(sub.loc[0.0, "rb"])
            if 1.0 in sub.index:
                row[f"{contrast}_rb_lam1"] = float(sub.loc[1.0, "rb"])
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    lines.append(summary.round(4).to_markdown(index=False))
    lines.append("")
    lines.append("## Interpretation cheat sheet")
    lines.append("")
    lines.append(
        "- **Band whose rb-peak sits near λ=0** — topology-carrying: "
        "the branching pattern of rest_post resembles task_test more than "
        "rest_pre's, irrespective of heights.")
    lines.append(
        "- **Band whose rb-peak sits near λ=1** — heights-carrying: the "
        "pair-height structure of rest_post resembles task_test more than "
        "rest_pre's. Closely related to H2c (pair-averaged cophenetic drift).")
    lines.append(
        "- **Band whose rb-peak sits at intermediate λ** — blend: both "
        "topology and heights contribute; a scalar that's genuinely more "
        "informative than either H2c or a pure-topology metric.")
    lines.append(
        "- **Band where rb stays near 0 everywhere** — KC sees no "
        "cross-phase effect in that band at scalar resolution.")
    lines.append("")
    lines.append("## Per-band full λ-stats (rb, n_positive, p)")
    lines.append("")
    for band in BRAIN_BANDS_NAMES:
        lines.append(f"### {band}")
        lines.append("")
        for contrast in ("H2a", "H2b"):
            sub = stats_df[(stats_df["band"] == band) &
                           (stats_df["contrast"] == contrast)].sort_values("lam")
            lines.append(f"#### {contrast} ({band})")
            lines.append("")
            lines.append(sub[["lam", "n", "n_positive", "median", "z", "p", "rb"]]
                          .round(4).to_markdown(index=False))
            lines.append("")
    out.write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Computing KC (m, M) vectors for all trees "
          f"(cache size: {len(PATIENTS)}×3×6 = {len(PATIENTS)*3*6})…")
    cache = load_vectors()
    m_norm, M_norm = build_norms(cache)
    print(f"  Cached (m, M) for {len(cache)} trees.")

    print(f"[2/5] KC distances over {len(LAMBDAS)} λ values × 3 phase pairs × 9 patients × 6 bands…")
    dist_df = compute_distances(cache, m_norm, M_norm)
    dist_csv = OUT_DIR / "stage4_kc_distances.csv"
    dist_df.to_csv(dist_csv, index=False)
    print(f"  {len(dist_df)} rows → {dist_csv}")

    print("[3/5] Building H2a/H2b contrasts at each λ…")
    contrast_df = build_contrasts(dist_df)
    con_csv = OUT_DIR / "stage4_kc_contrasts.csv"
    contrast_df.to_csv(con_csv, index=False)
    print(f"  {len(contrast_df)} rows → {con_csv}")

    print("[4/5] Cohort stats per (band, λ, contrast)…")
    stats_df = cohort_stats_per_band_lam(contrast_df)
    stats_csv = OUT_DIR / "stage4_kc_stats.csv"
    stats_df.to_csv(stats_csv, index=False)
    print(f"  {len(stats_df)} rows → {stats_csv}")

    print("[5/5] Figures + report…")
    plot_lambda_curves(stats_df, FIG_DIR / "stage4_kc_lambda_per_band.pdf")
    plot_distance_distributions(dist_df, FIG_DIR / "stage4_kc_distance_distributions.pdf")
    plot_per_patient_detail(contrast_df, FIG_DIR / "stage4_kc_per_patient.pdf")
    write_report(stats_df, OUT_DIR / "stage4_kc_exploration.md")

    # Peak-summary stdout
    print("\n" + "=" * 100)
    print("PER-BAND λ-EXTREMA SUMMARY (H2a contrast)")
    print("=" * 100)
    rows = []
    for band in BRAIN_BANDS_NAMES:
        sub = stats_df[(stats_df["band"] == band) & (stats_df["contrast"] == "H2a")]
        if sub.empty:
            continue
        idx = sub["rb"].astype(float).idxmax()
        best = sub.loc[idx]
        row = {
            "band": band,
            "rb@lam0": float(sub.loc[sub["lam"] == 0.0, "rb"].values[0]),
            "rb@lam1": float(sub.loc[sub["lam"] == 1.0, "rb"].values[0]),
            "rb_peak_lam": float(best["lam"]),
            "rb_peak": float(best["rb"]),
            "npos@peak": int(best["n_positive"]),
            "p@peak": float(best["p"]),
        }
        rows.append(row)
    summary = pd.DataFrame(rows)
    print(summary.round(4).to_string(index=False))
    print("=" * 100)


if __name__ == "__main__":
    main()
