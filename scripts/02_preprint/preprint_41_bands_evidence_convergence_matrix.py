#!/usr/bin/env python3
r"""Per-band evidence-convergence matrix — the one-glance synthesis panel.

Rows = the six bands. Left strip = the per-patient spread (who traces, n=10).
Verdict columns, grouped so the reader cannot misread them as one crossed test:

    [ per-patient trace ] | TRACE EXISTS          | COGNITIVE CONTENT (rho^coph) | ANATOMY
                          |  rho^coph   Grassmann  |  encoding      inference    |  home

  * TRACE EXISTS — two INDEPENDENT probes that the band's trace is real:
      rho^coph  : per-pair cophenetic gate (paired Wilcoxon, matched-strength).
      Grassmann : leading-subspace rotation (cluster-mass permutation + LOO).
    They are complementary existence probes, NOT a content decomposition —
    hence the group bracket. beta clears both; alpha is cophenetic-only.
  * COGNITIVE CONTENT — the four-phase cophenetic decomposition (audit_103):
      encoding  : T_learn      = rho(D_taskLearn-D_preA , D_post-D_preB)
      inference : T_infspec.e  = partial-rho(D_test-D_taskLearn , persist | encoding)
    Both are cophenetic-only (that is the point of the bracket). beta is the
    only band with an inference-specific component beyond encoding.
  * ANATOMY — the matched-strength localizer home (audit_83/_83c LOO).
    beta -> OFC is the only LOO-robust home; theta/high-gamma "localize" on a
    band with ~zero net trace (hatched = localization WITHOUT a cohort trace).

Colour = verdict tier {verified/locked, subset·directional, null, loc-without-
trace}. Project rule: NO numeric text inside the axes (region tags + glyphs are
categorical); the full numeric table is printed to stdout. Nothing recomputed
except the cohort Wilcoxon of the arc per-patient cache (printed for audit).

Sources (all locked, n=10, matched-strength R=200/1000):
    rho^coph gate   data/audit/matched_strength_surrogate_split_baseline/{cohort_summary,per_patient_per_band}.csv
    Grassmann       data/audit/grassmann_cluster_extent/cohort_summary.csv
    arc content     data/audit/consolidation_arc/arc_null_per_patient.csv
    localization    data/audit/localization_atlas/{per_band_taxonomy_verdict,carrier_loo}.csv
Synthesis: .agents/reports/2026-06-26_per-band-phenomenology-vision.md
  (supersedes loc_verdict.tier for alpha [cohort/no-home] and low_gamma
   [consistent -> strong-subset, PFC LOO-failed]; tiers derived here from the
   quantitative columns + that report, NOT from the stale `tier` column).

Output: data/preprint/figures/all_bands/fig_bands_evidence_convergence_matrix.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08",
          "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ALPHA = 0.05

# verdict tiers ---------------------------------------------------------------
TIER_COLOR = {
    "verified": "#1f7a1f",   # locked / clears the mandatory null
    "subset":   "#e67e22",   # positive-but-not-cohort / directional / LOO-fragile
    "null":     "#c9cdcf",   # tested, no trace
}
TIER_GLYPH = {"verified": "✓", "subset": "~", "null": "·"}
LOCWO_FACE, LOCWO_HATCH, LOCWO_EDGE = "#f4f4f4", "////", "#e67e22"  # localizes w/o trace
# per-patient strip
FLAG_COLOR = {"trace": "#1f7a1f", "null": "#dcdfe0", "anti": "#c0392b"}


def _cohort_wilcoxon_greater(obs: np.ndarray, ref: np.ndarray) -> float:
    """One-sided paired Wilcoxon, H1: median(obs-ref) > 0 (audit_103 convention)."""
    d = np.asarray(obs) - np.asarray(ref)
    if np.allclose(d, 0):
        return 1.0
    return float(wilcoxon(obs, ref, alternative="greater", zero_method="wilcox").pvalue)


def _loo_robust(loo: pd.DataFrame, band: str, carrier: str) -> bool:
    sub = loo[(loo.band == band) & (loo.carrier == carrier)]
    if sub.empty:
        return False
    return bool((sub.carrier_q <= ALPHA).all())


def main() -> None:
    A = ROOT / "data" / "audit"
    coph = pd.read_csv(A / "matched_strength_surrogate_split_baseline"
                       / "cohort_summary.csv").set_index("band")
    perpat = pd.read_csv(A / "matched_strength_surrogate_split_baseline"
                         / "per_patient_per_band.csv")
    grass = pd.read_csv(A / "grassmann_cluster_extent"
                        / "cohort_summary.csv").set_index("band")
    arc = pd.read_csv(A / "consolidation_arc" / "arc_null_per_patient.csv")
    loc = pd.read_csv(A / "localization_atlas"
                      / "per_band_taxonomy_verdict.csv").set_index("band")
    loo = pd.read_csv(A / "localization_atlas" / "carrier_loo.csv")

    # ---- per-patient trace flags (matched-strength individual clearance) -----
    pflag = {}
    prho = {}
    for b in BANDS:
        sub = arc  # placeholder to keep linters calm
        for pt in COHORT:
            row = perpat[(perpat.patient == pt) & (perpat.band == b)]
            if row.empty:
                pflag[(b, pt)] = "null"; prho[(b, pt)] = np.nan; continue
            p1 = float(row.obs_p_one_sided.iloc[0])
            prho[(b, pt)] = float(row.obs_rho.iloc[0])
            pflag[(b, pt)] = ("trace" if p1 < ALPHA
                              else "anti" if p1 > 1 - ALPHA else "null")

    # ---- arc cohort content verdicts (compute + audit-print) ----------------
    enc_p, inf_p = {}, {}
    for b in BANDS:
        s = arc[arc.band == b].set_index("patient").loc[COHORT]
        enc_p[b] = _cohort_wilcoxon_greater(s.T_learn_obs.values,
                                            s.T_learn_surr_p50.values)
        inf_p[b] = _cohort_wilcoxon_greater(s.T_infspec_pe_obs.values,
                                            s.T_infspec_pe_surr_p50.values)

    # ---- assemble per-band tier for each column -----------------------------
    def trace_tier(b):
        p = float(coph.loc[b, "paired_wilcoxon_p"]); m = float(coph.loc[b, "obs_median_rho"])
        return "verified" if (p < ALPHA and m > 0) else "subset" if m > 0 else "null"

    def grass_tier(b):
        pm = float(grass.loc[b, "cluster_p_cluster_mass"])
        lo = float(grass.loc[b, "cluster_p_mass_loo_max"])
        return "verified" if (pm < ALPHA and lo < ALPHA) else "subset" if pm < ALPHA else "null"

    def enc_tier(b):
        p = enc_p[b]
        return "verified" if p < ALPHA else "subset" if p < 0.10 else "null"

    def inf_tier(b):
        return "verified" if inf_p[b] < ALPHA else "null"

    def home_cell(b):
        carriers = str(loc.loc[b, "robust_carriers"]).strip()
        locwo = bool(loc.loc[b, "loc_without_cohort_trace"])
        tr = trace_tier(b)
        if carriers in ("", "—", "-", "nan"):
            return "null", "—", False
        top = carriers.split(",")[0].strip()
        tag = {"OFC": "OFC", "PFC": "PFC", "MTL": "MTL",
               "lateral_temporal": "MTL", "parietal": "par.",
               "occipital": "occ.", "insula": "ins."}.get(top, top[:4])
        if locwo:
            return "locwo", tag, True                      # theta, high_gamma
        if _loo_robust(loo, b, top) and tr == "verified":
            return "verified", tag, False                  # beta -> OFC
        return "subset", tag + "*", False                  # low_gamma -> PFC* (LOO-fail)

    COLS = ["rho", "grass", "enc", "inf", "home"]
    tiers = {b: {"rho": trace_tier(b), "grass": grass_tier(b),
                 "enc": enc_tier(b), "inf": inf_tier(b)} for b in BANDS}
    homes = {b: home_cell(b) for b in BANDS}

    # ---- stdout numeric table (project rule: numbers here, colour in figure) -
    print("\n=== per-band evidence-convergence (n=10, matched-strength) ===")
    print(f"{'band':<11}{'coph_p':>9}{'med':>7}{'n/10':>6}  "
          f"{'grass_p':>9}{'grLOO':>7}  {'enc_p':>8}{'inf_p':>8}  home")
    for b in BANDS:
        cp = float(coph.loc[b, "paired_wilcoxon_p"]); cm = float(coph.loc[b, "obs_median_rho"])
        cn = str(coph.loc[b, "n_above_surrogate"])
        gp = float(grass.loc[b, "cluster_p_cluster_mass"]); gl = float(grass.loc[b, "cluster_p_mass_loo_max"])
        ht, tag, hw = homes[b]
        note = " [loc w/o trace]" if hw else ""
        print(f"{b:<11}{cp:>9.4f}{cm:>+7.3f}{cn:>6}  {gp:>9.4f}{gl:>7.3f}  "
              f"{enc_p[b]:>8.4f}{inf_p[b]:>8.4f}  {tag} ({ht}){note}")
    print("\nper-patient trace flags (v=trace p<.05, a=anti p>.95, .=null):")
    print(f"{'band':<11}" + "".join(f"{p.split('_')[1]:>5}" for p in COHORT))
    for b in BANDS:
        gl = {"trace": "v", "anti": "a", "null": "."}
        print(f"{b:<11}" + "".join(f"{gl[pflag[(b,p)]]:>5}" for p in COHORT))
    print("\nbeta OFC LOO: all drops q<=0.05 ->",
          _loo_robust(loo, "beta", "OFC"),
          "| low_gamma PFC LOO:", _loo_robust(loo, "low_gamma", "PFC"))

    # ---- geometry -----------------------------------------------------------
    nb = len(BANDS)
    TW = 0.17                      # per-patient tick width
    x_strip0 = 0.0
    x_strip1 = TW * len(COHORT)    # 1.70
    xv = x_strip1 + 0.55           # verdict block start
    CW = 1.0
    col_x = {c: xv + i * CW for i, c in enumerate(COLS)}
    x_end = xv + len(COLS) * CW

    fig, ax = plt.subplots(figsize=(10.0, 5.4))

    def row_y(b):                  # delta on top
        return nb - 1 - BANDS.index(b)

    # per-patient strip
    for b in BANDS:
        y = row_y(b)
        for j, pt in enumerate(COHORT):
            ax.add_patch(Rectangle((x_strip0 + j * TW, y + 0.08), TW * 0.86, 0.84,
                                    facecolor=FLAG_COLOR[pflag[(b, pt)]],
                                    edgecolor="white", linewidth=0.4, zorder=2))

    # verdict cells
    for b in BANDS:
        y = row_y(b)
        for c in ["rho", "grass", "enc", "inf"]:
            t = tiers[b][c]
            ax.add_patch(Rectangle((col_x[c], y), CW, 1.0, facecolor=TIER_COLOR[t],
                                    edgecolor="white", linewidth=1.4, zorder=2))
            fg = "white" if t in ("verified", "subset") else "#555"
            ax.text(col_x[c] + CW / 2, y + 0.5, TIER_GLYPH[t], ha="center",
                    va="center", fontsize=12, color=fg, zorder=3)
        # home cell (region tag)
        ht, tag, hw = homes[b]
        if hw:
            ax.add_patch(Rectangle((col_x["home"], y), CW, 1.0, facecolor=LOCWO_FACE,
                                   edgecolor=LOCWO_EDGE, linewidth=1.4,
                                   hatch=LOCWO_HATCH, zorder=2))
            tcol = "#8a5a12"
        else:
            ax.add_patch(Rectangle((col_x["home"], y), CW, 1.0,
                                   facecolor=TIER_COLOR[ht], edgecolor="white",
                                   linewidth=1.4, zorder=2))
            tcol = "white" if ht in ("verified", "subset") else "#555"
        ax.text(col_x["home"] + CW / 2, y + 0.5, tag, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color=tcol, zorder=3)

    # row (band) labels
    for b in BANDS:
        ax.text(-0.18, row_y(b) + 0.5, BRAIN_BAND_TEX_DICT.get(b, b),
                ha="right", va="center", fontsize=12)

    # column headers
    heads = {"rho": r"$\rho^{\mathrm{coph}}$", "grass": "Grassmann",
             "enc": "encoding", "inf": "inference", "home": "home"}
    for c in COLS:
        ax.text(col_x[c] + CW / 2, nb + 0.06, heads[c], ha="center", va="bottom",
                fontsize=9.5)
    ax.text(x_strip1 / 2, nb + 0.06, "per-patient trace\n(n=10)", ha="center",
            va="bottom", fontsize=8.5, color="0.35")

    # group brackets
    def bracket(x0, x1, label, color):
        yb = nb + 0.92
        ax.plot([x0, x1], [yb, yb], color=color, lw=1.3, clip_on=False)
        ax.plot([x0, x0], [yb - 0.08, yb], color=color, lw=1.3, clip_on=False)
        ax.plot([x1, x1], [yb - 0.08, yb], color=color, lw=1.3, clip_on=False)
        ax.text((x0 + x1) / 2, yb + 0.07, label, ha="center", va="bottom",
                fontsize=7.8, color=color, fontweight="bold")
    bracket(col_x["rho"], col_x["grass"] + CW, "TRACE EXISTS", "#2c3e50")
    bracket(col_x["enc"], col_x["inf"] + CW, "COPHENETIC CONTENT", "#2c3e50")
    bracket(col_x["home"], col_x["home"] + CW, "ANATOMY", "#2c3e50")

    ax.set_xlim(-0.75, x_end + 0.15)
    ax.set_ylim(-0.15, nb + 1.35)
    ax.axis("off")
    ax.text(-0.72, nb + 1.2, "a", fontsize=13, fontweight="bold", va="top")

    # legend
    handles = [
        Patch(facecolor=TIER_COLOR["verified"], edgecolor="white", label="verified / locked"),
        Patch(facecolor=TIER_COLOR["subset"], edgecolor="white", label="subset · directional · LOO-fragile"),
        Patch(facecolor=TIER_COLOR["null"], edgecolor="white", label="tested, null"),
        Patch(facecolor=LOCWO_FACE, edgecolor=LOCWO_EDGE, hatch=LOCWO_HATCH, label="localizes without a trace"),
        Patch(facecolor=FLAG_COLOR["anti"], edgecolor="white", label="per-patient anti (strip)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=8.6)
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_evidence_convergence_matrix.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\n[preprint_41] -> {out}")


if __name__ == "__main__":
    main()
