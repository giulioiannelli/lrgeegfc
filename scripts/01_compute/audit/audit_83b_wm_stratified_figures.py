#!/usr/bin/env python3
"""White-matter-stratified trace figures (mirror of preprint_28/29/30).

Reads ONLY the audit_83 / audit_84 cohort + per-patient caches and renders:

  1. fig_wm_verdict_matrix_coph.pdf  — band × 6 configs, cophenetic ρ_split,
     colored by sensitivity_flag (config vs full baseline).
  2. fig_wm_verdict_matrix_raw.pdf   — same, raw |ImCoh| ρ_split^raw.
  3. fig_wm_exclude_forest.pdf       — per-patient full vs exclude_wm obs_stat
     (cophenetic) for the headline bands, with own-surrogate p50 reference.
  4. fig_wm_grassmann_kspan.pdf      — band × {exclude_wm, wm_only}, color =
     longest separated k-run (Grassmann T_G), with full's run as text-free ref.

No in-axes numeric text (project rule): numbers go to stdout. Figure-level
legends, PDF-only, transparent, use_lrg_style, no suptitle.

Inputs : data/audit/wm_stratified/{cophenetic_raw_cohort,cophenetic_raw_per_patient,
         grassmann_cohort}.csv
Output : data/audit/wm_stratified/figures/*.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

SRC = ROOT / "data" / "audit" / "wm_stratified"
OUT = SRC / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HEADLINE_BANDS = ["alpha", "beta", "low_gamma"]
CONFIG_ORDER = ["full", "gray_gray", "exclude_wm", "cross", "wm_wm", "wm_only"]
CONFIG_LABEL = {"full": "full", "gray_gray": "gray–gray",
                "exclude_wm": "excl-WM", "cross": "cross",
                "wm_wm": "WM–WM", "wm_only": "WM-only"}
FLAG_COLOR = {"baseline": "#2c5f8a", "persist": "#1f7a1f", "weaken": "#c0392b",
              "emerge": "#e67e22", "absent": "#95a5a6", "undefined": "#ffffff",
              "unknown": "#d9d9d9"}
FLAG_ORDER = ["baseline", "persist", "weaken", "emerge", "absent", "undefined"]
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]


# ---------------------------------------------------------------------------
# 1-2. Verdict matrix (per substrate)
# ---------------------------------------------------------------------------
def verdict_matrix(cohort: pd.DataFrame, substrate: str, out_name: str) -> None:
    cs = cohort[cohort.substrate == substrate]
    flag = {(r.band, r.config): str(r.sensitivity_flag) for r in cs.itertuples()}
    obsm = {(r.band, r.config): r.obs_median for r in cs.itertuples()}
    pval = {(r.band, r.config): r.paired_wilcoxon_p for r in cs.itertuples()}

    nb, nc = len(BANDS), len(CONFIG_ORDER)
    fig, ax = plt.subplots(figsize=(1.05 * nc + 1.6, 0.85 * nb + 1.2))

    print(f"\n=== {substrate} verdict matrix (obs_median / paired Wilcoxon p) ===")
    print(f"{'band':<11}" + "".join(f"{CONFIG_LABEL[c]:>13}" for c in CONFIG_ORDER))
    for bi, band in enumerate(BANDS):
        line = f"{band:<11}"
        for ci, cfg in enumerate(CONFIG_ORDER):
            f = flag.get((band, cfg), "undefined")
            color = FLAG_COLOR.get(f, "#d9d9d9")
            y = nb - 1 - bi
            if f == "undefined":
                ax.add_patch(plt.Rectangle((ci, y), 1, 1, facecolor="white",
                                           edgecolor="#bbbbbb", hatch="////",
                                           linewidth=0.5))
            else:
                ax.add_patch(plt.Rectangle((ci, y), 1, 1, facecolor=color,
                                           edgecolor="white", linewidth=1.2))
            om = obsm.get((band, cfg), np.nan)
            pv = pval.get((band, cfg), np.nan)
            cell = "n/d" if f == "undefined" else f"{om:+.2f}/p{pv:.3f}"
            line += f"{cell:>13}"
        print(line)

    ax.set_xlim(0, nc); ax.set_ylim(0, nb)
    ax.set_xticks(np.arange(nc) + 0.5)
    ax.set_xticklabels([CONFIG_LABEL[c] for c in CONFIG_ORDER], rotation=30,
                       ha="right")
    ax.set_yticks(np.arange(nb) + 0.5)
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in reversed(BANDS)])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    handles = [Patch(facecolor=FLAG_COLOR[f], edgecolor="white", label=f)
               if f != "undefined" else
               Patch(facecolor="white", edgecolor="#bbbbbb", hatch="////",
                     label="undefined")
               for f in FLAG_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=len(FLAG_ORDER),
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = OUT / out_name
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[audit_83b] -> {out}")


# ---------------------------------------------------------------------------
# 3. Forest: per-patient full vs exclude_wm obs (cophenetic), headline bands
# ---------------------------------------------------------------------------
def exclude_forest(per_pat: pd.DataFrame, out_name: str) -> None:
    d = per_pat[(per_pat.substrate == "cophenetic") & per_pat.defined]
    bands = [b for b in HEADLINE_BANDS if b in set(d.band)]
    if not bands:
        print("[audit_83b] forest: no headline bands present; skip")
        return
    fig, axes = plt.subplots(1, len(bands), figsize=(3.0 * len(bands) + 0.5, 4.2),
                             sharey=True)
    if len(bands) == 1:
        axes = [axes]
    ypos = {p: i for i, p in enumerate(reversed(COHORT))}
    for ax, band in zip(axes, bands):
        db = d[d.band == band]
        full = db[db.config == "full"].set_index("patient")
        excl = db[db.config == "exclude_wm"].set_index("patient")
        for p in COHORT:
            if p not in full.index or p not in excl.index:
                continue
            y = ypos[p]
            xf, xe = full.loc[p, "obs_stat"], excl.loc[p, "obs_stat"]
            ax.plot([xf, xe], [y, y], color="#cccccc", lw=1.0, zorder=1)
            ax.scatter([xf], [y], s=26, color="#2c5f8a", zorder=3,
                       label="full" if p == COHORT[0] else None)
            ax.scatter([xe], [y], s=26, color="#1f7a1f", zorder=3,
                       label="exclude_wm" if p == COHORT[0] else None)
            # own-surrogate p50 reference tick (gray-only)
            sp50 = excl.loc[p, "surr_p50"]
            ax.scatter([sp50], [y], s=10, marker="|", color="#c0392b", zorder=2)
        ax.axvline(0.0, color="black", lw=0.6, ls=":")
        ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band))
        ax.set_xlabel(r"$\rho_{\mathrm{split}}$ (positive = trace)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(list(ypos.values()))
    axes[0].set_yticklabels(list(reversed(COHORT)))
    handles = [plt.Line2D([], [], marker="o", ls="", color="#2c5f8a", label="full"),
               plt.Line2D([], [], marker="o", ls="", color="#1f7a1f", label="exclude_wm"),
               plt.Line2D([], [], marker="|", ls="", color="#c0392b",
                          label="exclude_wm surrogate median")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = OUT / out_name
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[audit_83b] -> {out}")


# ---------------------------------------------------------------------------
# 4. Grassmann k-span strip: band × {full, exclude_wm, wm_only}
# ---------------------------------------------------------------------------
def grassmann_kspan(gcohort: pd.DataFrame, out_name: str) -> None:
    cfgs = ["full", "exclude_wm", "wm_only"]
    label = {"full": "full", "exclude_wm": "excl-WM", "wm_only": "WM-only"}
    # Judge Grassmann on the cohort paired-Wilcoxon p (the audit_66/78 basis).
    # The strict cophenetic-calibrated "separated" flag adds a med_surr≈0
    # condition that the positive-offset chordal statistic fails even when it
    # clears matched-strength, so it undercounts the Grassmann trace.
    run = {}
    nk_sep = {}
    for (band, cfg), sub in gcohort.groupby(["band", "config"]):
        sep = (sub.sort_values("k").paired_wilcoxon_p < 0.05).values
        longest = cur = 0
        for v in sep:
            cur = cur + 1 if v else 0
            longest = max(longest, cur)
        run[(band, cfg)] = longest
        nk_sep[(band, cfg)] = int(sep.sum())

    nb, nc = len(BANDS), len(cfgs)
    vmax = max([1] + list(run.values()))
    fig, ax = plt.subplots(figsize=(1.2 * nc + 2.0, 0.85 * nb + 1.2))
    cmap = plt.get_cmap("viridis")
    print("\n=== Grassmann longest p<0.05 k-run (n_k with Wilcoxon p<0.05) ===")
    print(f"{'band':<11}" + "".join(f"{label[c]:>12}" for c in cfgs))
    for bi, band in enumerate(BANDS):
        line = f"{band:<11}"
        for ci, cfg in enumerate(cfgs):
            r = run.get((band, cfg), 0)
            y = nb - 1 - bi
            frac = r / vmax if vmax else 0.0
            ax.add_patch(plt.Rectangle((ci, y), 1, 1,
                                       facecolor=cmap(frac) if r > 0 else "#f0f0f0",
                                       edgecolor="white", linewidth=1.2))
            line += f"{str(r)+'('+str(nk_sep.get((band,cfg),0))+')':>12}"
        print(line)
    ax.set_xlim(0, nc); ax.set_ylim(0, nb)
    ax.set_xticks(np.arange(nc) + 0.5)
    ax.set_xticklabels([label[c] for c in cfgs], rotation=20, ha="right")
    ax.set_yticks(np.arange(nb) + 0.5)
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in reversed(BANDS)])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=0, vmax=vmax))
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("longest Wilcoxon p<0.05 k-run")
    fig.tight_layout()
    out = OUT / out_name
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"[audit_83b] -> {out}")


def main() -> None:
    coh = SRC / "cophenetic_raw_cohort.csv"
    pp = SRC / "cophenetic_raw_per_patient.csv"
    if not coh.exists() or not pp.exists():
        raise SystemExit(f"[audit_83b] missing {coh} / {pp}; run audit_83 first")
    cohort = pd.read_csv(coh)
    per_pat = pd.read_csv(pp)
    verdict_matrix(cohort, "cophenetic", "fig_wm_verdict_matrix_coph.pdf")
    verdict_matrix(cohort, "raw", "fig_wm_verdict_matrix_raw.pdf")
    exclude_forest(per_pat, "fig_wm_exclude_forest.pdf")

    gco = SRC / "grassmann_cohort.csv"
    if gco.exists():
        grassmann_kspan(pd.read_csv(gco), "fig_wm_grassmann_kspan.pdf")
    else:
        print(f"[audit_83b] {gco} absent; skip Grassmann strip (run audit_84)")


if __name__ == "__main__":
    main()
