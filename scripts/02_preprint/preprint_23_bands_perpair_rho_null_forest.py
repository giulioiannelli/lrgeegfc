#!/usr/bin/env python3
"""Per-pair cophenetic trace, shown the way the GATE actually tests it.

The joint-density renders (preprint_17/18) fail on the cophenetic substrate
because ``Δ^coph`` is ~90% tied: a 2D density piles the non-mover atom on the
rank-centre (= the diagonal) and fakes a trace in every band. Spearman ρ — the
§5.3 statistic of record — is tie-corrected, so it ignores that centre and is
carried by the pairs that actually move. This figure renders ρ directly, so the
ties never enter.

Per band (2×3, δ θ α / β γ_l γ_h), one row per patient:

  • light grey span  = the patient's matched-strength surrogate null, p5–p95
  • tick at centre   = surrogate median ρ (null centre)
  • coloured stem    = obs ρ − null-median (the paired difference the Wilcoxon
                       gate ranks); GREEN = trace direction (obs > null median),
                       RED = anti
  • dot at obs ρ     = filled if it clears its OWN null (obs > p95, or < p5 for
                       a strong anti); open if it sits inside the null band

The band is read on TWO axes, never collapsed to a single pass/fail (see
`.agents/reports/2026-06-25_cophenetic-gate-presence-vs-consistency.md`):

  • PRESENCE  — do more patients clear their OWN null than the 5 % chance
                baseline? Binomial(n, 0.05) on the filled-green count. Detects a
                real trace in *some* patients. Only θ fails this.
  • CONSISTENCY — does the cohort *centre* sit on the trace side? The locked
                one-sided paired Wilcoxon (`gate_p`). β, α pass; γ_l, δ, γ_h
                fail because they are directionally split, not empty.

The 2×2 gives three band verdicts: a cohort-wide trace = presence ✓ AND
consistency ✓ (β, α); present but split = presence ✓ but consistency ✗ — large
effects, mixed sign, subgroup-carried (γ_l: Pat_05/02/06 strongly pro,
Pat_14/15/07/13 anti; also δ, γ_h); absent = presence ✗ (θ only). A failed
consistency test is reported as SPLIT, never as "no trace" — that is the whole
point: γ_l clears its null in 5/10 patients (binom p=6e-5, same as α) yet the
one-sided Wilcoxon, which γ_l fails at p=0.116, would otherwise label it empty.

Row ordering is selectable (``--sort``):
  • ``rho``     — patients sorted by observed ρ within each panel (default; the
                  classic forest, best at top). Rows are anonymous.
  • ``patient`` — patients in fixed canonical order (Pat_02 … Pat_15, Pat_02 at
                  top) identical across all six panels, so the SAME patient sits
                  at the SAME row in every band — read down a column of panels to
                  follow one patient across bands. Patient labels on the left.
The sort type is appended to the output filename.

The band symbol is coloured by the two-axis verdict (``--verdict two-axis``,
default): green = cohort-wide trace, amber = present but split, grey = absent;
a glyph beside the letter depicts the cell (full disc / half disc / open ring).
``--verdict binary`` restores the legacy one-axis colouring (green if the locked
Wilcoxon gate passes, grey otherwise). Presence (binomial) and consistency
(Wilcoxon) p-values plus the per-patient clear-counts print to stdout (kept out
of the figure per the no-in-axes-text rule); they belong in the caption.

Input (locked, instant — no recompute): the §5.3 split-baseline matched-strength
per-patient table
``data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv``
(cophenetic substrate). Output:
``data/preprint/figures/all_bands/patients_forest_plot/
fig_bands_perpair_rho_null_forest_coph_sort-{rho,patient}.pdf``
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from scipy.stats import binomtest, wilcoxon

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

# Single typography knob — every text size below is set RELATIVE to this (via
# matplotlib's named sizes), so bumping BASE_FONTSIZE rescales the whole
# figure's text at once. No per-element fontsize is hardcoded.
BASE_FONTSIZE = 16
plt.rcParams.update({
    "font.size": BASE_FONTSIZE,
    "axes.labelsize": "large",      # axis labels
    "xtick.labelsize": "medium",    # ρ-axis numbers
    "ytick.labelsize": "medium",    # patient row labels (sort=patient)
    "legend.fontsize": "large",
})

PERPAT_CSV = (ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
             / "per_patient_per_band.csv")
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# canonical cohort order (n=10), used when --sort patient (Pat_02 at top)
CANON = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
         "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

C_PRO = "#1a7c3e"      # obs ρ above null median — trace direction
C_ANTI = "#c0392b"     # obs ρ below null median — anti
C_NULL = "#b9b9b9"     # matched-strength null span (p5–p95)
C_NULLMID = "#7a7a7a"  # null median tick

# Two-axis verdict palette (presence × consistency). See band_verdict().
C_COHORT = C_PRO       # present AND cohort-consistent — a cohort-wide trace
C_HETERO = "#cc8400"   # present but split — subgroup trace, no cohort direction
C_ABSENT = "0.55"      # presence test fails — the genuine null
ALPHA = 0.05


def gate_p(d: pd.DataFrame) -> float:
    """CONSISTENCY axis: paired one-sided Wilcoxon of (obs ρ − own surrogate
    median). Asks 'does the cohort *centre* sit on the trace side?'. This is the
    locked §5.3 C3 gate — but on its own it cannot tell a split cohort (large
    effects, mixed sign) from a true null, so it is paired with presence_p()."""
    diff = (d.obs_rho - d.surr_p50).values
    return float(wilcoxon(diff, alternative="greater").pvalue)


def presence_p(d: pd.DataFrame) -> tuple[int, int, float]:
    """PRESENCE axis: how many patients clear their OWN matched-strength null,
    vs the 5 % per-patient chance baseline. Under the global null each patient
    clears its one-sided p95 with prob 0.05 independently, so #clearing ~
    Binomial(n, 0.05). Returns (pos_clear, neg_clear, binomial p for pos_clear).
    A calibrated existence test — NOT an arbitrary patient-count threshold."""
    pos_clear = int((d.obs_rho > d.surr_p95).sum())
    neg_clear = int((d.obs_rho < d.surr_p5).sum())
    p = binomtest(pos_clear, len(d), ALPHA, alternative="greater").pvalue
    return pos_clear, neg_clear, float(p)


def band_verdict(d: pd.DataFrame) -> str:
    """Read the cophenetic trace on two axes instead of collapsing to pass/fail.

      cohort  — presence ✓ AND consistency ✓  (a cohort-wide trace: β, α)
      hetero  — presence ✓ but consistency ✗  (present but split / subgroup-
                carried, no shared cohort direction: γ_l, δ, γ_h)
      absent  — presence ✗                    (the genuine null: θ)

    'no trace' is reserved for `absent` only; a failed consistency test is
    reported as `hetero`, never as absence."""
    _, _, p_pres = presence_p(d)
    p_cons = gate_p(d)
    if p_pres >= ALPHA:
        return "absent"
    return "cohort" if p_cons < ALPHA else "hetero"


_VERDICT_STYLE = {  # colour, marker fillstyle (full / split-half / open)
    "cohort": (C_COHORT, "full"),
    "hetero": (C_HETERO, "right"),
    "absent": (C_ABSENT, "none"),
}


def plot_band(ax, d: pd.DataFrame, xlim, p_gate: float, sort: str,
              show_ylabels: bool, verdict_mode: str = "two-axis") -> None:
    if sort == "patient":
        order = [p for p in CANON if p in set(d.patient)]      # top → bottom
        d = d.set_index("patient").loc[order].reset_index()
    else:                                                      # sort == "rho"
        order = None
        d = d.sort_values("obs_rho").reset_index(drop=True)

    for i, row in d.iterrows():
        # matched-strength null band (p5–p95) + median tick
        ax.plot([row.surr_p5, row.surr_p95], [i, i], color=C_NULL, lw=5.5,
                solid_capstyle="round", alpha=0.6, zorder=1)
        ax.plot([row.surr_p50], [i], marker="|", ms=8, color=C_NULLMID,
                mew=1.3, zorder=2)
        # paired-difference stem from null centre to observed ρ
        pro = row.obs_rho > row.surr_p50
        col = C_PRO if pro else C_ANTI
        ax.plot([row.surr_p50, row.obs_rho], [i, i], color=col, lw=1.8, zorder=3)
        beyond = (row.obs_rho > row.surr_p95) or (row.obs_rho < row.surr_p5)
        ax.plot([row.obs_rho], [i], marker="o", ms=7.5, zorder=4,
                mfc=(col if beyond else "white"), mec=col, mew=1.5)

    ax.axvline(0.0, color="0.45", ls="--", lw=0.8, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("0.6")

    if sort == "patient":
        ax.invert_yaxis()                       # Pat_02 (i=0) at the top
        ax.spines["left"].set_color("0.6")
        if show_ylabels:
            ax.set_yticks(range(len(d)))
            ax.set_yticklabels([p.replace("Pat_", "Pat ") for p in order],
                               size="medium")
        else:
            ax.set_yticks([])
    else:
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)

    # band symbol + verdict glyph
    band = d.band.iloc[0]
    tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    if verdict_mode == "binary":
        # legacy one-axis colouring: green = locked gate passes, grey = fails.
        ax.text(0.97, 0.07, tex, transform=ax.transAxes, ha="right",
                va="bottom", size="xx-large", fontweight="bold",
                color=(C_PRO if p_gate < ALPHA else "0.45"))
        return
    # two-axis: letter coloured by the presence×consistency verdict, with a
    # glyph that depicts the cell — full disc = cohort-consistent, half disc =
    # split, open ring = absent.
    verdict = band_verdict(d)
    col, fill = _VERDICT_STYLE[verdict]
    ax.text(0.97, 0.07, tex, transform=ax.transAxes, ha="right", va="bottom",
            size="xx-large", fontweight="bold", color=col)
    ax.plot([0.80], [0.135], transform=ax.transAxes, clip_on=False, zorder=6,
            marker=MarkerStyle("o", fillstyle=fill), ms=15,
            mfc=col, mfcalt="white", mec=col, mew=1.7)


def main() -> Path:
    ap = argparse.ArgumentParser(
        description="Per-pair cophenetic trace forest (band-grouped).")
    ap.add_argument("--sort", choices=["rho", "patient"], default="rho",
                    help="row order within each band panel: 'rho' (default, "
                         "forest sorted by observed ρ) or 'patient' (fixed "
                         "canonical order, same patient at the same row across "
                         "all bands).")
    ap.add_argument("--verdict", choices=["two-axis", "binary"],
                    default="two-axis",
                    help="band-letter colouring. 'two-axis' (default) colours by "
                         "presence×consistency (green=cohort-wide, amber=present-"
                         "but-split, grey=absent) and draws a verdict glyph; "
                         "'binary' is the legacy green/grey on the locked "
                         "one-sided Wilcoxon gate alone.")
    args = ap.parse_args()

    df = pd.read_csv(PERPAT_CSV)
    lo = float(df.obs_rho.min()); hi = float(df.obs_rho.max())
    pad = 0.06 * (hi - lo)
    xlim = (lo - pad, hi + pad)

    left = 0.085 if args.sort == "patient" else 0.04
    fig = plt.figure(figsize=(15.5, 9.5))
    outer = fig.add_gridspec(2, 3, wspace=0.10, hspace=0.20,
                             left=left, right=0.985, top=0.93, bottom=0.20)
    print(f"Per-pair cophenetic trace — TWO AXES (presence × consistency); "
          f"sort={args.sort}, verdict={args.verdict}:")
    print(f"  {'band':11s} {'pos_clr':>7s} {'neg_clr':>7s} "
          f"{'binom_p(pres)':>13s} {'wilcox_p(cons)':>14s} {'verdict':>8s}")
    for idx, band in enumerate(BAND_ORDER):
        d = df[df.band == band]
        rr, cc = idx // 3, idx % 3
        ax = fig.add_subplot(outer[rr, cc])
        p = gate_p(d)
        plot_band(ax, d, xlim, p, args.sort, show_ylabels=(cc == 0),
                  verdict_mode=args.verdict)
        if rr == 1:
            ax.set_xlabel(r"$\rho^{\mathrm{coph}}$", color="0.20")
        pos_clr, neg_clr, p_pres = presence_p(d)
        label = {"cohort": "cohort", "hetero": "split", "absent": "absent"}
        print(f"  {band:11s} {pos_clr:7d} {neg_clr:7d} {p_pres:13.2e} "
              f"{p:14.4f} {label[band_verdict(d)]:>8s}")

    handles = [
        Line2D([0], [0], color=C_NULL, lw=9, solid_capstyle="round",
               alpha=0.6, label="matched-strength null (p5–p95)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=2.6, mfc=C_PRO, mec=C_PRO,
               ms=13, label="patient clears null (filled)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=2.6, mfc="white", mec=C_PRO,
               ms=13, label="patient within null (open)"),
        Line2D([0], [0], marker="o", color=C_ANTI, lw=2.6, mfc=C_ANTI,
               mec=C_ANTI, ms=13, label="anti"),
    ]
    if args.verdict == "two-axis":
        handles += [
            Line2D([0], [0], lw=0, marker=MarkerStyle("o", fillstyle="full"),
                   mfc=C_COHORT, mec=C_COHORT, ms=14,
                   label="band: cohort-wide trace"),
            Line2D([0], [0], lw=0, marker=MarkerStyle("o", fillstyle="right"),
                   mfc=C_HETERO, mfcalt="white", mec=C_HETERO, mew=1.7, ms=14,
                   label="band: present but split"),
            Line2D([0], [0], lw=0, marker=MarkerStyle("o", fillstyle="none"),
                   mec=C_ABSENT, mew=1.7, ms=14, label="band: absent"),
        ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01), ncol=4, frameon=False,
               handlelength=2.0, handletextpad=0.6, columnspacing=1.8)

    out_dir = (ROOT / "data" / "preprint" / "figures" / "all_bands"
               / "patients_forest_plot")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_bands_perpair_rho_null_forest_coph_sort-{args.sort}.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
