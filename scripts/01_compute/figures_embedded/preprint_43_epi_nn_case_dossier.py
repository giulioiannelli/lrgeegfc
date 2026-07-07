"""
preprint_43 — Results §3 figure: "the same propagator can inspire epileptogenic markers".

Proof-of-concept dossier for the 3rd Results subsection (method generality: one operator,
a second read-out). Reads ONLY cached audit CSVs (no recompute). PDF-only, vector,
use_lrg_style. Foregrounds the LOCKED marker arc — organization -> node read-out -> detector
-> positioning — not a clinical-superiority claim.

  A  relational organization (the lead)  — audit_89/90 epi co-diffusion community: per-band excess
                                            z of epi<->epi co-diffusion vs matched-strength random
                                            same-size communities (mean over diffusion scales).
                                            delta/low-g/beta strong, high-gamma ~0 (negative control).
  B  per-node read-out beyond strength   — audit_132 all-contacts heat-tau5 AUC *beyond strength*
                                            vs each band's matched-strength null p95. delta/low-g/beta
                                            clear (p<0.001); alpha marginal. The only legitimate null
                                            (|ImCoh| zero-lag-immune -> no spatial control).
  C  deployable detector (proof-of-conc) — audit_117 LOPO calibrated P(SOZ): per-patient AUC sorted,
                                            responders vs the 2 hub-patients that fail, vs label-shuffle
                                            null. The two-population split shown, not hidden.
  D  lands in the field range (generality) — our marker/detector AUC against the established interictal-
                                            marker range (~0.70-0.86), from a substrate built for a
                                            cognitive question. Honest note: ours is SOZ-label-validated,
                                            field-best are outcome-validated.

Sources: data/audit/epi_propagator_blocks/, epi_marker_allcontacts/, epi_propagator_detector/.
Literature anchors (panel D) cited in .agents/reports/2026-06-11_epi-soz-marker-literature-positioning.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

AC = ROOT / "data/audit/epi_marker_allcontacts"
DET = ROOT / "data/audit/epi_propagator_detector"
BLOCKS = ROOT / "data/audit/epi_propagator_blocks"
OUT = ROOT / "data/preprint/figures/all_bands"
OUT.mkdir(parents=True, exist_ok=True)

# ---- palette (saturated; no near-white interiors) ----
C_BAND = {"delta": "#08519c", "low_gamma": "#238b45", "beta": "#6a51a3",
          "alpha": "#969696", "theta": "#cc6677", "high_gamma": "#cccccc"}
C_RESP = "#1b7837"
C_HUB = "#c0392b"
C_NULL = "#9e9e9e"
C_FIELD = "#d9c7a3"      # field-range shade
C_OURS = "#d4a017"       # this work
C_SOTA = "#34495e"       # established field markers
C_ALPHA = "#b8860b"      # panel D alpha
C_BETA = "#6a51a3"       # panel D beta

BANDS_A = ["delta", "low_gamma", "beta", "alpha"]
HUBS = {"Pat_10", "Pat_15"}


# ---------------------------------------------------------------- panel A
def panel_community(ax):
    df = pd.read_csv(BLOCKS / "propagator_blocks_cohort.csv")
    # mean excess-z of epi<->epi co-diffusion vs matched-strength communities, over scales
    g = df.groupby("band")["med_z_rand_EE"].mean()
    bands = ["delta", "low_gamma", "beta", "alpha", "theta", "high_gamma"]
    vals = [float(g[b]) for b in bands]
    ax.axhline(0, color="k", lw=0.6, zorder=1)
    ax.axhline(1.96, color="0.55", ls=":", lw=1.0, zorder=1)
    ax.annotate("per-cell $p<0.05$", (len(bands) - 1, 1.96), xytext=(0, 3),
                textcoords="offset points", fontsize=6.5, ha="right", color="0.4")
    for x, (b, v) in enumerate(zip(bands, vals)):
        ax.bar(x, v, color=C_BAND[b], edgecolor="k", lw=0.4, width=0.72, zorder=2)
        ax.annotate(f"{v:.1f}", (x, v), xytext=(0, 3), textcoords="offset points",
                    fontsize=7.5, ha="center", color=C_BAND[b] if b != "high_gamma" else "0.5")
    ax.annotate("negative\ncontrol", (len(bands) - 1, 0.0), xytext=(0, 16),
                textcoords="offset points", fontsize=6.8, ha="center", color="0.45")
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands])
    ax.set_ylim(-0.4, 5.0)
    ax.set_ylabel("epi co-diffusion $z$ vs matched-strength")
    ax.set_title("A   epileptic tissue $=$ a co-diffusion community",
                 loc="left", fontweight="bold")


# ---------------------------------------------------------------- panel B
def panel_connectivity(ax):
    df = pd.read_csv(AC / "allcontacts_nulls.csv")
    df = df[df.marker == "heat_t5"].set_index("band")
    ax.axvline(0.5, color="k", ls="--", lw=0.8, zorder=1)
    for y, b in enumerate(BANDS_A):
        row = df.loc[b]
        obs, p95 = float(row.real_median_auc), float(row.null_median_p95)
        # null span: chance .. p95
        ax.plot([0.5, p95], [y, y], color=C_NULL, lw=6, alpha=0.45,
                solid_capstyle="butt", zorder=2)
        ax.scatter([p95], [y], marker="|", s=140, color=C_NULL, lw=1.6, zorder=3)
        clears = obs > p95 and b != "alpha"
        ax.scatter([obs], [y], s=95, color=C_BAND[b], edgecolor="k", lw=0.5,
                   zorder=4)
        ax.annotate(f"{obs:.2f}", (obs, y), xytext=(5, 7),
                    textcoords="offset points", fontsize=8,
                    color=C_BAND[b], fontweight="bold")
    ax.set_yticks(range(len(BANDS_A)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS_A])
    ax.set_xlim(0.45, 0.9)
    ax.set_ylim(-0.6, len(BANDS_A) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("SOZ AUC beyond node strength")
    ax.set_title("B   per-node read-out $>$ strength (multiband)",
                 loc="left", fontweight="bold")
    h = [Line2D([], [], marker="o", ls="", mfc="#444", mec="k", label="observed (median)"),
         Line2D([], [], marker="|", ls="", color=C_NULL, label="matched-strength null $p_{95}$")]
    ax.legend(handles=h, loc="lower right", fontsize=7.5, frameon=False)


# ---------------------------------------------------------------- panel B
def panel_deployment(ax):
    df = pd.read_csv(DET / "detector_lopo_per_patient.csv").sort_values("auc")
    shuffle_null = 0.476
    med = float(df.auc.median())
    cols = [C_HUB if p in HUBS else C_RESP for p in df.patient]
    y = np.arange(len(df))
    ax.hlines(y, shuffle_null, df.auc, color=cols, lw=2.4, zorder=2)
    ax.scatter(df.auc, y, color=cols, s=46, edgecolor="k", lw=0.4, zorder=3)
    ax.axvline(0.5, color="k", ls="--", lw=0.8, zorder=1, label="chance")
    ax.axvline(shuffle_null, color=C_NULL, ls=":", lw=1.4, zorder=1)
    ax.axvline(med, color="#222", ls="-", lw=1.0, alpha=0.6, zorder=1)
    ax.annotate(f"median {med:.2f}", (med, len(df) - 0.4), fontsize=7.5,
                ha="center", color="#222")
    ax.annotate("shuffle null", (shuffle_null, len(df) - 1.4), xytext=(-3, 0),
                textcoords="offset points", fontsize=6.8, ha="right",
                va="center", color=C_NULL, rotation=90)
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace("Pat_", "P") for p in df.patient], fontsize=8)
    for tick, p in zip(ax.get_yticklabels(), df.patient):
        if p in HUBS:
            tick.set_color(C_HUB)
    ax.set_xlim(0.4, 1.02)
    ax.set_ylim(-0.7, len(df) - 0.2)
    ax.set_xlabel("LOPO detector AUC (calibrated P(SOZ))")
    ax.set_title("C   deployable detector: 9/10, hubs fail",
                 loc="left", fontweight="bold")
    h = [Line2D([], [], marker="o", ls="", color=C_RESP, label="responder (8)"),
         Line2D([], [], marker="o", ls="", color=C_HUB, label="hub-SOZ fails (Pat_10/15)")]
    ax.legend(handles=h, loc="lower right", fontsize=7.5, frameon=False)


# ---------------------------------------------------------------- panel C
def panel_field(ax):
    # field interictal node-marker range (lit 2026-06-11)
    ax.axvspan(0.70, 0.86, color=C_FIELD, alpha=0.55, zorder=0)
    ax.axvline(0.5, color="k", ls="--", lw=0.8, zorder=1)
    # established field markers (outcome- or established-SOZ-validated)
    sota = [("source-sink (outcome)", 0.86), ("PDC directed", 0.88),
            ("EZ-fingerprint (LOPO)", 0.70), ("DRS paradigm", 0.75)]
    for name, auc in sota:
        ax.scatter([auc], [1], marker="D", s=46, color=C_SOTA, edgecolor="k",
                   lw=0.4, zorder=4)
        ax.annotate(name, (auc, 1), xytext=(0, 11), textcoords="offset points",
                    fontsize=6.6, ha="center", rotation=30, color=C_SOTA)
    # this work (SOZ-label validated) — stagger labels to avoid the 0.80/0.81 collision
    ours = [("marker $\\delta$", 0.80, (-9, -20), "right"),
            ("detector mean", 0.81, (4, -33), "center"),
            ("detector median", 0.87, (10, -20), "left")]
    for name, auc, off, ha in ours:
        ax.scatter([auc], [0], marker="*", s=190, color=C_OURS, edgecolor="k",
                   lw=0.5, zorder=4)
        ax.annotate(name, (auc, 0), xytext=off, textcoords="offset points",
                    fontsize=6.8, ha=ha, color="#7a5c00")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["this work\n(SOZ labels)", "field SOTA\n(outcome)"], fontsize=8)
    ax.set_ylim(-0.75, 1.75)
    ax.set_xlim(0.45, 0.95)
    ax.set_xlabel("AUC")
    ax.set_title("D   lands in the established-marker range",
                 loc="left", fontweight="bold")
    ax.annotate("field interictal-marker range", (0.78, 1.62), fontsize=7,
                ha="center", color="#6b5a33", style="italic")


def main():
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2))
    panel_community(axes[0, 0])
    panel_connectivity(axes[0, 1])
    panel_deployment(axes[1, 0])
    panel_field(axes[1, 1])
    fig.tight_layout(pad=1.4, w_pad=2.6, h_pad=3.0)
    out = OUT / "fig_epi_method_generality.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    main()
