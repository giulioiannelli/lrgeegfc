#!/usr/bin/env python3
r"""fig:arc_g — the cingulate carries different content in different rhythms (§2, R2.7).

CORE MESSAGE (the synthesis / multiplexing twist): the two cingulate results describe one
region holding different task content at different frequencies. In low-gamma the cingulate
carries the MEMORY (encoding) trace and not the inference component; in beta it carries the
INFERENCE component and not the memory trace. The same cortex, read at two frequencies,
holds two different parts of the task -- the premises the patient was shown surface in one
rhythm, the relations the patient reasoned out in another.

The two halves stand on different footing (drawn asymmetrically):
  - memory @ low-gamma : cohort-clearing absolute within-region trace (R2.5)      -> FIRM
  - inference @ beta    : full-length enrichment lead, but length-assisted (R2.6)  -> PROVISIONAL
So the picture is firm for memory and provisional for inference.

No new compute: this reads R2.5 (within-system absolute test, low-gamma column) and
R2.3/R2.6 (system-enrichment demeaned test, beta column). Both report BH q.

Reads : data/audit/inference_localization_rhosym/within_system_trace_rhosym_include.csv
        data/audit/inference_localization_rhosym/{encoding,inference_pe}_localization_rhosym_include.csv
Writes: data/reports/results_section2/fig_arc_g_cingulate_multiplexing.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BASE = ROOT / "data/audit/inference_localization_rhosym"
OUT = ROOT / "data/reports/results_section2/fig_arc_g_cingulate_multiplexing.pdf"

DROP = {"other", "non_anatomical"}
CING = "cingulate"
C_MEM, C_INF, C_NULL = "#2a9d5c", "#d1491c", "#e0e0e0"
Q_SIG = 0.05


def within_q(target: str) -> dict:
    d = pd.read_csv(BASE / "within_system_trace_rhosym_include.csv")
    r = d[(d.band == "low_gamma") & (d.target == target) & (d.system == CING)].iloc[0]
    return dict(q=float(r.bh_q), rho=float(r.median_obs_rho),
                n_pos=int(r.n_pos), K=int(r.K_implanted))


def enrich_q(tag: str) -> dict:
    d = pd.read_csv(BASE / f"{tag}_localization_rhosym_include.csv")
    d = d[(d.band == "beta") & (d.granularity == "system") &
          (~d.unit.isin(DROP))].copy()
    n = int(d.n_surr.iloc[0])
    up = np.clip(d.matched_strength_p.astype(float), 1.0 / (n + 1), 1.0)
    d["q_up"] = bh_fdr(up.tolist())
    row = d[d.unit == CING].iloc[0]
    return dict(q=float(row.q_up), ms_p=float(row.matched_strength_p))


def draw_cell(ax, cx, cy, cell):
    clears = cell["clears"]
    prov = cell.get("prov", False)
    if clears and not prov:
        face, ec, hatch, txt_c = C_MEM, C_MEM, None, "white"
    elif clears and prov:
        face, ec, hatch, txt_c = C_INF, C_INF, "////", "white"
    else:
        face, ec, hatch, txt_c = C_NULL, "#b6b6b6", None, "0.35"
    box = FancyBboxPatch((cx - 0.44, cy - 0.44), 0.88, 0.88,
                         boxstyle="round,pad=0.02,rounding_size=0.06",
                         linewidth=2.0 if clears else 1.2,
                         linestyle="--" if prov else "-",
                         facecolor=face, edgecolor=ec, alpha=0.92 if clears else 1.0,
                         hatch=hatch, zorder=2)
    ax.add_patch(box)
    head = r"$\checkmark$ clears" if clears else "n.s."
    ax.text(cx, cy + 0.26, head, ha="center", va="center", fontsize=15,
            color=txt_c, fontweight="bold", zorder=3)
    ax.text(cx, cy + 0.02, cell["qtxt"], ha="center", va="center", fontsize=13,
            color=txt_c, zorder=3)
    if cell.get("sub"):
        ax.text(cx, cy - 0.20, cell["sub"], ha="center", va="center", fontsize=10,
                color=txt_c, style="italic", zorder=3)
    if cell.get("conf"):
        ax.text(cx, cy - 0.35, cell["conf"], ha="center", va="center", fontsize=10,
                color=txt_c, fontweight="bold", zorder=3)


def main():
    mem_lg, inf_lg = within_q("encoding"), within_q("inference_pe")
    mem_b, inf_b = enrich_q("encoding"), enrich_q("inference_pe")

    # rows: memory (top), inference (bottom); cols: low-gamma (left), beta (right)
    cells = {
        ("mem", "lg"): dict(clears=mem_lg["q"] < Q_SIG,
                            qtxt=rf"$q = {mem_lg['q']:.3f}$",
                            sub=rf"{mem_lg['n_pos']}/{mem_lg['K']} patients",
                            conf="firm"),
        ("mem", "b"): dict(clears=mem_b["q"] < Q_SIG,
                           qtxt=rf"$q = {mem_b['q']:.2f}$"),
        ("inf", "lg"): dict(clears=inf_lg["q"] < Q_SIG,
                            qtxt=(rf"$q = {inf_lg['q']:.2f}$" if inf_lg["q"] < 0.995
                                  else r"$q > 0.9$")),
        ("inf", "b"): dict(clears=inf_b["q"] < Q_SIG, prov=True,
                           qtxt=rf"$q = {inf_b['q']:.3f}$",
                           sub="length-assisted", conf="provisional"),
    }
    colx = {"lg": 1.0, "b": 2.0}
    rowy = {"mem": 2.0, "inf": 1.0}

    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    for (rk, ck), cell in cells.items():
        draw_cell(ax, colx[ck], rowy[rk], cell)

    # column headers (frequency)
    for ck, lbl in (("lg", BRAIN_BAND_TEX_DICT["low_gamma"]),
                    ("b", BRAIN_BAND_TEX_DICT["beta"])):
        ax.text(colx[ck], 2.72, lbl, ha="center", va="center", fontsize=22)
    ax.text(1.5, 3.08, "read at frequency", ha="center", va="center",
            fontsize=12, color="0.4", style="italic")
    # row headers (task content)
    ax.text(0.30, rowy["mem"], "memory\n(encoding)", ha="center", va="center",
            fontsize=14, color=C_MEM, fontweight="bold", rotation=90)
    ax.text(0.30, rowy["inf"], "inference", ha="center", va="center",
            fontsize=14, color=C_INF, fontweight="bold", rotation=90)
    ax.text(-0.02, 1.5, "task content", ha="center", va="center", fontsize=12,
            color="0.4", style="italic", rotation=90)

    # region identity + reading
    ax.text(1.5, 0.28, "one region — the cingulate — its content switched by rhythm",
            ha="center", va="center", fontsize=11.5, color="0.25")

    ax.set_xlim(-0.15, 2.6)
    ax.set_ylim(0.05, 3.25)
    ax.set_aspect("equal")
    ax.axis("off")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print("fig:arc_g — cingulate multiplexing 2x2 (rho_sym)")
    print(f"  memory   low-gamma: q={mem_lg['q']:.3f} rho={mem_lg['rho']:+.3f} "
          f"{mem_lg['n_pos']}/{mem_lg['K']}  -> {'CLEARS (firm)' if mem_lg['q']<Q_SIG else 'null'}")
    print(f"  memory   beta     : q={mem_b['q']:.3f} (MS_p={mem_b['ms_p']:.3f})"
          f"  -> {'clears' if mem_b['q']<Q_SIG else 'null'}")
    print(f"  inference low-gamma: q={inf_lg['q']:.3f} rho={inf_lg['rho']:+.3f} "
          f"-> {'clears' if inf_lg['q']<Q_SIG else 'null'}")
    print(f"  inference beta     : q={inf_b['q']:.3f} (MS_p={inf_b['ms_p']:.3f})"
          f"  -> {'CLEARS (provisional, length-assisted)' if inf_b['q']<Q_SIG else 'null'}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
