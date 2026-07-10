#!/usr/bin/env python3
r"""fig_encoding_inference_dissociation — the beta encoding-vs-inference DOUBLE
DISSOCIATION under the rho_sym estimator (Results section 1, collaborator meeting).

rho_sym-native re-render; adapted from audit_110b_inference_localization_figure.py,
which is LEFT UNTOUCHED. This script only READS the already-computed rho_sym CSVs
(no surrogates are run here) and re-points the audit_110b layout at them.

SIGNIFICANCE CONVENTION (locked to the compute audits, NOT to fig_trace_b):
  Direction-separated BH-FDR over the system-granularity family, exactly as
    - audit_158 (encoding): q_up = BH(matched_strength_p) -> CARRIER;
                            q_lo = BH(matched_strength_p_lower) -> DEPLETED.
    - audit_160 (inference_pe): q_up = BH(matched_strength_p) -> LEAD only.
      Inference is ENRICHMENT-ONLY; the audit never tests inference depletion, so we
      never star a depleted inference cell (its lollipop position is still drawn, but
      it carries no significance claim).

CORE MESSAGE (band beta, systems level):
  - ENCODING (task_learn) over-accumulates the cophenetic trace in ORBITOFRONTAL
    cortex and DEPLETES prefrontal cortex. OFC is the only enriched carrier system
    (matched-strength BH q_up = 0.040 include / 0.010 exclude, canonical audit_158);
    PFC is the depleted one (q_lo ~ 0.010). VERIFIED.
  - INFERENCE-SPECIFIC (task_test beyond learn, encoding partialled out) over-
    accumulates in the CINGULATE (BH q_up = 0.030 include / 0.040 exclude at full
    length, canonical audit_160) — the sole enriched concentration. BUT this lead
    FAILS the length-matched control (cingulate BH q = 0.15 include / 0.30 exclude,
    audit_161): it is a DIRECTIONAL, LENGTH-ASSISTED lead, NOT an established result,
    and is annotated as such (hatched diamond + dashed ring + a dedicated legend
    entry). It is never shown at equal status to the encoding->OFC result.

A double dissociation: OFC lights up for encoding but not inference; cingulate lights
up for inference but not encoding; PFC is depleted in encoding.

FRAMING (locked): significance is marked UNIFORMLY (a star on every BH-q<0.05 cell in
a *tested* direction) — no threshold lines, no pre-shaded windows. Direction is left
(depleted) / right (enriched). The length-matched caveat is carried by a distinct
glyph on the single cingulate-inference cell, not by hiding its full-length star.

Reads:
  data/audit/inference_localization_rhosym/encoding_localization_rhosym_{include,exclude}.csv
  data/audit/inference_localization_rhosym/inference_pe_localization_rhosym_{include,exclude}.csv
  data/audit/inference_localization_rhosym/lenmatched_null_rhosym_R200_{include,exclude}.csv
Writes:
  data/preprint/figures/results_section1/fig_encoding_inference_dissociation.pdf

Usage:
    python fig_encoding_inference_dissociation.py
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

# reuse the exact pooled-coords helper from the current rho_sym brain figure
# (import, do NOT fork — no edit is made to that file)
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
try:
    from fig_trace_b_ofc_localization import pooled_coords_systems  # type: ignore
    from nilearn.plotting import plot_glass_brain
    _BRAIN_OK = True
except Exception as exc:  # pragma: no cover — brain panel is optional
    print(f"[warn] glass-brain machinery unavailable ({exc}); shipping panel a only")
    _BRAIN_OK = False

BASE = ROOT / "data" / "audit" / "inference_localization_rhosym"
OUT = ROOT / "data" / "preprint" / "figures" / "results_section1" / "fig_encoding_inference_dissociation.pdf"

# limbic/paralimbic core first (OFC + cingulate = the two dissociating hotspots)
SYS_ORDER = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal",
             "PFC", "parietal", "sensorimotor", "occipital"]
HOTSPOT = {"encoding": "OFC", "inference_pe": "cingulate"}

# contrast identity colours (encoding vs inference-specific)
C_ENC, C_INF = "#1f6fb2", "#d1491c"
CONTRAST_COL = {"encoding": C_ENC, "inference_pe": C_INF}
CONTRAST_LAB = {"encoding": "encoding (task_learn)",
                "inference_pe": "inference-specific ( | encoding)"}
# within-brain direction colours (enriched / depleted / n.s.)
C_ENR, C_DEP, C_NS = "#2a9d5c", "#d1791f", "#9a9a9a"
C_ENR_PALE = "#8fd0ab"          # length-assisted cingulate (distinct, not equal status)

CAP = 3.0                        # -log10 p cap (n_surr=1000 -> min p ~ 1/1001)
Q_SIG = 0.05


# --------------------------------------------------------------------------- #
# data (READ ONLY — no recompute)
# --------------------------------------------------------------------------- #
def load_sys(tag: str, inc: str) -> dict:
    """system -> dict(enr, p_up, p_lo, q_up, q_lo, signed) for band beta.

    DIRECTION-SEPARATED BH-FDR, matching the compute audits exactly: q_up is BH over
    matched_strength_p (enriched / CARRIER / LEAD), q_lo is BH over
    matched_strength_p_lower (depleted). Both families span the FULL set of
    system-granularity units (audit_158/audit_160 do not drop non-anatomical from the
    BH family; we simply do not *plot* those extra units). The star gating that uses
    these q's is applied per tag: encoding tests both directions, inference q_up only.
    """
    df = pd.read_csv(BASE / f"{tag}_localization_rhosym_{inc}.csv")
    d = df[df.granularity == "system"].copy()          # full BH family (audit convention)
    r = int(d.n_surr.iloc[0])
    p_up = np.clip(d.matched_strength_p.astype(float), 1.0 / (r + 1), 1.0)
    p_lo = np.clip(d.matched_strength_p_lower.astype(float), 1.0 / (r + 1), 1.0)
    q_up = np.asarray(bh_fdr(p_up.tolist()))
    q_lo = np.asarray(bh_fdr(p_lo.tolist()))
    enr = (d.M_obs > d.surr_median).to_numpy()
    p_dir = np.where(enr, p_up, p_lo)
    signed = np.where(enr, 1.0, -1.0) * np.minimum(-np.log10(p_dir), CAP)
    out = {}
    for u, e, pu, pl, qu, ql, sv in zip(d.unit, enr, p_up, p_lo, q_up, q_lo, signed):
        out[str(u)] = dict(enr=bool(e), p_up=float(pu), p_lo=float(pl),
                           q_up=float(qu), q_lo=float(ql), signed=float(sv))
    return out


def sig_sets(sysdir: dict, tag: str) -> tuple[set, set]:
    """(enriched-significant, depleted-significant) plotted systems, direction-separated
    BH, matching audit_158 (encoding: CARRIER q_up / DEPLETED q_lo) and audit_160
    (inference: LEAD q_up only — no depletion test)."""
    enr_sig, dep_sig = set(), set()
    for s in SYS_ORDER:
        r = sysdir.get(s)
        if r is None:
            continue
        if r["enr"] and r["q_up"] < Q_SIG:
            enr_sig.add(s)
        elif (not r["enr"]) and tag == "encoding" and r["q_lo"] < Q_SIG:
            dep_sig.add(s)
    return enr_sig, dep_sig


def lenmatched_cingulate_q(inc: str) -> float:
    """BH q for the cingulate under the LENGTH-MATCHED inference null (the caveat)."""
    df = pd.read_csv(BASE / f"lenmatched_null_rhosym_R200_{inc}.csv")
    row = df[df.unit == "cingulate"]
    return float(row.bh_q_system.iloc[0]) if not row.empty else float("nan")


# --------------------------------------------------------------------------- #
# panel a — systems-level diverging lollipop (the primary, self-sufficient panel)
# --------------------------------------------------------------------------- #
def draw_lollipop(ax, enc_inc, inf_inc, enc_exc, inf_exc):
    y0 = np.arange(len(SYS_ORDER))[::-1]           # OFC at top
    ypos = {s: y for s, y in zip(SYS_ORDER, y0)}
    dy = {"encoding": +0.20, "inference_pe": -0.20}
    prim = {"encoding": enc_inc, "inference_pe": inf_inc}
    excl = {"encoding": enc_exc, "inference_pe": inf_exc}
    sigs = {tag: sig_sets(prim[tag], tag) for tag in prim}

    ax.axvline(0, color="0.55", lw=0.9, zorder=1)  # enrichment / depletion divide

    for tag in ("encoding", "inference_pe"):
        col = CONTRAST_COL[tag]
        enr_sig, dep_sig = sigs[tag]
        for s in SYS_ORDER:
            rec = prim[tag].get(s)
            if rec is None:
                continue
            y = ypos[s] + dy[tag]
            x = rec["signed"]
            sig = (s in enr_sig) or (s in dep_sig)
            la = (tag == "inference_pe" and s == "cingulate")   # length-assisted cell

            ax.plot([0, x], [y, y], color=col, lw=1.7, alpha=0.9, zorder=2)

            # epilepsy-excluded replicate: small pale shadow marker (uniform, all cells)
            xe = excl[tag].get(s, {}).get("signed")
            if xe is not None:
                ax.scatter([xe], [y + 0.085], s=15, marker="o", facecolors="none",
                           edgecolors=col, linewidths=0.8, alpha=0.5, zorder=3)

            if la:
                # DISTINCT glyph: hatched diamond + dashed-look double ring.
                ax.scatter([x], [y], s=140, marker="D", facecolors=C_ENR_PALE,
                           edgecolors=col, linewidths=1.6, hatch="////", zorder=4)
                ax.scatter([x], [y], s=260, marker="D", facecolors="none",
                           edgecolors=col, linewidths=1.1, alpha=0.55, zorder=4)
            elif sig:
                ax.scatter([x], [y], s=95, marker="o", facecolors=col,
                           edgecolors="black", linewidths=1.2, zorder=4)
            else:
                ax.scatter([x], [y], s=70, marker="o", facecolors="white",
                           edgecolors=col, linewidths=1.4, zorder=4)

            if sig:  # uniform significance star on every BH-q<0.05 cell (tested direction)
                ax.annotate(r"$\ast$", (x, y), textcoords="offset points",
                            xytext=(7 if x >= 0 else -7, 3), ha="center",
                            fontsize=13, color=col, zorder=5)

    ax.set_yticks(y0)
    ax.set_yticklabels(SYS_ORDER)
    for lab in ax.get_yticklabels():        # emphasise the two dissociating hotspots
        if lab.get_text() in ("OFC", "cingulate"):
            lab.set_fontweight("bold")
    ax.set_ylim(y0.min() - 0.65, y0.max() + 0.65)
    ax.set_xlim(-CAP - 0.5, CAP + 0.5)
    ax.set_xlabel(r"$\leftarrow$ depleted      signed enrichment  "
                  r"$\mathrm{sign}(M_{\mathrm{obs}}-\tilde{M}_0)\cdot"
                  r"-\log_{10}p_{\mathrm{MS}}$      enriched $\rightarrow$")
    ax.spines[["top", "right"]].set_visible(False)


# --------------------------------------------------------------------------- #
# panels b/c — glass-brain confirmation (secondary; skipped if nilearn fiddly)
# --------------------------------------------------------------------------- #
def draw_brain(fig, rect, coords, systems, sysdir, tag):
    disp = plot_glass_brain(None, display_mode="lzr", axes=rect, figure=fig)
    winner = HOTSPOT[tag]

    def mask(names):
        return np.isin(systems, list(names))

    enr_sig, dep_sig = sig_sets(sysdir, tag)
    ns = set(SYS_ORDER) - enr_sig - dep_sig

    m_ns = mask(ns | (set(np.unique(systems)) - set(SYS_ORDER)))
    if m_ns.any():
        disp.add_markers(coords[m_ns], marker_color=C_NS, marker_size=6, alpha=0.20)
    m_dep = mask(dep_sig)
    if m_dep.any():
        disp.add_markers(coords[m_dep], marker_color=C_DEP, marker_size=26, alpha=0.9)
    m_enr = mask(enr_sig - {winner})
    if m_enr.any():
        disp.add_markers(coords[m_enr], marker_color=C_ENR, marker_size=26, alpha=0.9)
    m_win = mask({winner}) if winner in enr_sig else np.zeros(len(systems), bool)
    if np.any(m_win):
        if tag == "inference_pe":     # length-assisted: pale face + heavy ring (distinct)
            disp.add_markers(coords[m_win], marker_color=C_ENR_PALE, marker_size=70,
                             alpha=0.98, edgecolors=C_INF, linewidths=1.8)
        else:                         # encoding OFC: saturated + black ring
            disp.add_markers(coords[m_win], marker_color=C_ENR, marker_size=70,
                             alpha=0.98, edgecolors="black", linewidths=1.4)
    return disp


# --------------------------------------------------------------------------- #
def main():
    enc_inc, inf_inc = load_sys("encoding", "include"), load_sys("inference_pe", "include")
    enc_exc, inf_exc = load_sys("encoding", "exclude"), load_sys("inference_pe", "exclude")
    la_q = {"include": lenmatched_cingulate_q("include"),
            "exclude": lenmatched_cingulate_q("exclude")}

    # console provenance — direction-separated BH, matching audit_158 / audit_160
    print("== beta encoding-vs-inference dissociation (rho_sym); direction-separated BH ==")
    for tag, di, de in (("encoding", enc_inc, enc_exc),
                        ("inference_pe", inf_inc, inf_exc)):
        car = sorted([(u, r) for u, r in di.items()
                      if u in SYS_ORDER and r["enr"] and r["q_up"] < Q_SIG],
                     key=lambda kv: kv[1]["q_up"])
        for u, r in car:
            print(f"  {tag:13s} CARRIER/LEAD {u:16s} p_up={r['p_up']:.3f} "
                  f"q_up(inc)={r['q_up']:.3f} q_up(exc)={de.get(u, {}).get('q_up', float('nan')):.3f}")
        if tag == "encoding":
            dep = sorted([(u, r) for u, r in di.items()
                          if u in SYS_ORDER and (not r["enr"]) and r["q_lo"] < Q_SIG],
                         key=lambda kv: kv[1]["q_lo"])
            for u, r in dep:
                print(f"  {tag:13s} DEPLETED     {u:16s} p_lo={r['p_lo']:.3f} "
                      f"q_lo(inc)={r['q_lo']:.3f} q_lo(exc)={de.get(u, {}).get('q_lo', float('nan')):.3f}")
    print(f"  cingulate inference full-length q_up: include={inf_inc['cingulate']['q_up']:.3f}"
          f" exclude={inf_exc['cingulate']['q_up']:.3f}")
    print(f"  cingulate inference LENGTH-MATCHED q: include={la_q['include']:.3f}"
          f" exclude={la_q['exclude']:.3f}  (fails -> length-assisted lead)")

    fig = plt.figure(figsize=(13.6, 6.7))
    axa = fig.add_axes([0.075, 0.14, 0.45, 0.79])
    draw_lollipop(axa, enc_inc, inf_inc, enc_exc, inf_exc)
    fig.text(0.012, 0.965, r"$\mathbf{a}$", fontsize=17, va="top", fontweight="bold")

    if _BRAIN_OK:
        try:
            coords, systems = pooled_coords_systems()
            draw_brain(fig, (0.575, 0.545, 0.40, 0.40), coords, systems, enc_inc, "encoding")
            draw_brain(fig, (0.575, 0.115, 0.40, 0.40), coords, systems, inf_inc, "inference_pe")
            fig.text(0.565, 0.965, r"$\mathbf{b}$", fontsize=17, va="top",
                     fontweight="bold", color=C_ENC)
            fig.text(0.565, 0.535, r"$\mathbf{c}$", fontsize=17, va="top",
                     fontweight="bold", color=C_INF)
            print(f"  brain: pooled {len(coords)} contacts, "
                  f"{int((systems != 'non_anatomical').sum())} anatomical")
        except Exception as exc:                     # pragma: no cover
            print(f"[warn] brain panels failed ({exc}); shipping panel a only")

    # figure-level legend
    handles = [
        Line2D([], [], color=C_ENC, marker="o", ls="-", ms=8, label=CONTRAST_LAB["encoding"]),
        Line2D([], [], color=C_INF, marker="o", ls="-", ms=8, label=CONTRAST_LAB["inference_pe"]),
        Line2D([], [], color="0.35", marker="o", ls="", mfc="0.35", mec="black", mew=1.2,
               ms=9, label=r"matched-strength BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color="0.35", marker="o", ls="", mfc="white", mec="0.35", mew=1.4,
               ms=9, label=r"n.s. ($q\geq0.05$)"),
        Line2D([], [], color=C_INF, marker="D", ls="", mfc=C_ENR_PALE, mec=C_INF, mew=1.6,
               ms=11, label="inference lead — length-assisted (fails length-matched null)"),
        Line2D([], [], color=C_ENC, marker="o", ls="", mfc="none", mec=C_ENC, mew=0.9,
               ms=6, alpha=0.6, label="epilepsy-excluded replicate"),
    ]
    if _BRAIN_OK:
        handles += [
            Patch(facecolor=C_ENR, edgecolor="none", label="system enriched (brain)"),
            Patch(facecolor=C_DEP, edgecolor="none", label="system depleted (brain)"),
        ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.045),
               ncol=4, frameon=False, fontsize=8.8, handletextpad=0.5, columnspacing=1.3)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
