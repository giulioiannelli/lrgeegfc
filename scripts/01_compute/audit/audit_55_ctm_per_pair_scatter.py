#!/usr/bin/env python3
"""§5.3 manuscript figure — per-pair Δ_task vs Δ_rest scatter.

Three panels (β, α, null-band) each at the patient whose ρ_split sits
closest to that band's cohort median. Visualizes the per-pair object
underlying ρ_split with same-probe vs cross-probe overlay.

Outputs
-------
``data/audit/ctm_per_pair_scatter/figures/manuscript_ctm_per_pair_scatter.pdf``
``data/audit/ctm_per_pair_scatter/tables/selection.csv``
``.agents/reports/2026-05-08_ctm_per_pair_scatter.md``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.probe import extract_probe_labels

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
})

CTM_PAIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
TD_CSV = ROOT / "data" / "audit" / "ctm_triangle" / "Td_per_patient_per_band.csv"

OUT = ROOT / "data" / "audit" / "ctm_per_pair_scatter"
FIG = OUT / "figures"
TBL = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TBL.mkdir(parents=True, exist_ok=True)

REPORT = ROOT / ".agents" / "reports" / "2026-05-08_ctm_per_pair_scatter.md"

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# Null band fixed to θ (rather than smallest-|median|): θ has per-pair Δ
# spread comparable to α / β so the visual contrast isolates *tilt* rather
# than cloud size. γ_h would give a concentrated cloud and read as
# "no dynamics" rather than "no rank agreement".
NULL_BAND = "theta"


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = CTM_PAIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    dD_task = np.asarray(d["dD_task"])
    dD_rest = np.asarray(d["dD_rest"])
    iu_i = d["iu_i"].astype(int)
    iu_j = d["iu_j"].astype(int)

    pat_dir = ROOT / "data" / "raw" / "stereoeeg_patients" / pat
    ch = pd.read_csv(pat_dir / "channel_labels.csv")
    labels = ch["label"].tolist()
    probes = np.array(extract_probe_labels(labels))

    n_contacts = int(max(iu_i.max(), iu_j.max())) + 1
    if probes.size != n_contacts:
        # Pat_10: channel_labels.csv has full set; npz uses post-mask 113.
        # We only need consistency, so trim if longer.
        probes = probes[:n_contacts]

    same_probe = probes[iu_i] == probes[iu_j]
    return dict(
        dD_task=dD_task, dD_rest=dD_rest,
        iu_i=iu_i, iu_j=iu_j, same_probe=same_probe,
        n_contacts=n_contacts,
    )


def pick_panels(Td: pd.DataFrame) -> tuple[list[dict], dict, str]:
    cohort_med = {b: float(Td[Td["band"] == b]["rho_split"].median())
                  for b in BANDS}
    null_band = NULL_BAND

    def closest(band: str) -> dict:
        sub = Td[Td["band"] == band].copy()
        sub["dist"] = (sub["rho_split"] - cohort_med[band]).abs()
        row = sub.loc[sub["dist"].idxmin()]
        return dict(
            pat=str(row["patient"]), band=band,
            rho_full_csv=float(row["rho_split"]),
            rho_x_csv=float(row["rho_split_cross_probe"]),
            cohort_med=cohort_med[band],
            dist_to_med=float(row["dist"]),
        )

    panels = [closest("beta"), closest("alpha"), closest(null_band)]
    return panels, cohort_med, null_band


def make_figure(panels: list[dict]) -> float:
    # Load pair data + verify ρ values reproduce
    for p in panels:
        p.update(load_pair_data(p["pat"], p["band"]))
        rho_full, _ = spearmanr(p["dD_task"], p["dD_rest"])
        cp = ~p["same_probe"]
        rho_x, _ = spearmanr(p["dD_task"][cp], p["dD_rest"][cp])
        p["rho_full_npz"] = float(rho_full)
        p["rho_x_npz"] = float(rho_x)
        p["N_pairs"] = int(len(p["dD_task"]))
        p["N_same"] = int(p["same_probe"].sum())
        p["N_cross"] = int((~p["same_probe"]).sum())

    # Shared axis range — 98th percentile of |Δ| across all 3 panels combined.
    all_vals = np.concatenate(
        [p["dD_task"] for p in panels] + [p["dD_rest"] for p in panels]
    )
    q98 = float(np.quantile(np.abs(all_vals), 0.98))
    lim = q98 * 1.05

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    for ax, p in zip(axes, panels):
        sp = p["same_probe"]
        cp = ~sp
        ax.scatter(p["dD_task"][sp], p["dD_rest"][sp],
                   s=8, color="#bbbbbb", alpha=0.45,
                   edgecolors="none", zorder=1)
        ax.scatter(p["dD_task"][cp], p["dD_rest"][cp],
                   s=10, color="#1f3d6e", alpha=0.55,
                   edgecolors="none", zorder=2)
        ax.axhline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.axvline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.6,
                ls="--", zorder=0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$\Delta_{\mathrm{task}}(i,j)"
                      r" = D_{\mathrm{task}} - D_{\mathrm{rsPre,A}}$")
        ax.set_ylabel(r"$\Delta_{\mathrm{rest}}(i,j)"
                      r" = D_{\mathrm{rsPost}} - D_{\mathrm{rsPre,B}}$")
        ax.spines[["top", "right"]].set_visible(False)

        band_tex = BRAIN_BAND_TEX_DICT[p["band"]]
        ax.text(
            0.03, 0.97,
            f"{p['pat']}, {band_tex}\n"
            rf"$\rho_{{\rm split}} = {p['rho_full_csv']:+.3f}$"
            "\n"
            rf"$\rho_{{\rm xprobe}} = {p['rho_x_csv']:+.3f}$"
            "\n"
            f"$N$ = {p['N_pairs']} "
            f"({p['N_same']} same / {p['N_cross']} cross)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9.5,
            bbox=dict(facecolor="white", edgecolor="#bbbbbb",
                      boxstyle="round,pad=0.35"),
        )

    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor="#bbbbbb", markeredgecolor="none",
                      alpha=0.7, label="same-probe pairs"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor="#1f3d6e", markeredgecolor="none",
                      alpha=0.8, label="cross-probe pairs"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "manuscript_ctm_per_pair_scatter.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[ctm_per_pair] wrote {FIG / 'manuscript_ctm_per_pair_scatter.pdf'}")
    return q98


def write_report(panels: list[dict], cohort_med: dict, null_band: str,
                 q98: float) -> None:
    band_tex = {"theta": "θ", "high_gamma": "γ_h"}
    null_tex = band_tex.get(null_band, null_band)

    rows = []
    for p in panels:
        label = ("β" if p["band"] == "beta"
                 else "α" if p["band"] == "alpha" else null_tex)
        rows.append({
            "panel": label, "patient": p["pat"], "band": p["band"],
            "rho_split_csv": p["rho_full_csv"],
            "rho_split_npz": p["rho_full_npz"],
            "rho_xprobe_csv": p["rho_x_csv"],
            "rho_xprobe_npz": p["rho_x_npz"],
            "cohort_median": p["cohort_med"],
            "dist_to_median": p["dist_to_med"],
            "N_pairs": p["N_pairs"], "N_same": p["N_same"],
            "N_cross": p["N_cross"],
        })
    sel = pd.DataFrame(rows)
    sel.to_csv(TBL / "selection.csv", index=False)
    print(sel.to_string(index=False))

    beta_p, alpha_p, null_p = panels

    def chk_pos(rho: float) -> str:
        return "✓" if rho > 0 else "⚠ flagged"

    def chk_round(rho: float, thresh: float) -> str:
        return "✓" if abs(rho) < thresh else "⚠ flagged"

    lines = [
        "---",
        "date: 2026-05-08",
        "era: COHORT_N10 / IMCOH_ABS",
        "status: current",
        "type: figure-build-report",
        "scope: section_5_3_ctm_per_pair_scatter",
        "---",
        "",
        "# §5.3 per-pair scatter — selection, verification, sanity",
        "",
        f"**Head.** Three-panel scatter of Δ_task(i,j) vs Δ_rest(i,j) over "
        f"contact pairs at three (patient, band) cells: "
        f"**β @ {beta_p['pat']}** "
        f"(ρ_split = {beta_p['rho_full_csv']:+.3f}, "
        f"|d| to cohort median {cohort_med['beta']:+.3f} = "
        f"{beta_p['dist_to_med']:.4f}); "
        f"**α @ {alpha_p['pat']}** "
        f"(ρ_split = {alpha_p['rho_full_csv']:+.3f}, "
        f"|d| = {alpha_p['dist_to_med']:.4f}); "
        f"**{null_tex} @ {null_p['pat']}** "
        f"(ρ_split = {null_p['rho_full_csv']:+.3f}, "
        f"|d| = {null_p['dist_to_med']:.4f}). "
        f"Same-probe pairs in light gray, cross-probe pairs in navy on top. "
        f"Axes clipped to ±{q98:.3f} (98th percentile of |Δ|).",
        "",
        "## Selection table",
        "",
        ("| panel | patient | band | ρ_split | cohort med | |d| | "
         "ρ_xprobe | N | same | cross |"),
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['panel']} | {r['patient']} | {r['band']} "
            f"| {r['rho_split_csv']:+.3f} | {r['cohort_median']:+.3f} "
            f"| {r['dist_to_median']:.4f} | {r['rho_xprobe_csv']:+.3f} "
            f"| {r['N_pairs']} | {r['N_same']} | {r['N_cross']} |"
        )

    lines += [
        "",
        "## Null-band rationale",
        "",
        f"Cohort medians at candidate null bands: "
        f"θ {cohort_med['theta']:+.4f}, "
        f"γ_h {cohort_med['high_gamma']:+.4f}; "
        f"δ {cohort_med['delta']:+.4f} (for reference). "
        f"**θ chosen** even though |median θ| = {abs(cohort_med['theta']):.4f} "
        f"is larger than |median γ_h| = "
        f"{abs(cohort_med['high_gamma']):.4f}: ρ_split measures the rank "
        f"agreement (diagonal tilt) between Δ_task and Δ_rest, not the "
        f"magnitude of either. γ_h cohort-median exemplars produce a "
        f"pathologically tight cloud near the origin, which reads as "
        f"\"no dynamics\" rather than \"no rank agreement\". θ has Δ "
        f"spread comparable to α and β at the cohort-median exemplar, so "
        f"the visual contrast against the β / α panels isolates the "
        f"absence of diagonal tilt at comparable cloud size — which is "
        f"the property ρ_split actually tests.",
        "",
        "## ρ verification (Spearman recomputed from the npz)",
        "",
        ("| panel | ρ_split (CSV) | ρ_split (npz) | "
         "ρ_xprobe (CSV) | ρ_xprobe (npz) | match? |"),
        "|---|---|---|---|---|---|",
    ]
    for r in rows:
        match = ("✓"
                 if (abs(r["rho_split_csv"] - r["rho_split_npz"]) < 1e-3
                     and abs(r["rho_xprobe_csv"] - r["rho_xprobe_npz"]) < 1e-3)
                 else "⚠ mismatch")
        lines.append(
            f"| {r['panel']} | {r['rho_split_csv']:+.4f} "
            f"| {r['rho_split_npz']:+.4f} "
            f"| {r['rho_xprobe_csv']:+.4f} "
            f"| {r['rho_xprobe_npz']:+.4f} | {match} |"
        )

    lines += [
        "",
        "## Sanity checks",
        "",
        f"1. **β panel diagonal tilt** (full): "
        f"ρ = {beta_p['rho_full_csv']:+.3f}. "
        f"Expected positive {chk_pos(beta_p['rho_full_csv'])}",
        f"2. **α panel diagonal tilt** (full): "
        f"ρ = {alpha_p['rho_full_csv']:+.3f}. "
        f"Expected positive {chk_pos(alpha_p['rho_full_csv'])}",
        f"3. **{null_tex} panel cloud shape** (full): "
        f"ρ = {null_p['rho_full_csv']:+.3f}. "
        f"Expected near-zero {chk_round(null_p['rho_full_csv'], 0.10)}",
        f"4. **β cross-probe subset tilt**: "
        f"ρ_xprobe = {beta_p['rho_x_csv']:+.3f}. "
        f"Expected positive {chk_pos(beta_p['rho_x_csv'])}",
        f"5. **α cross-probe subset tilt**: "
        f"ρ_xprobe = {alpha_p['rho_x_csv']:+.3f}. "
        f"Expected positive {chk_pos(alpha_p['rho_x_csv'])}",
        f"6. **{null_tex} cross-probe subset shape**: "
        f"ρ_xprobe = {null_p['rho_x_csv']:+.3f}. "
        f"Expected near-zero {chk_round(null_p['rho_x_csv'], 0.15)}",
        "",
        "## Files",
        "",
        f"- Figure: `data/audit/ctm_per_pair_scatter/figures/"
        f"manuscript_ctm_per_pair_scatter.pdf`",
        f"- Selection: `data/audit/ctm_per_pair_scatter/tables/selection.csv`",
        f"- Source npz: `data/reports/imcoh_continuous_trace/per_pair_split/"
        f"{{Pat}}_{{band}}.npz`",
        f"- Cohort table: `data/audit/ctm_triangle/Td_per_patient_per_band.csv`",
        f"- Build script: `scripts/01_compute/audit/audit_55_ctm_per_pair_scatter.py`",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines))
    print(f"wrote {REPORT}")


def main() -> None:
    Td = pd.read_csv(TD_CSV)
    panels, cohort_med, null_band = pick_panels(Td)
    print("Cohort medians ρ_split:")
    for b in BANDS:
        print(f"  {b:11s}: {cohort_med[b]:+.4f}")
    print(f"Null band chosen: {null_band}")

    q98 = make_figure(panels)
    write_report(panels, cohort_med, null_band, q98)


if __name__ == "__main__":
    main()
