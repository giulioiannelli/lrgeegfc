#!/usr/bin/env python3
"""Critical validation of epileptic-node findings under ImCoh-LRG.

Memory-safe driver: runs per-patient computations in isolated subprocesses
(see _validation_one_patient.py), then merges CSVs and produces aggregate
statistics, a summary figure, and a report.

  V1. Cross-probe epi-epi edge strength
  V2. Cross-probe nearest-neighbor enrichment
  V3. Multiple comparison burden (BH-FDR)
  V4. Epileptic node spatial distribution
  V5. Logic chain (edge ratio ↔ ultrametric ratio)

Run:  python scripts/08_epileptic/epileptic_imcoh_validation.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io import load_epileptic_nodes

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)

OUT = FIGURES_ROOT / "epileptic_imcoh_validation"
OUT.mkdir(parents=True, exist_ok=True)

PAT_COLORS = {
    "Pat_02": "#e74c3c", "Pat_03": "#3498db", "Pat_05": "#9b59b6",
    "Pat_07": "#2ecc71", "Pat_08": "#f39c12",
}

LINES: list[str] = []


def R(s: str = "") -> None:
    print(s)
    LINES.append(s)


def load_ch(pat):
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not p.exists():
        p = SEEG_DATAPATH / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [str(l).strip('"').split(",")[0].strip().replace(" ", "")
            for l in df.iloc[:, 0]]


def probe(label):
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def benjamini_hochberg(pvals, alpha=0.05):
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=bool)
    order = np.argsort(pvals)
    sorted_p = np.asarray(pvals)[order]
    thresh = alpha * np.arange(1, n + 1) / n
    below = sorted_p <= thresh
    if not below.any():
        return np.zeros(n, dtype=bool)
    k_max = np.where(below)[0][-1]
    discoveries = np.zeros(n, dtype=bool)
    discoveries[order[: k_max + 1]] = True
    return discoveries


# ── Step 1: Run per-patient workers ───────────────────────────────────
R("=" * 90)
R("CRITICAL VALIDATION — EPILEPTIC ImCoh-LRG FINDINGS")
R("=" * 90)
R()
R("Running per-patient workers (isolated subprocesses for memory safety)")
R("-" * 90)

worker = Path(__file__).parent / "_validation_one_patient.py"
for pat in PATIENTS:
    print(f"  Running {pat}...", flush=True)
    r = subprocess.run(
        [sys.executable, str(worker), pat],
        capture_output=True, text=True, timeout=120,
    )
    if r.returncode != 0:
        R(f"    ERROR ({pat}): {r.stderr}")
    else:
        R(f"  {r.stdout.strip()}")

R()

# ── Step 2: Merge CSVs ────────────────────────────────────────────────
v1df = pd.concat([pd.read_csv(OUT / f"_v1_{p}.csv") for p in PATIENTS], ignore_index=True)
v2df = pd.concat([pd.read_csv(OUT / f"_v2_{p}.csv") for p in PATIENTS], ignore_index=True)
v5df_raw = pd.concat([pd.read_csv(OUT / f"_v5_{p}.csv") for p in PATIENTS], ignore_index=True)

v1df["bh_sig"] = benjamini_hochberg(v1df["p_value"].values)
v2df["bh_sig"] = benjamini_hochberg(v2df["p_value"].values)

# ── V4: Spatial distribution ──────────────────────────────────────────
R("V4: Epileptic node spatial distribution")
R("-" * 90)
R()

v4_data = {}
for pat in PATIENTS:
    ch = load_ch(pat)
    epi_all = load_epileptic_nodes(pat)
    epi_set = set(epi_all)
    epi_idx = [i for i, l in enumerate(ch) if l in epi_set]
    probes_arr = np.array([probe(l) for l in ch])

    n_epi_epi = 0
    n_sp = 0
    for ii in range(len(epi_idx)):
        for jj in range(ii + 1, len(epi_idx)):
            n_epi_epi += 1
            if probes_arr[epi_idx[ii]] == probes_arr[epi_idx[jj]]:
                n_sp += 1

    epi_probes = {}
    for i in epi_idx:
        pr = probes_arr[i]
        epi_probes[pr] = epi_probes.get(pr, 0) + 1

    v4_data[pat] = dict(
        N=len(ch), n_epi=len(epi_idx), n_epi_epi=n_epi_epi,
        n_sp=n_sp, n_cp=n_epi_epi - n_sp,
        frac_sp=n_sp / max(n_epi_epi, 1) * 100,
        epi_probes=epi_probes,
        total_probes=len(set(probes_arr)),
    )

R(f"{'Patient':<10s} {'N':>4s} {'N_epi':>6s} {'epi pairs':>10s} "
  f"{'same-probe':>11s} {'cross-probe':>12s} {'%same':>6s}")
R("-" * 70)
for pat in PATIENTS:
    d = v4_data[pat]
    R(f"  {pat:<8s} {d['N']:>4d} {d['n_epi']:>6d} {d['n_epi_epi']:>10d} "
      f"{d['n_sp']:>11d} {d['n_cp']:>12d} {d['frac_sp']:>5.1f}%")
R()

R("Epileptic probe composition:")
for pat in PATIENTS:
    d = v4_data[pat]
    comp = ", ".join(f"{pr}({n})" for pr, n in sorted(d["epi_probes"].items()))
    R(f"  {pat}: {len(d['epi_probes'])}/{d['total_probes']} probes: {comp}")
R()

# ── V1: Summary ────────────────────────────────────────────────────────
R("V1: Cross-probe epi-epi edge strength (same-probe pairs removed)")
R("-" * 90)
R()
R("Original finding (B7): epi-epi ImCoh edges 2-5× stronger in ALL 5 patients.")
R("Question: Does this survive after removing same-probe pairs?")
R()

n_sig_raw = int((v1df["p_value"] < 0.05).sum())
n_sig_bh = int(v1df["bh_sig"].sum())
R(f"Total conditions: {len(v1df)}")
R(f"Raw p < 0.05: {n_sig_raw}/{len(v1df)} ({n_sig_raw/len(v1df):.0%})")
R(f"BH-FDR q < 0.05: {n_sig_bh}/{len(v1df)} ({n_sig_bh/len(v1df):.0%})")
R(f"Mean cross-probe ratio: {v1df['ratio'].mean():.3f}")
R(f"Mean all-pairs ratio:   {v1df['ratio_all'].mean():.3f}")
R()

for pat in PATIENTS:
    sub = v1df[v1df["patient"] == pat]
    ns_raw = int((sub["p_value"] < 0.05).sum())
    ns_bh = int(sub["bh_sig"].sum())
    mr = sub["ratio"].mean()
    mr_all = sub["ratio_all"].mean()
    R(f"  {pat}: CP={mr:.3f} (all={mr_all:.3f}, Δ={mr/mr_all - 1:+.0%}), "
      f"raw={ns_raw}/{len(sub)}, BH={ns_bh}/{len(sub)}")
R()

R("V1 by band (cross-probe):")
for band in BANDS:
    sub = v1df[v1df["band"] == band]
    mr = sub["ratio"].mean()
    n_bh = int(sub["bh_sig"].sum())
    R(f"  {band:12s}: ratio={mr:.3f}, BH sig={n_bh}/{len(sub)}")
R()

v1df.to_csv(OUT / "V1_cross_probe_edge_strength.csv", index=False)
R(f"  → {OUT / 'V1_cross_probe_edge_strength.csv'}")
R()

# ── V2: Summary ────────────────────────────────────────────────────────
R("V2: Cross-probe nearest-neighbor enrichment")
R("-" * 90)
R()
R("Original finding (A3): NN of epileptic nodes is also epileptic 2.8-7.3× expected.")
R("Question: Does this survive when NN search excludes same-probe candidates?")
R()

n_sig_raw = int((v2df["p_value"] < 0.05).sum())
n_sig_bh = int(v2df["bh_sig"].sum())
R(f"Total conditions: {len(v2df)}")
R(f"Raw p < 0.05: {n_sig_raw}/{len(v2df)} ({n_sig_raw/len(v2df):.0%})")
R(f"BH-FDR q < 0.05: {n_sig_bh}/{len(v2df)} ({n_sig_bh/len(v2df):.0%})")
R(f"Mean enrichment: {v2df['enrichment'].mean():.2f}×")
R()

for pat in PATIENTS:
    sub = v2df[v2df["patient"] == pat]
    ns_raw = int((sub["p_value"] < 0.05).sum())
    ns_bh = int(sub["bh_sig"].sum())
    me = sub["enrichment"].mean()
    R(f"  {pat}: enrichment={me:.2f}×, raw={ns_raw}/{len(sub)}, BH={ns_bh}/{len(sub)}")
R()

R("V2 by band:")
for band in BANDS:
    sub = v2df[v2df["band"] == band]
    me = sub["enrichment"].mean()
    n_bh = int(sub["bh_sig"].sum())
    R(f"  {band:12s}: enrichment={me:.2f}×, BH sig={n_bh}/{len(sub)}")
R()

v2df.to_csv(OUT / "V2_cross_probe_nn_enrichment.csv", index=False)
R(f"  → {OUT / 'V2_cross_probe_nn_enrichment.csv'}")
R()

# ── V3: Multiple comparison ───────────────────────────────────────────
R("V3: Multiple comparison burden analysis")
R("-" * 90)
R()

analyses = {
    "V1 (cross-probe edges)": v1df,
    "V2 (cross-probe NN)": v2df,
}
for csv_name, key in [
    ("probe_matched_null.csv", "A1 (probe-matched null)"),
    ("nn_enrichment.csv", "A3 (original NN)"),
]:
    p_orig = FIGURES_ROOT / "epileptic_imcoh_deep" / csv_name
    if p_orig.exists():
        analyses[key] = pd.read_csv(p_orig)

R(f"{'Analysis':<30s} {'N':>5s} {'Raw':>6s} {'E[FP]':>7s} {'BH':>5s} {'Bonf':>5s}")
R("-" * 65)
mc_data = {}
for name, df in analyses.items():
    pvals = df["p_value"].values
    n = len(pvals)
    n_raw = int((pvals < 0.05).sum())
    expected_fp = n * 0.05
    bh = benjamini_hochberg(pvals)
    n_bh = int(bh.sum())
    n_bonf = int((pvals < 0.05 / n).sum())
    R(f"  {name:<28s} {n:>5d} {n_raw:>6d} {expected_fp:>7.1f} {n_bh:>5d} {n_bonf:>5d}")
    mc_data[name] = dict(n=n, raw=n_raw, expected_fp=expected_fp, bh=n_bh, bonf=n_bonf)
R()
R("Reading this: Raw >> E[FP] with BH retaining most → real signal.")
R()

# ── V5: Logic chain ───────────────────────────────────────────────────
R("V5: Logic chain — edge ratio ↔ ultrametric distance ratio")
R("-" * 90)
R()
R("Claim: strong epi-epi ImCoh → close in LRG dendrogram.")
R("Test: Spearman(edge_ratio, ultrametric_ratio) — expected negative.")
R()

v5_merged = v5df_raw.merge(v1df[["patient", "phase", "band", "ratio"]],
                           on=["patient", "phase", "band"])
rho_overall, p_overall = np.nan, np.nan
if len(v5_merged) > 5:
    rho_overall, p_overall = spearmanr(v5_merged["ratio"], v5_merged["um_ratio"])
    consistent = rho_overall < 0 and p_overall < 0.05
    R(f"Overall: ρ = {rho_overall:.3f}, p = {p_overall:.2e} (n = {len(v5_merged)})")
    R(f"  → {'CONSISTENT' if consistent else 'INCONSISTENT'}")
    R()
    R("Per patient:")
    for pat in PATIENTS:
        sub = v5_merged[v5_merged["patient"] == pat]
        if len(sub) < 5:
            continue
        rp, pp = spearmanr(sub["ratio"], sub["um_ratio"])
        R(f"  {pat}: ρ = {rp:.3f}, p = {pp:.3f} (n = {len(sub)})")
    R()

v5_merged.to_csv(OUT / "V5_logic_chain.csv", index=False)
R(f"  → {OUT / 'V5_logic_chain.csv'}")
R()

# ── Summary figure ────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# (0,0) V1 ratio by patient
ax = axes[0, 0]
for pi, pat in enumerate(PATIENTS):
    sub = v1df[v1df["patient"] == pat]
    ax.boxplot(sub["ratio"].values, positions=[pi], widths=0.6,
               boxprops=dict(color=PAT_COLORS[pat]),
               medianprops=dict(color=PAT_COLORS[pat]),
               whiskerprops=dict(color=PAT_COLORS[pat]),
               capprops=dict(color=PAT_COLORS[pat]),
               flierprops=dict(markeredgecolor=PAT_COLORS[pat], markersize=3))
ax.axhline(1.0, color="k", ls="--", lw=0.8)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels(PATIENTS, fontsize=9)
ax.set_ylabel("Edge ratio (epi-epi / non-non)")
ax.set_title("V1: Cross-probe edge strength", fontsize=10)

# (0,1) V2 enrichment by patient
ax = axes[0, 1]
for pi, pat in enumerate(PATIENTS):
    sub = v2df[v2df["patient"] == pat]
    if len(sub) == 0:
        continue
    ax.boxplot(sub["enrichment"].values, positions=[pi], widths=0.6,
               boxprops=dict(color=PAT_COLORS[pat]),
               medianprops=dict(color=PAT_COLORS[pat]),
               whiskerprops=dict(color=PAT_COLORS[pat]),
               capprops=dict(color=PAT_COLORS[pat]),
               flierprops=dict(markeredgecolor=PAT_COLORS[pat], markersize=3))
ax.axhline(1.0, color="k", ls="--", lw=0.8)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels(PATIENTS, fontsize=9)
ax.set_ylabel("NN enrichment (obs / expected)")
ax.set_title("V2: Cross-probe NN enrichment", fontsize=10)

# (0,2) V4 same-probe %
ax = axes[0, 2]
sp_fracs = [v4_data[p]["frac_sp"] for p in PATIENTS]
ax.bar(range(len(PATIENTS)), sp_fracs,
       color=[PAT_COLORS[p] for p in PATIENTS], alpha=0.7)
for i, v in enumerate(sp_fracs):
    ax.text(i, v + 1, f"{v:.0f}%", ha="center", fontsize=9)
ax.set_xticks(range(len(PATIENTS)))
ax.set_xticklabels(PATIENTS, fontsize=9)
ax.set_ylabel("% epi-epi pairs on same probe")
ax.set_title("V4: Spatial bias potential", fontsize=10)
ax.set_ylim(0, 100)

# (1,0) all-pairs vs cross-probe
ax = axes[1, 0]
for pat in PATIENTS:
    sub = v1df[v1df["patient"] == pat]
    ax.scatter(sub["ratio_all"], sub["ratio"], s=12, alpha=0.4,
               color=PAT_COLORS[pat], label=pat)
lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])]
ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
ax.set_xlabel("All-pairs edge ratio")
ax.set_ylabel("Cross-probe edge ratio")
ax.set_title("V1: Effect of removing same-probe", fontsize=10)
ax.legend(fontsize=7, markerscale=2)

# (1,1) V5 logic chain
ax = axes[1, 1]
if len(v5_merged) > 5:
    for pat in PATIENTS:
        sub = v5_merged[v5_merged["patient"] == pat]
        ax.scatter(sub["ratio"], sub["um_ratio"], s=12, alpha=0.4,
                   color=PAT_COLORS[pat], label=pat)
    ax.set_xlabel("Cross-probe edge ratio (epi/non)")
    ax.set_ylabel("Ultrametric distance ratio (epi/non)")
    ax.set_title(f"V5: Logic chain (ρ={rho_overall:.2f}, p={p_overall:.1e})",
                 fontsize=10)
    ax.legend(fontsize=7, markerscale=2)
    ax.axhline(1.0, color="gray", ls=":", lw=0.5)
    ax.axvline(1.0, color="gray", ls=":", lw=0.5)

# (1,2) V3 multiple comparison
ax = axes[1, 2]
labels = list(mc_data.keys())
x = np.arange(len(labels))
w = 0.25
ax.bar(x - w, [mc_data[l]["raw"] for l in labels], w,
       label="Raw p<.05", alpha=0.7, color="#3498db")
ax.bar(x, [mc_data[l]["bh"] for l in labels], w,
       label="BH-FDR", alpha=0.7, color="#2ecc71")
ax.bar(x + w, [mc_data[l]["bonf"] for l in labels], w,
       label="Bonferroni", alpha=0.7, color="#e74c3c")
for i, l in enumerate(labels):
    ax.plot([i - 1.5 * w, i + 1.5 * w],
            [mc_data[l]["expected_fp"]] * 2, "k--", lw=1.5)
ax.set_xticks(x)
ax.set_xticklabels([l.split("(")[0].strip() for l in labels],
                   fontsize=8, rotation=15)
ax.set_ylabel("# significant")
ax.set_title("V3: Multiple comparison correction", fontsize=10)
ax.legend(fontsize=7)

fig.suptitle("Critical Validation — Epileptic ImCoh-LRG Findings\n"
             "All cross-probe tests remove same-probe pairs",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT / "validation_summary.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "validation_summary.png", bbox_inches="tight", dpi=200)
plt.close(fig)
R(f"  → {OUT / 'validation_summary.pdf'}")
R()

# ── Verdict ────────────────────────────────────────────────────────────
R("=" * 90)
R("VALIDATION SUMMARY")
R("=" * 90)
R()

v1_survives = v1df["ratio"].mean() > 1.0 and int(v1df["bh_sig"].sum()) > 0
v2_survives = v2df["enrichment"].mean() > 1.0 and int(v2df["bh_sig"].sum()) > 0
v5_consistent = (len(v5_merged) > 5 and not np.isnan(rho_overall)
                 and rho_overall < 0 and p_overall < 0.05)

R(f"V1 (cross-probe edges):    {'SURVIVES' if v1_survives else 'FAILS'}")
R(f"    Mean ratio = {v1df['ratio'].mean():.3f}, "
  f"BH-FDR = {int(v1df['bh_sig'].sum())}/{len(v1df)}")
R()
R(f"V2 (cross-probe NN):       {'SURVIVES' if v2_survives else 'FAILS'}")
R(f"    Mean enrichment = {v2df['enrichment'].mean():.2f}×, "
  f"BH-FDR = {int(v2df['bh_sig'].sum())}/{len(v2df)}")
R()
R(f"V3 (multi. comparison):    see table above")
R(f"V4 (spatial bias):         see table above")
R()
R(f"V5 (logic chain):          {'CONSISTENT' if v5_consistent else 'INCONSISTENT'}")
if not np.isnan(rho_overall):
    R(f"    ρ = {rho_overall:.3f}, p = {p_overall:.2e}")
R()

# Save report
report_path = OUT / "VALIDATION_REPORT.md"
with open(report_path, "w") as f:
    f.write("\n".join(LINES))
R(f"Full report: {report_path}")

# Clean up per-patient intermediate CSVs
for p in PATIENTS:
    for v in ("v1", "v2", "v5"):
        (OUT / f"_{v}_{p}.csv").unlink(missing_ok=True)
