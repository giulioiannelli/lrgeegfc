#!/usr/bin/env python3
"""Do epileptic nodes participate in H2 beta reorganization?

Three-part test:
  1. Pair-class decomposed disagreement
     At each height h, compute partition disagreement between phase pairs
     restricted to each pair class (epi-epi, non-non, cross, cross-probe-only).
     H2a_class(h) = d_class(rsPre, rsPost) − d_class(taskTest, rsPost).
     Large H2a_non-non + small H2a_epi-epi = epi tissue insulated from
     reorganization.

  2. Ablation (separate subprocess per patient, memory-isolated)
     Zero out epi-epi ImCoh edges, recompute LRG per phase, recompute H2a.
     If H2a beta magnitude is unchanged → epi tissue doesn't drive it.
     If H2a beta magnitude drops → epi tissue contributes to the reorganization.

  3. Summary
     Per-band H2a decomposition figure + ablation delta table.

Run:  python scripts/08_epileptic/epileptic_h2_stratified.py
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
from scipy.cluster.hierarchy import fcluster

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)

OUT = FIGURES_ROOT / "epileptic_h2_stratified"
OUT.mkdir(parents=True, exist_ok=True)

H_GRID = np.linspace(0.01, 0.99, 60)  # 60 heights, tractable

PAT_COLORS = {
    "Pat_02": "#e74c3c", "Pat_03": "#3498db", "Pat_05": "#9b59b6",
    "Pat_07": "#2ecc71", "Pat_08": "#f39c12",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}


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


def pairwise_same_matrix(labels):
    """NxN boolean: same[i,j] = True iff labels[i]==labels[j]."""
    labels = np.asarray(labels)
    return labels[:, None] == labels[None, :]


def stratified_disagreement(labels1, labels2, pair_mask):
    """Fraction of pairs (from pair_mask, upper-triangle) that disagree on
    cluster membership between two partitions."""
    same1 = pairwise_same_matrix(labels1)
    same2 = pairwise_same_matrix(labels2)
    ut = np.triu(np.ones(same1.shape, dtype=bool), k=1) & pair_mask
    if not ut.any():
        return np.nan
    disagree = (same1 != same2) & ut
    return disagree.sum() / ut.sum()


def compute_h2a_stratified(linkages_by_phase, pair_masks):
    """Return dict {class: H2a(h) array} for the given (patient, band)."""
    results = {cls: np.full(len(H_GRID), np.nan) for cls in pair_masks}
    if "rsPre" not in linkages_by_phase or "rsPost" not in linkages_by_phase \
            or "taskTest" not in linkages_by_phase:
        return results

    for ih, h in enumerate(H_GRID):
        lab_pre = fcluster(linkages_by_phase["rsPre"], t=h, criterion="distance")
        lab_post = fcluster(linkages_by_phase["rsPost"], t=h, criterion="distance")
        lab_tt = fcluster(linkages_by_phase["taskTest"], t=h, criterion="distance")

        for cls, mask in pair_masks.items():
            d_prepost = stratified_disagreement(lab_pre, lab_post, mask)
            d_ttpost = stratified_disagreement(lab_tt, lab_post, mask)
            if not (np.isnan(d_prepost) or np.isnan(d_ttpost)):
                results[cls][ih] = d_prepost - d_ttpost
    return results


# ── Part 1: stratified H2a on cached LRG ──────────────────────────────
print("=" * 90)
print("PART 1: Pair-class-stratified H2a from existing LRG cache (|ImCoh|)")
print("=" * 90)

h2a_rows = []  # long format: patient, band, h, class, value
for pat in PATIENTS:
    ch = load_ch(pat)
    N = len(ch)
    epi_set = set(load_epileptic_nodes(pat))
    epi_mask = np.array([l in epi_set for l in ch])
    probes_arr = np.array([probe(l) for l in ch])

    # NxN masks for pair classes (upper triangle applied in disagreement fn)
    epi_epi = epi_mask[:, None] & epi_mask[None, :]
    non_non = ~epi_mask[:, None] & ~epi_mask[None, :]
    cross = epi_mask[:, None] ^ epi_mask[None, :]
    same_probe = probes_arr[:, None] == probes_arr[None, :]
    cross_probe = ~same_probe
    all_pairs = np.ones_like(epi_epi, dtype=bool)

    pair_masks = {
        "all": all_pairs,
        "epi-epi": epi_epi,
        "non-non": non_non,
        "cross": cross,
        "epi-epi_CP": epi_epi & cross_probe,
        "non-non_CP": non_non & cross_probe,
    }

    for band in BANDS:
        linkages_by_phase = {}
        for phase in PHASES:
            lrg = load_lrg_result(pat, phase, band, "imcoh_abs")
            if lrg is None or lrg.n_nodes != N:
                continue
            linkages_by_phase[phase] = lrg.linkage_matrix

        h2a = compute_h2a_stratified(linkages_by_phase, pair_masks)
        for cls, curve in h2a.items():
            for ih, h in enumerate(H_GRID):
                if not np.isnan(curve[ih]):
                    h2a_rows.append(dict(
                        patient=pat, band=band, h=h, pair_class=cls,
                        h2a=curve[ih],
                    ))
    print(f"  {pat} done ({len(h2a_rows)} rows so far)")

h2a_df = pd.DataFrame(h2a_rows)
h2a_df.to_csv(OUT / "h2a_stratified.csv", index=False)
print(f"  → {OUT / 'h2a_stratified.csv'}")
print()

# ── Part 2: ablation subprocesses ─────────────────────────────────────
print("=" * 90)
print("PART 2: Ablation — zero out epi-epi edges and recompute LRG+H2a")
print("=" * 90)
print("Running per-patient subprocesses (memory safety).")
print()

worker = Path(__file__).parent / "_h2_ablation_one_patient.py"
for pat in PATIENTS:
    print(f"  Running ablation for {pat}...", flush=True)
    r = subprocess.run(
        [sys.executable, str(worker), pat],
        capture_output=True, text=True, timeout=1800,
    )
    if r.returncode != 0:
        print(f"    ERROR ({pat}): {r.stderr[-500:]}")
    else:
        print(f"  {r.stdout.strip()[-300:]}")
print()

# Merge ablation CSVs
ablation_rows = []
for pat in PATIENTS:
    p = OUT / f"_ablation_{pat}.csv"
    if p.exists():
        ablation_rows.append(pd.read_csv(p))
if ablation_rows:
    ab_df = pd.concat(ablation_rows, ignore_index=True)
    ab_df.to_csv(OUT / "h2a_ablation.csv", index=False)
    print(f"  → {OUT / 'h2a_ablation.csv'}")
else:
    ab_df = pd.DataFrame()
print()

# ── Part 3: figures + summary ─────────────────────────────────────────
print("=" * 90)
print("PART 3: Summary and figures")
print("=" * 90)
print()

# H2a curves per band: mean across patients for each pair class
CLASSES = ["non-non_CP", "epi-epi_CP", "cross", "all"]
CLASS_LABELS = {
    "non-non_CP": "non–non (cross-probe)",
    "epi-epi_CP": "epi–epi (cross-probe)",
    "cross": "epi–non (all)",
    "all": "all pairs",
}
CLASS_COLORS = {
    "non-non_CP": "#2ecc71",
    "epi-epi_CP": "#e74c3c",
    "cross": "#9b59b6",
    "all": "#7f8c8d",
}

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
for ax, band in zip(axes.flat, BANDS):
    for cls in CLASSES:
        sub = h2a_df[(h2a_df["band"] == band) & (h2a_df["pair_class"] == cls)]
        if len(sub) == 0:
            continue
        # Mean and SEM across patients at each h
        agg = sub.groupby("h")["h2a"].agg(["mean", "std", "count"]).reset_index()
        sem = agg["std"] / np.sqrt(agg["count"])
        ax.plot(agg["h"], agg["mean"], color=CLASS_COLORS[cls],
                label=CLASS_LABELS[cls], lw=1.5)
        ax.fill_between(agg["h"], agg["mean"] - sem, agg["mean"] + sem,
                        color=CLASS_COLORS[cls], alpha=0.2)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_title(f"{band}", fontsize=11)
    if band in ("low_gamma", "high_gamma", "beta"):
        ax.set_xlabel("height h")
    if band in ("delta", "beta"):
        ax.set_ylabel("H2a (disagreement)")

axes[0, 0].legend(fontsize=8, loc="upper right")
fig.suptitle("H2a = d(rsPre,rsPost) − d(taskTest,rsPost) — stratified by pair class\n"
             "Positive = rsPost more similar to taskTest than to rsPre on that pair class\n"
             "Mean ± SEM across 5 patients",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT / "h2a_stratified_by_band.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "h2a_stratified_by_band.png", bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"  → {OUT / 'h2a_stratified_by_band.pdf'}")

# Integrated H2a per band per class (area under positive part)
summary = (
    h2a_df.groupby(["band", "pair_class", "patient"])["h2a"]
    .apply(lambda x: np.trapz(np.clip(x, 0, None), H_GRID[:len(x)]))
    .reset_index(name="auc_pos")
)
summary_mean = summary.groupby(["band", "pair_class"])["auc_pos"].mean().reset_index()

print()
print("INTEGRATED H2a (positive area under curve) per band × pair class:")
print(f"{'band':<12s}  " + "  ".join([f"{CLASS_LABELS[c]:<22s}" for c in CLASSES]))
print("-" * (12 + 4 + 24 * 4))
for band in BANDS:
    row = [f"{band:<12s}"]
    for cls in CLASSES:
        val = summary_mean[
            (summary_mean["band"] == band) & (summary_mean["pair_class"] == cls)
        ]["auc_pos"]
        row.append(f"{val.iloc[0]:>10.4f}            " if len(val) > 0 else " " * 24)
    print("  ".join(row))
print()

# Beta detail: is epi-epi_CP H2a smaller than non-non_CP H2a?
beta_ee = summary[(summary["band"] == "beta") & (summary["pair_class"] == "epi-epi_CP")]["auc_pos"].values
beta_nn = summary[(summary["band"] == "beta") & (summary["pair_class"] == "non-non_CP")]["auc_pos"].values
print(f"Beta H2a_positive AUC:")
print(f"  non-non cross-probe: mean {beta_nn.mean():.4f}, per-pat = {beta_nn.round(4)}")
print(f"  epi-epi cross-probe: mean {beta_ee.mean():.4f}, per-pat = {beta_ee.round(4)}")
if len(beta_nn) == len(beta_ee) and len(beta_nn) >= 3:
    from scipy.stats import wilcoxon
    try:
        _, p = wilcoxon(beta_nn, beta_ee, alternative="greater")
        print(f"  Wilcoxon non-non > epi-epi: p = {p:.4f}")
    except Exception as e:
        print(f"  Wilcoxon failed: {e}")
print()

# Ablation summary
if len(ab_df) > 0:
    print("ABLATION: H2a_positive AUC, original vs ablated (epi-epi zeroed):")
    print(f"{'band':<12s}  {'orig mean':>10s}  {'ablated mean':>13s}  {'Δ%':>7s}")
    for band in BANDS:
        sub = ab_df[ab_df["band"] == band]
        if len(sub) == 0:
            continue
        o = sub["h2a_auc_original"].mean()
        a = sub["h2a_auc_ablated"].mean()
        dp = (a - o) / abs(o) * 100 if abs(o) > 1e-10 else np.nan
        print(f"  {band:<10s}  {o:>10.4f}  {a:>13.4f}  {dp:>+6.1f}%")

# Clean up per-patient ablation CSVs
for p in PATIENTS:
    (OUT / f"_ablation_{p}.csv").unlink(missing_ok=True)

print()
print("Done.")
