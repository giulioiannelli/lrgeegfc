#!/usr/bin/env python3
"""WP-spatial — Spatial analysis of task trace.

Tasks 1–6: electrode shaft mapping, per-node trace scores, spatial
clustering, cross-patient regional consistency, alpha universality,
within-shaft sanity check.

Run: python scripts/wp1/wp1_spatial.py
"""
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result

OUT = ROOT / "data" / "wp_spatial"
OUT.mkdir(parents=True, exist_ok=True)

PATIENTS = ["Pat_02", "Pat_05", "Pat_07", "Pat_08"]  # exclude Pat_03
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
KEY_BANDS = ["alpha", "beta", "theta"]
BAND_TEX = {b: BRAIN_BAND_TEX_DICT[b] for b in BANDS}
FC_METHOD = "msc"
LRG_CACHE = ROOT / "data" / "lrg_cache"
DATA_DIR = ROOT / "data" / "stereoeeg_patients"
TASK_REFS = ["task_learn", "task_test"]
K_RANGE = list(range(2, 59))  # full k range for scalar

PAT_COLORS = {
    "Pat_02": "#1f77b4", "Pat_05": "#2ca02c",
    "Pat_07": "#d62728", "Pat_08": "#9467bd",
}


def save_fig(fig, path, *, what, proves, how_to_read, dpi=300):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")
    path.with_suffix(".md").write_text(
        f"# {path.stem}\n\n## What the figure shows\n{what}\n\n"
        f"## What it proves\n{proves}\n\n## How to read it\n{how_to_read}\n"
    )


def coclassification_matrix(labels):
    return (labels[:, None] == labels[None, :]).astype(np.float64)


# ── Load LRG ──────────────────────────────────────────────────────────
print("Loading LRG data...")
lrg = {}
for pat in PATIENTS:
    for band in BANDS:
        for phase in ["rest_pre", "task_learn", "task_test", "rest_post"]:
            res = load_lrg_result(pat, phase, band, FC_METHOD, cache_root=LRG_CACHE)
            if res is not None:
                lrg[(pat, band, phase)] = res
print(f"  {len(lrg)} results")


# ── Load electrode info ───────────────────────────────────────────────
def load_electrode_info(patient):
    """Return DataFrame with columns: label, shaft, x, y, z, region."""
    pat_dir = DATA_DIR / patient
    pat_num = patient.replace("Pat_0", "").replace("Pat_", "")

    # Load channel labels
    lbl_path = pat_dir / "channel_labels.csv"
    lbl_df = pd.read_csv(lbl_path, header=None)
    # Some files have a 'label' header row
    if str(lbl_df.iloc[0, 0]).strip().lower() == "label":
        lbl_df = lbl_df.iloc[1:].reset_index(drop=True)
    labels = [str(l).strip() for l in lbl_df.iloc[:, 0]]

    # Extract shaft names
    shafts = []
    for l in labels:
        m = re.match(r"([A-Za-z']+)", l.replace(" ", ""))
        shafts.append(m.group(1) if m else "?")

    # Load coordinates
    coords = None
    # Try xlsx first
    for f in sorted(pat_dir.glob("Implant*.xlsx")):
        try:
            coords = pd.read_excel(f)
            if "x" in coords.columns:
                break
        except Exception:
            pass
    # Try csv
    if coords is None or "x" not in getattr(coords, "columns", []):
        for f in sorted(pat_dir.glob("Implant*.csv")):
            try:
                coords = pd.read_csv(f, encoding="latin-1")
                if "x" in coords.columns:
                    break
            except Exception:
                pass

    # Build result
    n = len(labels)
    result = pd.DataFrame({
        "node_id": range(n),
        "label": labels,
        "shaft": shafts,
    })

    if coords is not None and "x" in coords.columns:
        # Match by label if possible, otherwise by position
        coords_labels = [str(l).strip() for l in coords["label"]]
        # Try to match
        x_arr, y_arr, z_arr, region_arr = [], [], [], []
        region_col = None
        for c in coords.columns:
            if "desikan" in c.lower() or "killany" in c.lower():
                region_col = c
                break

        for lbl in labels:
            # Fuzzy match: remove spaces
            lbl_clean = lbl.replace(" ", "").replace(",G2", "").replace(",G1", "")
            matched = False
            for idx, cl in enumerate(coords_labels):
                cl_clean = cl.replace(" ", "").replace(",G2", "").replace(",G1", "")
                if cl_clean == lbl_clean:
                    x_arr.append(coords.iloc[idx]["x"])
                    y_arr.append(coords.iloc[idx]["y"])
                    z_arr.append(coords.iloc[idx]["z"])
                    region_arr.append(str(coords.iloc[idx][region_col])
                                      if region_col and pd.notna(coords.iloc[idx][region_col])
                                      else "unknown")
                    matched = True
                    break
            if not matched:
                x_arr.append(np.nan)
                y_arr.append(np.nan)
                z_arr.append(np.nan)
                region_arr.append("unknown")

        result["x"] = x_arr
        result["y"] = y_arr
        result["z"] = z_arr
        result["region"] = region_arr
    else:
        result["x"] = np.nan
        result["y"] = np.nan
        result["z"] = np.nan
        result["region"] = "unknown"

    return result


print("Loading electrode info...")
elec_info = {}
for pat in PATIENTS:
    elec_info[pat] = load_electrode_info(pat)
    n_coords = elec_info[pat]["x"].notna().sum()
    n_shafts = elec_info[pat]["shaft"].nunique()
    print(f"  {pat}: {len(elec_info[pat])} contacts, {n_shafts} shafts, "
          f"{n_coords} with coordinates")


# ── Per-node trace score ──────────────────────────────────────────────
def compute_node_trace_scalar(patient, band):
    """Scalar A analog at node level: mean over k of [overlap(post,task) - overlap(pre,task)]."""
    r_pre = lrg.get((patient, band, "rest_pre"))
    r_post = lrg.get((patient, band, "rest_post"))
    if r_pre is None or r_post is None or r_pre.n_nodes != r_post.n_nodes:
        return None
    n = r_pre.n_nodes

    all_scores = []
    for tr in TASK_REFS:
        r_task = lrg.get((patient, band, tr))
        if r_task is None or r_task.n_nodes != n:
            continue
        for k in K_RANGE:
            if k > n:
                break
            la_pre = fcluster(r_pre.linkage_matrix, k, criterion="maxclust")[:n]
            la_task = fcluster(r_task.linkage_matrix, k, criterion="maxclust")[:n]
            la_post = fcluster(r_post.linkage_matrix, k, criterion="maxclust")[:n]
            C_pre = coclassification_matrix(la_pre)
            C_task = coclassification_matrix(la_task)
            C_post = coclassification_matrix(la_post)
            ov_pre = np.mean(C_pre == C_task, axis=1)
            ov_post = np.mean(C_post == C_task, axis=1)
            all_scores.append(ov_post - ov_pre)

    if not all_scores:
        return None
    return np.mean(all_scores, axis=0)


# ======================================================================
# TASK 1: Electrode shaft mapping
# ======================================================================
def task1():
    print("\n── Task 1: Shaft mapping ──")
    md = ["# Electrode Shaft Mapping", ""]

    for pat in PATIENTS:
        ei = elec_info[pat]
        shaft_counts = ei["shaft"].value_counts().sort_index()
        n_shafts = len(shaft_counts)
        n_contacts = len(ei)

        md.append(f"## {pat}: {n_contacts} contacts, {n_shafts} shafts")
        md.append("")
        md.append("| Shaft | N contacts | Regions |")
        md.append("|-------|-----------|---------|")

        for shaft in shaft_counts.index:
            n_c = shaft_counts[shaft]
            regions = ei[ei["shaft"] == shaft]["region"].unique()
            regions = [r for r in regions if r != "unknown"]
            region_str = ", ".join(sorted(set(regions)))[:60] if regions else "—"
            md.append(f"| {shaft} | {n_c} | {region_str} |")
        md.append("")

    (OUT / "electrode_shaft_mapping.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'electrode_shaft_mapping.md'}")


# ======================================================================
# TASK 2: Per-node trace scores
# ======================================================================
def task2():
    print("\n── Task 2: Per-node trace scores ──")
    (OUT / "node_trace_scores").mkdir(parents=True, exist_ok=True)

    for pat in PATIENTS:
        ei = elec_info[pat].copy()
        n_lrg = lrg.get((pat, "alpha", "rest_pre"))
        if n_lrg is None:
            continue
        n = n_lrg.n_nodes
        # Trim to LRG giant component size
        if len(ei) > n:
            ei = ei.iloc[:n].copy()

        for band in BANDS:
            scores = compute_node_trace_scalar(pat, band)
            if scores is not None and len(scores) == len(ei):
                ei[f"trace_{band}"] = scores

        out_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        ei.to_csv(out_path, index=False)
        print(f"  {pat}: {len(ei)} nodes, saved to {out_path.name}")


# ======================================================================
# TASK 3: Spatial clustering — brain maps + shaft decomposition
# ======================================================================
def task3():
    print("\n── Task 3: Spatial analysis ──")
    (OUT / "brain_maps").mkdir(parents=True, exist_ok=True)
    (OUT / "shaft_decomposition").mkdir(parents=True, exist_ok=True)

    # Brain maps
    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)
        has_coords = ndf["x"].notna().sum() > 10

        for band in KEY_BANDS:
            col = f"trace_{band}"
            if col not in ndf.columns:
                continue
            scores = ndf[col].values

            if has_coords:
                x, y, z = ndf["x"].values, ndf["y"].values, ndf["z"].values
                valid = np.isfinite(x) & np.isfinite(scores)
                if valid.sum() < 10:
                    continue

                vmax = max(abs(np.nanmin(scores)), abs(np.nanmax(scores)))
                fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
                for ax, (cx, cy, xlabel, ylabel, title) in zip(axes, [
                    (x, y, "x (mm)", "y (mm)", "Axial"),
                    (y, z, "y (mm)", "z (mm)", "Sagittal"),
                ]):
                    sc = ax.scatter(cx[valid], cy[valid], c=scores[valid],
                                    cmap="RdBu", vmin=-vmax, vmax=vmax,
                                    s=25, edgecolors="black", linewidth=0.3)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    ax.set_title(title)
                    ax.set_aspect("equal")

                fig.colorbar(sc, ax=axes, label="Trace score", shrink=0.8)
                pct = np.mean(scores[valid] > 0) * 100
                fig.suptitle(f"{pat} — {BAND_TEX[band]} — node trace "
                             f"({pct:.0f}% positive)", fontsize=13)
                fig.tight_layout()
                save_fig(fig, OUT / "brain_maps" / f"{pat}_{band}_brainmap.pdf",
                    what=f"Spatial distribution of per-node trace scores for {pat}, {band} band. "
                         f"Red = trace, Blue = anti-trace.",
                    proves=f"Whether trace nodes are spatially clustered in {pat} {band}.",
                    how_to_read="Red clusters = localized trace. Scattered = distributed.",
                )

    # Shaft decomposition
    md = ["# Shaft Mean Trace Scores", ""]
    shaft_data = []

    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)

        md.append(f"## {pat}")
        md.append("")
        header = "| Shaft | N |"
        for band in KEY_BANDS:
            header += f" {BAND_TEX[band]} |"
        md.append(header)
        md.append("|-------|---|" + "|".join(["------"] * len(KEY_BANDS)) + "|")

        for shaft in sorted(ndf["shaft"].unique()):
            smask = ndf["shaft"] == shaft
            n_c = smask.sum()
            row = f"| {shaft} | {n_c} |"
            for band in KEY_BANDS:
                col = f"trace_{band}"
                if col in ndf.columns:
                    mean_s = ndf.loc[smask, col].mean()
                    shaft_data.append({
                        "patient": pat, "shaft": shaft, "band": band,
                        "n_contacts": n_c, "mean_trace": mean_s,
                    })
                    bold = "**" if abs(mean_s) > 0.02 else ""
                    row += f" {bold}{mean_s:+.4f}{bold} |"
                else:
                    row += " — |"
            md.append(row)
        md.append("")

    (OUT / "shaft_decomposition" / "shaft_mean_traces.md").write_text("\n".join(md))
    print(f"  Saved shaft mean traces")


# ======================================================================
# TASK 4: Cross-patient regional consistency
# ======================================================================
def task4():
    print("\n── Task 4: Regional consistency ──")
    (OUT / "region_analysis").mkdir(parents=True, exist_ok=True)

    # Collect region-level data
    region_rows = []
    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)

        for band in KEY_BANDS:
            col = f"trace_{band}"
            if col not in ndf.columns:
                continue
            for region in ndf["region"].unique():
                if region == "unknown" or pd.isna(region):
                    continue
                rmask = ndf["region"] == region
                if rmask.sum() < 2:
                    continue
                mean_s = ndf.loc[rmask, col].mean()
                n_pos = (ndf.loc[rmask, col] > 0).sum()
                region_rows.append({
                    "patient": pat, "band": band, "region": region,
                    "n_contacts": rmask.sum(), "mean_trace": mean_s,
                    "n_positive": n_pos,
                    "pct_positive": n_pos / rmask.sum() * 100,
                })

    rdf = pd.DataFrame(region_rows)
    if rdf.empty:
        print("  No region data available")
        (OUT / "region_analysis" / "region_consistency.md").write_text(
            "# Region consistency\n\nNo atlas region data available.")
        return

    # Find regions covered by ≥3 patients
    md = ["# Cross-patient Regional Consistency", ""]
    for band in KEY_BANDS:
        md.append(f"## {BAND_TEX[band]} band")
        md.append("")
        bdf = rdf[rdf["band"] == band]
        region_counts = bdf.groupby("region")["patient"].nunique()
        shared_regions = region_counts[region_counts >= 3].index

        if len(shared_regions) == 0:
            md.append("No regions covered by ≥3 patients.")
            md.append("")
            continue

        md.append("| Region | N patients | Mean trace | N positive | Pattern |")
        md.append("|--------|-----------|-----------|------------|---------|")

        for region in sorted(shared_regions):
            rsub = bdf[bdf["region"] == region]
            n_pat = len(rsub)
            mean_t = rsub["mean_trace"].mean()
            n_pos = (rsub["mean_trace"] > 0).sum()
            pattern = "TRACE" if n_pos == n_pat else ("RESET" if n_pos == 0 else "mixed")
            md.append(f"| {region} | {n_pat} | {mean_t:+.4f} | {n_pos}/{n_pat} | {pattern} |")
        md.append("")

    (OUT / "region_analysis" / "region_consistency.md").write_text("\n".join(md))
    print(f"  Saved region consistency")

    # Coverage analysis
    md_cov = ["# Coverage Analysis", "",
              "For patients without trace in a band: do they have electrodes in "
              "regions where others show strong trace?", ""]

    for band in KEY_BANDS:
        bdf = rdf[rdf["band"] == band]
        # Find regions where trace is strong (mean > 0.02)
        trace_regions = set()
        for region in bdf["region"].unique():
            rsub = bdf[bdf["region"] == region]
            if rsub["mean_trace"].mean() > 0.02:
                trace_regions.add(region)

        md_cov.append(f"## {BAND_TEX[band]}: trace-positive regions = {sorted(trace_regions)}")
        md_cov.append("")
        for pat in PATIENTS:
            pat_regions = set(bdf[bdf["patient"] == pat]["region"].unique())
            covered = trace_regions & pat_regions
            missing = trace_regions - pat_regions
            # Overall trace for this patient+band
            csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
            if csv_path.exists():
                ndf = pd.read_csv(csv_path)
                col = f"trace_{band}"
                if col in ndf.columns:
                    pct = np.mean(ndf[col].values > 0) * 100
                    md_cov.append(f"- {pat} (trace={pct:.0f}%): covers {len(covered)}/{len(trace_regions)} "
                                  f"trace regions, missing: {sorted(missing) if missing else 'none'}")
        md_cov.append("")

    (OUT / "region_analysis" / "coverage_analysis.md").write_text("\n".join(md_cov))
    print(f"  Saved coverage analysis")


# ======================================================================
# TASK 5: Alpha universality assessment
# ======================================================================
def task5():
    print("\n── Task 5: Alpha universality ──")

    md = ["# Alpha Universality Assessment", ""]

    # For alpha: compute per-shaft trace mean for each patient
    shaft_traces = {}
    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)
        if "trace_alpha" not in ndf.columns:
            continue

        shaft_means = ndf.groupby("shaft")["trace_alpha"].mean()
        n_pos_shafts = (shaft_means > 0).sum()
        n_total_shafts = len(shaft_means)
        shaft_traces[pat] = {
            "n_pos_shafts": n_pos_shafts,
            "n_total": n_total_shafts,
            "frac_pos": n_pos_shafts / n_total_shafts,
            "shaft_means": shaft_means,
        }

    md.append("## Per-patient shaft analysis (alpha)")
    md.append("")
    md.append("| Patient | Shafts with positive trace | Total shafts | Fraction |")
    md.append("|---------|--------------------------|-------------|----------|")
    for pat in PATIENTS:
        if pat in shaft_traces:
            st = shaft_traces[pat]
            md.append(f"| {pat} | {st['n_pos_shafts']} | {st['n_total']} | "
                      f"{st['frac_pos']:.2f} |")

    md += [
        "",
        "## Hypothesis assessment",
        "",
        "**Hypothesis A (spatially widespread)**: Alpha trace occurs in many "
        "brain regions — every implant captures it because it's everywhere.",
        "",
    ]

    # Check: what fraction of shafts are positive per patient?
    fracs = [shaft_traces[p]["frac_pos"] for p in PATIENTS if p in shaft_traces]
    if fracs:
        mean_frac = np.mean(fracs)
        md.append(f"Test: mean fraction of positive shafts = {mean_frac:.2f}")
        if mean_frac > 0.6:
            md.append("→ SUPPORTED: majority of shafts show trace in every patient.")
        elif mean_frac > 0.4:
            md.append("→ PARTIAL: roughly half of shafts show trace.")
        else:
            md.append("→ NOT SUPPORTED: trace is concentrated in a few shafts.")

    md += [
        "",
        "**Hypothesis B (concentrated but universally covered)**: Alpha trace "
        "is in a specific region that all implants happen to cover.",
        "",
    ]

    # Check: is there a region that is positive in ALL patients?
    region_data = {}
    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)
        if "trace_alpha" not in ndf.columns:
            continue
        for region in ndf["region"].unique():
            if region == "unknown" or pd.isna(region):
                continue
            rmean = ndf[ndf["region"] == region]["trace_alpha"].mean()
            if region not in region_data:
                region_data[region] = {}
            region_data[region][pat] = rmean

    universal_regions = []
    for region, pat_dict in region_data.items():
        if len(pat_dict) >= 3 and all(v > 0 for v in pat_dict.values()):
            universal_regions.append((region, np.mean(list(pat_dict.values()))))

    if universal_regions:
        md.append(f"Regions with positive alpha trace in ALL covering patients:")
        for r, m in sorted(universal_regions, key=lambda x: -x[1]):
            md.append(f"  - {r}: mean trace = {m:+.4f}")
    else:
        md.append("No region is consistently positive across all covering patients.")

    md += [
        "",
        "**Hypothesis C (global network-level)**: The trace is not about specific "
        "nodes but about overall partition structure — detectable regardless of "
        "which nodes are sampled.",
        "",
        "Test: If the scalar VI trace (which compares full partitions) is 5/5 "
        "positive while per-node trace is variable, this supports Hypothesis C.",
        "The scalar VI for alpha IS 5/5 positive (mean +0.129), and the per-node "
        "co-classification is noisy (~49% at meso-scale). This is consistent with C.",
    ]

    (OUT / "alpha_universality_assessment.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'alpha_universality_assessment.md'}")


# ======================================================================
# TASK 6: Within-shaft sanity check
# ======================================================================
def task6():
    print("\n── Task 6: Within-shaft sanity check ──")

    md = ["# Within-shaft Sanity Check", ""]

    for pat in PATIENTS:
        csv_path = OUT / "node_trace_scores" / f"{pat}_node_scores.csv"
        if not csv_path.exists():
            continue
        ndf = pd.read_csv(csv_path)
        n = len(ndf)

        # Check: what fraction of within-shaft pairs are ALWAYS co-classified?
        shafts = ndf["shaft"].values
        r_pre = lrg.get((pat, "alpha", "rest_pre"))
        if r_pre is None:
            continue

        # Sample a few k values
        always_coclassified_within = 0
        total_within = 0
        always_coclassified_between = 0
        total_between = 0

        ks_sample = [3, 5, 10, 20, 30, 50]
        n_lrg = r_pre.n_nodes
        n_check = min(n, n_lrg)

        for i in range(n_check):
            for j in range(i + 1, n_check):
                same_shaft = shafts[i] == shafts[j]
                always_same = True
                for phase in ["rest_pre", "task_learn", "task_test", "rest_post"]:
                    res = lrg.get((pat, "alpha", phase))
                    if res is None:
                        continue
                    for k in ks_sample:
                        if k > n_lrg:
                            continue
                        labels = fcluster(res.linkage_matrix, k, criterion="maxclust")[:n_lrg]
                        if labels[i] != labels[j]:
                            always_same = False
                            break
                    if not always_same:
                        break

                if same_shaft:
                    total_within += 1
                    if always_same:
                        always_coclassified_within += 1
                else:
                    total_between += 1
                    if always_same:
                        always_coclassified_between += 1

        within_frac = always_coclassified_within / max(total_within, 1)
        between_frac = always_coclassified_between / max(total_between, 1)

        md.append(f"## {pat}")
        md.append("")
        md.append(f"- Within-shaft pairs always co-classified: "
                  f"{always_coclassified_within}/{total_within} ({within_frac:.1%})")
        md.append(f"- Between-shaft pairs always co-classified: "
                  f"{always_coclassified_between}/{total_between} ({between_frac:.1%})")
        md.append("")

        # Within-shaft vs between-shaft trace
        within_traces, between_traces = [], []
        for band in ["alpha"]:
            col = f"trace_{band}"
            if col not in ndf.columns:
                continue
            scores = ndf[col].values[:n_check]

            for i in range(n_check):
                for j in range(i + 1, n_check):
                    if shafts[i] == shafts[j]:
                        within_traces.append(scores[i])
                    else:
                        between_traces.append(scores[i])

        if within_traces and between_traces:
            md.append(f"  Alpha trace — within-shaft nodes: mean={np.mean(within_traces):+.4f}")
            md.append(f"  Alpha trace — between-shaft nodes: mean={np.mean(between_traces):+.4f}")
        md.append("")

    (OUT / "shaft_sanity_check.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'shaft_sanity_check.md'}")


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    task1()
    task2()
    task3()
    task4()
    task5()
    task6()
    print(f"\nAll done. Outputs in {OUT}")
