#!/usr/bin/env python3
"""Audit Step 11 — Consensus-subtree CBR (variant D).

Per scope: .agents/guides/task-persistence-investigation/2026-04-25_consensus-subtree.md

Strictest CBR variant. A leafset L is CON-trace iff:
  L ∈ F_taskT_w  AND  L ∈ F_rsPost_w  AND  L ∉ F_rsPre,
where F_φ_w is the size-windowed leafset family of phase φ. Exact
set equality, no thresholds.

Outputs:
  data/audit/per_patient_hierarchy_consensus/{Pat_XX}/{band}_hierarchy-crossphase.pdf
  data/audit/per_patient_hierarchy_consensus/leaf_assignment.csv
  data/audit/per_patient_hierarchy_consensus/consensus_leafsets.csv
  data/audit/per_patient_hierarchy_consensus/CENSUS.md
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, PATTERN_PRIORITY,
    SIZE_MIN, SIZE_MAX,
    nodes_with_leaves, _load_Z,
    plot_dendrogram_with_strip,
)


OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_consensus"


def leafset_family(Z: np.ndarray, size_min: int | None = None,
                    size_max: int | None = None
                    ) -> dict[frozenset[int], dict]:
    """Return {leafset -> {row, h, h_rel, size}} for every internal node."""
    out: dict[frozenset[int], dict] = {}
    if Z is None:
        return out
    dmax = float(Z[-1, 2])
    for v in nodes_with_leaves(Z):
        s = v["size"]
        if size_min is not None and s < size_min:
            continue
        if size_max is not None and s > size_max:
            continue
        out[v["leaves"]] = {
            "row": v["row"], "h": v["h"],
            "h_rel": v["h"] / dmax if dmax > 0 else 0.0,
            "size": s,
        }
    return out


def classify_CON(in_taskT_w: bool, in_rsPost_w: bool,
                  in_rsPre: bool) -> str | None:
    if in_taskT_w and in_rsPost_w and not in_rsPre:
        return "trace"
    if in_taskT_w and in_rsPost_w and in_rsPre:
        return "persist"
    if in_rsPre and not (in_taskT_w or in_rsPost_w):
        return "reset"
    if in_rsPost_w and not in_taskT_w and not in_rsPre:
        return "rearrange"
    return None


def consensus_for_patient_band(pat: str, band: str
                                  ) -> tuple[np.ndarray | None, list[dict],
                                               list[dict]]:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    if Zs["rest_post"] is None or Zs["rest_pre"] is None or Zs["task_learn"] is None:
        return None, [], []
    if Zs["task_test"] is None:
        return None, [], []

    F_rpre = leafset_family(Zs["rest_pre"])  # full, no size window
    F_tt_w = leafset_family(Zs["task_test"], SIZE_MIN, SIZE_MAX)
    F_rpost_w = leafset_family(Zs["rest_post"], SIZE_MIN, SIZE_MAX)

    n = Zs["rest_post"].shape[0] + 1
    leaf_patterns: list[set[str]] = [set() for _ in range(n)]
    consensus_rows: list[dict] = []

    # Iterate every distinct windowed leafset (rsPost ∪ taskT)
    candidate_leafsets = set(F_tt_w.keys()) | set(F_rpost_w.keys())
    for L in candidate_leafsets:
        in_tt = L in F_tt_w
        in_rp = L in F_rpost_w
        in_pre = L in F_rpre
        label = classify_CON(in_tt, in_rp, in_pre)
        if label is None:
            continue
        consensus_rows.append({
            "patient": pat, "band": band,
            "leafset_size": len(L),
            "leafset_hash": hex(hash(L) & 0xFFFFFFFFFFFFFFFF)[2:].rjust(16, "0"),
            "h_rel_in_taskT": F_tt_w[L]["h_rel"] if in_tt else float("nan"),
            "h_rel_in_rsPost": F_rpost_w[L]["h_rel"] if in_rp else float("nan"),
            "h_rel_in_rsPre": F_rpre[L]["h_rel"] if in_pre else float("nan"),
            "label": label,
        })
        # Project to leaves only via rsPost-supported labels
        if in_rp:
            for lid in L:
                leaf_patterns[lid].add(label)

    out = np.empty(n, dtype=object)
    for lid in range(n):
        chosen = "rearrange"
        for p in PATTERN_PRIORITY:
            if p in leaf_patterns[lid]:
                chosen = p
                break
        out[lid] = chosen

    leaf_rows = [{"patient": pat, "band": band, "leaf_id": int(lid),
                   "category": cat} for lid, cat in enumerate(out)]
    return out, leaf_rows, consensus_rows


def figure_for(pat: str, band: str, cats: np.ndarray) -> Path:
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
    n_phases = len(PHASE_ORDER)
    fig = plt.figure(figsize=(4.5 * n_phases, 5.4))
    outer = fig.add_gridspec(1, n_phases, wspace=0.18, top=0.84, bottom=0.05)
    for ip, phase in enumerate(PHASE_ORDER):
        cell = outer[0, ip].subgridspec(2, 1, height_ratios=[8, 1], hspace=0.04)
        ax_dend = fig.add_subplot(cell[0])
        ax_strip = fig.add_subplot(cell[1])
        Z = Zs[phase]
        if Z is None:
            ax_dend.text(0.5, 0.5, "missing", ha="center", va="center",
                          transform=ax_dend.transAxes, color="#999")
            ax_dend.set_xticks([]); ax_dend.set_yticks([])
            ax_strip.axis("off")
            ax_dend.set_title(PHASE_TEX[phase], fontsize=11)
        else:
            plot_dendrogram_with_strip(ax_dend, ax_strip, Z, cats,
                                         title=PHASE_TEX[phase])
        if ip == 0:
            ax_dend.set_ylabel(BRAIN_BAND_TEX_DICT.get(band, band), fontsize=12)

    handles = [Patch(color=PATTERN_COLOR[p], label=p.upper())
                for p in PATTERN_PRIORITY]
    fig.legend(handles=handles, loc="upper center", ncol=4,
                bbox_to_anchor=(0.5, 0.98), fontsize=10, frameon=False,
                title=f"{pat} — {BRAIN_BAND_TEX_DICT.get(band, band)} — "
                      f"Consensus-subtree (exact equality, size ∈ [{SIZE_MIN}, {SIZE_MAX}])",
                title_fontsize=11)

    out_dir = OUT_DIR / pat
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_hierarchy-crossphase.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def write_census_md(consensus: list[dict]) -> None:
    """Markdown census: one row per (patient, band)."""
    by_pb: dict[tuple[str, str], dict[str, int]] = {}
    for r in consensus:
        key = (r["patient"], r["band"])
        by_pb.setdefault(key, {p: 0 for p in PATTERN_PRIORITY})
        by_pb[key][r["label"]] += 1

    lines = ["# Consensus-subtree census (CON, exact-equality CBR)\n",
              f"Size window: [{SIZE_MIN}, {SIZE_MAX}].",
              f"Generated by `scripts/01_compute/audit/audit_11_consensus_subtree.py`.\n",
              "| patient | band | trace | persist | reset | rearrange |",
              "|:--------|:-----|------:|--------:|------:|----------:|"]
    totals = {p: 0 for p in PATTERN_PRIORITY}
    for pat in PATIENTS_4PHASE:
        for band in BAND_ORDER:
            t = by_pb.get((pat, band))
            if t is None:
                continue
            lines.append(f"| {pat} | {band} | {t['trace']} | {t['persist']} | "
                          f"{t['reset']} | {t['rearrange']} |")
            for p in PATTERN_PRIORITY:
                totals[p] += t[p]
    lines.append(f"| **TOTAL** | | **{totals['trace']}** | "
                  f"**{totals['persist']}** | **{totals['reset']}** | "
                  f"**{totals['rearrange']}** |")
    (OUT_DIR / "CENSUS.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_leaves: list[dict] = []
    all_consensus: list[dict] = []
    skipped: list[str] = []
    n_written = 0

    for pat in PATIENTS_4PHASE:
        print(f"{pat} ...", flush=True)
        for band in BAND_ORDER:
            cats, leaves, consensus = consensus_for_patient_band(pat, band)
            if cats is None:
                skipped.append(f"{pat}/{band}")
                continue
            all_leaves.extend(leaves)
            all_consensus.extend(consensus)
            out = figure_for(pat, band, cats)
            n_written += 1
            print(f"  wrote {out.relative_to(ROOT)}")

    leaf_csv = OUT_DIR / "leaf_assignment.csv"
    with leaf_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient", "band", "leaf_id", "category"])
        w.writeheader(); w.writerows(all_leaves)

    cons_csv = OUT_DIR / "consensus_leafsets.csv"
    with cons_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "patient", "band", "leafset_size", "leafset_hash",
            "h_rel_in_taskT", "h_rel_in_rsPost", "h_rel_in_rsPre", "label",
        ])
        w.writeheader(); w.writerows(all_consensus)

    write_census_md(all_consensus)

    if skipped:
        (OUT_DIR / "missing.txt").write_text("\n".join(skipped) + "\n")

    n_trace = sum(1 for r in all_consensus if r["label"] == "trace")
    print(f"\n{n_written} figures, {len(skipped)} skipped")
    print(f"Total CON-trace consensus leafsets across cohort: {n_trace}")


if __name__ == "__main__":
    main()
