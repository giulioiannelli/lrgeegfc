#!/usr/bin/env python3
"""WP0 Task 4 — Multiscale exploration (h-multiscale and tau-multiscale).

Task 4a: Vary dendrogram cut level h at fixed tau = tau' = 1/lambda_max.
Task 4b: Vary diffusion time tau with optimal h(tau) at each tau.
Task 4c: Added-value assessment comparing approaches.

Produces (under data/wp0_metric_exploration/task4_multiscale/):
  h_multiscale/
    ari_profiles_all_patients.pdf + .md
    ari_matrices_examples.pdf + .md
  tau_multiscale/
    tau_profiles_all_patients.pdf + .md
    tau_sensitive_bands.pdf + .md
  added_value_assessment.md

Run: python scripts/wp0/task4_multiscale.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "wp0"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

from _common import (
    ALL_PAIRS,
    BANDS,
    FC_METHOD,
    LRG_CACHE,
    MSC_CACHE,
    OUT_ROOT,
    PATIENTS,
    PHASES,
    classify_pair,
    load_all_lrg,
    load_all_fc,
    save_fig,
)
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, DEFAULT_NPERSEG

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER]
OUT_H = OUT_ROOT / "task4_multiscale" / "h_multiscale"
OUT_TAU = OUT_ROOT / "task4_multiscale" / "tau_multiscale"

PAIR_SHORT = {
    (a, b): f"{a[:2]}{a[-4:]}-{b[:2]}{b[-4:]}" for a, b in ALL_PAIRS
}

# Phase pair type colors
PAIR_COLORS = {
    "within": "#2196F3",
    "cross": "#F44336",
}


# ══════════════════════════════════════════════════════════════════════════
# Task 4a — h-multiscale: ARI profiles across dendrogram cut levels
# ══════════════════════════════════════════════════════════════════════════
def task_4a(lrg_data: dict):
    print("\n── Task 4a: h-multiscale ──")
    OUT_H.mkdir(parents=True, exist_ok=True)

    # For each (patient, band, pair): compute ARI at n_clusters = 2, 3, ..., max_n
    N_MAX_FRAC = 0.5  # go up to N/2 clusters
    all_profiles = []  # List of dicts for later analysis

    for pat in PATIENTS:
        for band in BAND_ORDER:
            # Get all phase results for this (patient, band)
            results = {}
            for phase in PHASES:
                key = (pat, phase, band)
                if key in lrg_data:
                    results[phase] = lrg_data[key]

            if len(results) < 2:
                continue

            # Determine max n_clusters (limited by smallest giant component)
            n_nodes_list = [r.n_nodes for r in results.values()]
            max_n = int(min(n_nodes_list) * N_MAX_FRAC)
            max_n = max(max_n, 3)
            n_range = list(range(2, max_n + 1))

            for pa, pb in ALL_PAIRS:
                if pa not in results or pb not in results:
                    continue
                res_a, res_b = results[pa], results[pb]

                # Different n_nodes? Skip (can't compare partitions)
                if res_a.n_nodes != res_b.n_nodes:
                    continue

                Za, Zb = res_a.linkage_matrix, res_b.linkage_matrix
                ari_profile = []
                for n in n_range:
                    la = fcluster(Za, n, criterion="maxclust")
                    lb = fcluster(Zb, n, criterion="maxclust")
                    ari_profile.append(adjusted_rand_score(la, lb))

                all_profiles.append({
                    "patient": pat,
                    "band": band,
                    "phase_a": pa,
                    "phase_b": pb,
                    "pair_type": classify_pair(pa, pb),
                    "n_range": n_range,
                    "ari_profile": ari_profile,
                })

    print(f"  Computed {len(all_profiles)} ARI profiles")

    # -- Figure 1: ARI profiles faceted by band --
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for idx, band in enumerate(BAND_ORDER):
        ax = axes[idx]
        band_profiles = [p for p in all_profiles if p["band"] == band]

        for prof in band_profiles:
            pair_type = prof["pair_type"]
            color = PAIR_COLORS[pair_type]
            alpha = 0.3 if pair_type == "cross" else 0.6
            lw = 1.5 if pair_type == "within" else 0.8
            ax.plot(prof["n_range"], prof["ari_profile"],
                    color=color, alpha=alpha, lw=lw)

        ax.set_xlabel("n clusters", fontsize=9)
        ax.set_ylabel("ARI", fontsize=9)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(0, color="gray", ls=":", alpha=0.3)

    # Legend
    from matplotlib.lines import Line2D
    axes[0].legend(handles=[
        Line2D([0], [0], color=PAIR_COLORS["within"], lw=2, label="Within (rest-rest, task-task)"),
        Line2D([0], [0], color=PAIR_COLORS["cross"], lw=1, label="Cross (rest-task)"),
    ], fontsize=8, loc="upper right")

    fig.suptitle("Task 4a: ARI profiles across dendrogram cut levels (h-multiscale)", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT_H / "ari_profiles_all_patients.pdf",
        what="ARI(n) curves for all (patient, band, phase-pair) combinations, faceted by frequency band. "
             "X-axis = number of clusters (dendrogram cut level h). Y-axis = Adjusted Rand Index. "
             "Blue = within-type pairs (rest-rest or task-task), Red = cross-type pairs (rest-task). "
             "Each line = one (patient, pair) combination.",
        proves="Whether ARI reveals structure beyond the single Ψ-optimal cut. "
               "If within-type ARI stays high across many n while cross-type ARI drops, "
               "the reorganization is robust across hierarchical levels. "
               "If curves converge at certain n, that scale is not informative.",
        how_to_read="Look for SEPARATION between blue and red curves: larger gaps = "
                    "the metric discriminates within vs cross at that hierarchical level. "
                    "If separation only exists at specific n ranges, those are the informative scales.",
    )

    # -- Figure 2: ARI cross-resolution matrices for examples --
    # Select 2 patients × 2 bands × 1 cross pair as examples
    example_configs = [
        (PATIENTS[0], "theta", "rest_pre", "task_learn"),
        (PATIENTS[0], "alpha", "rest_pre", "task_learn"),
        (PATIENTS[1], "theta", "rest_pre", "task_learn"),
        (PATIENTS[1], "alpha", "rest_pre", "task_learn"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    for ax, (pat, band, pa, pb) in zip(axes, example_configs):
        key_a = (pat, pa, band)
        key_b = (pat, pb, band)
        if key_a not in lrg_data or key_b not in lrg_data:
            ax.set_visible(False)
            continue

        res_a, res_b = lrg_data[key_a], lrg_data[key_b]
        if res_a.n_nodes != res_b.n_nodes:
            ax.set_visible(False)
            continue

        Za, Zb = res_a.linkage_matrix, res_b.linkage_matrix
        max_n = min(int(res_a.n_nodes * N_MAX_FRAC), 30)
        n_range = list(range(2, max_n + 1))
        nn = len(n_range)

        ari_matrix = np.zeros((nn, nn))
        for i, ni in enumerate(n_range):
            la = fcluster(Za, ni, criterion="maxclust")
            for j, nj in enumerate(n_range):
                lb = fcluster(Zb, nj, criterion="maxclust")
                ari_matrix[i, j] = adjusted_rand_score(la, lb)

        im = ax.imshow(ari_matrix, cmap="viridis", aspect="auto", vmin=-0.05, vmax=1)
        # Sparse tick labels
        tick_step = max(1, nn // 6)
        ticks = list(range(0, nn, tick_step))
        ax.set_xticks(ticks)
        ax.set_xticklabels([n_range[t] for t in ticks], fontsize=7)
        ax.set_yticks(ticks)
        ax.set_yticklabels([n_range[t] for t in ticks], fontsize=7)
        ax.set_xlabel(f"n clusters ({pb})", fontsize=8)
        ax.set_ylabel(f"n clusters ({pa})", fontsize=8)
        ax.set_title(f"{pat} {BRAIN_BAND_TEX_DICT[band]}", fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.7)

    fig.suptitle("Task 4a: Cross-resolution ARI matrices (examples)", fontsize=13, y=1.02)
    fig.tight_layout()
    save_fig(fig, OUT_H / "ari_matrices_examples.pdf",
        what="ARI(n, n') cross-resolution matrices for 4 example (patient, band) configs. "
             "Rows = cut levels for phase φ (rest_pre), Columns = cut levels for phase φ' (task_learn). "
             "Color = ARI between the two partitions at those cut levels.",
        proves="Whether phases agree at different hierarchical levels (off-diagonal structure). "
               "Diagonal = same-resolution comparison. Off-diagonal = phase φ at coarser resolution "
               "matches phase φ' at finer resolution.",
        how_to_read="Bright diagonal = phases agree when cut at the same level. "
                    "Off-diagonal bright regions = phases match at SHIFTED hierarchical levels. "
                    "Uniform dark = no agreement at any level combination.",
    )

    return all_profiles


# ══════════════════════════════════════════════════════════════════════════
# Task 4b — tau-multiscale: diffusion time sweep
# ══════════════════════════════════════════════════════════════════════════
def task_4b(fc_data: dict, lrg_data: dict):
    print("\n── Task 4b: tau-multiscale ──")
    OUT_TAU.mkdir(parents=True, exist_ok=True)

    import networkx as nx
    from lrgsglib.utils.lrg.infocomm import lapl_dists, extract_ultrametric_matrix
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    from lrgsglib.core import get_giant_component
    from lrgsglib.utils.lrg.clustering import compute_optimal_threshold
    from scipy.cluster.hierarchy import linkage as scipy_linkage

    N_TAU = 20  # number of tau values to sweep

    # Cache file for tau sweep results
    tau_cache_dir = OUT_ROOT / "task4_multiscale" / "tau_cache"
    tau_cache_dir.mkdir(parents=True, exist_ok=True)
    tau_csv = tau_cache_dir / "tau_sweep_results.csv"

    if tau_csv.exists():
        print(f"  Loading cached tau sweep from {tau_csv}")
        tau_df = pd.read_csv(tau_csv)
    else:
        print("  Computing tau sweep (this may take a moment)...")
        t0 = time.time()

        # Step 1: For each (patient, phase, band), compute graph properties
        graph_data = {}  # (pat, phase, band) -> (L, w, n_nodes, giant)
        for pat in PATIENTS:
            for band in BAND_ORDER:
                for phase in PHASES:
                    key = (pat, phase, band)
                    W = fc_data.get(key)
                    if W is None:
                        continue
                    graph = nx.from_numpy_array(W)
                    giant = get_giant_component(graph)
                    L, w = get_graph_lspectrum(giant)
                    graph_data[key] = (L, w, giant.number_of_nodes())

        print(f"    Graphs prepared: {len(graph_data)}")

        # Step 2: For each (patient, band), determine shared tau range and sweep
        tau_rows = []
        for pat in PATIENTS:
            for band in BAND_ORDER:
                # Collect tau' for all phases
                tau_primes = []
                for phase in PHASES:
                    key = (pat, phase, band)
                    if key not in graph_data:
                        continue
                    L, w, n = graph_data[key]
                    w_pos = w[w > 1e-10]
                    if len(w_pos) == 0:
                        continue
                    tau_primes.append(1.0 / np.max(w_pos))

                if len(tau_primes) < 2:
                    continue

                # Use geometric mean of tau' values, sweep 2 decades around it
                tau_center = np.exp(np.mean(np.log(tau_primes)))
                tau_grid = np.logspace(
                    np.log10(tau_center / 10),
                    np.log10(tau_center * 100),
                    N_TAU,
                )

                # Compute ultrametric at each (phase, tau)
                ultra_by_phase_tau = {}  # (phase, tau_idx) -> (U_condensed, Z, labels)
                for phase in PHASES:
                    key = (pat, phase, band)
                    if key not in graph_data:
                        continue
                    L, w, n = graph_data[key]

                    for ti, tau in enumerate(tau_grid):
                        try:
                            dists = lapl_dists(L, tau=tau)
                            Z = scipy_linkage(dists, method="average")
                            threshold, *_ = compute_optimal_threshold(Z)
                            labels = fcluster(Z, threshold, criterion="distance")
                            ultra_by_phase_tau[(phase, ti)] = (dists, Z, labels)
                        except Exception:
                            pass

                # Compare phase pairs at each tau
                for pa, pb in ALL_PAIRS:
                    for ti, tau in enumerate(tau_grid):
                        ka = (pa, ti)
                        kb = (pb, ti)
                        if ka not in ultra_by_phase_tau or kb not in ultra_by_phase_tau:
                            continue

                        dists_a, Za, la = ultra_by_phase_tau[ka]
                        dists_b, Zb, lb = ultra_by_phase_tau[kb]

                        # ARI at optimal partition
                        n_a = graph_data[(pat, pa, band)][2]
                        n_b = graph_data[(pat, pb, band)][2]
                        if n_a == n_b:
                            ari = adjusted_rand_score(la, lb)
                        else:
                            ari = np.nan

                        # Scaled Frobenius of ultrametric vectors (if same size)
                        if len(dists_a) == len(dists_b):
                            log_a = np.log(np.clip(dists_a, 1e-15, None))
                            log_b = np.log(np.clip(dists_b, 1e-15, None))
                            na_norm = np.linalg.norm(log_a)
                            nb_norm = np.linalg.norm(log_b)
                            if na_norm > 1e-15 and nb_norm > 1e-15:
                                scaled_dist = np.linalg.norm(
                                    log_a / na_norm - log_b / nb_norm
                                )
                            else:
                                scaled_dist = np.nan
                        else:
                            scaled_dist = np.nan

                        tau_rows.append({
                            "patient": pat,
                            "band": band,
                            "phase_a": pa,
                            "phase_b": pb,
                            "pair_type": classify_pair(pa, pb),
                            "tau_idx": ti,
                            "tau": tau,
                            "log10_tau": np.log10(tau),
                            "tau_over_tau_prime": tau / tau_center,
                            "ari": ari,
                            "scaled_distance": scaled_dist,
                        })

                print(f"    {pat} {band}: {len(tau_grid)} tau values, "
                      f"tau' ~ {tau_center:.2e}")

        tau_df = pd.DataFrame(tau_rows)
        tau_df.to_csv(tau_csv, index=False)
        elapsed = time.time() - t0
        print(f"  Tau sweep complete: {len(tau_df)} rows, {elapsed:.1f}s")
        print(f"  Cached to {tau_csv}")

    # -- Figure 1: distance and ARI profiles vs log10(tau) --
    fig, axes = plt.subplots(2, 6, figsize=(24, 8), sharex=True)

    for bi, band in enumerate(BAND_ORDER):
        ax_dist = axes[0, bi]
        ax_ari = axes[1, bi]
        sub = tau_df[tau_df["band"] == band]

        for (pat, pa, pb), grp in sub.groupby(["patient", "phase_a", "phase_b"]):
            grp = grp.sort_values("log10_tau")
            ptype = classify_pair(pa, pb)
            color = PAIR_COLORS[ptype]
            alpha = 0.3 if ptype == "cross" else 0.6
            lw = 1.0

            ax_dist.plot(grp["log10_tau"], grp["scaled_distance"],
                         color=color, alpha=alpha, lw=lw)
            ax_ari.plot(grp["log10_tau"], grp["ari"],
                        color=color, alpha=alpha, lw=lw)

        ax_dist.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax_dist.set_ylabel("Scaled distance" if bi == 0 else "")
        ax_ari.set_ylabel("ARI" if bi == 0 else "")
        ax_ari.set_xlabel("log₁₀(τ)")
        ax_ari.set_ylim(-0.05, 1.05)

    from matplotlib.lines import Line2D
    axes[0, 0].legend(handles=[
        Line2D([0], [0], color=PAIR_COLORS["within"], lw=2, label="Within"),
        Line2D([0], [0], color=PAIR_COLORS["cross"], lw=1, label="Cross"),
    ], fontsize=8, loc="upper right")

    fig.suptitle("Task 4b: τ-multiscale profiles — distance and ARI vs diffusion time", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT_TAU / "tau_profiles_all_patients.pdf",
        what="Two rows of 6 panels (one per band). Top row: scaled ultrametric distance vs log₁₀(τ). "
             "Bottom row: ARI at Ψ-optimal cut vs log₁₀(τ). Blue = within-type pairs, Red = cross-type. "
             "Each line = one (patient, pair).",
        proves="Whether reorganization is concentrated at specific τ or spread across diffusion scales. "
               "If curves separate (within vs cross) only at certain τ, that scale carries the reorganization signal.",
        how_to_read="Top: within-type pairs should have LOWER distance (more similar) than cross-type. "
                    "Bottom: within-type should have HIGHER ARI. "
                    "Look for τ values where the gap is largest — those are the most informative scales. "
                    "If gap is constant across τ, the single τ' is sufficient.",
    )

    # -- Figure 2: tau sensitivity by band --
    # For each band: compute the gap (mean_cross - mean_within) as function of tau
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, metric_col, ylabel in [
        (axes[0], "scaled_distance", "Gap: mean(cross) − mean(within) [distance]"),
        (axes[1], "ari", "Gap: mean(within) − mean(cross) [ARI]"),
    ]:
        for band in BAND_ORDER:
            sub = tau_df[tau_df["band"] == band].dropna(subset=[metric_col])
            if sub.empty:
                continue

            gap_by_tau = []
            for ti in sorted(sub["tau_idx"].unique()):
                tsub = sub[sub["tau_idx"] == ti]
                within = tsub[tsub["pair_type"] == "within"][metric_col].mean()
                cross = tsub[tsub["pair_type"] == "cross"][metric_col].mean()

                if metric_col == "scaled_distance":
                    gap = cross - within  # cross should be larger (more different)
                else:
                    gap = within - cross  # within should have higher ARI

                log_tau = tsub["log10_tau"].iloc[0]
                gap_by_tau.append((log_tau, gap))

            if gap_by_tau:
                taus, gaps = zip(*gap_by_tau)
                ax.plot(taus, gaps, label=BRAIN_BAND_TEX_DICT[band], lw=2)

        ax.set_xlabel("log₁₀(τ)")
        ax.set_ylabel(ylabel)
        ax.axhline(0, color="gray", ls=":", alpha=0.4)
        ax.legend(fontsize=8)

    fig.suptitle("Task 4b: Within-vs-cross gap across τ — band comparison", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT_TAU / "tau_sensitive_bands.pdf",
        what="Two panels showing the within-vs-cross gap as a function of log₁₀(τ), one curve per band. "
             "Left: gap in scaled distance (positive = cross pairs more distant, which supports reorganization). "
             "Right: gap in ARI (positive = within pairs more similar).",
        proves="Whether some bands show τ-dependent reorganization structure — i.e., reorganization that "
               "emerges or vanishes at certain diffusion scales. If curves are flat, τ' is sufficient. "
               "If some bands peak at different τ, the multiscale approach reveals band-specific organization.",
        how_to_read="Higher curve = stronger within-vs-cross separation at that τ. "
                    "If all bands peak at the same τ, a single scale is enough. "
                    "If bands peak at different τ, τ-multiscale adds value and different bands "
                    "organize at different spatial scales.",
    )

    return tau_df


# ══════════════════════════════════════════════════════════════════════════
# Task 4c — Added value assessment
# ══════════════════════════════════════════════════════════════════════════
def task_4c(h_profiles: list, tau_df: pd.DataFrame):
    print("\n── Task 4c: Added value assessment ──")

    md_lines = [
        "# Task 4c — Added Value Assessment", "",
        "Comparison of discriminative power across approaches:", "",
        "1. **Single-point** — one τ (= τ'), one h (= Ψ-optimal)",
        "2. **h-multiscale** — ARI across dendrogram cut levels at fixed τ'",
        "3. **τ-multiscale** — distance/ARI across diffusion times with optimal h(τ)",
        "4. **Raw FC** — no LRG processing at all", "",
        "---", "",
    ]

    # h-multiscale: compute gap between within and cross ARI at each n
    md_lines.append("## h-multiscale summary")
    md_lines.append("")

    for band in BAND_ORDER:
        band_profiles = [p for p in h_profiles if p["band"] == band]
        if not band_profiles:
            continue

        # Find n where within-cross gap is maximized
        # Collect all profiles by pair type
        within_arrs = []
        cross_arrs = []
        min_len = min(len(p["ari_profile"]) for p in band_profiles)
        n_range = band_profiles[0]["n_range"][:min_len]

        for p in band_profiles:
            arr = p["ari_profile"][:min_len]
            if p["pair_type"] == "within":
                within_arrs.append(arr)
            else:
                cross_arrs.append(arr)

        if not within_arrs or not cross_arrs:
            continue

        within_mean = np.mean(within_arrs, axis=0)
        cross_mean = np.mean(cross_arrs, axis=0)
        gap = within_mean - cross_mean
        best_n_idx = np.argmax(gap)
        best_n = n_range[best_n_idx]
        best_gap = gap[best_n_idx]

        md_lines.append(
            f"- **{band}**: max gap = {best_gap:.3f} at n = {best_n} clusters "
            f"(within ARI = {within_mean[best_n_idx]:.3f}, "
            f"cross ARI = {cross_mean[best_n_idx]:.3f})"
        )

    md_lines += ["", "## τ-multiscale summary", ""]

    for band in BAND_ORDER:
        sub = tau_df[(tau_df["band"] == band) & tau_df["ari"].notna()]
        if sub.empty:
            continue

        # Gap at each tau
        best_gap = -np.inf
        best_tau = None
        for ti in sorted(sub["tau_idx"].unique()):
            tsub = sub[sub["tau_idx"] == ti]
            within_ari = tsub[tsub["pair_type"] == "within"]["ari"].mean()
            cross_ari = tsub[tsub["pair_type"] == "cross"]["ari"].mean()
            gap = within_ari - cross_ari
            if gap > best_gap:
                best_gap = gap
                best_tau = tsub["log10_tau"].iloc[0]

        md_lines.append(
            f"- **{band}**: max within-cross ARI gap = {best_gap:.3f} "
            f"at log₁₀(τ) = {best_tau:.2f}"
        )

    md_lines += [
        "", "## Comparison with single-point and raw FC", "",
        "See task3_analysis results for single-point and raw FC comparisons.", "",
        "Key questions answered:", "",
        "1. **Does h-multiscale add value?** Check if the best-n gap exceeds the "
        "Ψ-optimal gap (n* from single-point). If yes, the optimal cut is not "
        "always the most discriminative.",
        "",
        "2. **Does τ-multiscale add value?** Check if the best-τ gap exceeds the "
        "gap at τ'. If some bands peak at τ ≠ τ', τ-multiscale reveals "
        "band-specific organization scales.",
        "",
        "3. **Does LRG add value over raw FC?** Compare ARI/distance gap with "
        "raw FC rank-distance gap from task3. If LRG achieves better within-vs-cross "
        "separation, it reveals hidden structure.",
    ]

    out_path = OUT_ROOT / "task4_multiscale" / "added_value_assessment.md"
    out_path.write_text("\n".join(md_lines))
    print(f"  Saved: {out_path}")


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.time()

    print("Loading data...")
    lrg_data = load_all_lrg()
    fc_data = load_all_fc()

    h_profiles = task_4a(lrg_data)
    tau_df = task_4b(fc_data, lrg_data)
    task_4c(h_profiles, tau_df)

    elapsed = time.time() - t0
    print(f"\nTask 4 complete ({elapsed:.1f}s). Outputs in {OUT_ROOT / 'task4_multiscale'}")


if __name__ == "__main__":
    main()
