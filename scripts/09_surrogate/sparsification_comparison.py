#!/usr/bin/env python
"""Compare MSC sparsification methods across patient/band/phase triplets.

Two-phase design to control memory:
  Phase 1 (--compute): runs each sparsification method in a SEPARATE subprocess
           so surrogate arrays are freed between methods.
  Phase 2 (--plot):    loads cached .npy matrices and generates all figures.

Usage:
  python scripts/sparsification_comparison.py --compute   # heavy: surrogates
  python scripts/sparsification_comparison.py --plot       # lightweight: figures only
  python scripts/sparsification_comparison.py              # both phases sequentially
"""

from __future__ import annotations

import argparse
import gc
import subprocess
import sys
import textwrap
import time
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.sparse.csgraph import laplacian

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

TRIPLETS = [
    ("Pat_02", "alpha", "rest_pre"),
    ("Pat_03", "beta",  "rest_pre"),
    ("Pat_05", "theta", "rest_pre"),
]

# Methods that don't need surrogates first, so they're fast and low-memory
METHODS_NO_SURR = ["none", "disparity", "ecm"]
# Methods that need surrogates — expensive
METHODS_SURR = ["soft", "fdr", "hybrid"]
METHODS = METHODS_NO_SURR + METHODS_SURR

METHOD_LABELS = {
    "none":      "Dense (none)",
    "soft":      "Soft surrogate",
    "fdr":       "FDR (q=0.05)",
    "disparity": "Disparity filter",
    "hybrid":    "Hybrid (surr+disp)",
    "ecm":       "ECM",
}

METHOD_COLORS = {
    "none":      "#888888",
    "soft":      "#4C72B0",
    "fdr":       "#DD8452",
    "disparity": "#55A868",
    "hybrid":    "#C44E52",
    "ecm":       "#8172B3",
}

from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH

NPERSEG = 4096
N_SURROGATES = 200
OUTPUT_DIR = FIGURES_ROOT / "sparsification_comparison"


# -----------------------------------------------------------------------
# Phase 1: Compute — each method in its own subprocess
# -----------------------------------------------------------------------

_COMPUTE_WORKER_SIMPLE = textwrap.dedent("""\
import sys, gc
sys.path.insert(0, "src")
from lrg_eegfc.workflow.msc import compute_msc_matrix

result = compute_msc_matrix(
    "{patient}", "{phase}", "{band}",
    sparsify="{method}",
    n_surrogates={n_surr},
    nperseg={nperseg},
    n_workers=2,
    verbose=True,
)
print(f"DONE  edges_nonzero={{(result.adjacency_matrix > 0).sum() // 2}}  "
      f"mean_msc={{result.mean_msc:.6f}}")
del result
gc.collect()
""")

# Compute all surrogate-based methods in ONE subprocess to avoid
# regenerating surrogates 3 times.
_COMPUTE_WORKER_SURR_ALL = textwrap.dedent("""\
import sys, gc
import numpy as np
from pathlib import Path
sys.path.insert(0, "src")

from lrg_eegfc.config.const import BRAIN_BANDS, DEFAULT_NPERSEG
from lrg_eegfc.utils.io import load_patient_dataset_robust
from lrg_eegfc.utils.fc.msc import (
    surrogate_msc_null, soft_sparsify_surrogate,
    fdr_sparsify_surrogate, hybrid_sparsify,
)
from lrg_eegfc.workflow.msc import (
    get_msc_cache_path, load_msc_matrix, DEFAULT_MSC_CACHE_ROOT,
)

patient, phase, band = "{patient}", "{phase}", "{band}"
nperseg = {nperseg}
n_surrogates = {n_surr}
cache_root = DEFAULT_MSC_CACHE_ROOT

# Check what's already cached
methods_todo = []
for method in ["soft", "fdr", "hybrid"]:
    kwargs = dict(patient=patient, phase=phase, band=band,
                  sparsify=method, nperseg=nperseg, n_surrogates=n_surrogates)
    if method == "fdr":
        kwargs["fdr_q"] = 0.05
    if method == "hybrid":
        kwargs["disparity_alpha"] = 0.05
    cached = load_msc_matrix(**kwargs)
    if cached is not None:
        nz = (cached > 0).sum() // 2
        print(f"[{{method}}] already cached — {{nz}} edges")
    else:
        methods_todo.append(method)

if not methods_todo:
    print("All surrogate methods already cached!")
    sys.exit(0)

print(f"Need to compute: {{methods_todo}}")

# Load dense MSC (always cached at this point)
dense_msc = load_msc_matrix(patient, phase, band, sparsify="none", nperseg=nperseg)
assert dense_msc is not None, "Dense MSC must be cached first"
print(f"Dense MSC loaded: {{dense_msc.shape}}")

# Load timeseries for surrogates
dataset = load_patient_dataset_robust(patient, SEEG_DATAPATH, phases=[phase])
recording = dataset[phase]
data = recording.timeseries
fs = float(recording.parameters.get('fs', 2048.0))
print(f"Data loaded: {{data.shape}}, fs={{fs}} Hz")

# Compute surrogates ONCE
bands_dict = {{band: BRAIN_BANDS[band]}}
print(f"Computing {{n_surrogates}} surrogates...")
W_null = surrogate_msc_null(data, fs, bands_dict, n_surrogates,
                            nperseg=nperseg, n_workers=2)

# Free timeseries
del data, dataset, recording
gc.collect()

W_null_band = W_null[band]
print(f"Surrogates computed: {{W_null_band.shape}}")

# Apply each method
for method in methods_todo:
    print(f"Applying {{method}}...")
    if method == "soft":
        A = soft_sparsify_surrogate(dense_msc, W_null_band)
        cache_path = get_msc_cache_path(patient, phase, band, cache_root,
            sparsify="soft", n_surrogates=n_surrogates, nperseg=nperseg)
    elif method == "fdr":
        A = fdr_sparsify_surrogate(dense_msc, W_null_band, q=0.05)
        cache_path = get_msc_cache_path(patient, phase, band, cache_root,
            sparsify="fdr", n_surrogates=n_surrogates, nperseg=nperseg, fdr_q=0.05)
    elif method == "hybrid":
        A = hybrid_sparsify(dense_msc, W_null_band, alpha=0.05)
        cache_path = get_msc_cache_path(patient, phase, band, cache_root,
            sparsify="hybrid", n_surrogates=n_surrogates, nperseg=nperseg,
            disparity_alpha=0.05)

    np.fill_diagonal(A, 0.0)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_path, A)
    nz = (A > 0).sum() // 2
    print(f"  [{{method}}] saved {{cache_path.name}} — {{nz}} edges")
    del A

del W_null, W_null_band
gc.collect()
print("DONE all surrogate methods")
""")


def run_subprocess(code: str, label: str, timeout: int = 1800):
    """Run a Python code snippet in a child process."""
    print(f"  [{label}] launching subprocess …")
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, timeout=timeout,
    )
    elapsed = time.time() - t0
    for line in proc.stdout.strip().splitlines():
        print(f"    {line}")
    if proc.returncode != 0:
        print(f"  [{label}] FAILED (rc={proc.returncode}) in {elapsed:.0f}s")
        for line in proc.stderr.strip().splitlines()[-15:]:
            print(f"    ERR: {line}")
        return False
    print(f"  [{label}] OK in {elapsed:.0f}s")
    return True


def phase_compute():
    """Phase 1: compute all sparsified matrices (subprocess per method)."""
    for patient, band, phase in TRIPLETS:
        print(f"\n{'='*60}")
        print(f"COMPUTING: {patient} / {band} / {phase}")
        print(f"{'='*60}")

        # Non-surrogate methods: each in its own subprocess (fast)
        for method in METHODS_NO_SURR:
            code = _COMPUTE_WORKER_SIMPLE.format(
                patient=patient, phase=phase, band=band,
                method=method, n_surr=0, nperseg=NPERSEG,
            )
            run_subprocess(code, method)

        # Surrogate methods: ONE subprocess computes surrogates once,
        # then applies soft/fdr/hybrid from the same null.
        code = _COMPUTE_WORKER_SURR_ALL.format(
            patient=patient, phase=phase, band=band,
            n_surr=N_SURROGATES, nperseg=NPERSEG,
        )
        run_subprocess(code, "soft+fdr+hybrid", timeout=1800)


# -----------------------------------------------------------------------
# Phase 2: Load from cache + metrics + plots
# -----------------------------------------------------------------------

def load_cached_matrix(patient, phase, band, method) -> np.ndarray | None:
    """Load a sparsified MSC matrix from cache."""
    from lrg_eegfc.workflow.msc import load_msc_matrix
    kwargs = dict(
        patient=patient, phase=phase, band=band,
        sparsify=method, nperseg=NPERSEG,
    )
    if method in ("soft", "fdr", "hybrid"):
        kwargs["n_surrogates"] = N_SURROGATES
    if method == "fdr":
        kwargs["fdr_q"] = 0.05
    if method in ("disparity", "hybrid"):
        kwargs["disparity_alpha"] = 0.05
    if method == "ecm":
        kwargs["ecm_alpha"] = 0.05
        kwargs["ecm_n_ensemble"] = 100
        kwargs["ecm_weight_scale"] = 1000
    return load_msc_matrix(**kwargs)


def compute_graph_metrics(A: np.ndarray) -> dict:
    """Compute graph metrics from adjacency matrix."""
    N = A.shape[0]

    # Build topology graph (nonzero edges only)
    G = nx.Graph()
    G.add_nodes_from(range(N))
    triu_i, triu_j = np.triu_indices(N, k=1)
    for idx in range(len(triu_i)):
        w = A[triu_i[idx], triu_j[idx]]
        if w > 0:
            G.add_edge(int(triu_i[idx]), int(triu_j[idx]), weight=float(w))

    n_edges = G.number_of_edges()
    max_edges = N * (N - 1) / 2
    density = n_edges / max_edges if max_edges > 0 else 0
    degrees = np.array([G.degree(n) for n in range(N)])

    weights_upper = A[triu_i, triu_j]
    nonzero_weights = weights_upper[weights_upper > 0]
    mean_weight = float(nonzero_weights.mean()) if len(nonzero_weights) > 0 else 0.0

    n_components = nx.number_connected_components(G)
    largest_cc = max(len(c) for c in nx.connected_components(G)) if n_components > 0 else 0
    isolated_nodes = int((degrees == 0).sum())

    try:
        clustering = nx.average_clustering(G, weight="weight")
    except Exception:
        clustering = 0.0
    try:
        assortativity = nx.degree_assortativity_coefficient(G)
    except Exception:
        assortativity = float("nan")
    try:
        transitivity = nx.transitivity(G)
    except Exception:
        transitivity = 0.0

    # Laplacian spectrum
    L = laplacian(A, normed=False)
    if hasattr(L, "toarray"):
        L = L.toarray()
    eigenvalues = np.sort(np.linalg.eigvalsh(L))
    spectral_gap = float(eigenvalues[1] - eigenvalues[0]) if N > 1 else 0.0
    algebraic_connectivity = float(eigenvalues[1]) if N > 1 else 0.0

    return {
        "n_nodes": N,
        "n_edges": n_edges,
        "density": density,
        "mean_degree": float(degrees.mean()),
        "mean_weight": mean_weight,
        "n_components": n_components,
        "largest_cc": largest_cc,
        "isolated_nodes": isolated_nodes,
        "clustering": clustering,
        "assortativity": assortativity,
        "transitivity": transitivity,
        "spectral_gap": spectral_gap,
        "algebraic_connectivity": algebraic_connectivity,
        "eigenvalues": eigenvalues,
        "degrees": degrees,
    }


# -----------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------

def plot_adjacency_comparison(results, triplet_key, output_dir):
    methods_avail = [m for m in METHODS if m in results]
    n = len(methods_avail)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows))
    if nrows * ncols == 1:
        axes = np.array([[axes]])
    axes = np.atleast_2d(axes)
    fig.suptitle(f"Adjacency Matrices — {triplet_key}", fontsize=16, fontweight="bold", y=0.98)

    for idx, method in enumerate(methods_avail):
        ax = axes[idx // ncols, idx % ncols]
        A = results[method]["matrix"]
        m = results[method]["metrics"]
        im = ax.imshow(A, cmap="viridis", aspect="equal", interpolation="none")
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(
            f"{METHOD_LABELS[method]}\n"
            f"density={m['density']:.3f}, edges={m['n_edges']}, comp={m['n_components']}",
            fontsize=10,
        )
        ax.set_xlabel("Channel")
        ax.set_ylabel("Channel")

    # hide unused
    for idx in range(n, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = output_dir / f"adj_matrices_{triplet_key}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_network_comparison(results, triplet_key, output_dir):
    methods_avail = [m for m in METHODS if m in results]
    n = len(methods_avail)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    # Layout from dense graph
    A_dense = results.get("none", next(iter(results.values())))["matrix"]
    G_layout = nx.from_numpy_array(A_dense)
    pos = nx.spring_layout(G_layout, seed=42, k=0.5, iterations=100)

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 7 * nrows))
    if nrows * ncols == 1:
        axes = np.array([[axes]])
    axes = np.atleast_2d(axes)
    fig.suptitle(f"Networks — {triplet_key}", fontsize=16, fontweight="bold", y=0.98)

    for idx, method in enumerate(methods_avail):
        ax = axes[idx // ncols, idx % ncols]
        A = results[method]["matrix"]
        m = results[method]["metrics"]
        N = A.shape[0]

        G = nx.Graph()
        G.add_nodes_from(range(N))
        ti, tj = np.triu_indices(N, k=1)
        for k in range(len(ti)):
            w = A[ti[k], tj[k]]
            if w > 0:
                G.add_edge(int(ti[k]), int(tj[k]), weight=float(w))

        edges = list(G.edges(data=True))
        if edges:
            weights = [d["weight"] for _, _, d in edges]
            max_w = max(weights) if weights else 1
            widths = [0.1 + 2.0 * w / max_w for w in weights]
            nx.draw_networkx_edges(G, pos, ax=ax, width=widths, alpha=0.25,
                                   edge_color=METHOD_COLORS[method])

        degs = np.array([G.degree(n_) for n_ in range(N)])
        node_c = degs / max(degs.max(), 1)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=35, node_color=node_c,
                               cmap=plt.cm.YlOrRd, vmin=0, vmax=1,
                               edgecolors="black", linewidths=0.3)
        ax.set_title(
            f"{METHOD_LABELS[method]}\nedges={m['n_edges']}, <k>={m['mean_degree']:.1f}, "
            f"CC={m['clustering']:.3f}", fontsize=10)
        ax.axis("off")

    for idx in range(n, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = output_dir / f"networks_{triplet_key}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_metrics_comparison(results, triplet_key, output_dir):
    methods_avail = [m for m in METHODS if m in results]
    metric_keys = [
        ("density", "Density"),
        ("n_edges", "Edges"),
        ("mean_degree", "Mean Degree"),
        ("mean_weight", "Mean Weight"),
        ("n_components", "Components"),
        ("isolated_nodes", "Isolated Nodes"),
        ("clustering", "Clustering (wtd)"),
        ("assortativity", "Assortativity"),
        ("transitivity", "Transitivity"),
        ("spectral_gap", "Spectral Gap (λ₂-λ₁)"),
        ("algebraic_connectivity", "Algebraic Conn. (λ₂)"),
    ]

    fig, axes = plt.subplots(3, 4, figsize=(22, 14))
    fig.suptitle(f"Graph Metrics — {triplet_key}", fontsize=16, fontweight="bold", y=0.98)

    for i, (key, label) in enumerate(metric_keys):
        ax = axes[i // 4, i % 4]
        vals = [results[m]["metrics"][key] for m in methods_avail]
        cols = [METHOD_COLORS[m] for m in methods_avail]
        xlabels = [METHOD_LABELS[m].split("(")[0].strip() for m in methods_avail]
        bars = ax.bar(range(len(vals)), vals, color=cols, edgecolor="black", linewidth=0.5)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=8)
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        for bar, v in zip(bars, vals):
            txt = f"{v:.4f}" if isinstance(v, float) and abs(v) < 10 else f"{v}"
            ax.text(bar.get_x() + bar.get_width() / 2, max(bar.get_height(), 0),
                    txt, ha="center", va="bottom", fontsize=7)
    axes[2, 3].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = output_dir / f"metrics_{triplet_key}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_laplacian_spectrum(results, triplet_key, output_dir):
    methods_avail = [m for m in METHODS if m in results]
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    fig.suptitle(f"Laplacian Spectrum — {triplet_key}", fontsize=16, fontweight="bold", y=1.02)

    ax = axes[0]
    for method in methods_avail:
        eigs = results[method]["metrics"]["eigenvalues"]
        ax.plot(eigs, label=METHOD_LABELS[method], color=METHOD_COLORS[method], lw=1.5)
    ax.set_xlabel("Index"); ax.set_ylabel("λ")
    ax.set_title("Full Spectrum", fontweight="bold"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    for method in methods_avail:
        eigs = results[method]["metrics"]["eigenvalues"]
        n_show = min(30, len(eigs))
        ax.plot(range(n_show), eigs[:n_show], "o-", label=METHOD_LABELS[method],
                color=METHOD_COLORS[method], lw=1.5, ms=4)
    ax.set_xlabel("Index"); ax.set_ylabel("λ")
    ax.set_title("First 30 (Gap Region)", fontweight="bold"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[2]
    for method in methods_avail:
        eigs = results[method]["metrics"]["eigenvalues"]
        eigs_pos = eigs[eigs > 1e-10]
        if len(eigs_pos) > 0:
            ax.semilogy(range(len(eigs_pos)), eigs_pos, label=METHOD_LABELS[method],
                        color=METHOD_COLORS[method], lw=1.5)
    ax.set_xlabel("Index (nonzero)"); ax.set_ylabel("λ (log)")
    ax.set_title("Log Scale", fontweight="bold"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    plt.tight_layout()
    out = output_dir / f"laplacian_spectrum_{triplet_key}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_degree_distributions(results, triplet_key, output_dir):
    methods_avail = [m for m in METHODS if m in results]
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    fig.suptitle(f"Degree Distributions — {triplet_key}", fontsize=16, fontweight="bold", y=1.02)

    ax = axes[0]
    for method in methods_avail:
        degs = results[method]["metrics"]["degrees"]
        if degs.max() > 0:
            bins = np.arange(0, degs.max() + 2) - 0.5
            ax.hist(degs, bins=bins, alpha=0.35, label=METHOD_LABELS[method],
                    color=METHOD_COLORS[method], edgecolor=METHOD_COLORS[method], lw=0.8)
    ax.set_xlabel("Degree k"); ax.set_ylabel("Count")
    ax.set_title("Degree Histogram", fontweight="bold"); ax.legend(fontsize=7); ax.grid(alpha=0.3)

    ax = axes[1]
    for method in methods_avail:
        degs = results[method]["metrics"]["degrees"]
        sd = np.sort(degs)
        ccdf = 1.0 - np.arange(1, len(sd) + 1) / len(sd)
        ax.plot(sd, ccdf, label=METHOD_LABELS[method], color=METHOD_COLORS[method], lw=1.5)
    ax.set_xlabel("Degree k"); ax.set_ylabel("P(K > k)")
    ax.set_title("Degree CCDF", fontweight="bold"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[2]
    for method in methods_avail:
        A = results[method]["matrix"]
        w = A[np.triu_indices(A.shape[0], k=1)]
        nz = w[w > 0]
        if len(nz) > 0:
            ax.hist(nz, bins=50, alpha=0.35, density=True,
                    label=f"{METHOD_LABELS[method]} (n={len(nz)})",
                    color=METHOD_COLORS[method], edgecolor=METHOD_COLORS[method], lw=0.8)
    ax.set_xlabel("Edge weight"); ax.set_ylabel("Density")
    ax.set_title("Weight Distribution (nonzero)", fontweight="bold"); ax.legend(fontsize=7); ax.grid(alpha=0.3)

    plt.tight_layout()
    out = output_dir / f"degree_dist_{triplet_key}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def print_summary_table(all_results):
    header = (
        f"{'Triplet':<25} {'Method':<18} {'Edges':>6} {'Density':>8} "
        f"{'<k>':>6} {'<w>':>8} {'Comp':>5} {'Isol':>5} "
        f"{'CC':>7} {'Assort':>7} {'l2':>8} {'Gap':>8}"
    )
    print("\n" + "=" * len(header))
    print("SPARSIFICATION COMPARISON SUMMARY")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for triplet_key, results in all_results.items():
        for method in METHODS:
            if method not in results:
                continue
            m = results[method]["metrics"]
            assort = f"{m['assortativity']:>7.4f}" if not np.isnan(m['assortativity']) else "     NA"
            print(
                f"{triplet_key:<25} {METHOD_LABELS[method]:<18} "
                f"{m['n_edges']:>6} {m['density']:>8.4f} "
                f"{m['mean_degree']:>6.1f} {m['mean_weight']:>8.5f} "
                f"{m['n_components']:>5} {m['isolated_nodes']:>5} "
                f"{m['clustering']:>7.4f} {assort} "
                f"{m['algebraic_connectivity']:>8.4f} {m['spectral_gap']:>8.4f}"
            )
        print("-" * len(header))


def phase_plot():
    """Phase 2: load cached matrices and generate all figures."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = {}

    for patient, band, phase in TRIPLETS:
        triplet_key = f"{patient}_{band}_{phase}"
        print(f"\n--- Loading: {triplet_key} ---")
        results = {}

        for method in METHODS:
            A = load_cached_matrix(patient, phase, band, method)
            if A is None:
                print(f"  [{method}] not cached — skipping")
                continue
            metrics = compute_graph_metrics(A)
            results[method] = {"matrix": A, "metrics": metrics}
            print(f"  [{method}] loaded {A.shape} — density={metrics['density']:.4f}, "
                  f"comp={metrics['n_components']}, l2={metrics['algebraic_connectivity']:.4f}")

        if len(results) < 2:
            print(f"  Not enough methods cached, skipping figures")
            continue

        all_results[triplet_key] = results

        print(f"  Generating figures …")
        plot_adjacency_comparison(results, triplet_key, OUTPUT_DIR)
        plot_network_comparison(results, triplet_key, OUTPUT_DIR)
        plot_metrics_comparison(results, triplet_key, OUTPUT_DIR)
        plot_laplacian_spectrum(results, triplet_key, OUTPUT_DIR)
        plot_degree_distributions(results, triplet_key, OUTPUT_DIR)

    if all_results:
        print_summary_table(all_results)

    print(f"\nAll figures saved to: {OUTPUT_DIR}")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="MSC sparsification comparison")
    parser.add_argument("--compute", action="store_true", help="Phase 1: compute (subprocess per method)")
    parser.add_argument("--plot", action="store_true", help="Phase 2: load cache + plot")
    args = parser.parse_args()

    # Default: both phases
    if not args.compute and not args.plot:
        args.compute = True
        args.plot = True

    if args.compute:
        phase_compute()
    if args.plot:
        phase_plot()


if __name__ == "__main__":
    main()
