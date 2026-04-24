#!/usr/bin/env python3
"""Diagnostic — noise sensitivity of the H2c / H2-RAW measures.

Scientific question (user-raised): do our Spearman-on-D measures
reflect actual tree topology, or are they arithmetic on D-values that
could be shifted by tiny numerical perturbations in the FC matrix?

Method: for each patient × band × phase, add Gaussian noise of scale
σ ∈ {0.5%, 1%, 2%, 5%} of the FC matrix's typical value to the raw
FC matrix and recompute LRG. Compare:

  (A) ρ_noise = Spearman(D_original, D_perturbed)
      — how much does noise alone drag down the D-matrix correlation?
  (B) ρ_cross_phase = Spearman(D^rpost_original, D^ttest_original)
      — the actual H2a-raw cross-phase correlation.
  (C) unique_bipartitions(Z_original, Z_perturbed)
      — how much does noise disrupt tree TOPOLOGY (count of internal
      bipartitions that exist in original but not in perturbed, out of
      N−1 possible)?

If ρ_noise at small σ stays near 1.0 and cross-phase ρ is notably
lower, our measures are topology-locked and the signal is real.
If ρ_noise drops quickly to the cross-phase level, our measures are
arithmetic-fragile — the cross-phase "signal" could be noise-level.

Runs on 3 patients × 2 bands for speed, each σ with 5 realizations.
No persistent cache (keeps perturbations isolated).

Writes: data/reports/imcoh_vi/diag_noise_sensitivity.{md,csv,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.hierarchy import leaders

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import compute_lrg_analysis, load_lrg_result


PATIENTS = ("Pat_02", "Pat_06", "Pat_08")
BANDS    = ("theta", "beta")   # one ergodic-candidate, one trace-candidate
SIGMAS   = (0.005, 0.01, 0.02, 0.05)
N_REPS   = 5


def _upper_tri(D: np.ndarray) -> np.ndarray:
    D = np.asarray(D)
    return D if D.ndim == 1 else D[np.triu_indices(D.shape[0], k=1)]


def _bipartitions(Z: np.ndarray, n: int) -> set[frozenset[int]]:
    """Return the set of non-trivial bipartitions induced by Z.

    Each internal node of Z defines a subset of leaves; that subset is
    a bipartition of the leaf set. Size of set = N - 1 (internal nodes).
    """
    members: dict[int, set[int]] = {i: {i} for i in range(n)}
    out: set[frozenset[int]] = set()
    for i, (a, b, *_) in enumerate(Z):
        a, b = int(a), int(b)
        node_id = n + i
        members[node_id] = members[a] | members[b]
        out.add(frozenset(members[node_id]))
    return out


def _perturbed_lrg(A: np.ndarray, sigma: float, seed: int,
                    pat: str, phase: str, band: str):
    """Recompute LRG on A + Gaussian noise scaled to sigma * std(A)."""
    rng = np.random.default_rng(seed)
    noise = rng.normal(0.0, sigma * float(A.std()), size=A.shape)
    noise = (noise + noise.T) / 2.0
    np.fill_diagonal(noise, 0.0)
    A_pert = np.maximum(A + noise, 0.0)  # keep non-negative for Laplacian
    # Use a synthetic phase label so we don't contaminate canonical caches.
    synth = f"{phase}__noise{sigma:g}_{seed}"
    r = compute_lrg_analysis(
        A_pert, pat, synth, band, fc_method="imcoh_abs",
        cache_root=IMCOH_LRG_CACHE / "_noise_probe",
        use_cache=False, overwrite_cache=True,
    )
    return r


def main() -> None:
    rows = []
    for pat in PATIENTS:
        for band in BANDS:
            # Reference: canonical LRG at rest_pre and task_test
            r_pre  = load_lrg_result(pat, "rest_pre",  band, "imcoh_abs", IMCOH_LRG_CACHE)
            r_tt   = load_lrg_result(pat, "task_test", band, "imcoh_abs", IMCOH_LRG_CACHE)
            r_post = load_lrg_result(pat, "rest_post", band, "imcoh_abs", IMCOH_LRG_CACHE)
            if r_pre is None or r_tt is None or r_post is None:
                continue
            D_pre   = _upper_tri(r_pre.ultrametric_matrix)
            D_tt    = _upper_tri(r_tt.ultrametric_matrix)
            D_post  = _upper_tri(r_post.ultrametric_matrix)
            # Reference cross-phase ρ
            rho_cross, _ = stats.spearmanr(D_post, D_tt)
            rho_base, _  = stats.spearmanr(D_post, D_pre)
            N = r_pre.linkage_matrix.shape[0] + 1
            bps_pre = _bipartitions(np.asarray(r_pre.linkage_matrix), N)

            # Perturb the rest_pre FC matrix and measure how much
            # Spearman + bipartitions shift under noise.
            A_pre = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs")
            for sigma in SIGMAS:
                for seed in range(N_REPS):
                    r_p = _perturbed_lrg(A_pre, sigma, seed, pat, "rest_pre", band)
                    D_p = _upper_tri(r_p.ultrametric_matrix)
                    rho_noise, _ = stats.spearmanr(D_pre, D_p)
                    bps_p = _bipartitions(np.asarray(r_p.linkage_matrix), N)
                    overlap = len(bps_pre & bps_p) / (N - 1)
                    rows.append({
                        "patient": pat, "band": band,
                        "sigma": sigma, "seed": seed,
                        "rho_noise": float(rho_noise),
                        "bip_overlap": float(overlap),
                        "rho_cross_ref": float(rho_cross),
                        "rho_base_ref": float(rho_base),
                    })
            print(f"  {pat}/{band}: rho_cross_phase={rho_cross:+.3f}  "
                  f"rho_base_phase={rho_base:+.3f}")

    df = pd.DataFrame(rows)
    out_dir = REPORTS_ROOT / "imcoh_vi"
    df.to_csv(out_dir / "diag_noise_sensitivity.csv", index=False)

    # ── figure ───────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), dpi=150)
    for ax, metric, ylabel, title in [
        (axes[0], "rho_noise",
         "Spearman(D_orig, D_noisy) — robustness of the ρ measure",
         "ρ under FC perturbation"),
        (axes[1], "bip_overlap",
         "fraction of tree bipartitions preserved",
         "Topology under FC perturbation"),
    ]:
        for band in BANDS:
            sub = df[df["band"] == band]
            means = sub.groupby("sigma")[metric].mean()
            stds = sub.groupby("sigma")[metric].std()
            ax.errorbar(means.index, means.values, yerr=stds.values,
                        marker="o", capsize=4, lw=1.8,
                        label=BRAIN_BAND_TEX_DICT[band])
        # Reference horizontal lines: typical cross-phase ρ values
        cross_mean = df["rho_cross_ref"].mean()
        ax.axhline(cross_mean, color="crimson", lw=1.0, linestyle="--",
                    label=f"cross-phase ρ reference ≈ {cross_mean:.2f}")
        ax.set_xscale("log")
        ax.set_xlabel("noise scale σ (fraction of std(FC))", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(0.0, 1.05)
        ax.legend(frameon=False, fontsize=9, loc="lower left")
        ax.grid(alpha=0.3)

    fig.suptitle(
        "Noise-sensitivity probe: does our ρ measure (left) survive FC "
        "perturbation better than tree topology (right)?",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out_png = out_dir / "figures" / "diag_noise_sensitivity.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out_png}")
    plt.close(fig)


if __name__ == "__main__":
    main()
