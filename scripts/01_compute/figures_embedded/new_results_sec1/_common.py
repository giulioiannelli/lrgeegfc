"""Shared backbone + diffusion-scale layer for the ``new_results_sec1`` figures.

Every Results-§1 figure rebuilt on the **mst@0.20 multiscale backbone** imports from
here so the sparsifier, the reporting scale ``s`` and the output directory are defined
once. The dense single-scale (``tau=1/lambda_max``) generators under
``figures_embedded/`` are left untouched — these wrappers are the drop-in replacement
for their ``lrg_linkage`` / ``coph_square`` calls.

Pipeline (matches ``scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py``):
    dense imcoh_abs FC  ->  mst_union_top_fraction(W, 0.20)  ->  laplacian_eig
                        ->  cophenetic / linkage at scale s = tau * lambda_max.

Nothing manuscript-specific leaks into the library; only these script-side wrappers
know about mst@0.20 and s=5.6.
"""
from __future__ import annotations

import sys
import numpy as np
from numpy.typing import NDArray
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

# canonical cohort + 4-phase loader (dense imcoh_abs FC; A,B split-half rest_pre)
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase  # noqa: E402

from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction  # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import (  # noqa: E402
    laplacian_eig,
    cophenetic_at_scale,
    linkage_at_scale,
    rho_sym,
)
from lrg_eegfc.visuals.styles import use_lrg_style, band_color  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

# --------------------------------------------------------------------------- #
# backbone + scale constants (single source; keep in lock-step with scripts 12/13)
# --------------------------------------------------------------------------- #
FRAC = 0.20                                        # mst_union_top_fraction density
SGRID = np.logspace(0.0, np.log10(180.0), 16)      # s = tau*lambda_max, hardcoded grid
S_REPORT = float(SGRID[5])                         # 5.646 — honest joint alpha/beta scale
I_REPORT = 5
PHASES = ("A", "B", "task_test", "rest_post")

# band -> TeX (canonical source; falls back to the sparsified-arc convention)
try:
    from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT as BTeX  # noqa: E402
except Exception:                                                    # pragma: no cover
    BTeX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
            "beta": r"$\beta$", "low_gamma": r"$\gamma_{\mathrm{low}}$",
            "high_gamma": r"$\gamma_{\mathrm{high}}$"}

FIGDIR = ROOT / "data" / "preprint" / "figures" / "new_results_sec1"
DRAFTDIR = FIGDIR / "_drafts"

# axis label conventions (match the sparsified-arc figures)
XLAB_S = r"$s=\tau\lambda_{\max}$"
YLAB_RHO = r"$\rho_{\mathrm{sym}}^{\mathrm{coph}}$"

# data roots (CSV only)
SA = ROOT / "data" / "sparsified_arc"
MS = SA / "ms_mst020"
RECOV = SA / "sparsification_recovery"
ENCINF = SA / "enc_inf_arc_mst020"
LADDER = SA / "controls_ladder_apples"   # apples-to-apples: ALL descriptors on mst@0.20 (2026-07-13)
EPI = SA / "epi_arc_mst020"
REACH = SA / "spatial_reach_mst020"      # diffusion physical reach R(s) vs scale (2026-07-13)


# --------------------------------------------------------------------------- #
# controls-ladder shared layer  (fig_controls_ladder + fig_descriptor_bloom)
# reads data/sparsified_arc/controls_ladder_apples/per_cell.csv (all descriptors
# on the mst@0.20 backbone). Single source for the representativeness gate + the
# dot-matrix glyphs so both figures are pixel-identical.
# --------------------------------------------------------------------------- #
FUNCS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
FUNC_LABEL = {
    "T_test":       r"trace  ($T_{\mathrm{test}}$)",
    "T_learn":      r"encoding  ($T_{\mathrm{learn}}$)",
    "T_infspec":    r"inference  ($T_{\mathrm{inf}}$)",
    "T_infspec_pe": r"inference $\mid$ encoding  ($T_{\mathrm{inf}\mid e}$)",
}
LADDER_DESCS = [
    ("raw_fc",     "raw FC edges"),
    ("strength",   "node strength"),
    ("clustering", "clustering coef."),
    ("geodesic",   "graph geodesic"),
    ("resistance", "resistance (spectral)"),
    ("coph_meso",  "cophenetic (multiscale)"),
]
LADDER_HIER = {"coph_meso"}
P_SIG = 0.05
DEGEN_P = 0.999          # gate p >= this => degenerate by construction
MIN_ABOVE = 3            # >=3/10 patients must clear their own null (not 1-2-driven)


def _cohort_p(diff: NDArray) -> float:
    """One-sided paired Wilcoxon of (obs - surrogate p50) over patients."""
    diff = np.asarray(diff, float)
    diff = diff[~np.isnan(diff)]
    if diff.size == 0 or np.all(diff == 0):
        return 1.0
    return float(wilcoxon(diff, alternative="greater")[1])


def representativeness_gate(df, desc: str, band: str, functional: str,
                           min_above: int = MIN_ABOVE) -> dict:
    """Representativeness verdict for one (descriptor, band, functional) cell.

    A cell is a REPRESENTATIVE clear only if it passes a gate with no target band:
    cohort gate p<0.05 AND survives leave-one-out (worst p<0.05) AND >= ``min_above``
    of 10 patients clear their OWN matched-strength null. Returns
    ``level in {rep, fragile, null, degenerate}`` plus the numbers.
    """
    F = functional
    g = df[(df.descriptor == desc) & (df.band == band)]
    if g.empty:
        return dict(level="null", n_above=0, cohort_p=np.nan, loo_worst=np.nan)
    diff = (g[f"{F}_obs"] - g[f"{F}_surr_p50"]).values
    p = _cohort_p(diff)
    n_above = int((g[f"{F}_p"] < P_SIG).sum())
    if p >= DEGEN_P:
        return dict(level="degenerate", n_above=n_above, cohort_p=p, loo_worst=p)
    loo_worst = max(_cohort_p(np.delete(diff, k)) for k in range(len(diff)))
    if p < P_SIG and loo_worst < P_SIG and n_above >= min_above:
        level = "rep"
    elif p < P_SIG:
        level = "fragile"
    else:
        level = "null"
    return dict(level=level, n_above=n_above, cohort_p=p, loo_worst=loo_worst)


def gate_area(n_above: int) -> float:
    """Dot area (pt^2) encoding per-patient support n_above in 0..10."""
    return 26.0 + 42.0 * float(n_above)


def plot_gate_cell(ax, j: float, y: float, band: str, res: dict) -> None:
    """Draw one ladder cell: fill=representative, ring=sig-not-rep, open=null, dash=degen."""
    lvl, n = res["level"], res["n_above"]
    col = band_color(band)
    if lvl == "degenerate":
        ax.scatter(j, y, marker="_", s=190, color="0.62", lw=2.2, zorder=3)
    elif lvl == "rep":
        ax.scatter(j, y, s=gate_area(n), facecolor=col, edgecolor="white", lw=0.7, zorder=4)
    elif lvl == "fragile":
        ax.scatter(j, y, s=gate_area(n), facecolor="none", edgecolor=col, lw=1.9, zorder=3)
    else:
        ax.scatter(j, y, s=gate_area(n), facecolor="none", edgecolor="0.6", lw=1.1, zorder=2)


def eig_backbone(W: NDArray, frac: float = FRAC) -> tuple[NDArray, NDArray]:
    """``mst_union_top_fraction(W, frac)`` -> combinatorial-Laplacian eigenpair."""
    return laplacian_eig(mst_union_top_fraction(np.asarray(W, float), frac))


def tree_at_scale(W: NDArray, s: float = S_REPORT, frac: float = FRAC) -> NDArray:
    """UPGMA linkage ``Z`` of the mst@``frac`` diffusion tree at scale ``s``.

    Drop-in for the dense ``lrg_linkage(W)`` (which is fixed at ``tau=1/lambda_max``);
    here the tree is addressed by the dimensionless scale ``s = tau*lambda_max``.
    """
    ev, V = eig_backbone(W, frac)
    return linkage_at_scale(ev, V, s)


def coph_square_at_scale(W: NDArray, s: float = S_REPORT, frac: float = FRAC) -> NDArray:
    """Full ``N x N`` cophenetic distance matrix of the mst@``frac`` tree at scale ``s``.

    Drop-in for the dense ``coph_square(W)``.
    """
    ev, V = eig_backbone(W, frac)
    return squareform(cophenetic_at_scale(ev, V, s))


def eig_by_phase(pat: str, band: str, frac: float = FRAC) -> dict:
    """``{phase: (ev, V)}`` on the mst@``frac`` backbone for the 4-phase trace."""
    return {ph: eig_backbone(load_phase(pat, ph, band), frac) for ph in PHASES}


def s_index(s: float) -> int:
    """Nearest grid index of ``SGRID`` to ``s`` (for reading the CSV sweeps)."""
    return int(np.argmin(np.abs(SGRID - s)))


def log_height_linkage(Z: NDArray, margin: float = 0.8) -> NDArray:
    """Display transform: remap linkage merge heights to ``log10`` space.

    mst@0.20 cophenetic heights ``T = 1/rho`` are log-distributed across many
    decades (a near-disconnected node on the sparse backbone gives one enormous
    merge), so a *linear*-radius circular dendrogram
    (``plot_circular_dendrogram``, radius ``= r_outer - h/h_max * ...``) crushes
    ~95-99% of merges onto the rim and only the top 2-3 splits show. Remapping
    ``Z[:, 2] -> log10(h)`` spreads the hierarchy evenly over the radius — the same
    thing ``ax.set_yscale("log")`` does for a rectangular dendrogram. Topology
    (``Z[:, 0:2]``) and merge ORDER are untouched (``log10`` is monotone); only the
    radial spacing changes. Pair with ``plot_circular_dendrogram`` (its default
    ``h_max = Z[-1, 2]`` then normalises by the transformed root).
    """
    Z = np.asarray(Z, float).copy()
    h = Z[:, 2]
    hmin = float(h[h > 0].min())
    Z[:, 2] = np.log10(np.maximum(h, hmin)) - np.log10(hmin * margin)
    return Z
