"""Shared helpers for the epileptic-node-stratified trace audits.

Private sibling module (cf. ``_fc_split_half.py``) imported by
``audit_77_epi_stratified_cophenetic`` and
``audit_78_epi_stratified_grassmann`` so the configuration vocabulary,
node/pair masking, surrogate-cache namespacing, per-cell RNG, and cohort
verdict logic live in exactly one place and cannot drift between the two
substrates.

Two families of "epi view" on the cross-phase trace:

- **Subgraph configs** rebuild the LRG on a node submatrix
  ``W[np.ix_(keep, keep)]``: ``full`` (all nodes), ``exclude_epi`` (non-epi
  nodes only), ``epi_only`` (epileptic nodes only). The matched-strength
  surrogate is regenerated on the submatrix (its own cache directory so an
  ``N``-mismatched ensemble never clobbers the canonical one).

- **Pair-class configs** keep the *full-graph* LRG cophenetic distances and
  restrict the *pair set*: ``nonepi_nonepi`` (both endpoints non-epi),
  ``cross`` (epi↔non-epi interface), ``epi_epi`` (both endpoints epi). The
  surrogate is the full-graph matched-strength ensemble with the per-pair
  statistic recomputed on the restricted pair set (the ``audit_71``
  ``pair_keep = keep[iu_i] & keep[iu_j]`` mechanism). Pair-class is
  cophenetic-only — the Grassmann subspace statistic lives on the whole
  eigenbasis and has no pair-restriction analogue.

Epi-node identity comes from the canonical ``build_epi_masks`` (red-font
implant-XLSX labels joined to FC-index-ordered channel labels) — never a
local copy (see ``audit_epi_mask_pattern`` memory; the ``np.isin(int, str)``
silent no-op bug).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.surrogate import surrogate_cache_path


# ---------------------------------------------------------------------------
# Vocabulary (single source of truth)
# ---------------------------------------------------------------------------
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PHASES_4 = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
PHASES_3 = ("rest_pre_A", "task_test", "rest_post")

SUBGRAPH_CONFIGS = ("full", "exclude_epi", "epi_only")
PAIRCLASS_CONFIGS = ("nonepi_nonepi", "cross", "epi_epi")
ALL_CONFIGS = SUBGRAPH_CONFIGS + PAIRCLASS_CONFIGS

CONFIG_CLASS = {**{c: "subgraph" for c in SUBGRAPH_CONFIGS},
                **{c: "pairclass" for c in PAIRCLASS_CONFIGS}}

# Honest hard floors (surfaced downstream as defined=False, never silent).
MIN_NODES = 3      # an LRG on < 3 leaves is degenerate
MIN_PAIRS = 10     # a Spearman over < 10 pairs is uninterpretable

# Default seeds. full / pair-class reuse the canonical seed-20260511 ensemble
# (audit_66 cache); exclude_epi reuses audit_67's seed-20260514 epiX cache;
# epi_only is a fresh ensemble keyed on the canonical seed.
SEED_CANONICAL = 20260511
SEED_EXCLUDE_EPI = 20260514

N_SURROGATES = 200
SWAP_FACTOR = 20


# ---------------------------------------------------------------------------
# Epi mask
# ---------------------------------------------------------------------------
def epi_mask_for(pat: str, n_expected: int | None = None) -> np.ndarray:
    """FC-index-ordered boolean epi mask via the canonical ``build_epi_masks``.

    Raises if its length disagrees with the FC matrix size (the alignment
    contract audit_67 also enforces).
    """
    mask = np.asarray(build_epi_masks(pat).epi_mask, dtype=bool)
    if n_expected is not None and mask.size != n_expected:
        raise ValueError(
            f"{pat}: epi_mask length {mask.size} != FC N {n_expected}"
        )
    return mask


# ---------------------------------------------------------------------------
# Node / pair masks per config
# ---------------------------------------------------------------------------
def node_mask_for_config(epi_mask: np.ndarray, config: str,
                         min_nodes: int = MIN_NODES) -> np.ndarray | None:
    """Boolean node-keep mask for a SUBGRAPH config.

    Returns ``None`` when the config is structurally undefined for this
    patient (e.g. ``epi_only`` with 0 epi nodes, or below ``min_nodes``).
    """
    if config == "full":
        return np.ones_like(epi_mask, dtype=bool)
    if config == "exclude_epi":
        keep = ~epi_mask
    elif config == "epi_only":
        keep = epi_mask.copy()
    else:
        raise ValueError(f"not a subgraph config: {config}")
    return keep if int(keep.sum()) >= min_nodes else None


def upper_tri_indices(n: int) -> tuple[np.ndarray, np.ndarray]:
    """``(iu_i, iu_j)`` in the order ``scipy.spatial.distance.squareform``
    emits — so a condensed cophenetic vector aligns elementwise with these.
    Matches the order produced by ``lrg_ultrametric_condensed`` /
    ``cophenetic_condensed_from_eigs``.
    """
    return np.triu_indices(n, k=1)


def pair_keep_for_config(epi_mask: np.ndarray, config: str,
                         min_pairs: int = MIN_PAIRS) -> np.ndarray | None:
    """Boolean pair-keep mask (condensed upper-tri order) for a PAIRCLASS
    config. ``None`` if the class has fewer than ``min_pairs`` pairs.
    """
    n = epi_mask.size
    iu_i, iu_j = upper_tri_indices(n)
    ei, ej = epi_mask[iu_i], epi_mask[iu_j]
    if config == "epi_epi":
        pk = ei & ej
    elif config == "nonepi_nonepi":
        pk = (~ei) & (~ej)
    elif config == "cross":
        pk = ei ^ ej            # exactly one endpoint epileptic
    else:
        raise ValueError(f"not a pair-class config: {config}")
    return pk if int(pk.sum()) >= min_pairs else None


# ---------------------------------------------------------------------------
# Surrogate-eig cache path per config (no clobber across configs)
# ---------------------------------------------------------------------------
SURR_DIR_EXCLUDE = CACHE_ROOT / "matched_strength_surrogate_epi_excluded_lrg"
SURR_DIR_EPIONLY = CACHE_ROOT / "matched_strength_surrogate_epi_only_lrg"


def config_seed(config: str) -> int:
    """Seed pinned per config so cache hits are deterministic and the
    existing canonical / epiX ensembles are reused."""
    return SEED_EXCLUDE_EPI if config == "exclude_epi" else SEED_CANONICAL


def surr_eig_path(config: str, pat: str, band: str, phase: str,
                  n_surr: int = N_SURROGATES, swap_factor: int = SWAP_FACTOR,
                  ) -> Path:
    """Cache path for the matched-strength eigendecomposition ensemble that a
    config's null needs.

    - ``full`` and all pair-class configs use the canonical full-graph cache
      (``surrogate_cache_path``, seed 20260511) — pair-class restricts that
      same ensemble, so it never needs its own files.
    - ``exclude_epi`` reuses audit_67's ``epiX`` cache (seed 20260514).
    - ``epi_only`` gets a fresh ``epiONLY`` cache (seed 20260511).
    """
    seed = config_seed(config)
    if config == "full" or config in PAIRCLASS_CONFIGS:
        return surrogate_cache_path(pat, band, phase, n_surr, swap_factor, seed)
    if config == "exclude_epi":
        return (SURR_DIR_EXCLUDE / pat
                / f"{band}_{phase}_epiX_R{n_surr}_swap{swap_factor}"
                  f"_seed{seed}_imcoh_abs.npz")
    if config == "epi_only":
        return (SURR_DIR_EPIONLY / pat
                / f"{band}_{phase}_epiONLY_R{n_surr}_swap{swap_factor}"
                  f"_seed{seed}_imcoh_abs.npz")
    raise ValueError(config)


# ---------------------------------------------------------------------------
# Deterministic per-cell RNG (cache reproducibility regardless of loop order)
# ---------------------------------------------------------------------------
_PHASE_IDX = {p: i for i, p in enumerate(PHASES_4)}


def cell_rng(pat: str, band: str, phase: str, config: str) -> np.random.Generator:
    """A reproducible generator seeded by ``(seed, config, patient, band,
    phase)`` so a freshly-computed ensemble is identical on re-run no matter
    the iteration order. (Pre-existing cache files are loaded by shape, so
    this only governs first-compute of new cells.)
    """
    pat_n = int(pat.split("_")[-1])
    ss = np.random.SeedSequence([
        config_seed(config), ALL_CONFIGS.index(config), pat_n,
        ALL_BANDS.index(band), _PHASE_IDX[phase],
    ])
    return np.random.default_rng(ss)


# ---------------------------------------------------------------------------
# Cohort verdict + sensitivity flag (shared by both substrates)
# ---------------------------------------------------------------------------
def cohort_verdict(obs: np.ndarray, surr_med: np.ndarray,
                   n_above: int) -> dict:
    """Cohort one-sided paired Wilcoxon (alternative='greater') of
    ``obs - surr_med`` plus the audit_63/66/67 verdict thresholds.

    ``n_above`` = count of patients with own-surrogate upper-tail
    ``obs_p_one_sided < 0.05`` (passed in by the caller, which has it).
    """
    from scipy.stats import wilcoxon

    obs = np.asarray(obs, dtype=float)
    surr_med = np.asarray(surr_med, dtype=float)
    m = np.isfinite(obs) & np.isfinite(surr_med)
    obs, surr_med = obs[m], surr_med[m]
    n = int(obs.size)
    if n < 3:
        return dict(wilcoxon_z=np.nan, wilcoxon_p=np.nan, verdict="undefined",
                    med_obs=np.nan, med_surr=np.nan, n_defined=n)
    try:
        wz, wp = wilcoxon(obs - surr_med, alternative="greater")
        wz, wp = float(wz), float(wp)
    except Exception:
        wz, wp = float("nan"), float("nan")
    med_obs = float(np.median(obs))
    med_surr = float(np.median(surr_med))
    if (np.isfinite(wp) and wp < 0.05
            and abs(med_surr) < 0.05 * max(1.0, abs(med_obs))
            and n_above >= 8):
        verdict = "separated"
    elif med_obs > 0 and med_surr > 0.5 * med_obs:
        verdict = "also_positive"
    else:
        verdict = "intermediate"
    return dict(wilcoxon_z=wz, wilcoxon_p=wp, verdict=verdict,
                med_obs=med_obs, med_surr=med_surr, n_defined=n)


def sensitivity_flag(full_p: float, cfg_p: float,
                     cfg_defined: bool = True) -> str:
    """Config verdict vs the ``full`` baseline (audit_67 ``_flag``):
    persist / weaken / emerge / absent, plus ``undefined`` / ``unknown``.
    """
    if not cfg_defined:
        return "undefined"
    try:
        fp, cp = float(full_p), float(cfg_p)
    except (TypeError, ValueError):
        return "unknown"
    if not (np.isfinite(fp) and np.isfinite(cp)):
        return "unknown"
    if fp < 0.05 and cp < 0.05:
        return "persist"
    if fp < 0.05 and cp >= 0.05:
        return "weaken"
    if fp >= 0.05 and cp < 0.05:
        return "emerge"
    return "absent"
