"""The analysed graph: one entry point, one locked substrate.

Every multiscale result in this project is a function of three choices made
*before* any operator runs -- the connectivity estimator's magnitude transform,
the sparsifier that takes the dense matrix to the graph the diffusion actually
sees, and the phase set the cross-phase statistic is defined on. Historically
each analysis script re-made those choices inline, which is how a project ends
up with results computed on subtly different graphs.

This module holds them once. Any analysis wanting "the graph" calls
:func:`canonical_graph`; anything wanting the diffusion spectrum calls
:func:`canonical_eig`; anything cross-phase calls :func:`canonical_phase_eigs`
and hands the result straight to
:func:`lrg_eegfc.utils.fc.heat_multiscale.cross_phase_functionals_over_scales`.

The substrate itself is the frozen :class:`Substrate` in :data:`CANONICAL`.
Overriding it is possible per call (robustness sweeps need to), but the default
is a single object so two analyses cannot silently disagree.

Knob-integrated readouts
------------------------
When a sparsifier's density is not fixed by construction, reporting one
fraction is reporting a choice. :func:`canonical_graph_ensemble` yields the
graph at every fraction of an invariance plateau, so a caller can integrate a
statistic over the knob and report the plateau median with the knob uncertainty
attached, instead of a single tuned value. Use it whenever
``CANONICAL.is_parameter_free`` is ``False``.

Phase set
---------
The phase set is **data**, never a hardcoded tuple. The transitive-inference
paradigm has five phases -- two split-half baseline arms plus ``task_learn``
(premises presented), ``task_test`` (novel non-adjacent pairs judged) and
``rest_post`` -- and separating *learning* a structure from *applying* it is
part of the science, so ``task_learn`` is a first-class phase everywhere.
Callers that only need the standard trace pass a four-phase set and get the
same code path; see :data:`CROSS_PHASE_ROLES`.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Iterator, Mapping, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from ..config.paths import IMCOH_HALVES_CACHE
from ..utils.fc.backbone import backbone_structure, select_backbone
from ..utils.fc.heat_multiscale import (
    CROSS_PHASE_ROLES,
    laplacian_eig,
)
from .fc import load_fc_matrix

__all__ = [
    "Substrate",
    "CANONICAL",
    "CANONICAL_PHASES",
    "CANONICAL_SCALES",
    "canonical_scale_grid",
    "canonical_graph",
    "canonical_graph_ensemble",
    "canonical_eig",
    "canonical_phase_graphs",
    "canonical_phase_eigs",
    "canonical_structure",
]


# --------------------------------------------------------------------------- #
# the locked substrate
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Substrate:
    """The three choices that define the graph an analysis runs on.

    Parameters
    ----------
    transform
        ImCoh magnitude transform, ``"abs"`` (``⟨|Im C|⟩_f``) or ``"sq"``
        (``⟨(Im C)²⟩_f``). Selected by the W0-A pre-registered rule R1.
    backbone
        Sparsifier family, as understood by
        :func:`lrg_eegfc.utils.fc.backbone.select_backbone`: ``"mst"``,
        ``"thresh"``, ``"tmfg"``, ``"pmfg"``, ``"perc"``, ``"disparity"``,
        ``"dense"``.
    frac
        Edge fraction, for the density-parameterised families only.
    disparity_alpha
        Significance level, for the disparity filter only.
    plateau
        The interval of ``frac`` over which the verdict was shown invariant
        (``(lo, hi)``), or ``None`` if the backbone has no free parameter.
        Reported alongside any number produced on this substrate.
    plateau_fracs
        The explicit fractions inside ``plateau`` that
        :func:`canonical_graph_ensemble` iterates, so a knob-integrated readout
        is reproducible rather than depending on a grid built at call time.
    """

    transform: str = "abs"
    backbone: str = "mst"
    frac: float = 0.20
    disparity_alpha: float = 0.20
    plateau: Optional[tuple[float, float]] = None
    plateau_fracs: tuple[float, ...] = ()
    notes: str = ""

    @property
    def is_parameter_free(self) -> bool:
        """``True`` when the backbone's density is fixed by construction."""
        return self.backbone in ("tmfg", "pmfg", "perc", "dense")

    @property
    def fc_method(self) -> str:
        """The ``fc_method`` string :func:`load_fc_matrix` expects."""
        return f"imcoh_{self.transform}"

    def label(self) -> str:
        """Short, filesystem-safe identifier of this substrate."""
        if self.backbone == "disparity":
            return f"imcoh_{self.transform}__disparity_a{self.disparity_alpha:g}"
        if self.is_parameter_free:
            return f"imcoh_{self.transform}__{self.backbone}"
        return f"imcoh_{self.transform}__{self.backbone}{self.frac:g}"

    def with_(self, **overrides) -> "Substrate":
        """A copy with fields replaced -- for robustness sweeps."""
        return replace(self, **overrides)


#: The substrate every lane inherits. Set by the W0-A contract; see
#: ``.agents/preprint/locked/PIPELINE_CONTRACT.md`` for the evidence and the
#: pre-registered rule that chose it. Change this in ONE place or not at all.
CANONICAL = Substrate(
    transform="abs",
    backbone="mst",
    frac=0.20,
    plateau=None,
    plateau_fracs=(),
    notes="provisional -- populated by the W0-A contract",
)

#: Default cross-phase phase set (five phases; ``A``/``B`` are the split halves
#: of ``rest_pre``). Ordered as roles are consumed, not alphabetically.
CANONICAL_PHASES: tuple[str, ...] = (
    "A", "B", "task_learn", "task_test", "rest_post",
)

#: Dimensionless diffusion scales ``s = tau * lambda_max``. Fixed grid so every
#: lane reports on the same axis; per-scale always, never a best-scale collapse.
CANONICAL_SCALES: NDArray = np.logspace(0.0, np.log10(180.0), 16)


def canonical_scale_grid() -> NDArray:
    """A copy of :data:`CANONICAL_SCALES` (callers must not mutate the module's)."""
    return CANONICAL_SCALES.copy()


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def _dense_fc(patient: str, phase: str, band: str, substrate: Substrate,
              halves_cache: Optional[Path]) -> NDArray:
    """Dense, cleaned FC for one phase under the substrate's transform.

    ``A`` / ``B`` are the split halves of ``rest_pre`` and come from
    :data:`lrg_eegfc.config.paths.IMCOH_HALVES_CACHE`; every other phase goes
    through :func:`load_fc_matrix`. A phase named ``"<phase>_A"`` /
    ``"<phase>_B"`` addresses the halves of that phase.
    """
    root = IMCOH_HALVES_CACHE if halves_cache is None else Path(halves_cache)
    if phase in ("A", "B"):
        src, half = "rest_pre", phase
    elif phase.endswith(("_A", "_B")):
        src, half = phase[:-2], phase[-1]
    else:
        src = half = None

    if half is not None:
        path = root / patient / f"{band}_{src}_{half}_imcoh_{substrate.transform}.npy"
        if not path.exists():
            raise FileNotFoundError(
                f"no split-half FC at {path}. Build it with "
                f"scripts/01_compute/paper_final/w0a_01b_cache_halves_both_transforms.py"
            )
        W = np.load(path)
    else:
        W = load_fc_matrix(patient, phase, band, fc_method=substrate.fc_method)
        if W is None:
            raise FileNotFoundError(
                f"no cached {substrate.fc_method} FC for {patient}/{phase}/{band}"
            )
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def canonical_graph(patient: str, phase: str, band: str, *,
                    substrate: Optional[Substrate] = None,
                    halves_cache: Optional[Path] = None,
                    dense: bool = False,
                    **overrides) -> NDArray:
    """The analysed graph for one (patient, phase, band).

    Dense FC under the substrate's transform, sparsified by the substrate's
    backbone. This is the single function every downstream lane calls, so no
    two analyses can run on different graphs by accident.

    Parameters
    ----------
    phase
        Any cached phase name; ``"A"`` / ``"B"`` are the split halves of
        ``rest_pre``, and ``"<phase>_A"`` / ``"<phase>_B"`` the halves of any
        other phase.
    substrate
        Override the locked :data:`CANONICAL` wholesale (robustness sweeps).
    dense
        Return the dense matrix without sparsifying -- the explicit escape
        hatch for a raw-FC baseline, so "raw" is a visible argument rather than
        a different code path.
    **overrides
        Field-level overrides of the substrate, e.g. ``frac=0.10``.

    Returns
    -------
    ndarray
        Symmetric, non-negative, zero-diagonal ``(N, N)`` adjacency.
    """
    sub = (substrate or CANONICAL)
    if overrides:
        sub = sub.with_(**overrides)
    W = _dense_fc(patient, phase, band, sub, halves_cache)
    if dense or sub.backbone == "dense":
        return W
    return select_backbone(W, sub.backbone, frac=sub.frac,
                           disparity_alpha=sub.disparity_alpha)


def canonical_graph_ensemble(patient: str, phase: str, band: str, *,
                             substrate: Optional[Substrate] = None,
                             halves_cache: Optional[Path] = None,
                             **overrides) -> Iterator[tuple[float, NDArray]]:
    """Yield ``(frac, graph)`` across the substrate's invariance plateau.

    The knob-integrated alternative to reporting one fraction: run the
    statistic at every yielded graph and report the plateau median with the
    spread, so the number carries its knob uncertainty. For a parameter-free
    backbone this yields the single graph with ``frac = nan``.
    """
    sub = (substrate or CANONICAL)
    if overrides:
        sub = sub.with_(**overrides)
    if sub.is_parameter_free or not sub.plateau_fracs:
        yield (float("nan"),
               canonical_graph(patient, phase, band, substrate=sub,
                               halves_cache=halves_cache))
        return
    W = _dense_fc(patient, phase, band, sub, halves_cache)
    for f in sub.plateau_fracs:
        yield float(f), select_backbone(W, sub.backbone, frac=float(f),
                                        disparity_alpha=sub.disparity_alpha)


def canonical_eig(patient: str, phase: str, band: str, **kwargs
                  ) -> tuple[NDArray, NDArray]:
    """Combinatorial-Laplacian eigendecomposition of :func:`canonical_graph`.

    ``L = D - W`` always (never a normalised variant). Returns
    ``(eigenvalues, eigenvectors)`` ready for the heat-kernel helpers.
    """
    return laplacian_eig(canonical_graph(patient, phase, band, **kwargs))


def canonical_phase_graphs(patient: str, band: str,
                           phases: Optional[Sequence[str]] = None,
                           **kwargs) -> dict[str, NDArray]:
    """``{phase: graph}`` over a phase set -- the set is data, not hardcoded.

    Defaults to :data:`CANONICAL_PHASES` (five phases). Raises if the phases do
    not share a node count, since a cross-phase per-pair statistic is undefined
    across different node sets.
    """
    phs = list(CANONICAL_PHASES if phases is None else phases)
    out = {ph: canonical_graph(patient, ph, band, **kwargs) for ph in phs}
    sizes = {ph: W.shape[0] for ph, W in out.items()}
    if len(set(sizes.values())) > 1:
        raise ValueError(f"phases disagree on node count for {patient}/{band}: {sizes}")
    return out


def canonical_phase_eigs(patient: str, band: str,
                         phases: Optional[Sequence[str]] = None,
                         **kwargs) -> dict[str, tuple[NDArray, NDArray]]:
    """``{phase: (eigenvalues, eigenvectors)}``, ready for the cross-phase helpers.

    Feed the result straight to
    :func:`lrg_eegfc.utils.fc.heat_multiscale.cross_phase_functionals_over_scales`,
    which reads the roles in :data:`CROSS_PHASE_ROLES` off whatever phases are
    present -- so dropping ``task_learn`` returns the standard trace alone with
    no other change.
    """
    return {ph: laplacian_eig(W)
            for ph, W in canonical_phase_graphs(patient, band, phases, **kwargs).items()}


def canonical_structure(patient: str, phase: str, band: str, **kwargs) -> dict:
    """Structural covariates of the analysed graph, plus what it kept.

    :func:`lrg_eegfc.utils.fc.backbone.backbone_structure` on the sparsified
    graph, plus ``weight_fraction`` (share of the dense graph's total edge
    weight retained) and the substrate label -- so any figure or table can
    state the graph it was computed on without re-deriving it.
    """
    sub = kwargs.pop("substrate", None) or CANONICAL
    field_overrides = {k: v for k, v in kwargs.items()
                       if k in ("frac", "backbone", "transform", "disparity_alpha")}
    if field_overrides:
        sub = sub.with_(**field_overrides)
    halves = kwargs.get("halves_cache")
    W = _dense_fc(patient, phase, band, sub, halves)
    B = canonical_graph(patient, phase, band, substrate=sub, halves_cache=halves)
    st = backbone_structure(B)
    st.update(substrate=sub.label(), patient=patient, phase=phase, band=band,
              weight_fraction=float(B.sum() / W.sum()) if W.sum() > 0 else float("nan"))
    return st
