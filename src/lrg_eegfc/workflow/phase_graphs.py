"""Cross-phase graph substrate: one injectable entry point for a four-phase cell.

The cross-phase analyses in this project all read the same four graphs per
(patient, band):

``A``, ``B``
    the two split halves of ``rest_pre`` -- the within-baseline reference pair
    that makes the cross-phase estimator split-half-referenced;
``task_test``, ``rest_post``
    the task episode and the rest that follows it.

``A`` / ``B`` come from the split-half FC cache; the two whole phases come from
the unified :func:`~lrg_eegfc.workflow.fc.load_fc_matrix`. That asymmetry was
previously re-implemented at the top of every script (and imported *from* a
script by six others), so a change to the substrate -- the FC transform, the
sparsifier, the clipping convention -- had six places to miss.

:func:`phase_graphs` is that single place. The substrate is an argument:
``fc_method`` selects the connectivity transform and ``backbone`` /
``backbone_kwargs`` select the sparsifier, so a downstream analysis states which
substrate it ran on instead of hardcoding one.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import NDArray

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix

__all__ = [
    "CROSS_PHASE_PHASES",
    "HALVES_CACHE",
    "load_phase_fc",
    "phase_graphs",
]

#: The four graphs of a cross-phase cell, in canonical order.
CROSS_PHASE_PHASES: tuple[str, ...] = ("A", "B", "task_test", "rest_post")

#: Split-half FC cache (``{band}_{phase}_{half}_{fc_method}.npy`` per patient).
HALVES_CACHE = CACHE_ROOT / "imcoh_halves_fc"


def load_phase_fc(
    patient: str,
    phase: str,
    band: str,
    *,
    fc_method: str = "imcoh_abs",
    source_phase: str = "rest_pre",
) -> NDArray:
    """Dense FC matrix for one phase of a cross-phase cell.

    ``phase`` in ``{"A", "B"}`` loads the split half of ``source_phase`` from
    :data:`HALVES_CACHE`; any other value is forwarded to
    :func:`~lrg_eegfc.workflow.fc.load_fc_matrix`. Either way the result is
    symmetrised, zero-diagonalled and clipped to ``[0, 1]`` -- the conventions
    the Laplacian downstream requires.
    """
    if phase in ("A", "B"):
        W = np.load(HALVES_CACHE / patient /
                    f"{band}_{source_phase}_{phase}_{fc_method}.npy")
    else:
        W = load_fc_matrix(patient, phase, band, fc_method=fc_method)
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def phase_graphs(
    patient: str,
    band: str,
    *,
    fc_method: str = "imcoh_abs",
    phases: tuple[str, ...] = CROSS_PHASE_PHASES,
    backbone: Optional[str] = None,
    backbone_kwargs: Optional[dict] = None,
) -> dict[str, NDArray]:
    """All four graphs of a cross-phase cell, optionally sparsified.

    Parameters
    ----------
    backbone
        ``None`` (default) returns the dense weighted FC matrices. A name
        accepted by :func:`~lrg_eegfc.utils.fc.backbone.select_backbone`
        (``"mst020"``, ``"tmfg"``, ``"pmfg"``, ``"disparity"``, ...) returns the
        sparsified graphs, applying the *same* filter to every phase so the four
        graphs remain comparable.
    backbone_kwargs
        Forwarded to ``select_backbone`` (e.g. ``{"frac": 0.20}``).

    Raises ``ValueError`` if the phases disagree on node count -- a mismatch
    means the cell cannot be compared pairwise and must not be silently padded.
    """
    Ws = {ph: load_phase_fc(patient, ph, band, fc_method=fc_method)
          for ph in phases}
    sizes = {ph: W.shape[0] for ph, W in Ws.items()}
    if len(set(sizes.values())) != 1:
        raise ValueError(f"{patient}/{band}: node-count mismatch across phases {sizes}")
    if backbone is None:
        return Ws
    from lrg_eegfc.utils.fc.backbone import select_backbone

    kw = dict(backbone_kwargs or {})
    return {ph: select_backbone(W, backbone, **kw) for ph, W in Ws.items()}
