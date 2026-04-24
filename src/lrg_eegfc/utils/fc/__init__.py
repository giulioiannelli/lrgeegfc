"""Functional connectivity utilities (correlation + coherence).

The coherence family (MSC, ImCoh, future variants) lives under the
:mod:`coherence` subpackage.  The older :mod:`msc` subpackage remains as
a thin compatibility shim and still forwards ``compute_msc_welch`` /
``band_average_msc`` / ``coherence_fc_pipeline`` to the new backend.
"""
from .corr import *  # noqa: F401,F403
from .msc import *  # noqa: F401,F403
from . import coherence  # noqa: F401

__all__ = []
