"""``lrg-eegfc plot`` subcommand group.

Phase 4-B split (2026-05-29): the 17 plot subcommands previously
lived in a single 1025-LOC ``plot.py``. The module is now a package
with one submodule per command family:

- ``_fc.py``     -- corr, msc, msc-grid, msc-all-patients, msc-validation, cleaning, comparison
- ``_lrg.py``    -- lrg, lrg-phase-grid, lrg-video
- ``_reorg.py``  -- reorganization, reorg-metrics, reorg-summary, metric-correlation, metastable
- ``_misc.py``   -- cross, time-windows

Each command loads from cache and delegates to a visuals module function.
No computation is performed here.
"""
from __future__ import annotations

import click


@click.group()
def plot() -> None:
    """Generate visualizations from cached results."""


# Side-effect imports: each submodule's ``@plot.command(...)`` decorators
# register on the ``plot`` group above. Order doesn't matter -- click
# treats commands as a set on the group.
from . import _fc      # noqa: F401, E402
from . import _lrg     # noqa: F401, E402
from . import _reorg   # noqa: F401, E402
from . import _misc    # noqa: F401, E402

__all__ = ["plot"]
