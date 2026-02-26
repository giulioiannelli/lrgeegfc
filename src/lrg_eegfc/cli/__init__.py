"""Command-line interface for :mod:`lrg_eegfc`.

The public entry point is :data:`app`, a Click group registered as the
``lrg-eegfc`` console script in *pyproject.toml*.
"""

from ._app import app

__all__ = ["app"]
