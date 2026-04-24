"""Backwards compatible wrapper around :mod:`lrg_eegfc.cli`."""

from __future__ import annotations

from lrg_eegfc.cli._legacy import main


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
