"""Figure-level layout helpers.

Promote the repeated ``fig.legend(...)`` and ``fig.text(provenance)``
boilerplate from individual scripts into a single source of truth.

See ``.agents/guides/05_plotting/legends.md`` and
``.agents/guides/05_plotting/multi-axis-figures.md`` for the rationale.
"""
from __future__ import annotations

from typing import Iterable, Literal, Sequence

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.legend import Legend

__all__ = [
    "figure_legend",
    "add_provenance_footer",
]

_Where = Literal["bottom", "top", "right"]


def figure_legend(
    fig: Figure,
    handles: Sequence,
    *,
    where: _Where = "bottom",
    fontsize: int = 8,
    pad: float = 0.04,
    ncol: int | None = None,
    frameon: bool = False,
) -> Legend:
    """Place a figure-level legend in the canonical position.

    Use whenever a legend describes content shared across more than one
    axis. Never let such a legend live inside one panel — it breaks the
    visual symmetry of multi-axis grids.

    Parameters
    ----------
    fig
        The matplotlib Figure to attach the legend to.
    handles
        Iterable of artists / proxy artists with ``label`` attributes.
        Typical: ``[Line2D(..., label=...), ...]``.
    where
        ``"bottom"`` (default), ``"top"``, or ``"right"``.
    fontsize
        Default 8. Use 6-7 for very dense grids.
    pad
        Distance from the figure edge in figure-fraction units.
        Default 0.04. The caller is expected to reserve the same strip
        in ``fig.tight_layout(rect=...)``.
    ncol
        Override the auto-chosen number of columns. By default,
        horizontal legends use ``ncol = len(handles)`` (single row);
        right-side vertical uses ``ncol = 1``.
    frameon
        Default ``False`` — borderless looks cleaner.

    Returns
    -------
    matplotlib.legend.Legend
        The created legend, in case the caller wants to tweak it.

    Notes
    -----
    After calling this, prefer:
    ``fig.tight_layout(rect=[0, pad, 1, 1])`` for ``where="bottom"``,
    or analogous reservations for ``"top"``/``"right"``.
    """
    if where == "bottom":
        kw = dict(loc="lower center", bbox_to_anchor=(0.5, -pad))
        default_ncol = len(handles)
    elif where == "top":
        kw = dict(loc="upper center", bbox_to_anchor=(0.5, 1 + pad))
        default_ncol = len(handles)
    elif where == "right":
        kw = dict(loc="center right", bbox_to_anchor=(1 + pad, 0.5))
        default_ncol = 1
    else:
        raise ValueError(f"unknown 'where': {where!r}; "
                         "use 'bottom', 'top', or 'right'")

    return fig.legend(
        handles=list(handles),
        ncol=ncol if ncol is not None else default_ncol,
        fontsize=fontsize,
        frameon=frameon,
        **kw,
    )


def add_provenance_footer(
    fig: Figure,
    label: str,
    *,
    fontsize: int = 6,
    color: str = "gray",
) -> None:
    """Add a small provenance footer in the bottom-right corner.

    Convention adopted from ``audit_25/26/28``. Helps trace stray PDFs
    back to the producing script and page.
    """
    fig.text(0.99, 0.01, label, ha="right",
             fontsize=fontsize, color=color)
