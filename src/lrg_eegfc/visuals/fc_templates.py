"""FC adjacency-matrix templates — compact, options-driven helpers.

Library home for the figure-class introduced under
``.agents/guides/05_plotting/fc_templates/``.  The shipped templates
(``single_adjacency``, future ``row_per_phase``, mosaics) are thin
scripts that call into these helpers.

Style invariants (see ``05_plotting/`` rules):

- Full vector PDFs.  Never call ``set_rasterized(True)``.
- No ``fig.suptitle``; no watermark by default.
- Math axis labels (``$i$``, ``$j$``) by default.
- ``imshow_colorbar_caxdivider`` for the colorbar.
"""
from __future__ import annotations

import re
from typing import List, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrgsglib.plotlib.colorbars import imshow_colorbar_caxdivider


__all__ = [
    "plot_fc_adjacency",
    "plot_fc_adjacency_row",
    "plot_fc_adjacency_grid",
    "fc_method_colorbar_label",
    "probe_groups",
]


TickMode = Literal["generic", "index", "chnames"]


_PROBE_RE = re.compile(r"^([A-Za-z]+)")

# Fractional widening of the right-most panel in shared-scale row mode
# so that ``imshow_colorbar_caxdivider`` (size="4%", pad=0.06) carves
# its cax + pad out of the *extra* width without shrinking the imshow.
# Empirical: 4% size + ~2% pad-equivalent on a 2.4-inch panel.
_ROW_CB_FRAC = 0.06

# Templates use ``{band}`` as the trailing-subscript placeholder.  When
# ``band`` is provided, it is replaced by the LaTeX band glyph
# (`\beta`, `\alpha`, …) from ``BRAIN_BAND_TEX_DICT``; when absent, it
# falls back to ``f`` ("frequency-averaged").
_FC_METHOD_CB_LABEL = {
    "corr":      r"$\mathrm{{corr}}_{{ij}}^{{({band})}}$",
    "msc":       r"$\langle \mathrm{{msc}}_{{ij}} \rangle_{{{band}}}$",
    "imcoh_abs": r"$\langle |\mathrm{{imcoh}}_{{ij}}| \rangle_{{{band}}}$",
    "imcoh_sq":  r"$\langle |\mathrm{{imcoh}}_{{ij}}|^{{2}} \rangle_{{{band}}}$",
}


def fc_method_colorbar_label(
    fc_method: Optional[str],
    band: Optional[str] = None,
) -> str:
    """Canonical LaTeX colorbar label for an FC method (and band), or ''.

    Format is ``<|imcoh_{ij}|>_{<band>}`` and analogues.  When
    ``band`` is a recognised key (``alpha``, ``beta``, …), the band
    glyph from ``BRAIN_BAND_TEX_DICT`` is stripped of its ``$`` and
    inserted as the subscript; otherwise the literal band string is
    used (or ``f`` if ``band is None``).
    """
    template = _FC_METHOD_CB_LABEL.get(fc_method or "", "")
    if not template:
        return ""
    if band is None:
        glyph = "f"
    else:
        tex = BRAIN_BAND_TEX_DICT.get(band, band)
        glyph = tex.strip("$") if isinstance(tex, str) else str(tex)
    return template.format(band=glyph)


def probe_groups(
    channel_labels: Sequence[str],
) -> List[Tuple[str, int, int]]:
    """Return ``[(prefix, start, end_inclusive), ...]`` for runs of
    consecutive labels that share the same leading-letter prefix.

    The contact label may contain a comma-separated annotation
    (``"A 1,G2"``); only the substring before the first comma is
    parsed.  The probe prefix is the leading letters before any
    digit/space.
    """
    groups: List[Tuple[str, int, int]] = []
    if not channel_labels:
        return groups
    current = None
    start = 0
    for i, raw in enumerate(channel_labels):
        head = str(raw).split(",", 1)[0].strip()
        m = _PROBE_RE.match(head)
        prefix = m.group(1) if m else head
        if current is None:
            current, start = prefix, i
        elif prefix != current:
            groups.append((current, start, i - 1))
            current, start = prefix, i
    groups.append((current, start, len(channel_labels) - 1))
    return groups


def _is_first_col(ax) -> bool:
    """Return True when the axis sits in the first column of its gridspec.

    Used to suppress the y-axis label and y-tick labels on non-leftmost
    panels in row/grid layouts.  Falls back to True for non-grid axes.
    """
    try:
        return ax.get_subplotspec().is_first_col()
    except Exception:
        return True


def _is_last_row(ax) -> bool:
    """Return True when the axis sits in the last row of its gridspec.

    Used to suppress the x-axis label and x-tick labels on non-bottom
    rows in grid layouts.  Falls back to True for non-grid axes
    (single panels, single rows where every panel is "bottom").
    """
    try:
        return ax.get_subplotspec().is_last_row()
    except Exception:
        return True


def _apply_tick_labels(
    ax,
    *,
    mode: TickMode,
    channel_labels: Optional[Sequence[str]],
    is_first_col: Optional[bool] = None,
    is_last_row: Optional[bool] = None,
) -> None:
    """Apply tick policy.  ``is_first_col`` / ``is_last_row`` override
    the gridspec auto-detection — required when the gridspec contains
    extra rows/cols dedicated to colorbars (those would shift the
    "last row" past the last panel row otherwise).
    """
    if mode == "generic":
        ax.set_xticks([])
        ax.set_yticks([])
        return
    is_leftmost = _is_first_col(ax) if is_first_col is None else is_first_col
    is_bottom = _is_last_row(ax) if is_last_row is None else is_last_row
    if mode == "index":
        if is_bottom:
            ax.set_xlabel(r"$i$", fontsize=10)
        else:
            ax.tick_params(axis="x", which="both",
                           bottom=False, labelbottom=False)
        if is_leftmost:
            ax.set_ylabel(r"$j$", fontsize=10)
        else:
            ax.tick_params(axis="y", which="both",
                           left=False, labelleft=False)
        ax.tick_params(labelsize=6)
        return
    if mode == "chnames":
        if channel_labels is None:
            raise ValueError(
                "tick_labels='chnames' requires channel_labels=<sequence>."
            )
        groups = probe_groups(channel_labels)
        major = [(s + e) / 2 for _, s, e in groups]
        labels = [p for p, _, _ in groups]
        boundaries = [s - 0.5 for _, s, _ in groups[1:]]

        # X ticks/labels only on the bottom row of the gridspec.
        if is_bottom:
            ax.set_xticks(major)
            ax.set_xticklabels(labels, fontsize=7)
            ax.set_xticks(boundaries, minor=True)
            ax.set_xlabel("probe", fontsize=8)
        else:
            ax.set_xticks([])
            ax.tick_params(axis="x", which="both",
                           bottom=False, labelbottom=False)

        # Y ticks/labels only on the leftmost col.  Explicit
        # ``set_yticklabels`` overrides ``sharey``, so we must gate.
        if is_leftmost:
            ax.set_yticks(major)
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_yticks(boundaries, minor=True)
            ax.set_ylabel("probe", fontsize=8)
        else:
            ax.set_yticks([])
            ax.tick_params(axis="y", which="both",
                           left=False, labelleft=False)

        ax.tick_params(which="major", length=0)
        ax.tick_params(which="minor", length=2.5, color="0.5", width=0.6)
        return
    raise ValueError(f"unknown tick_labels mode {mode!r}")


def plot_fc_adjacency(
    M: np.ndarray,
    *,
    ax=None,
    cmap: str = "magma",
    vmin: Optional[float] = 0.0,
    vmax: Optional[float] = None,
    zero_diagonal: bool = True,
    log_scale: bool = False,
    tick_labels: TickMode = "generic",
    channel_labels: Optional[Sequence[str]] = None,
    colorbar: bool = True,
    colorbar_label: Optional[str] = None,
    fc_method: Optional[str] = None,
    band: Optional[str] = None,
    figsize: Tuple[float, float] = (3.4, 3.0),
    is_first_col: Optional[bool] = None,
    is_last_row: Optional[bool] = None,
):
    """Plot a single FC adjacency matrix.

    Compact, options-driven.  Returns ``(fig, ax, im)`` for further
    customisation by the caller (titles, watermarks, etc. are NOT
    added here — keep the helper minimal).

    Parameters
    ----------
    M : (N, N) array
        FC matrix.  Diagonal is zeroed in-place on a copy when
        ``zero_diagonal=True`` (default).
    ax : Axes, optional
        Existing axis to draw into.  If ``None``, a new
        ``(figsize)`` figure is created.
    cmap, vmin, vmax : standard ``imshow`` controls.
        ``vmax=None`` → ``np.nanmax(M)`` after diagonal zeroing.
    zero_diagonal : bool, default True
        Force the diagonal to zero before plotting (FC self-coherence
        is uninformative and dominates the colour scale otherwise).
    log_scale : bool, default False
        Use a logarithmic colour scale (``matplotlib.colors.LogNorm``).
        Non-positive entries are masked to NaN and rendered as the
        colormap's "bad" colour.  Useful when most off-diagonal
        values cluster near zero so the linear scale looks black.
    tick_labels : {"generic", "index", "chnames"}, default "generic"
        - ``"generic"``: no ticks; axes labelled as math ``$i$``,
          ``$j$``.
        - ``"index"``: matplotlib's auto numeric ticks.
        - ``"chnames"``: one tick per probe family at the family
          midpoint, with light minor ticks at probe boundaries;
          requires ``channel_labels``.
    channel_labels : sequence of str, optional
        Required for ``tick_labels="chnames"``.  Each entry is a
        contact label like ``"A 1,G2"``; the leading letters define
        the probe family.
    colorbar : bool, default True
        Add a right-side colorbar via ``imshow_colorbar_caxdivider``.
    colorbar_label : str, optional
        Explicit label.  Overrides the auto-detected ``fc_method``
        label.
    fc_method : str, optional
        If given (one of ``corr | msc | imcoh_abs | imcoh_sq``),
        sets a canonical LaTeX colorbar label, parameterised by
        ``band`` (e.g. ``<|imcoh_{ij}|>_{β}``).
    band : str, optional
        Band key (``alpha``, ``beta``, …).  Subscripts the colorbar
        label with the LaTeX glyph from ``BRAIN_BAND_TEX_DICT``.
        Falls back to ``f`` when omitted.
    figsize : (w, h), default ``(3.4, 3.0)``
        Used only when ``ax is None``.
    """
    M_view = np.asarray(M, dtype=float)
    M_plot = M_view.copy() if zero_diagonal else M_view.copy()
    if zero_diagonal:
        np.fill_diagonal(M_plot, 0.0)
    if vmax is None:
        vmax = float(np.nanmax(M_plot)) if np.any(np.isfinite(M_plot)) else 1.0

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # interpolation="nearest" — FC matrices are discrete (i, j) cells;
    # antialiasing blurs probe-block boundaries.
    if log_scale:
        M_plot = np.where(M_plot > 0, M_plot, np.nan)
        positive = M_plot[np.isfinite(M_plot)]
        log_vmin = float(positive.min()) if positive.size else 1e-6
        if vmin is not None and vmin > 0:
            log_vmin = max(log_vmin, vmin)
        norm = LogNorm(vmin=log_vmin, vmax=vmax)
        im = ax.imshow(M_plot, cmap=cmap, norm=norm,
                       aspect="equal", interpolation="nearest")
    else:
        im = ax.imshow(M_plot, cmap=cmap, vmin=vmin, vmax=vmax,
                       aspect="equal", interpolation="nearest")

    _apply_tick_labels(
        ax, mode=tick_labels, channel_labels=channel_labels,
        is_first_col=is_first_col, is_last_row=is_last_row,
    )

    if colorbar:
        _, _, clb = imshow_colorbar_caxdivider(im, ax, size="4%", pad=0.05)
        if colorbar_label is None:
            colorbar_label = fc_method_colorbar_label(fc_method, band=band)
        if colorbar_label:
            clb.set_label(colorbar_label, fontsize=8)
        clb.ax.tick_params(labelsize=6)

    return fig, ax, im


def plot_fc_adjacency_row(
    matrices: Sequence[np.ndarray],
    *,
    titles: Optional[Sequence[str]] = None,
    fc_method: Optional[str] = None,
    band: Optional[str] = None,
    bands: Optional[Sequence[str]] = None,
    tick_labels: TickMode = "generic",
    channel_labels: Optional[Sequence[str]] = None,
    log_scale: bool = False,
    shared_scale: bool = True,
    cmap: str = "magma",
    panel_size: float = 2.4,
    fig=None,
    axes=None,
):
    """Plot a row of FC adjacency matrices.

    Two layout modes:

    - ``shared_scale=True`` (default): one ``(vmin, vmax)`` across the
      whole row + a single colorbar at the right of the rightmost
      axis.  Use this when the matrices are directly comparable
      (e.g. one patient × one band × four phases).
    - ``shared_scale=False``: each panel has its own colour scale and
      its own colorbar.  Use this when the matrices have different
      natural amplitudes (e.g. one patient × one phase × six bands —
      γ matrices typically smaller than α/β).

    Parameters
    ----------
    matrices : sequence of (N, N) array
        One per panel, in display order.
    titles : sequence of str, optional
        One per panel (TeX OK).  Placed via ``ax.set_title``.
    fc_method, band, bands : str / sequence of str, optional
        ``band`` for a row that shares a band; ``bands`` for a row
        whose panels iterate over bands (one per panel).  Either
        feeds the canonical colorbar label.
    tick_labels, channel_labels : as in :func:`plot_fc_adjacency`.
        Tick labels are emitted on every panel; ``sharey`` suppresses
        repeated y-tick labels in the row.
    log_scale, cmap : as in :func:`plot_fc_adjacency`.
    shared_scale : bool, default True.
    panel_size : float, default 2.4
        Width and height of each panel in inches.
    fig, axes : optional matplotlib figure / axes array
        If both are given, the row is drawn into them; otherwise a
        new ``(1, n)`` figure is created with ``sharey=True``.

    Returns
    -------
    fig, axes
    """
    n = len(matrices)
    if n == 0:
        raise ValueError("matrices is empty")

    if fig is None or axes is None:
        # In shared-scale mode the right-most panel hosts the single
        # colorbar via ``imshow_colorbar_caxdivider`` (which carves
        # ``cb_size + cb_pad`` out of its host axis).  We pre-pad the
        # last gridspec column by exactly that fraction so the imshow
        # portion of the last panel matches the others after the cax
        # steal.  In per-panel mode every panel surrenders the same
        # fraction, so widths stay uniform with equal ratios.
        if shared_scale:
            width_ratios = [1.0] * (n - 1) + [1.0 + _ROW_CB_FRAC]
            fig_w = panel_size * (n - 1 + (1.0 + _ROW_CB_FRAC)) + 0.6
        else:
            width_ratios = [1.0] * n
            fig_w = panel_size * n + 0.6
        fig, axarr = plt.subplots(
            1, n,
            figsize=(fig_w, panel_size + 0.4),
            sharey=True, squeeze=False,
            gridspec_kw={"width_ratios": width_ratios},
        )
        axes = axarr[0]

    if shared_scale:
        flat = np.concatenate(
            [np.asarray(m, dtype=float)[np.isfinite(m)].ravel()
             for m in matrices]
        )
        row_vmax = float(np.nanmax(flat)) if flat.size else 1.0
    else:
        row_vmax = None

    last_idx = n - 1
    last_im = None
    for i, M in enumerate(matrices):
        cur_band = bands[i] if bands is not None else band
        if shared_scale:
            cb_label_i: Optional[str] = ""
            panel_colorbar = False
        else:
            # Per-panel mode: every cb carries its own band-specific
            # label (the glyph for that panel's band; if all panels
            # share a band, the same fixed glyph repeats).  This was
            # previously only the rightmost cb with the generic ``_f``
            # glyph — which is wrong when the row spans bands and a
            # reader wants to know which scale they're looking at.
            cb_label_i = fc_method_colorbar_label(fc_method, band=cur_band)
            panel_colorbar = True
        _, _, im = plot_fc_adjacency(
            M,
            ax=axes[i],
            cmap=cmap,
            vmin=0.0,
            vmax=row_vmax if shared_scale else None,
            log_scale=log_scale,
            tick_labels=tick_labels,
            channel_labels=channel_labels,
            colorbar=panel_colorbar,
            colorbar_label=cb_label_i,
            fc_method=fc_method,
            band=cur_band,
        )
        last_im = im
        if titles is not None:
            axes[i].set_title(titles[i], fontsize=10)

    if shared_scale and last_im is not None:
        # Single shared colorbar via ``imshow_colorbar_caxdivider``.
        # The cax + pad eats from the right-most panel's host axis;
        # the gridspec ``width_ratios`` above pre-padded that panel
        # by ``_ROW_CB_FRAC`` so the imshow stays the same width as
        # the others.  Label uses the row's ``band`` glyph when the
        # row shares a band; otherwise no label (panels disagree).
        _, _, clb = imshow_colorbar_caxdivider(
            last_im, axes[-1], size="4%", pad=0.06
        )
        cb_label = (
            fc_method_colorbar_label(fc_method, band=band)
            if bands is None else ""
        )
        if cb_label:
            clb.set_label(cb_label, fontsize=8)
        clb.ax.tick_params(labelsize=6)

    return fig, axes


_ColorbarMode = Literal["shared", "per_row", "per_col"]


def plot_fc_adjacency_grid(
    matrices: Sequence[Sequence[Optional[np.ndarray]]],
    *,
    row_titles: Optional[Sequence[str]] = None,
    col_titles: Optional[Sequence[str]] = None,
    fc_method: Optional[str] = None,
    band: Optional[str] = None,
    bands: Optional[Sequence[str]] = None,
    row_bands: Optional[Sequence[str]] = None,
    col_bands: Optional[Sequence[str]] = None,
    tick_labels: TickMode = "chnames",
    channel_labels: Optional[Sequence[str]] = None,
    channel_labels_per_row: Optional[Sequence[Sequence[str]]] = None,
    log_scale: bool = True,
    colorbar_mode: _ColorbarMode = "shared",
    cmap: str = "magma",
    panel_size: float = 2.0,
    wspace: float = 0.05,
    hspace: float = 0.08,
):
    """Plot a grid (mosaic) of FC adjacency matrices.

    ``matrices[i][j]`` is the FC matrix at row ``i``, column ``j``.
    ``None`` is allowed for empty cells (the panel is hidden).

    Three colorbar policies, selected per the row/col semantics:

    - ``"shared"`` — single full-height cb on the right; one global
      ``vmax``.  Use when every panel is directly comparable
      (e.g. patients × phases at a fixed band).
    - ``"per_row"`` — n cbs stacked on the right, one per row.  Each
      row's cb represents that row's ``vmax``.  Use when each row has
      its own natural amplitude (e.g. bands × phases at a fixed
      patient — γ rows would be crushed by a global scale).
    - ``"per_col"`` — n horizontal cbs along the bottom, one per
      column.  Each column's cb represents that column's ``vmax``.
      Use when each column has its own natural amplitude (e.g.
      patients × bands at a fixed phase — band columns disagree on
      amplitude across the cohort).

    The right (``shared`` / ``per_row``) or bottom (``per_col``) cb
    column / row is allocated as an extra cell in the gridspec so
    panels keep their full size — no ``make_axes_locatable`` shrink.

    ``tick_labels`` follows the same gating as
    :func:`plot_fc_adjacency_row`: in ``"chnames"`` / ``"index"``
    modes only the leftmost column gets y-ticks, only the bottom row
    gets x-ticks.  Pass ``channel_labels_per_row`` when rows contain
    different patients (each row uses its own probe layout).

    ``row_titles`` are emitted as rotated bold text to the left of
    the leftmost column; ``col_titles`` go on the top row's
    ``set_title``.  TeX OK in both.
    """
    n_rows = len(matrices)
    if n_rows == 0:
        raise ValueError("matrices is empty")
    n_cols = max(len(r) for r in matrices)
    if colorbar_mode not in ("shared", "per_row", "per_col"):
        raise ValueError(f"unknown colorbar_mode {colorbar_mode!r}")

    # ``bands`` is the legacy alias for ``col_bands`` (length = n_cols).
    if col_bands is None and bands is not None:
        col_bands = bands

    def _band_for_panel(i: int, j: int) -> Optional[str]:
        if col_bands is not None:
            return col_bands[j]
        if row_bands is not None:
            return row_bands[i]
        return band

    def _band_for_row(i: int) -> Optional[str]:
        if row_bands is not None:
            return row_bands[i]
        return band

    def _band_for_col(j: int) -> Optional[str]:
        if col_bands is not None:
            return col_bands[j]
        return band

    def _flatten(M):
        if M is None:
            return np.array([])
        a = np.asarray(M, dtype=float)
        return a[np.isfinite(a)].ravel()

    if colorbar_mode == "shared":
        flat = np.concatenate([
            _flatten(matrices[i][j])
            for i in range(n_rows)
            for j in range(min(n_cols, len(matrices[i])))
        ])
        gvmax = float(np.nanmax(flat)) if flat.size else 1.0
        get_vmax = lambda i, j: gvmax  # noqa: E731
    elif colorbar_mode == "per_row":
        row_vmax: List[float] = []
        for i in range(n_rows):
            flat = np.concatenate([
                _flatten(matrices[i][j])
                for j in range(min(n_cols, len(matrices[i])))
            ])
            row_vmax.append(float(np.nanmax(flat)) if flat.size else 1.0)
        get_vmax = lambda i, j: row_vmax[i]  # noqa: E731
    else:  # per_col
        col_vmax: List[float] = []
        for j in range(n_cols):
            flat = np.concatenate([
                _flatten(matrices[i][j])
                for i in range(n_rows)
                if j < len(matrices[i])
            ])
            col_vmax.append(float(np.nanmax(flat)) if flat.size else 1.0)
        get_vmax = lambda i, j: col_vmax[j]  # noqa: E731

    # Gridspec: extra row or col reserved for the colorbar(s).
    cb_frac = 0.06
    if colorbar_mode in ("shared", "per_row"):
        gs_rows, gs_cols = n_rows, n_cols + 1
        width_ratios = [1.0] * n_cols + [cb_frac]
        height_ratios = [1.0] * n_rows
    else:  # per_col
        gs_rows, gs_cols = n_rows + 1, n_cols
        width_ratios = [1.0] * n_cols
        height_ratios = [1.0] * n_rows + [cb_frac]

    fig_w = panel_size * sum(width_ratios) + 1.0  # row-title margin
    fig_h = panel_size * sum(height_ratios) + 0.7  # col-title margin
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        gs_rows, gs_cols,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
        wspace=wspace, hspace=hspace,
    )

    axarr: List[List] = []
    for i in range(n_rows):
        row_axes = []
        for j in range(n_cols):
            row_axes.append(fig.add_subplot(gs[i, j]))
        axarr.append(row_axes)

    last_im_per_row: List = [None] * n_rows
    last_im_per_col: List = [None] * n_cols
    last_im_global = None
    for i in range(n_rows):
        for j in range(n_cols):
            ax = axarr[i][j]
            M = matrices[i][j] if j < len(matrices[i]) else None
            if M is None:
                ax.set_visible(False)
                continue
            cl = (
                channel_labels_per_row[i]
                if channel_labels_per_row is not None
                else channel_labels
            )
            cur_band = _band_for_panel(i, j)
            # Panel-position overrides — the gridspec's last row/last
            # col may be a colorbar axis (per_col / per_row mode), so
            # ``ax.get_subplotspec().is_last_row()`` would lie about the
            # last *panel* row.  Pass the panel-relative position
            # explicitly.
            _, _, im = plot_fc_adjacency(
                M,
                ax=ax, cmap=cmap,
                vmin=0.0, vmax=get_vmax(i, j),
                log_scale=log_scale,
                tick_labels=tick_labels,
                channel_labels=cl,
                colorbar=False,
                fc_method=fc_method,
                band=cur_band,
                is_first_col=(j == 0),
                is_last_row=(i == n_rows - 1),
            )
            last_im_per_row[i] = im
            last_im_per_col[j] = im
            last_im_global = im

    if col_titles:
        for j, ct in enumerate(col_titles):
            if j < n_cols:
                axarr[0][j].set_title(ct, fontsize=11)

    if row_titles:
        # Offset the row title far enough left to clear any y-axis
        # tick labels and "probe"/"$j$" axis label on the leftmost panel.
        if tick_labels == "chnames":
            offset_frac = -0.45
        elif tick_labels == "index":
            offset_frac = -0.30
        else:
            offset_frac = -0.10
        for i, rt in enumerate(row_titles):
            if i < n_rows:
                ax = axarr[i][0]
                ax.text(
                    offset_frac, 0.5, rt, transform=ax.transAxes,
                    rotation=90, fontsize=11, fontweight="bold",
                    ha="center", va="center",
                )

    # Colorbars — placed in the dedicated gridspec cells via
    # ``plt.colorbar(im, cax=...)``.  Labelling rule:
    #   - shared mode: one full-height cb, single label using the
    #     fixed ``band`` if all panels share a band, else ``_f``.
    #   - per_row mode: each row's cb gets its own label (the row's
    #     band-specific glyph if rows iterate over bands; otherwise
    #     the fixed band glyph; otherwise ``_f``).
    #   - per_col mode: each col's cb gets its own label (analogous).
    if colorbar_mode == "shared" and last_im_global is not None:
        cb_ax = fig.add_subplot(gs[:, -1])
        clb = plt.colorbar(last_im_global, cax=cb_ax)
        cb_label = fc_method_colorbar_label(fc_method, band=band)
        if cb_label:
            clb.set_label(cb_label, fontsize=8)
        clb.ax.tick_params(labelsize=6)
    elif colorbar_mode == "per_row":
        for i in range(n_rows):
            if last_im_per_row[i] is None:
                continue
            cb_ax = fig.add_subplot(gs[i, -1])
            clb = plt.colorbar(last_im_per_row[i], cax=cb_ax)
            clb.ax.tick_params(labelsize=6)
            lbl = fc_method_colorbar_label(fc_method, band=_band_for_row(i))
            if lbl:
                clb.set_label(lbl, fontsize=8)
    elif colorbar_mode == "per_col":
        for j in range(n_cols):
            if last_im_per_col[j] is None:
                continue
            cb_ax = fig.add_subplot(gs[-1, j])
            clb = plt.colorbar(
                last_im_per_col[j], cax=cb_ax, orientation="horizontal",
            )
            clb.ax.tick_params(labelsize=6)
            lbl = fc_method_colorbar_label(fc_method, band=_band_for_col(j))
            if lbl:
                clb.set_label(lbl, fontsize=8)

    return fig, axarr
