#!/usr/bin/env python3
r"""Compound assembly of Results §1 figures (fig:trace1, fig:trace2) directly in
matplotlib, replacing the LaTeX minipage assembly.

Each panel is drawn by its OWN script's ``draw(subfig, panel_label=...)`` function
(single source of truth — no forked panel code); this module only tiles those panels
into ``fig.subfigures`` whose width/height ratios follow the panels' native aspect
ratios, so the assembly reads as a proportional mosaic.

Panel g (θ|α|β 3D brains) is Plotly and has no matplotlib representation, so it is
embedded as its rendered high-res PNG via ``imshow`` — the one raster tile; every
other tile stays fully vector, exactly as it is standalone.

The compound canvas is sized at the printed figure width (``PAGE_W_IN``) so the
absolute point sizes the panels set land on the page unscaled — LaTeX then only
fits width to ``\linewidth`` (aspect preserved). Panel letters are owned here, so a
tile's source file name need not match its panel letter.

Writes:
  data/preprint/figures/results_section1/fig_trace1.pdf
  data/preprint/figures/results_section1/fig_trace2.pdf
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# sibling panel scripts (figures_embedded is not a package) importable by name
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lrg_eegfc.utils.scripting import setup_script_env  # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style  # noqa: E402

import matplotlib.text as mtext  # noqa: E402
from PIL import Image  # noqa: E402

import fig_trace_a_band_forest as p_violin  # writes fig_trace_a_band_violin.pdf  # noqa: E402
import fig_trace_b_ofc_localization as p_ofc  # noqa: E402
import fig_trace_c_diffmaps as p_diffmaps  # noqa: E402
import fig_trace_d_tissue_class as p_tissue  # noqa: E402
import fig_trace_e_reinstatement as p_reinst  # noqa: E402
import fig_trace_f_soz_divergence as p_soz  # noqa: E402
import fig_trace_perpatient_bands_forest as p_perpat  # noqa: E402

ROOT = setup_script_env()
use_lrg_style()

OUTDIR = ROOT / "data/preprint/figures/results_section1"
G_PNG = OUTDIR / "fig_trace_g_pairglow_3d.png"
G_PNG_VERTICAL = OUTDIR / "fig_trace_g_pairglow_3d_vertical.png"

# native standalone figure WIDTHS (inches) — used to shrink each panel's absolute
# fonts by (tile_width / native_width) so a tile reads like a scaled-down standalone
NATIVE_W = dict(violin=8.2, ofc=10.2, diffmaps=18.12, tissue=7.6, reinst=12.2,
                soz=9.6, perpat=13.2)


def _scale_fonts(sfig, factor):
    """Multiply every text point-size inside a SubFigure (preserves internal hierarchy)."""
    for t in sfig.findobj(mtext.Text):
        t.set_fontsize(t.get_fontsize() * factor)


def _scale_axes_fonts(axes, factor):
    """Same, for a list of flat Axes (fig2 has no SubFigures — nilearn forbids them)."""
    seen = set()
    for ax in axes:
        for t in ax.findobj(mtext.Text):
            if id(t) not in seen:
                seen.add(id(t))
                t.set_fontsize(t.get_fontsize() * factor)


def _drop_legends(sfig):
    for lg in list(sfig.legends):
        lg.remove()


def _remove_letter(sfig, ch):
    """Delete a panel's own bold letter (letters are re-added uniformly by the assembler)."""
    tgt = r"$\mathbf{%s}$" % ch
    for t in sfig.findobj(mtext.Text):
        if t.get_text() == tgt:
            t.remove()
            return True
    return False


LETTER_SIZE = 12


def _letter(sfig, ch):
    """One uniform bold panel letter at the tile's top-left corner."""
    sfig.text(0.008, 0.992, ch, fontsize=LETTER_SIZE, fontweight="bold",
              va="top", ha="left")

# printed width of a figure* (both-column span); aspect is what actually ships
PAGE_W_IN = 7.2

# native aspect ratios (w/h) measured from the standalone panel PDFs
AR = dict(g=2.757, violin=1.122, diffmaps=3.009, reinst=1.637,
          ofc=2.465, tissue=1.468, soz=1.269)

def _imshow_raster(subfig, png_path):
    """Full-bleed raster tile (the Plotly 3D panel g; its own letter/titles are baked in)."""
    ax = subfig.add_axes([0, 0, 1, 1])
    ax.imshow(plt.imread(str(png_path)))
    ax.axis("off")
    ax.set_facecolor("none")


def _row1_split(ar_left, ar_right):
    """Width fraction for the left tile so both tiles share one row height."""
    # w_l/ar_l == (1-w_l)/ar_r  ->  w_l = ar_l / (ar_l + ar_r)
    return ar_left / (ar_left + ar_right)


def _optimize_pdf(path: Path) -> None:
    """Losslessly shrink a matplotlib vector PDF via Ghostscript re-emission.

    matplotlib emits the panel-b diffmaps as ~170k individual ``pcolormesh`` quads
    (fully vector, per the no-raster rule), which it compresses poorly; Ghostscript
    re-encodes the content streams far more compactly while PRESERVING every raster
    at full resolution (no downsampling, Flate/lossless) and keeping all vector
    artists vector. Typically ~halves fig_trace1.pdf (4.3 -> 2.3 MB) with no visible
    change. No-op (keeps the original) if Ghostscript is unavailable or yields no gain.
    """
    gs = shutil.which("gs")
    if gs is None:
        print(f"  [optimize] ghostscript not found; kept {path.name} unoptimized")
        return
    tmp = path.with_suffix(".opt.pdf")
    before = path.stat().st_size
    cmd = [gs, "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.5", "-dPDFSETTINGS=/prepress",
           "-dAutoFilterColorImages=false", "-dColorImageFilter=/FlateEncode",
           "-dAutoFilterGrayImages=false", "-dGrayImageFilter=/FlateEncode",
           "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false",
           "-dDownsampleMonoImages=false", "-dNOPAUSE", "-dBATCH", "-dQUIET",
           f"-sOutputFile={tmp}", str(path)]
    try:
        subprocess.run(cmd, check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"  [optimize] ghostscript failed ({exc}); kept {path.name} unoptimized")
        tmp.unlink(missing_ok=True)
        return
    if tmp.exists() and tmp.stat().st_size < before:
        tmp.replace(path)
        print(f"  [optimize] {path.name}: {before/1e6:.2f} -> {path.stat().st_size/1e6:.2f} MB (lossless gs)")
    else:
        tmp.unlink(missing_ok=True)
        print(f"  [optimize] {path.name}: no gain, kept original")


def build_fig1():
    r"""Mosaic (super-grid 2x3; width_ratios 3:5:3, height_ratios 2:3):
          A A A  B B B B B  D D D
          A A A  B B B B B  D D D
          A A A  C C C C C  E E E
          A A A  C C C C C  E E E
          A A A  C C C C C  E E E
        A = g brain spectrum (vertical raster)     -> panel a
        B = cophenetic fingerprint maps (diffmaps) -> panel b
        C = post-task reinstatement                -> panel c
        D = per-patient held-trace violin          -> panel d
        E = per-patient / per-band forest          -> panel e
    """
    gw, gh = Image.open(G_PNG_VERTICAL).size
    # row-0 (b diffmaps / d violin) vs row-1 (c reinstatement / e forest). row-0 shrunk and
    # row-1 grown so c's taller tile fills the old white band between b and c (b's square maps
    # need less height than the old 2:3 split gave them).
    HR = [1.6, 3.2]
    # column A holds the vertical brain raster; its width is set so the tile aspect equals
    # the raster aspect exactly (wA/H == gw/gh) -> no side whitespace, whatever g's shape.
    WR = [(gw / gh) * sum(HR), 4.5, 2.24]
    W = 13.0
    H = W * sum(HR) / sum(WR)                    # square grid units -> true tile aspects
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(2, 3, width_ratios=WR, height_ratios=HR,
                          left=0.0, right=1.0, top=1.0, bottom=0.0,
                          wspace=0.0, hspace=0.0)
    sfA = fig.add_subfigure(gs[0:2, 0])      # left, full height  -> a
    sfB = fig.add_subfigure(gs[0, 1])        # top-center         -> b
    sfC = fig.add_subfigure(gs[1, 1])        # center             -> c
    sfD = fig.add_subfigure(gs[0, 2])        # top-right          -> d
    sfE = fig.add_subfigure(gs[1, 2])        # bottom-right       -> e

    wB = WR[1] / sum(WR) * W                    # center-column tile width (in)
    wD = WR[2] / sum(WR) * W                    # right-column tile width (in)

    _imshow_raster(sfA, G_PNG_VERTICAL);                                          _letter(sfA, "a")
    p_diffmaps.draw(sfB); _remove_letter(sfB, "c"); _scale_fonts(sfB, wB / NATIVE_W["diffmaps"]); _letter(sfB, "b")
    p_reinst.draw(sfC);   _remove_letter(sfC, "d"); _scale_fonts(sfC, wB / NATIVE_W["reinst"]);   _letter(sfC, "c")
    # d, e author their own text/markers at tile size (dense panels — a faithful shrink would
    # make labels illegible), so they are NOT font-rescaled here.
    p_violin.draw(sfD);   _remove_letter(sfD, "b"); _drop_legends(sfD);           _letter(sfD, "d")
    p_perpat.draw(sfE);                                                           _letter(sfE, "e")

    out = OUTDIR / "fig_trace1.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    _optimize_pdf(out)
    return out


def build_fig2():
    r"""Flat compound — a nilearn glass brain FORBIDS SubFigures (nilearn resolves
    ``figure=`` to the ROOT figure and lays the brain in root coordinates, so every
    SubFigure frame is ignored and panels sprawl / brains stack). All tiles are therefore
    flat ``fig.add_axes`` rects in FIGURE coordinates:

        a  OFC localization glass brains          (top band, full width)
        b  tissue-class carrier forest            (bottom-left)
        c  SOZ alpha/beta divergence              (bottom-right: 2 scatters + cohort forest)

    b and c are authored for wide standalones, so their fonts are scaled down to the tile.
    The band palette + clears/does-not-clear key is one shared legend along the bottom.
    """
    from matplotlib.lines import Line2D
    from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
    from lrg_eegfc.visuals.styles import band_color

    W, H = 7.2, 5.7
    fig = plt.figure(figsize=(W, H))

    # ---------- a : glass brains (top band, full width) ----------
    p_ofc.draw_brain(fig, rect=(0.03, 0.665, 0.955, 0.335))
    fig.text(0.008, 0.992, "a", fontsize=LETTER_SIZE, fontweight="bold", va="top", ha="left")
    fig.legend(handles=p_ofc.legend_handles(compact=True), loc="upper center",
               bbox_to_anchor=(0.5, 0.648), ncol=3, frameon=False, fontsize=8.4,
               handletextpad=0.4, columnspacing=1.3)

    # ---------- b : tissue-class carrier forest (bottom-left) ----------
    # the two-line "non-SOZ <-> non-SOZ" row label lets the plot grow leftward
    axB = fig.add_axes([0.175, 0.150, 0.275, 0.400])
    p_tissue.build_panel(axB, draw_letter=False)
    _scale_axes_fonts([axB], 0.60)
    fig.text(0.010, 0.565, "b", fontsize=LETTER_SIZE, fontweight="bold", va="top", ha="left")

    # ---------- c : SOZ alpha/beta divergence (bottom-right) ----------
    axC0 = fig.add_axes([0.575, 0.295, 0.185, 0.205])          # alpha scatter
    axC1 = fig.add_axes([0.800, 0.295, 0.185, 0.205])          # beta scatter
    axCf = fig.add_axes([0.595, 0.100, 0.390, 0.145])          # cohort forest (spans)
    p_soz.scatter_panel(axC0, "alpha", tag=None,
                        xlabel=r"rest$\rightarrow$task (rank)",
                        ylabel=r"rest$\rightarrow$post (rank)")
    p_soz.scatter_panel(axC1, "beta", tag=None,
                        xlabel=r"rest$\rightarrow$task (rank)", ylabel="")
    p_soz.forest_panel(axCf)                                   # symbols sit right of the row labels
    _scale_axes_fonts([axC0, axC1, axCf], 0.60)
    fig.text(0.490, 0.565, "c", fontsize=LETTER_SIZE, fontweight="bold", va="top", ha="left")
    soz_h = [
        Line2D([0], [0], marker="o", ls="", mfc=p_soz.C_SOZ, mec="white", ms=7, label="SOZ–SOZ"),
        Line2D([0], [0], marker="s", ls="", mfc=p_soz.C_SOZ, mec=p_soz.C_SOZ, ms=8,
               alpha=0.32, label=r"SOZ extent (1$\sigma$)"),
        Line2D([0], [0], marker="o", ls="", mfc=p_soz.C_SEA, mec="none", ms=7, label="non-SOZ"),
        Line2D([0], [0], ls="--", color=p_soz.C_DIAG, lw=1.8, label="held diagonal"),
    ]
    fig.legend(handles=soz_h, loc="upper center", bbox_to_anchor=(0.778, 0.578),
               ncol=2, frameon=False, fontsize=7.6, handletextpad=0.5, columnspacing=1.2)

    # ---------- shared band + style legend (very bottom, serves b and c) ----------
    band_h = [Line2D([0], [0], marker="o", ls="", mfc=band_color(b), mec="black", ms=8,
                     label=BRAIN_BAND_TEX_DICT.get(b, b)) for b in p_tissue.BANDS_D]
    style_h = [
        Line2D([0], [0], marker="o", ls="", mfc="0.4", mec="black", ms=9,
               label="clears matched-strength null"),
        Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.4", mew=1.4, ms=8,
               label="does not clear"),
    ]
    fig.legend(handles=band_h + style_h, loc="lower center", bbox_to_anchor=(0.5, 0.008),
               ncol=8, frameon=False, fontsize=8, handletextpad=0.4, columnspacing=1.0)

    out = OUTDIR / "fig_trace2.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    # fig2 is ~180 KB (no vector-quad bloat); left un-gs'd (optimizer is fig1-only).
    return out


def main():
    o1 = build_fig1()
    o2 = build_fig2()
    print(f"wrote {o1}\nwrote {o2}")


if __name__ == "__main__":
    main()
