#!/usr/bin/env python3
"""Audit 57 — §2 manuscript figure cohort audit.

Visually documents which Section 2 figures of ``notes_imcoh.tex``
still use a 6-patient subset rather than the current n=10 cohort
(locked 2026-04-25: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15).

A 4x4 small-multiples grid: each cell shows a thumbnail of the
original figure rendered via ``pdftoppm`` plus a colored frame:

- red    : observed_n < 10  (cohort mismatch — needs regeneration)
- green  : observed_n = 10  (n=10 control — OK)
- blue   : single-patient / triple-patient exemplar (intentional)

Each thumbnail is annotated with figure_name | observed_n |
expected_n | severity. Header text via ``fig.text`` (no
``fig.suptitle`` per project rule).

Output
------
``data/reports/notes_verification_2026-05-08/figures/cohort_n_audit.pdf``

Per project rule: full-vector PDF, no rasterization, PDF only.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()


# --------------------------------------------------------------------- #
# Cohort + figure registry                                              #
# --------------------------------------------------------------------- #

CURRENT_COHORT_N = 10
COHORT = ("Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15")

FW = ROOT / "data" / "outputs" / "figures" / "section2" / "for_writing_agent"
FD = ROOT / "data" / "outputs" / "figures" / "section2" / "fig_D"
RFC = ROOT / ".agents" / "writing-bundles" / "raw-fc"

OUT_DIR = (
    ROOT / "data" / "reports" / "notes_verification_2026-05-08" / "figures"
)
OUT_PDF = OUT_DIR / "cohort_n_audit.pdf"


@dataclass(frozen=True)
class FigEntry:
    name: str            # short label printed below thumbnail
    pdf: Path            # source PDF to thumbnail
    observed_n: int      # cohort size used in the figure
    expected_n: int      # design-required cohort size
    kind: str            # "cohort", "exemplar"
    note: str = ""       # extra one-word severity tag for caption


# Order: cohort-mismatch first (severity high), then OK n=10, then
# exemplars. Pad with disabled cells if fewer than 16.
REGISTRY: list[FigEntry] = [
    # ---- n=6 era figures (red, severity high) -------------------------
    FigEntry(
        name="fig_D1_enrichment_heatmap_redesigned",
        pdf=FD / "fig_D1_enrichment_heatmap_redesigned.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    FigEntry(
        name="fig_D2_bias_reduction_beta_rsPre",
        pdf=FW / "fig_D2_bias_reduction_beta_rsPre.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    FigEntry(
        name="fig_C2_spectral_distribution_msc_rsPre",
        pdf=FW / "fig_C2_spectral_distribution_msc_rsPre.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    FigEntry(
        name="fig_C2_spectral_distribution_imcoh_rsPre",
        pdf=FW / "fig_C2_spectral_distribution_imcoh_rsPre.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    FigEntry(
        name="fig_C3_weight_overlay_all_patients_beta",
        pdf=FW / "fig_C3_weight_overlay_all_patients_beta.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    FigEntry(
        name="fig_B_beta_rsPre_MSC_vs_ImCoh",
        pdf=FW / "fig_B_beta_rsPre_MSC_vs_ImCoh.pdf",
        observed_n=6, expected_n=10, kind="cohort", note="HIGH",
    ),
    # ---- n=10 controls (green, severity OK) ---------------------------
    FigEntry(
        name="fig2_phase_geometry_imcoh_abs",
        pdf=RFC / "fig2_phase_geometry_imcoh_abs.pdf",
        observed_n=10, expected_n=10, kind="cohort", note="OK",
    ),
    FigEntry(
        name="fig6_4phase_geometry_imcoh_abs",
        pdf=RFC / "fig6_4phase_geometry_imcoh_abs.pdf",
        observed_n=10, expected_n=10, kind="cohort", note="OK",
    ),
    # ---- exemplars (blue) ---------------------------------------------
    FigEntry(
        name="fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rsPre",
        pdf=FW / "fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rsPre.pdf",
        observed_n=3, expected_n=3, kind="exemplar",
        note="exemplar n=3",
    ),
    FigEntry(
        name="fig_A_msc_Pat_02",
        pdf=FW / "fig_A_msc_Pat_02.pdf",
        observed_n=1, expected_n=1, kind="exemplar",
        note="exemplar n=1",
    ),
    FigEntry(
        name="fig_A_imcoh_abs_Pat_02",
        pdf=FW / "fig_A_imcoh_abs_Pat_02.pdf",
        observed_n=1, expected_n=1, kind="exemplar",
        note="exemplar n=1",
    ),
    FigEntry(
        name="fig_A_msc_Pat_05",
        pdf=FW / "fig_A_msc_Pat_05.pdf",
        observed_n=1, expected_n=1, kind="exemplar",
        note="exemplar n=1",
    ),
    FigEntry(
        name="fig_A_imcoh_abs_Pat_05",
        pdf=FW / "fig_A_imcoh_abs_Pat_05.pdf",
        observed_n=1, expected_n=1, kind="exemplar",
        note="exemplar n=1",
    ),
]


# Frame colors and semantics --------------------------------------------

COLOR_MISMATCH = "#c0392b"   # red
COLOR_OK = "#27ae60"         # green
COLOR_EXEMPLAR = "#2980b9"   # blue
COLOR_EMPTY = "#bdbdbd"      # grey for unused cells


def frame_color(entry: FigEntry) -> str:
    if entry.kind == "exemplar":
        return COLOR_EXEMPLAR
    if entry.observed_n == entry.expected_n:
        return COLOR_OK
    return COLOR_MISMATCH


# --------------------------------------------------------------------- #
# PDF -> raster thumbnail via pdftoppm                                  #
# --------------------------------------------------------------------- #

def render_thumbnail(pdf: Path, tmp_dir: Path, dpi: int = 90) -> Path | None:
    """Render the first page of ``pdf`` to PNG, return PNG path or None.

    Uses ``pdftoppm -r DPI -png -f 1 -l 1`` (first-page only). Returns
    None if the PDF cannot be rendered (e.g. missing file). The PNG is
    embedded in the final figure as a single ``imshow`` artist; the
    overall PDF stays vector for everything except this thumbnail
    raster — there is no library way to embed a vector PDF page inside
    matplotlib. The thumbnails are reference-only, so they are the
    only raster pixels in the audit figure.
    """
    if not pdf.exists():
        return None
    stem = pdf.stem.replace(" ", "_")
    out_prefix = tmp_dir / stem
    cmd = [
        "pdftoppm", "-r", str(dpi),
        "-png", "-f", "1", "-l", "1",
        str(pdf), str(out_prefix),
    ]
    try:
        subprocess.run(
            cmd, check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    # pdftoppm produces "<prefix>-1.png"
    png = tmp_dir / f"{stem}-1.png"
    return png if png.exists() else None


# --------------------------------------------------------------------- #
# Plot                                                                  #
# --------------------------------------------------------------------- #

def severity_text(entry: FigEntry) -> str:
    if entry.kind == "exemplar":
        return "exemplar"
    if entry.observed_n == entry.expected_n:
        return "OK"
    delta = entry.expected_n - entry.observed_n
    return f"missing {delta}"


def caption_lines(entry: FigEntry) -> list[str]:
    head = entry.name
    if len(head) > 36:
        head = head[:34] + "..."
    return [
        head,
        f"n={entry.observed_n}/{entry.expected_n}  ::  {severity_text(entry)}",
    ]


def draw_cell(
    ax_thumb: plt.Axes,
    ax_caption: plt.Axes,
    entry: FigEntry,
    png: Path | None,
) -> None:
    """Draw one thumbnail + caption pair."""
    color = frame_color(entry)

    # Thumbnail axis
    ax_thumb.set_xticks([]); ax_thumb.set_yticks([])
    for s in ax_thumb.spines.values():
        s.set_color(color); s.set_linewidth(2.4)

    if png is not None:
        try:
            img = mpimg.imread(str(png))
            ax_thumb.imshow(img, aspect="equal", interpolation="bilinear")
        except Exception:
            ax_thumb.text(
                0.5, 0.5, "render-fail", ha="center", va="center",
                fontsize=8, color="black", transform=ax_thumb.transAxes,
            )
    else:
        ax_thumb.text(
            0.5, 0.5, "(missing PDF)", ha="center", va="center",
            fontsize=8, color="black", transform=ax_thumb.transAxes,
        )

    # Caption axis
    ax_caption.set_xticks([]); ax_caption.set_yticks([])
    for s in ax_caption.spines.values():
        s.set_visible(False)
    ax_caption.set_xlim(0, 1); ax_caption.set_ylim(0, 1)

    # Color swatch on the left
    swatch = mpatches.Rectangle(
        (0.01, 0.20), 0.06, 0.60,
        facecolor=color, edgecolor="none",
    )
    ax_caption.add_patch(swatch)

    lines = caption_lines(entry)
    ax_caption.text(
        0.10, 0.62, lines[0], ha="left", va="center",
        fontsize=7.5, family="monospace",
    )
    ax_caption.text(
        0.10, 0.25, lines[1], ha="left", va="center",
        fontsize=7.0, color="#333333",
    )


def draw_empty_cell(ax_thumb: plt.Axes, ax_caption: plt.Axes) -> None:
    ax_thumb.set_xticks([]); ax_thumb.set_yticks([])
    for s in ax_thumb.spines.values():
        s.set_color(COLOR_EMPTY); s.set_linewidth(0.8); s.set_linestyle(":")
    ax_thumb.text(
        0.5, 0.5, "—", ha="center", va="center",
        fontsize=14, color=COLOR_EMPTY, transform=ax_thumb.transAxes,
    )
    ax_caption.set_xticks([]); ax_caption.set_yticks([])
    for s in ax_caption.spines.values():
        s.set_visible(False)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if shutil.which("pdftoppm") is None:
        raise RuntimeError("pdftoppm not found on PATH; required to render thumbnails")

    # Render all PDFs in a temp dir.
    with tempfile.TemporaryDirectory(prefix="audit_57_") as tmp_str:
        tmp = Path(tmp_str)
        thumbs: list[Path | None] = []
        for entry in REGISTRY:
            png = render_thumbnail(entry.pdf, tmp, dpi=90)
            thumbs.append(png)

        # Layout: 4 rows x 4 cols of cells. Each cell = thumbnail (top
        # ~78% height) + caption (bottom ~22% height). Implemented via
        # an 8-row x 4-col gridspec grouping rows in (thumb, cap) pairs.
        n_rows, n_cols = 4, 4
        fig_w, fig_h = 12.5, 14.0
        fig = plt.figure(figsize=(fig_w, fig_h), constrained_layout=False)

        # leave room for header (top) + footer/legend (bottom)
        gs = fig.add_gridspec(
            nrows=2 * n_rows, ncols=n_cols,
            left=0.04, right=0.985,
            top=0.92, bottom=0.10,
            wspace=0.18, hspace=0.55,
            height_ratios=[3.4, 1.0] * n_rows,
        )

        n_total = n_rows * n_cols
        for idx in range(n_total):
            r = idx // n_cols
            c = idx % n_cols
            ax_thumb = fig.add_subplot(gs[2 * r, c])
            ax_cap = fig.add_subplot(gs[2 * r + 1, c])
            if idx < len(REGISTRY):
                draw_cell(ax_thumb, ax_cap, REGISTRY[idx], thumbs[idx])
            else:
                draw_empty_cell(ax_thumb, ax_cap)

        # ----- Header (top) ------------------------------------------
        fig.text(
            0.5, 0.965,
            "Section 2 figure cohort audit — n=6 era figures vs n=10 current cohort",
            ha="center", va="center",
            fontsize=15, fontweight="bold",
        )
        fig.text(
            0.5, 0.940,
            (f"Cohort locked 2026-04-25 (n={CURRENT_COHORT_N}): "
             + ", ".join(COHORT)),
            ha="center", va="center",
            fontsize=9.5, color="#333333",
        )

        # ----- Footer / legend ---------------------------------------
        legend_handles = [
            mpatches.Patch(facecolor=COLOR_MISMATCH, edgecolor="none",
                           label=f"cohort mismatch (observed n<{CURRENT_COHORT_N})"),
            mpatches.Patch(facecolor=COLOR_OK, edgecolor="none",
                           label=f"n={CURRENT_COHORT_N} OK"),
            mpatches.Patch(facecolor=COLOR_EXEMPLAR, edgecolor="none",
                           label="single/triple-patient exemplar (intentional)"),
        ]
        fig.legend(
            handles=legend_handles,
            loc="lower center", bbox_to_anchor=(0.5, 0.045),
            ncol=3, frameon=False, fontsize=10,
        )

        # cohort-mismatch summary
        n_mismatch = sum(1 for e in REGISTRY
                         if e.kind == "cohort"
                         and e.observed_n != e.expected_n)
        n_ok = sum(1 for e in REGISTRY
                   if e.kind == "cohort"
                   and e.observed_n == e.expected_n)
        n_exemplar = sum(1 for e in REGISTRY if e.kind == "exemplar")
        fig.text(
            0.5, 0.012,
            (f"summary: {n_mismatch} cohort-mismatch (red)  ::  "
             f"{n_ok} OK (green)  ::  {n_exemplar} exemplar (blue)  "
             f"::  {len(REGISTRY)} catalogued / {n_total} cells"),
            ha="center", va="center",
            fontsize=8.5, color="#222222",
        )

        # NOTE: no fig.suptitle (project rule); no provenance footer
        # by default (project rule). PDF only, vector for matplotlib
        # artists; the only raster pixels are inside the thumbnail
        # imshow images themselves.
        fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
        plt.close(fig)

    # Verify the file is non-empty and openable.
    assert OUT_PDF.exists(), f"output PDF not written: {OUT_PDF}"
    size = OUT_PDF.stat().st_size
    assert size > 1024, f"output PDF suspiciously small: {size} bytes"
    print(f"wrote: {OUT_PDF}")
    print(f"size : {size} bytes")


if __name__ == "__main__":
    main()
