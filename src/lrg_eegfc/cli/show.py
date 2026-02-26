"""``lrg-eegfc show`` subcommands.

Query cached results and print summary statistics to the terminal.
No figures are generated — this is for quick numerical inspection.
"""

from __future__ import annotations

from pathlib import Path

import click

from ._common import (
    CliReporter,
    band_phase_options,
    fc_method_option,
    resolve_bands,
    resolve_phases,
    single_patient_option,
    verbose_option,
)


@click.group()
def show() -> None:
    """Query cached results (print stats, no figures)."""


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _mat_stats(mat) -> dict:
    """Compute summary statistics for a matrix."""
    import numpy as np

    n = mat.shape[0]
    # Mask out diagonal for off-diagonal stats
    mask = ~np.eye(n, dtype=bool)
    off = mat[mask]

    nonzero = np.count_nonzero(off)
    total = off.size

    return {
        "shape": f"{n} x {n}",
        "mean": f"{off.mean():.4f}",
        "std": f"{off.std():.4f}",
        "min": f"{off.min():.4f}",
        "max": f"{off.max():.4f}",
        "median": f"{float(np.median(off)):.4f}",
        "density": f"{nonzero / total:.2%}" if total else "N/A",
        "nonzero_edges": f"{nonzero // 2}",  # symmetric, so /2
    }


def _print_table(rows: list[tuple[str, str]], indent: int = 2) -> None:
    """Print a two-column table."""
    if not rows:
        return
    w = max(len(r[0]) for r in rows) + 1
    pad = " " * indent
    for label, value in rows:
        click.echo(f"{pad}{label:<{w}}  {value}")


# ---------------------------------------------------------------------------
# show corr
# ---------------------------------------------------------------------------

@show.command()
@single_patient_option()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True)
@verbose_option()
@click.pass_context
def corr(ctx, patient, band, bands, phase, phases, cache_root, verbose):
    """Show correlation matrix statistics from cache."""
    from lrg_eegfc.workflow.corr import load_corr_matrix

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    for b in band_list:
        for p in phase_list:
            mat = load_corr_matrix(patient, p, b, cache_root=cache_root)
            if mat is None:
                rpt.skip(f"{b}/{p}: not cached")
                continue

            click.echo(f"\n{patient}  {b}  {p}  (correlation)")
            click.echo("-" * 50)
            _print_table(list(_mat_stats(mat).items()))


# ---------------------------------------------------------------------------
# show msc
# ---------------------------------------------------------------------------

@show.command()
@single_patient_option()
@band_phase_options()
@click.option("--sparsify", default="none", show_default=True)
@click.option("--n-surrogates", type=int, default=0, show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/msc_cache"), show_default=True)
@verbose_option()
@click.pass_context
def msc(ctx, patient, band, bands, phase, phases, sparsify, n_surrogates,
        nperseg, cache_root, verbose):
    """Show MSC matrix statistics from cache."""
    from lrg_eegfc.workflow.msc import load_msc_matrix

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    for b in band_list:
        for p in phase_list:
            mat = load_msc_matrix(
                patient, p, b,
                sparsify=sparsify,
                n_surrogates=n_surrogates,
                nperseg=nperseg,
                cache_root=cache_root,
            )
            if mat is None:
                rpt.skip(f"{b}/{p}: not cached")
                continue

            click.echo(f"\n{patient}  {b}  {p}  (MSC, sparsify={sparsify})")
            click.echo("-" * 50)
            _print_table(list(_mat_stats(mat).items()))


# ---------------------------------------------------------------------------
# show lrg
# ---------------------------------------------------------------------------

@show.command()
@single_patient_option()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/lrg_cache"), show_default=True)
@verbose_option()
@click.pass_context
def lrg(ctx, patient, band, bands, phase, phases, fc_method, cache_root,
        verbose):
    """Show LRG analysis results from cache."""
    import numpy as np
    from lrg_eegfc.workflow.lrg import load_lrg_result
    from scipy.spatial.distance import squareform

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    for b in band_list:
        for p in phase_list:
            r = load_lrg_result(
                patient, p, b,
                fc_method=fc_method,
                cache_root=cache_root,
            )
            if r is None:
                rpt.skip(f"{b}/{p}: not cached")
                continue

            click.echo(f"\n{patient}  {b}  {p}  (LRG, fc={fc_method})")
            click.echo("-" * 50)

            # Core info
            rows = [
                ("n_nodes", str(r.n_nodes)),
                ("optimal_threshold", f"{r.optimal_threshold:.4f}"),
            ]

            # Entropy curve summary
            if r.entropy_tau is not None and len(r.entropy_tau) > 0:
                rows.append(("entropy_tau_range",
                             f"[{r.entropy_tau[0]:.2e}, {r.entropy_tau[-1]:.2e}]"))
                rows.append(("entropy_steps", str(len(r.entropy_tau))))
            if r.entropy_1_minus_S is not None and len(r.entropy_1_minus_S) > 0:
                rows.append(("1-S range",
                             f"[{r.entropy_1_minus_S.min():.4f}, "
                             f"{r.entropy_1_minus_S.max():.4f}]"))
            if r.entropy_C is not None and len(r.entropy_C) > 0:
                rows.append(("C range",
                             f"[{r.entropy_C.min():.4f}, "
                             f"{r.entropy_C.max():.4f}]"))
                rows.append(("C peak",
                             f"{r.entropy_C.max():.4f} at tau="
                             f"{r.entropy_tau[np.argmax(r.entropy_C)]:.2e}"))

            # Ultrametric distance stats
            if r.ultrametric_matrix is not None:
                um = r.ultrametric_matrix
                if um.ndim == 1:
                    # condensed form
                    rows.append(("ultrametric_mean", f"{um.mean():.4f}"))
                    rows.append(("ultrametric_std", f"{um.std():.4f}"))
                    rows.append(("ultrametric_range",
                                 f"[{um.min():.4f}, {um.max():.4f}]"))
                else:
                    mask = ~np.eye(um.shape[0], dtype=bool)
                    off = um[mask]
                    rows.append(("ultrametric_mean", f"{off.mean():.4f}"))
                    rows.append(("ultrametric_std", f"{off.std():.4f}"))

            # Linkage info
            if r.linkage_matrix is not None:
                n_merges = r.linkage_matrix.shape[0]
                rows.append(("linkage_merges", str(n_merges)))
                rows.append(("max_merge_dist",
                             f"{r.linkage_matrix[:, 2].max():.4f}"))

            _print_table(rows)


# ---------------------------------------------------------------------------
# show cleaned
# ---------------------------------------------------------------------------

@show.command("cleaned")
@single_patient_option()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=Path("data/corr_cache"), show_default=True)
@verbose_option()
@click.pass_context
def cleaned(ctx, patient, band, bands, phase, phases, cache_root, verbose):
    """Show Marchenko-Pastur cleaning statistics from cache."""
    from lrg_eegfc.workflow.cleaning import load_cleaned_corr_matrix

    rpt = CliReporter.from_context(ctx, verbose=verbose)
    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)

    for b in band_list:
        for p in phase_list:
            try:
                mat, meta = load_cleaned_corr_matrix(
                    patient, p, b,
                    cache_root=cache_root,
                    load_metadata=True,
                )
            except (FileNotFoundError, TypeError):
                rpt.skip(f"{b}/{p}: not cached")
                continue

            if mat is None:
                rpt.skip(f"{b}/{p}: not cached")
                continue

            click.echo(f"\n{patient}  {b}  {p}  (cleaned correlation)")
            click.echo("-" * 50)

            # Matrix stats
            stats = _mat_stats(mat)
            _print_table(list(stats.items()))

            # Cleaning metadata
            if meta:
                click.echo()
                click.echo("  MP cleaning metadata:")
                mp_rows = []
                if "threshold" in meta:
                    mp_rows.append(("  threshold", f"{meta['threshold']:.6f}"))
                if "lambda_min" in meta:
                    mp_rows.append(("  lambda_min", f"{meta['lambda_min']:.4f}"))
                if "lambda_max" in meta:
                    mp_rows.append(("  lambda_max", f"{meta['lambda_max']:.4f}"))
                if "n_signal_components" in meta:
                    mp_rows.append(("  signal_components", str(meta["n_signal_components"])))
                if "signal_eigenvalues" in meta:
                    import numpy as np
                    ev = meta["signal_eigenvalues"]
                    mp_rows.append(("  signal_eig_range",
                                    f"[{ev.min():.4f}, {ev.max():.4f}]"))
                _print_table(mp_rows, indent=0)
