"""``lrg-eegfc show`` subcommands.

Query cached results and print summary statistics to the terminal.
No figures are generated — this is for quick numerical inspection.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click

from lrg_eegfc.config.paths import CORR_CACHE, MSC_CACHE, LRG_CACHE

from ._common import (
    band_phase_options,
    fc_method_option,
    patient_options,
    resolve_bands,
    resolve_patients,
    resolve_phases,
    single_patient_option,
    verbose_option,
)


@click.group()
def show() -> None:
    """Query cached results (print stats, no figures)."""


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

def _show_options() -> callable:
    """Stack ``--patient``, ``--patients``, and ``--all`` on a show command."""
    def decorator(f):
        f = click.option(
            "--all", "show_all", is_flag=True,
            help="Show results for all discovered patients.",
        )(f)
        f = single_patient_option(required=False)(f)
        f = patient_options()(f)
        return f
    return decorator


def _resolve_show_patients(
    ctx: click.Context,
    patients: tuple[str, ...],
    patient: Optional[str],
    show_all: bool,
) -> Optional[list[str]]:
    """Resolve patient list for show commands.

    Returns *None* (and prints ``--help``) when nothing is specified.
    """
    if patient:
        return [patient]
    if patients:
        return list(patients)
    if show_all:
        return resolve_patients(())
    # Nothing specified → show help
    click.echo(ctx.get_help())
    return None


# -- compact one-liners --

def _mat_compact(mat) -> str:
    """One-line summary for a matrix."""
    import numpy as np

    n = mat.shape[0]
    mask = ~np.eye(n, dtype=bool)
    off = mat[mask]
    nonzero = int(np.count_nonzero(off))
    total = off.size
    density = f"{nonzero / total:.0%}" if total else "?"
    return f"{n}x{n}  mean={off.mean():.3f}  std={off.std():.3f}  density={density}"


def _lrg_compact(r) -> str:
    """One-line summary for an LRG result."""
    import numpy as np

    parts = [f"n={r.n_nodes}", f"tau*={r.optimal_threshold:.4f}"]
    if r.entropy_C is not None and len(r.entropy_C) > 0:
        parts.append(f"C_peak={r.entropy_C.max():.3f}")
    if r.linkage_matrix is not None:
        parts.append(f"merges={r.linkage_matrix.shape[0]}")
    return "  ".join(parts)


def _cleaned_compact(mat, meta) -> str:
    """One-line summary for a cleaned correlation matrix."""
    import numpy as np

    n = mat.shape[0]
    mask = ~np.eye(n, dtype=bool)
    off = mat[mask]
    nonzero = int(np.count_nonzero(off))
    total = off.size
    density = f"{nonzero / total:.0%}" if total else "?"
    parts = [f"{n}x{n}", f"mean={off.mean():.3f}", f"density={density}"]
    if meta and "n_signal_components" in meta:
        parts.append(f"signal_comps={meta['n_signal_components']}")
    return "  ".join(parts)


# -- verbose detail tables --

def _mat_full(mat) -> list[tuple[str, str]]:
    """Detailed stats as table rows."""
    import numpy as np

    n = mat.shape[0]
    mask = ~np.eye(n, dtype=bool)
    off = mat[mask]
    nonzero = int(np.count_nonzero(off))
    total = off.size
    return [
        ("shape", f"{n} x {n}"),
        ("mean", f"{off.mean():.4f}"),
        ("std", f"{off.std():.4f}"),
        ("min", f"{off.min():.4f}"),
        ("max", f"{off.max():.4f}"),
        ("median", f"{float(np.median(off)):.4f}"),
        ("density", f"{nonzero / total:.2%}" if total else "N/A"),
        ("edges", f"{nonzero // 2}"),
    ]


def _lrg_full(r) -> list[tuple[str, str]]:
    """Detailed LRG stats as table rows."""
    import numpy as np

    rows: list[tuple[str, str]] = [
        ("n_nodes", str(r.n_nodes)),
        ("optimal_threshold", f"{r.optimal_threshold:.4f}"),
    ]
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
                      f"[{r.entropy_C.min():.4f}, {r.entropy_C.max():.4f}]"))
        rows.append(("C peak",
                      f"{r.entropy_C.max():.4f} at tau="
                      f"{r.entropy_tau[np.argmax(r.entropy_C)]:.2e}"))
    if r.ultrametric_matrix is not None:
        um = r.ultrametric_matrix
        if um.ndim == 1:
            rows.append(("ultrametric_mean", f"{um.mean():.4f}"))
            rows.append(("ultrametric_std", f"{um.std():.4f}"))
            rows.append(("ultrametric_range",
                          f"[{um.min():.4f}, {um.max():.4f}]"))
        else:
            mask = ~np.eye(um.shape[0], dtype=bool)
            off = um[mask]
            rows.append(("ultrametric_mean", f"{off.mean():.4f}"))
            rows.append(("ultrametric_std", f"{off.std():.4f}"))
    if r.linkage_matrix is not None:
        n_merges = r.linkage_matrix.shape[0]
        rows.append(("linkage_merges", str(n_merges)))
        rows.append(("max_merge_dist",
                      f"{r.linkage_matrix[:, 2].max():.4f}"))
    return rows


def _cleaned_full(mat, meta) -> list[tuple[str, str]]:
    """Detailed cleaned correlation stats as table rows."""
    import numpy as np

    rows = _mat_full(mat)
    if meta:
        rows.append(("", ""))
        rows.append(("-- MP cleaning --", ""))
        if "threshold" in meta:
            rows.append(("threshold", f"{meta['threshold']:.6f}"))
        if "lambda_min" in meta:
            rows.append(("lambda_min", f"{meta['lambda_min']:.4f}"))
        if "lambda_max" in meta:
            rows.append(("lambda_max", f"{meta['lambda_max']:.4f}"))
        if "n_signal_components" in meta:
            rows.append(("signal_components", str(meta["n_signal_components"])))
        if "signal_eigenvalues" in meta:
            ev = meta["signal_eigenvalues"]
            rows.append(("signal_eig_range",
                          f"[{ev.min():.4f}, {ev.max():.4f}]"))
    return rows


# -- printing --

def _print_table(rows: list[tuple[str, str]], indent: int = 2) -> None:
    """Print a two-column table."""
    if not rows:
        return
    w = max(len(r[0]) for r in rows) + 1
    pad = " " * indent
    for label, value in rows:
        click.echo(f"{pad}{label:<{w}}  {value}")


def _footer(found: int, hint_cmd: str, example_pat: str) -> None:
    """Print entry count or a no-results hint."""
    if found == 0:
        click.echo(f"\nNo cached results found.")
        click.echo(f"  Hint: lrg-eegfc {hint_cmd} --patients {example_pat}")
    else:
        click.echo(f"\n{found} entries shown.")


# ---------------------------------------------------------------------------
# show corr
# ---------------------------------------------------------------------------

@show.command()
@_show_options()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=CORR_CACHE, show_default=True)
@verbose_option()
@click.pass_context
def corr(ctx, patients, patient, show_all, band, bands, phase, phases,
         cache_root, verbose):
    """Show correlation matrix statistics from cache."""
    pat_list = _resolve_show_patients(ctx, patients, patient, show_all)
    if pat_list is None:
        return

    from lrg_eegfc.workflow.corr import load_corr_matrix

    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    found = 0

    for pat in pat_list:
        pat_printed = False
        for p in phase_list:
            phase_printed = False
            for b in band_list:
                mat = load_corr_matrix(pat, p, b, cache_root=cache_root)
                if mat is None:
                    continue
                if not pat_printed:
                    click.echo(f"\n{pat}")
                    pat_printed = True
                if not phase_printed:
                    click.echo(f"  {p}")
                    phase_printed = True
                found += 1
                if verbose:
                    click.echo(f"    {b}")
                    _print_table(_mat_full(mat), indent=6)
                else:
                    click.echo(f"    {b:14s} {_mat_compact(mat)}")

    _footer(found, "compute corr", pat_list[0] if pat_list else "Pat_02")


# ---------------------------------------------------------------------------
# show msc
# ---------------------------------------------------------------------------

@show.command()
@_show_options()
@band_phase_options()
@click.option("--sparsify", default="none", show_default=True)
@click.option("--n-surrogates", type=int, default=0, show_default=True)
@click.option("--nperseg", type=int, default=4096, show_default=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=MSC_CACHE, show_default=True)
@verbose_option()
@click.pass_context
def msc(ctx, patients, patient, show_all, band, bands, phase, phases,
        sparsify, n_surrogates, nperseg, cache_root, verbose):
    """Show MSC matrix statistics from cache."""
    pat_list = _resolve_show_patients(ctx, patients, patient, show_all)
    if pat_list is None:
        return

    from lrg_eegfc.workflow.msc import load_msc_matrix

    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    found = 0

    for pat in pat_list:
        pat_printed = False
        for p in phase_list:
            phase_printed = False
            for b in band_list:
                mat = load_msc_matrix(
                    pat, p, b,
                    sparsify=sparsify,
                    n_surrogates=n_surrogates,
                    nperseg=nperseg,
                    cache_root=cache_root,
                )
                if mat is None:
                    continue
                if not pat_printed:
                    click.echo(f"\n{pat}")
                    pat_printed = True
                if not phase_printed:
                    click.echo(f"  {p}")
                    phase_printed = True
                found += 1
                if verbose:
                    click.echo(f"    {b}  (sparsify={sparsify})")
                    _print_table(_mat_full(mat), indent=6)
                else:
                    click.echo(f"    {b:14s} {_mat_compact(mat)}")

    _footer(found, "compute msc", pat_list[0] if pat_list else "Pat_02")


# ---------------------------------------------------------------------------
# show lrg
# ---------------------------------------------------------------------------

@show.command()
@_show_options()
@band_phase_options()
@fc_method_option(required=True)
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=LRG_CACHE, show_default=True)
@verbose_option()
@click.pass_context
def lrg(ctx, patients, patient, show_all, band, bands, phase, phases,
        fc_method, cache_root, verbose):
    """Show LRG analysis results from cache."""
    pat_list = _resolve_show_patients(ctx, patients, patient, show_all)
    if pat_list is None:
        return

    from lrg_eegfc.workflow.lrg import load_lrg_result

    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    found = 0

    for pat in pat_list:
        pat_printed = False
        for p in phase_list:
            phase_printed = False
            for b in band_list:
                r = load_lrg_result(
                    pat, p, b,
                    fc_method=fc_method,
                    cache_root=cache_root,
                )
                if r is None:
                    continue
                if not pat_printed:
                    click.echo(f"\n{pat}")
                    pat_printed = True
                if not phase_printed:
                    click.echo(f"  {p}")
                    phase_printed = True
                found += 1
                if verbose:
                    click.echo(f"    {b}  (fc={fc_method})")
                    _print_table(_lrg_full(r), indent=6)
                else:
                    click.echo(f"    {b:14s} {_lrg_compact(r)}")

    _footer(
        found,
        f"compute lrg --fc-method {fc_method}",
        pat_list[0] if pat_list else "Pat_02",
    )


# ---------------------------------------------------------------------------
# show cleaned
# ---------------------------------------------------------------------------

@show.command("cleaned")
@_show_options()
@band_phase_options()
@click.option("--cache-root", type=click.Path(path_type=Path),
              default=CORR_CACHE, show_default=True)
@verbose_option()
@click.pass_context
def cleaned(ctx, patients, patient, show_all, band, bands, phase, phases,
            cache_root, verbose):
    """Show Marchenko-Pastur cleaning statistics from cache."""
    pat_list = _resolve_show_patients(ctx, patients, patient, show_all)
    if pat_list is None:
        return

    from lrg_eegfc.workflow.cleaning import load_cleaned_corr_matrix

    band_list = resolve_bands(bands, band)
    phase_list = resolve_phases(phases, phase)
    found = 0

    for pat in pat_list:
        pat_printed = False
        for p in phase_list:
            phase_printed = False
            for b in band_list:
                try:
                    mat, meta = load_cleaned_corr_matrix(
                        pat, p, b,
                        cache_root=cache_root,
                        load_metadata=True,
                    )
                except (FileNotFoundError, TypeError):
                    continue
                if mat is None:
                    continue
                if not pat_printed:
                    click.echo(f"\n{pat}")
                    pat_printed = True
                if not phase_printed:
                    click.echo(f"  {p}")
                    phase_printed = True
                found += 1
                if verbose:
                    click.echo(f"    {b}")
                    _print_table(_cleaned_full(mat, meta), indent=6)
                else:
                    click.echo(f"    {b:14s} {_cleaned_compact(mat, meta)}")

    _footer(found, "compute clean", pat_list[0] if pat_list else "Pat_02")
