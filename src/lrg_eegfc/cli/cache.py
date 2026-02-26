"""``lrg-eegfc cache`` subcommands."""

from __future__ import annotations

from pathlib import Path

import click

from ._common import patient_options, resolve_patients


# ---------------------------------------------------------------------------
# Cache root directories (mirrors workflow defaults)
# ---------------------------------------------------------------------------

_CACHE_ROOTS = {
    "corr":    Path("data/corr_cache"),
    "msc":     Path("data/msc_cache"),
    "lrg":     Path("data/lrg_cache"),
    "cleaned": Path("data/cleaned_corr_cache"),
}

_ALL_TYPES = list(_CACHE_ROOTS.keys())


def _human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TB"


# ---------------------------------------------------------------------------
# Group
# ---------------------------------------------------------------------------

@click.group()
def cache() -> None:
    """Manage cached computation results."""


# ---------------------------------------------------------------------------
# list
# ---------------------------------------------------------------------------

@cache.command("list")
@click.option(
    "--cache-type", "-t",
    type=click.Choice(_ALL_TYPES + ["all"]),
    default="all",
    show_default=True,
    help="Type of cache to inspect.",
)
@patient_options()
def list_cache(cache_type: str, patients: tuple[str, ...]) -> None:
    """List cached files with sizes and counts."""
    types = _ALL_TYPES if cache_type == "all" else [cache_type]

    for ctype in types:
        root = _CACHE_ROOTS[ctype]
        if not root.exists():
            click.echo(f"{ctype:8s}  {root}/  [not found]")
            continue

        # Optionally filter by patient
        if patients:
            pat_dirs = [root / p for p in patients if (root / p).exists()]
        else:
            pat_dirs = sorted(p for p in root.iterdir() if p.is_dir())

        total_files = 0
        total_bytes = 0
        for pdir in pat_dirs:
            files = list(pdir.glob("*.npy")) + list(pdir.glob("*.npz"))
            n = len(files)
            s = sum(f.stat().st_size for f in files)
            total_files += n
            total_bytes += s
            click.echo(f"  {ctype:8s}  {pdir.name:8s}  {n:4d} files  {_human_size(s)}")

        if not pat_dirs:
            click.echo(f"  {ctype:8s}  (empty)")
        else:
            click.echo(
                f"  {ctype:8s}  TOTAL     {total_files:4d} files  "
                f"{_human_size(total_bytes)}"
            )
        click.echo()


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@cache.command()
def info() -> None:
    """Print cache naming conventions."""
    click.echo(
        "Cache naming conventions\n"
        "========================\n\n"
        "Correlation:\n"
        "  {root}/Pat_XX/{band}_{phase}_corr_ftype-{abs|none}_zdiag-{True|False}.npy\n\n"
        "MSC (dense):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-none_nperseg-{N}.npy\n\n"
        "MSC (soft-sparsified):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-soft_nsurr-{N}_nperseg-{N}.npy\n\n"
        "MSC (FDR):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-fdr_nsurr-{N}_q-{q}_nperseg-{N}.npy\n\n"
        "MSC (disparity):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-disparity_alpha-{a}_nperseg-{N}.npy\n\n"
        "MSC (hybrid):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-hybrid_nsurr-{N}_alpha-{a}_nperseg-{N}.npy\n\n"
        "MSC (ECM):\n"
        "  {root}/Pat_XX/{band}_{phase}_msc_sparsify-ecm_alpha-{a}_nens-{N}_wscale-{W}_nperseg-{N}.npy\n\n"
        "LRG:\n"
        "  {root}/Pat_XX/{band}_{phase}_lrg_{corr|msc}.npz\n\n"
        "Dev mode (--filter-time) uses *_dev cache roots."
    )


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------

@cache.command()
@click.option(
    "--cache-type", "-t",
    type=click.Choice(_ALL_TYPES + ["all"]),
    default="all",
    show_default=True,
)
@patient_options()
def verify(cache_type: str, patients: tuple[str, ...]) -> None:
    """Verify cache integrity (files loadable, correct shapes)."""
    import numpy as np

    types = _ALL_TYPES if cache_type == "all" else [cache_type]
    ok_count = 0
    bad_count = 0

    for ctype in types:
        root = _CACHE_ROOTS[ctype]
        if not root.exists():
            continue

        if patients:
            pat_dirs = [root / p for p in patients if (root / p).exists()]
        else:
            pat_dirs = sorted(p for p in root.iterdir() if p.is_dir())

        for pdir in pat_dirs:
            for f in sorted(pdir.glob("*.np[yz]")):
                try:
                    if f.suffix == ".npy":
                        arr = np.load(f)
                        if not np.all(np.isfinite(arr)):
                            click.echo(f"  [WARN] {f}  contains non-finite values")
                    else:
                        np.load(f, allow_pickle=False)
                    ok_count += 1
                except Exception as exc:
                    click.echo(f"  [BAD]  {f}  {exc}")
                    bad_count += 1

    click.echo(f"\nVerified: {ok_count} ok, {bad_count} corrupt.")


# ---------------------------------------------------------------------------
# clean
# ---------------------------------------------------------------------------

@cache.command("clean")
@click.option("--dry-run", is_flag=True, help="Show what would be removed.")
@click.option("--dev-only", is_flag=True, help="Only remove _dev caches.")
@click.confirmation_option(prompt="Remove cache files?")
def clean_cache(dry_run: bool, dev_only: bool) -> None:
    """Remove stale or development cache files."""
    import shutil

    targets: list[Path] = []

    if dev_only:
        for root in _CACHE_ROOTS.values():
            dev_root = root.parent / (root.name + "_dev")
            if dev_root.exists():
                targets.append(dev_root)
    else:
        for root in _CACHE_ROOTS.values():
            if root.exists():
                targets.append(root)
            dev_root = root.parent / (root.name + "_dev")
            if dev_root.exists():
                targets.append(dev_root)

    if not targets:
        click.echo("Nothing to clean.")
        return

    total = 0
    for t in targets:
        size = sum(f.stat().st_size for f in t.rglob("*") if f.is_file())
        label = "would remove" if dry_run else "removing"
        click.echo(f"  {label}: {t}  ({_human_size(size)})")
        if not dry_run:
            shutil.rmtree(t)
        total += size

    click.echo(f"\n{'Would free' if dry_run else 'Freed'}: {_human_size(total)}")
