"""Root CLI application with lazy-loaded subcommand groups.

The LazyGroup pattern ensures that ``lrg-eegfc --help`` responds instantly
without importing heavyweight dependencies (numpy, scipy, matplotlib).
Each subcommand module is imported only when the user actually invokes it.
"""

from __future__ import annotations

import importlib
from typing import Dict, Optional

import click


# Short descriptions shown in --help without importing the modules.
_LAZY_HELP: Dict[str, str] = {
    "bundle":  "Build publication bundles.",
    "cache":   "Manage cached computation results.",
    "compute": "Compute FC matrices and analyses.",
    "config":  "Display project configuration.",
    "data":    "Inspect and manage patient data.",
    "plot":    "Generate visualizations from cache.",
    "show":    "Query cached results (print stats, no figures).",
}


class LazyGroup(click.Group):
    """Click group that defers subcommand imports until invocation.

    The ``--help`` listing uses pre-registered short descriptions so that
    no subcommand module is loaded just to display help.
    """

    def __init__(
        self,
        *args,
        lazy_subcommands: Optional[Dict[str, str]] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self._lazy_subcommands = lazy_subcommands or {}

    def list_commands(self, ctx: click.Context) -> list[str]:
        base = super().list_commands(ctx)
        lazy = sorted(self._lazy_subcommands.keys())
        return base + lazy

    def get_command(
        self, ctx: click.Context, cmd_name: str
    ) -> Optional[click.BaseCommand]:
        if cmd_name in self._lazy_subcommands:
            return self._import_command(cmd_name)
        return super().get_command(ctx, cmd_name)

    def format_commands(
        self, ctx: click.Context, formatter: click.HelpFormatter
    ) -> None:
        """Override to avoid importing lazy modules just for --help."""
        commands = []
        for subcommand in self.list_commands(ctx):
            # Use static help text for lazy commands
            if subcommand in self._lazy_subcommands:
                help_text = _LAZY_HELP.get(subcommand, "")
                commands.append((subcommand, help_text))
            else:
                cmd = self.get_command(ctx, subcommand)
                if cmd is None or cmd.hidden:
                    continue
                help_text = cmd.get_short_help_str(limit=150)
                commands.append((subcommand, help_text))

        if commands:
            with formatter.section("Commands"):
                formatter.write_dl(commands)

    def _import_command(self, cmd_name: str) -> click.BaseCommand:
        module_path, attr_name = self._lazy_subcommands[cmd_name].rsplit(":", 1)
        mod = importlib.import_module(module_path)
        return getattr(mod, attr_name)


@click.group(
    cls=LazyGroup,
    context_settings={"max_content_width": 120},
    lazy_subcommands={
        "compute": "lrg_eegfc.cli.compute:compute",
        "plot": "lrg_eegfc.cli.plot:plot",
        "show": "lrg_eegfc.cli.show:show",
        "data": "lrg_eegfc.cli.data:data",
        "cache": "lrg_eegfc.cli.cache:cache",
        "config": "lrg_eegfc.cli.config_cmd:config",
        "bundle": "lrg_eegfc.cli.bundle:bundle",
    },
)
@click.version_option(package_name="lrg-eegfc")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output.")
@click.pass_context
def app(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """LRG EEG Functional Connectivity Analysis Toolkit.

    Analyse stereo-EEG recordings through the Laplacian Renormalization Group
    framework: compute functional connectivity, run hierarchical analysis,
    and generate publication-ready figures.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
