"""Top-level command-line interface for KaryoScope.

This module wires the eight v0.1 subcommands into a single ``karyoscope``
entry point. Each subcommand lives in its own module under
``karyoscope.commands`` so that they can grow independently.

Logging vs. program output
==========================

Two channels coexist deliberately:

* **Program output** — what a command is doing for the user. Emitted via
  ``click.echo``. Always visible regardless of verbosity. Examples:
  "Installing KS_human_CHM13_v2 v2.0.0...", a listing table, a help
  message.
* **Logging / diagnostics** — behind-the-scenes information for
  developers and power users debugging issues. Emitted via the standard
  ``logging`` module. Hidden by default; opt in with ``-v`` or ``-vv``.
  Examples: which URL was fetched, cache hits, SHA-256 verifications,
  subprocess invocations.

The verbosity flags only affect the logging channel. Program output is
unaffected; if you want to suppress that, use command-specific flags
like ``--quiet`` on ``download``.
"""

from __future__ import annotations

import importlib
import logging
import sys

import click

from karyoscope import diskspace
from karyoscope._version import __version__

#: Subcommand name -> ``module:attribute`` holding its ``click.Command``.
#:
#: Imported on demand rather than up front. Registering the subcommands used to
#: mean importing every command module, and one of them (``download``)
#: pulls in ``requests`` -- ~190 ms of the ~290 ms it took to import this
#: module, paid by every invocation including ``karyoscope --version`` and
#: every ``annotate`` run, for an HTTP library they never touch.
#:
#: ``--help`` still imports everything, because Click asks each command for its
#: short help to build the listing. That is the one path where the old cost
#: remains; it is not a regression, and avoiding it would mean duplicating each
#: command's summary here, where it could drift from the docstring it mirrors.
_LAZY_COMMANDS: dict[str, str] = {
    "download": "karyoscope.commands.download:cmd",
    "register": "karyoscope.commands.register:cmd",
    "prep-bed": "karyoscope.commands.prep_bed:cmd",
    "build": "karyoscope.commands.build:cmd",
    "annotate": "karyoscope.commands.annotate:cmd",
    "scaffold": "karyoscope.commands.scaffold:cmd",
    "remap-bed": "karyoscope.commands.remap_bed:cmd",
    "bin": "karyoscope.commands.bin_cmd:cmd",
    "centromeres": "karyoscope.commands.centromeres:cmd",
    "karyotype": "karyoscope.commands.karyotype:cmd",
    "info": "karyoscope.commands.info:cmd",
    "version": "karyoscope.commands.version:cmd",
}

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "max_content_width": 100,
}


class _KaryoscopeGroup(click.Group):
    """Command group that turns a full disk into a message, not a traceback.

    Commands check free space before they start, but an estimate can be
    beaten: a shared filesystem fills up underneath a long run, a quota
    lands mid-write, or a database's registry entry understates its size.
    When that happens the failure is a bare ``OSError: [Errno 28]`` from
    whichever write happened to be unlucky — a stack trace through
    ``tarfile`` or a BED writer that tells the user nothing about their
    actual problem.

    Catching it here covers every subcommand, including ones that don't
    (yet) run their own space check.
    """

    def list_commands(self, ctx: click.Context) -> list[str]:
        """Every subcommand, without importing any of them."""
        return sorted({*super().list_commands(ctx), *_LAZY_COMMANDS})

    def get_command(self, ctx: click.Context, name: str) -> click.Command | None:
        """Resolve one subcommand, importing only its module."""
        cmd = super().get_command(ctx, name)
        if cmd is not None:
            return cmd
        target = _LAZY_COMMANDS.get(name)
        if target is None:
            return None
        module_name, _, attribute = target.partition(":")
        module = importlib.import_module(module_name)
        return getattr(module, attribute)

    def invoke(self, ctx: click.Context) -> object:
        try:
            return super().invoke(ctx)
        except OSError as exc:
            if not diskspace.is_enospc(exc):
                raise
            command = ctx.invoked_subcommand or "karyoscope"
            raise click.ClickException(
                str(diskspace.enospc_error(exc, what=f"running `karyoscope {command}`"))
            ) from exc


#: Maps a verbosity integer (negative for quiet, 0 default, positive for verbose)
#: to a stdlib logging level. Anything beyond ``2`` is clamped to DEBUG.
_VERBOSITY_TO_LEVEL = {
    -1: logging.ERROR,
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def _configure_logging(verbosity: int) -> None:
    """Install a stderr log handler with a format suited to a CLI tool.

    Format depends on verbosity:

    * Default (WARNING) / quiet (ERROR): no timestamp, no module name.
      Just ``LEVEL: message``. The few lines you see at this level are
      mostly user-actionable; a timestamp would be noise.
    * ``-v`` (INFO): ``HH:MM:SS LEVEL: message``. Timestamps make
      per-step wall times extractable from the log via diffing
      consecutive lines -- useful for one-off benchmarks without
      pulling in a dedicated profiler.
    * ``-vv`` (DEBUG): same timestamp plus the module name in
      brackets, so subsystem-level debugging is easier to trace.

    Replaces any pre-existing handlers on the root logger so that repeated
    calls (e.g., from tests) don't compound output.
    """
    level = _VERBOSITY_TO_LEVEL.get(verbosity, logging.DEBUG if verbosity > 0 else logging.ERROR)
    handler = logging.StreamHandler(stream=sys.stderr)
    if level <= logging.DEBUG:
        fmt = "%(asctime)s %(levelname)s [%(name)s] %(message)s"
    elif level <= logging.INFO:
        fmt = "%(asctime)s %(levelname)s: %(message)s"
    else:
        fmt = "%(levelname)s: %(message)s"
    handler.setFormatter(logging.Formatter(fmt, datefmt="%H:%M:%S"))
    root = logging.getLogger()
    root.handlers = [handler]
    root.setLevel(level)


@click.group(cls=_KaryoscopeGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(
    __version__,
    "-V",
    "--version",
    message="karyoscope %(version)s",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase logging verbosity. Repeat for more (-v=info, -vv=debug).",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Decrease logging verbosity to errors only. Conflicts with -v.",
)
@click.pass_context
def main(ctx: click.Context, verbose: int, quiet: bool) -> None:
    """KaryoScope: rapid, alignment-free sequence annotation for the pangenome era.

    Run a subcommand with ``--help`` to see its options, e.g.:

    \b
        karyoscope annotate --help

    For more information, see https://github.com/barthel-lab/KaryoScope.
    """
    if quiet and verbose:
        raise click.UsageError("--quiet and --verbose cannot be combined.")
    verbosity = -1 if quiet else verbose
    _configure_logging(verbosity)
    # --quiet has to reach the *program output* channel too, not just
    # logging: the milestone lines the long-running commands print (see
    # karyoscope.progress) go to stdout, so lowering the log level alone
    # would leave no way to get a silent run. Stashed on the context
    # rather than threaded through every command signature.
    ctx.ensure_object(dict)["quiet"] = quiet


# Register subcommands. The order here determines the order in `--help`.


if __name__ == "__main__":
    main()
