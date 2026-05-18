"""``karyoscope bin`` — aggregate base-pair annotations into larger bins.

This command is not yet implemented in v0.1.0.dev.

(The module is named ``bin_cmd`` rather than ``bin`` to avoid shadowing
Python's built-in :func:`bin`.)
"""

from __future__ import annotations

import click


@click.command(
    help="Aggregate base-pair annotation BEDs into larger bins.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will port bin_features.py logic."""
    raise click.ClickException("`karyoscope bin` is not yet implemented in this development build.")
