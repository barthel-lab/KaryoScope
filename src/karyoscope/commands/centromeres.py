"""``karyoscope centromeres`` — extract centromere coordinates from feature BEDs.

This command is not yet implemented in v0.1.0.dev.
"""

from __future__ import annotations

import click


@click.command(
    help="Extract centromere coordinates from feature annotation BEDs.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will port get_centromeres.py logic."""
    raise click.ClickException(
        "`karyoscope centromeres` is not yet implemented in this development build."
    )
