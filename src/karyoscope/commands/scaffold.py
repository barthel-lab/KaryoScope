"""``karyoscope scaffold`` — order, orient, and rename assembly contigs.

This command is not yet implemented in v0.1.0.dev.
"""

from __future__ import annotations

import click


@click.command(
    help="Order, orient, and rename assembly contigs into canonical scaffolds.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will port scaffold_stats and scaffold_features logic."""
    raise click.ClickException(
        "`karyoscope scaffold` is not yet implemented in this development build."
    )
