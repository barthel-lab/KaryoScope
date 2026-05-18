"""``karyoscope download`` — acquire pre-built databases from the registry.

This command is not yet implemented in v0.1.0.dev. Track progress at
https://github.com/barthel-lab/KaryoScope/issues
"""

from __future__ import annotations

import click


@click.command(
    help="Acquire pre-built databases from the KaryoScope registry.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will fetch and install databases from karyoscope-registry."""
    raise click.ClickException(
        "`karyoscope download` is not yet implemented in this development "
        "build. See https://github.com/barthel-lab/KaryoScope for status."
    )
