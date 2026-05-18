"""``karyoscope karyotype`` — render karyotype visualizations.

This command is not yet implemented in v0.1.0.dev.
"""

from __future__ import annotations

import click


@click.command(
    help="Render karyotype visualizations (SVG/PNG/PDF) from feature annotations.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will port KaryoScope_assembly.py logic."""
    raise click.ClickException(
        "`karyoscope karyotype` is not yet implemented in this development build."
    )
