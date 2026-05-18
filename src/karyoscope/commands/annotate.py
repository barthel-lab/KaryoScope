"""``karyoscope annotate`` — annotate sequences with k-mer features.

This command is not yet implemented in v0.1.0.dev.
"""

from __future__ import annotations

import click


@click.command(
    help="Annotate sequences (FASTA/FASTQ/BAM) with k-mer features.",
    no_args_is_help=False,
)
def cmd() -> None:
    """Stub. Will wrap get_featureIDs and smooth_features."""
    raise click.ClickException(
        "`karyoscope annotate` is not yet implemented in this development build."
    )
