"""``karyoscope annotate`` — annotate a FASTA against a KaryoScope database.

Produces one BED per feature set declared in the database. Each BED has
the schema ``seq_name TAB start TAB end TAB feature_name``, with feature
ids translated to their per-set names using the database's
``features.tsv``. Two sentinels:

* ``novel`` — k-mer not in the KMC index.
* ``Unknown`` — k-mer in the index but its feature id has no row in
  features.tsv. Usually a database / index mismatch, not a user error.

The output BEDs are bgzipped by default; pass ``--no-bgzip`` to write
plain ``.bed`` files. The intermediate combined BED produced by the
C++ helper is deleted by default; pass ``--keep-intermediates`` to keep
it for debugging.

Smoothing is not applied here — that's a separate stage of the pipeline
(``karyoscope smooth``, coming in 5c). The output filename includes
``.presmoothed.`` to make this explicit.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.core.annotate import annotate
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.exceptions import (
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    ManifestError,
)

logger = logging.getLogger(__name__)


@click.command(
    help="Annotate a FASTA against a KaryoScope database, producing per-feature-set BEDs.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Input FASTA file (plain or gzipped).",
)
@click.option(
    "--outdir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory to write output BEDs into. Default: same directory as --input.",
)
@click.option(
    "--db",
    "db_id",
    type=str,
    default=None,
    help="Database id to use (e.g., KS_human_CHM13_v2). "
    "Default: the unique installed database if there's exactly one.",
)
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override the database root directory (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
@click.option(
    "--feature-set",
    "feature_sets_arg",
    multiple=True,
    help="Restrict output to this feature set. Repeatable. "
    "Default: all feature sets declared in the database's manifest.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=0,
    show_default=True,
    help="Threads to use for k-mer querying. 0 means auto-detect.",
)
@click.option(
    "--keep-intermediates",
    is_flag=True,
    default=False,
    help="Keep the combined .featureIDs.bed from the C++ step (useful for debugging).",
)
@click.option(
    "--bgzip/--no-bgzip",
    default=True,
    show_default=True,
    help="bgzip the per-feature-set output BEDs. Pass --no-bgzip to write plain .bed files.",
)
def cmd(
    input_path: Path,
    output_dir: Path | None,
    db_id: str | None,
    db_root_arg: Path | None,
    feature_sets_arg: tuple[str, ...],
    threads: int,
    keep_intermediates: bool,
    bgzip: bool,
) -> None:
    """Run the annotate pipeline for ``--input``.

    \b
    Examples:
        # Annotate using the only installed database
        karyoscope annotate -i my_assembly.fa.gz -o results/

        # Pick a specific database and only the chromosome set
        karyoscope annotate -i reads.fa.gz --db KS_human_CHM13_v2 \\
                            --feature-set chromosome
    """
    db_root = paths.ensure_db_root(db_root_arg)
    if output_dir is None:
        output_dir = input_path.parent
    feature_sets = list(feature_sets_arg) if feature_sets_arg else None

    try:
        result = annotate(
            input_path=input_path,
            output_dir=output_dir,
            db_root=db_root,
            db_id=db_id,
            feature_sets=feature_sets,
            threads=threads,
            keep_intermediates=keep_intermediates,
            bgzip=bgzip,
        )
    except (
        DatabaseNotFoundError,
        DatabaseLayoutError,
        ManifestError,
        ToolNotFoundError,
        ExternalToolError,
        KaryoscopeError,
    ) as e:
        raise click.ClickException(str(e)) from e

    click.echo("Wrote:")
    for fs, path in result.output_paths.items():
        click.echo(f"  {fs}: {path}")
    if result.combined_intermediate is not None:
        click.echo(f"  (intermediate: {result.combined_intermediate})")
