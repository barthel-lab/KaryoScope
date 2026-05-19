"""``karyoscope bin`` — aggregate base-pair annotations into larger bins.

(The module is named ``bin_cmd`` rather than ``bin`` to avoid shadowing
Python's built-in :func:`bin`.)

Two ways to invoke it:

* **Database-mode**: ``--db DBID --feature-set FS`` resolves the
  database's ``hierarchy.tsv``, extracts the leaf set for the given
  feature set, and uses it to prioritise leaf features per bin (the
  recommended path; matches what ``annotate`` produced).
* **Bare**: omit both ``--db`` and ``--feature-set`` to skip leaf
  prioritisation entirely. The selection rule reduces to "feature
  with the largest within-bin overlap, with novel deprioritised".

Either path takes a sorted BED on ``-i`` and writes a binned BED to
``-o`` (or stdin / stdout via ``-``).

This is intentionally a small, single-purpose command: it is the
building block consumed by :mod:`karyoscope.commands.scaffold` (which
calls the same :func:`karyoscope.core.bin.bin_features` in-process to
produce the 1 Mb binned BEDs the orientation logic needs).
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.core.annotate import resolve_database
from karyoscope.core.bin import bin_features, leaves_for
from karyoscope.core.io.hierarchy import parse_hierarchy
from karyoscope.exceptions import (
    BinError,
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    ManifestError,
)
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


@click.command(
    help="Aggregate base-pair annotation BEDs into larger bins.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "input_path",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="Input BED file (sorted by chrom then start). Use '-' for stdin. .gz supported.",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="Output BED file. Use '-' for stdout. Output is gzipped iff the path ends in .gz.",
)
@click.option(
    "--bin-size",
    "-b",
    type=int,
    required=True,
    help="Bin size in base pairs (e.g. 1000000 for 1 Mb).",
)
@click.option(
    "--db",
    "db_id",
    type=str,
    default=None,
    help="Database id whose hierarchy.tsv defines the leaf-feature set. "
    "Default: the unique installed database if exactly one is installed and --feature-set is given.",
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
    "feature_set",
    type=str,
    default=None,
    help="Feature set to use for leaf prioritisation. Required when --db is given (or implied).",
)
def cmd(
    input_path: Path,
    output_path: Path,
    bin_size: int,
    db_id: str | None,
    db_root_arg: Path | None,
    feature_set: str | None,
) -> None:
    """Run the binner.

    \b
    Examples:
        # Bin a smoothed BED at 1 Mb, using the chromosome leaf set
        karyoscope bin -i sample.KS_v2.chromosome.smoothed.bed.gz \\
                       -o sample.KS_v2.chromosome.smoothed.binned1Mb.bed.gz \\
                       -b 1000000 --feature-set chromosome

        # Bare binning, no leaf prioritisation, stdin to stdout
        gunzip -c sample.bed.gz | karyoscope bin -i - -o - -b 100000
    """
    if bin_size < 1:
        raise click.UsageError("--bin-size must be a positive integer")

    leaf_set: set[str] | None = None
    if feature_set is not None:
        try:
            db_root = paths.ensure_db_root(db_root_arg)
            resolved_id, db_dir = resolve_database(db_root, db_id)
            manifest = validate_database_layout(db_dir)
            if feature_set not in manifest.feature_sets:
                raise click.UsageError(
                    f"feature set {feature_set!r} is not declared in database "
                    f"{resolved_id!r} (declared: {sorted(manifest.feature_sets)})"
                )
            hierarchy = parse_hierarchy(db_dir / manifest.hierarchy)
            leaf_set = leaves_for(hierarchy, feature_set)
            if not leaf_set:
                logger.warning(
                    "feature set %r has no leaf features in %s; falling back to "
                    "no leaf prioritisation",
                    feature_set,
                    manifest.hierarchy,
                )
                leaf_set = None
        except (
            DatabaseNotFoundError,
            DatabaseLayoutError,
            ManifestError,
            KaryoscopeError,
        ) as e:
            raise click.ClickException(str(e)) from e
    elif db_id is not None:
        raise click.UsageError("--db requires --feature-set; pass both or neither")

    try:
        bin_features(
            input_path=input_path,
            output_path=output_path,
            bin_size=bin_size,
            leaf_set=leaf_set,
        )
    except BinError as e:
        raise click.ClickException(str(e)) from e

    if str(output_path) != "-":
        click.echo(f"Wrote {output_path}")
