"""``karyoscope remap-bed`` -- apply an existing scaffold map to a foreign BED.

``scaffold`` builds a ``scaffold_map.tsv`` and rewrites the BEDs it annotated
in the same run. This command does the *standalone* version: take a BED that
was annotated separately -- possibly against a **different database** than the
one used to derive the map -- and rewrite it into the scaffolded coordinate
system using a previously-built map.

The motivating case is a two-database workflow: scaffold/lay out an assembly
with a roles-bearing database (e.g. ``KS_human_CHM13_v2``), annotate a feature
set from a plot-only database (e.g. a cytoband database) on the *original*
FASTA, then ``remap-bed`` that annotation into scaffold coordinates so it lines
up with the region BEDs. It is the file-producing counterpart of the in-memory
remap ``karyotype --scaffold-db`` performs while rendering.

The heavy lifting (and the BED/map compatibility checks) lives in
:func:`karyoscope.core.scaffold.remap_bed_with_map`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.core.scaffold import remap_bed_with_map
from karyoscope.exceptions import KaryoscopeError, ScaffoldError
from karyoscope.installed import resolve_database
from karyoscope.manifest import validate_database_layout

logger = logging.getLogger(__name__)


@click.command(
    help="Apply an existing scaffold map to a separately-annotated BED.",
    no_args_is_help=True,
)
@click.option(
    "--bed",
    "-b",
    "bed_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Annotation BED in original (unscaffolded) contig coordinates. .gz supported.",
)
@click.option(
    "--map",
    "-m",
    "map_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="A scaffold_map.tsv written by `karyoscope scaffold` (or `karyotype`).",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="Output BED in scaffolded coordinates. bgzipped iff the path ends in .gz.",
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Promote advisory checks (filename-stem mismatch; map contigs with no "
    "records in the BED) from warnings to hard errors.",
)
@click.option(
    "--db",
    "db_id",
    default=None,
    help="Database the BED was annotated against. With --feature-set, its manifest "
    "says how labels change on a flipped (reverse-complemented) contig -- required "
    "for a strand feature set, whose plus/minus must swap.",
)
@click.option(
    "--feature-set",
    "feature_set",
    default=None,
    help="Feature set the BED holds (looked up in --db's manifest feature_set_options).",
)
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override the database root directory (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
def cmd(
    bed_path: Path,
    map_path: Path,
    output_path: Path,
    strict: bool,
    db_id: str | None,
    feature_set: str | None,
    db_root_arg: Path | None,
) -> None:
    """Rewrite ``--bed`` into scaffold coordinates using ``--map``.

    \b
    Example (two-database: layout from KS_human_CHM13_v2, plot cytoband):
        karyoscope scaffold -i hap1=hap1.fa.gz --db KS_human_CHM13_v2 --mode bed
        karyoscope annotate -i hap1.fa.gz --db KS_human_CHM13_cytoband
        karyoscope remap-bed \\
            -b hap1.KS_human_CHM13_cytoband.cytoband.smoothed.bed.gz \\
            -m hap1.KS_human_CHM13_v2.scaffold_map.tsv \\
            -o hap1.cytoband.scaffolded.bed.gz
    """
    flip_labels: dict[str, str] | None = None
    if (db_id is None) != (feature_set is None):
        raise click.UsageError("--db and --feature-set go together")
    try:
        if db_id is not None:
            _, db_dir = resolve_database(paths.ensure_db_root(db_root_arg), db_id)
            manifest = validate_database_layout(db_dir)
            if feature_set not in manifest.feature_sets:
                raise click.UsageError(
                    f"feature set {feature_set!r} is not in {db_id!r}: {manifest.feature_sets!r}"
                )
            flip_labels = manifest.flip_labels.get(feature_set)
        stats = remap_bed_with_map(
            bed_path, output_path, map_path, strict=strict, flip_labels=flip_labels
        )
    except ScaffoldError as e:
        raise click.ClickException(str(e)) from e
    except KaryoscopeError as e:
        raise click.ClickException(str(e)) from e

    click.echo(
        f"Wrote {output_path} "
        f"({stats.mapped_contigs}/{stats.bed_contigs} contigs placed, "
        f"{stats.dropped_contigs} dropped as unscaffolded)"
    )
