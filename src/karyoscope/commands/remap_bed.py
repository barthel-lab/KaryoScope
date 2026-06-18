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

from karyoscope.core.scaffold import remap_bed_with_map
from karyoscope.exceptions import KaryoscopeError, ScaffoldError

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
    help="Output BED in scaffolded coordinates. Gzipped iff the path ends in .gz.",
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Promote advisory checks (filename-stem mismatch; map contigs with no "
    "records in the BED) from warnings to hard errors.",
)
def cmd(bed_path: Path, map_path: Path, output_path: Path, strict: bool) -> None:
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
    try:
        stats = remap_bed_with_map(bed_path, output_path, map_path, strict=strict)
    except ScaffoldError as e:
        raise click.ClickException(str(e)) from e
    except KaryoscopeError as e:
        raise click.ClickException(str(e)) from e

    click.echo(
        f"Wrote {output_path} "
        f"({stats.mapped_contigs}/{stats.bed_contigs} contigs placed, "
        f"{stats.dropped_contigs} dropped as unscaffolded)"
    )
