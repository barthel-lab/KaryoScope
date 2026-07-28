"""``karyoscope annotate`` — annotate a FASTA against a KaryoScope database.

Produces one or two BEDs per feature set declared in the database:

* ``<input>.<dbid>.<feature_set>.presmoothed.bed[.gz]`` — every k-mer
  translated through ``features.tsv`` to its name in the feature set,
  with adjacent same-name intervals merged. This is the "raw" output
  of the annotation.
* ``<input>.<dbid>.<feature_set>.smoothed.bed[.gz]`` — additionally
  runs the hierarchy-aware smoothing pass, which promotes noisy short
  intervals (typically ``novel`` runs flanked by real features) up to
  the lowest common ancestor of the flankers in the database's
  ``hierarchy.tsv``.

Default behaviour produces both outputs. ``--no-smooth`` skips the
smoothing pass and produces only the presmoothed BED.
``--no-keep-presmoothed`` skips writing the presmoothed BED and
produces only the smoothed one. Combining both flags is an error
(would produce no output).

Two sentinel feature names that may appear in the output:

* ``novel`` — k-mer not in the KMC index.
* ``categorized`` — only seen in the smoothed output, indicates an
  interval was promoted all the way up to the hierarchy's root.

The intermediate combined BED produced by the C++ helper is deleted by
default; pass ``--keep-intermediates`` to keep it for debugging.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope import progress as _progress
from karyoscope.commands._options import resolve_resource_check_flag
from karyoscope.core.annotate import annotate_batch
from karyoscope.exceptions import (
    KaryoscopeError,
)

logger = logging.getLogger(__name__)


@click.command(
    help="Annotate a FASTA against a KaryoScope database, producing per-feature-set BEDs.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "input_paths",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    multiple=True,
    help="Input sequence file. Accepts FASTA (.fasta/.fa/.fna, plain "
    "or .gz), FASTQ (.fastq/.fq, plain or .gz), or BAM (.bam). BAM "
    "inputs are piped through `samtools fasta` (requires samtools on "
    "PATH); no intermediate file is written. See --preserve-order for "
    "how the smoothing implementation adapts to the input type. "
    "Repeatable: pass -i several times to annotate multiple inputs in "
    "one run; on the HKS backend the index is loaded once per feature "
    "set for the whole cohort instead of once per input.",
)
@click.option(
    "--outdir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory to write output BEDs into. Default: same directory as "
    "--input (required when more than one -i is given). Output filenames are "
    "prefixed by each input's basename, so multiple inputs never collide.",
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
    help="Threads for both k-mer querying and smoothing. 0 means auto-detect.",
)
@click.option(
    "--k",
    "k",
    type=int,
    default=None,
    help="Query k-mer length. Defaults to the database's k. Only a variable-k "
    "index (built with `karyoscope build --variable-k`) accepts a value other "
    "than its k; use it for a k-sweep. Outputs are tagged .k<k> so runs into "
    "one directory don't collide.",
)
@click.option(
    "--smooth/--no-smooth",
    default=True,
    show_default=True,
    help="Produce the hierarchy-smoothed BED in addition to the presmoothed BED.",
)
@click.option(
    "--keep-presmoothed/--no-keep-presmoothed",
    default=True,
    show_default=True,
    help="Keep the presmoothed BED. Pass --no-keep-presmoothed to write only the smoothed output.",
)
@click.option(
    "--keep-intermediates",
    is_flag=True,
    default=False,
    help="Keep the combined .featureIDs.bed from the C++ step (useful for debugging).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Regenerate the combined intermediate even if a complete one already "
    "exists. By default a rerun reuses a verified combined BED left by a "
    "previous (e.g. OOM-killed) run and skips the get_featureIDs step, "
    "resuming straight into smoothing. A partial file from a killed run is "
    "never reused regardless of this flag.",
)
@click.option(
    "--bgzip/--no-bgzip",
    default=True,
    show_default=True,
    help="bgzip the per-feature-set output BEDs. Pass --no-bgzip to write plain .bed files. "
    "Note this shrinks the final output but not the peak disk usage: compression runs "
    "after every BED has been written in full.",
)
@click.option(
    "--no-resource-check",
    is_flag=True,
    default=False,
    help="Skip the up-front disk and memory checks. Disk: the output footprint "
    "is estimated from the input size (roughly 0.8 GB per feature set per Gbp), "
    "so pass this if that estimate is wrong for your input -- it is calibrated "
    "on human data and can be well off for other organisms. Memory: HKS "
    "databases need the index resident (~10 GB for a human database), read from "
    "the index files themselves; skipping that check risks the run being killed "
    "by the OS instead.",
)
@click.option(
    "--no-space-check",
    is_flag=True,
    default=False,
    hidden=True,
    help="Deprecated alias for --no-resource-check.",
)
@click.option(
    "--preserve-order/--no-preserve-order",
    "preserve_input_order",
    default=True,
    show_default=True,
    help="Write output BEDs with sequences in the same order as the input. "
    "The implementation strategy is chosen automatically based on the "
    "input file extension: FASTA inputs use per-sequence temp files "
    "(safe for whole-chromosome chunks); FASTQ/BAM inputs use streaming "
    "ordered dispatch (no temp files, scales to millions of reads). "
    "Pass --no-preserve-order for the fastest path when order doesn't "
    "matter downstream -- typically read data where you'll aggregate "
    "the results anyway, or long-read FASTA where the per-sequence "
    "temp files would scale poorly.",
)
def cmd(
    input_paths: tuple[Path, ...],
    output_dir: Path | None,
    db_id: str | None,
    db_root_arg: Path | None,
    feature_sets_arg: tuple[str, ...],
    threads: int,
    k: int | None,
    smooth: bool,
    keep_presmoothed: bool,
    keep_intermediates: bool,
    bgzip: bool,
    no_resource_check: bool,
    no_space_check: bool,
    preserve_input_order: bool,
    force: bool,
) -> None:
    """Run the annotate pipeline for one or more ``--input`` files.

    \b
    Examples:
        # Default: produce both presmoothed and smoothed BEDs
        karyoscope annotate -i my_assembly.fa.gz -o results/

        # Several inputs in one run (HKS: one index load per feature set)
        karyoscope annotate -i a.fa.gz -i b.fa.gz -i c.fa.gz -o results/

        # Only the presmoothed BED (skip smoothing)
        karyoscope annotate -i reads.fa.gz --no-smooth

        # Only the smoothed BED (don't keep presmoothed on disk)
        karyoscope annotate -i reads.fa.gz --no-keep-presmoothed

        # Pick a specific database and only the chromosome set
        karyoscope annotate -i reads.fa.gz --db KS_human_CHM13_v2 \\
                            --feature-set chromosome
    """
    skip_checks = resolve_resource_check_flag(no_resource_check, no_space_check, command="annotate")
    db_root = paths.ensure_db_root(db_root_arg)
    inputs = list(input_paths)
    if output_dir is None:
        if len(inputs) > 1:
            raise click.ClickException(
                "--outdir is required when annotating more than one input "
                "(multiple inputs share one output directory)."
            )
        output_dir = inputs[0].parent
    feature_sets = list(feature_sets_arg) if feature_sets_arg else None

    common = dict(
        output_dir=output_dir,
        db_root=db_root,
        db_id=db_id,
        feature_sets=feature_sets,
        threads=threads,
        k=k,
        smooth=smooth,
        keep_presmoothed=keep_presmoothed,
        keep_intermediates=keep_intermediates,
        bgzip=bgzip,
        preserve_input_order=preserve_input_order,
        force=force,
        check_space=not skip_checks,
        progress=_progress.from_context(),
    )

    try:
        # annotate_batch delegates a single input to the battle-tested
        # single-input annotate() and batches when there are several.
        results = annotate_batch(input_paths=inputs, **common)
    except KaryoscopeError as e:
        raise click.ClickException(str(e)) from e

    multi = len(results) > 1
    for input_path, result in results.items():
        if multi:
            click.echo(f"Wrote ({input_path.name}):")
        else:
            click.echo("Wrote:")
        for fs, path in result.presmoothed_paths.items():
            click.echo(f"  {fs} (presmoothed): {path}")
        for fs, path in result.smoothed_paths.items():
            click.echo(f"  {fs} (smoothed):    {path}")
        if result.combined_intermediate is not None:
            click.echo(f"  (intermediate: {result.combined_intermediate})")
