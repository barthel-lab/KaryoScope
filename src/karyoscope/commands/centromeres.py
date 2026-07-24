"""``karyoscope centromeres`` -- extract per-contig centromere coordinates.

Consumes the scaffolded binned region BEDs that ``karyoscope scaffold``
produces (or auto-derives them) and walks each contig's bins to find
the min/max coordinates of bins classified as centromeric. Output is
a 3-column BED (contig, start, end) per input -- coordinate-only, no
feature column, so the filename has no ``<feature_set>`` segment.

Like ``karyoscope scaffold``, this command auto-derives its
prerequisites: missing scaffolded BED triggers ``scaffold_run``
(which itself cascades through annotate, seqtk telo, and bin);
missing binned scaffolded BEDs trigger ``bin_features`` in-process.
Pass ``--no-auto`` to turn missing inputs into hard errors.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.commands.scaffold import _parse_named_path, _split_comma
from karyoscope.core.centromeres import (
    DEFAULT_COARSE_BIN_SIZE,
    DEFAULT_FINE_BIN_SIZE,
    centromeres_run,
)
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.core.scaffold import DEFAULT_HUMAN_ACROCENTRICS
from karyoscope.core.scaffold_run import InputSpec
from karyoscope.exceptions import (
    CentromereError,
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    ManifestError,
    ScaffoldError,
)

logger = logging.getLogger(__name__)


@click.command(
    help="Extract per-contig centromere coordinates from a genome assembly.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "inputs_raw",
    multiple=True,
    required=True,
    help="FASTA-format genome assembly. Repeat per haplotype. Form: "
    "'NAME=PATH' or bare 'PATH' to auto-infer the label from the filename "
    "stem. Read-level inputs (FASTQ / BAM) are rejected -- centromere "
    "detection needs scaffolded assembly contigs, which reads don't "
    "provide. Use `karyoscope annotate` for read-level annotation.",
)
@click.option(
    "--telo",
    "telo_raw",
    multiple=True,
    help="Optional precomputed seqtk telo output. Form: 'NAME=PATH'. "
    "Without this, the cascade runs seqtk telo internally.",
)
@click.option(
    "--split-haps",
    "split_haps_regex",
    type=str,
    default=None,
    help="Optional regex applied per contig name; capture group 1 is the hap label.",
)
@click.option(
    "--db",
    "db_id",
    type=str,
    default=None,
    help="Database id. Default: the unique installed database if exactly one is installed.",
)
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override the database root directory (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
@click.option(
    "--centromere-feature-set",
    "centromere_feature_set",
    type=str,
    default=None,
    help="Override which feature set drives centromere detection. Default: "
    "manifest.roles.centromere_detection, with chain fallback to "
    "roles.region_assignment and then the literal 'region'.",
)
@click.option(
    "--coarse-bin-size",
    type=int,
    default=DEFAULT_COARSE_BIN_SIZE,
    show_default=True,
    help="Bin size (bp) for the coarse pass. 1 Mb matches the manuscript benchmarks.",
)
@click.option(
    "--fine-bin-size",
    type=int,
    default=DEFAULT_FINE_BIN_SIZE,
    show_default=True,
    help="Bin size (bp) for the optional fine-refinement pass. Pass 0 to disable refinement.",
)
@click.option(
    "--min-scaffold-length",
    type=int,
    default=5_000_000,
    show_default=True,
    help="Drop contigs shorter than this (no telomere) during the scaffold step.",
)
@click.option(
    "--telo-motif",
    default=None,
    help="Telomere repeat motif for `seqtk telo` (its -m). Default: seqtk's CCCTAA "
    "(vertebrate TTAGGG). Plants (e.g. Arabidopsis) use CCCTAAA (TTTAGGG).",
)
@click.option(
    "--acrocentric",
    "acrocentrics_raw",
    multiple=True,
    help="Chromosome name to treat as acrocentric during scaffold's flip decision. "
    "Repeatable; accepts comma-separated lists. Default: human acrocentrics with a warning.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=0,
    show_default=True,
    help="Threads for auto-run annotate invocations.",
)
@click.option(
    "--bgzip/--no-bgzip",
    default=True,
    show_default=True,
    help="bgzip the output centromere BED.",
)
@click.option(
    "--auto/--no-auto",
    default=True,
    show_default=True,
    help="Auto-derive missing inputs (scaffold cascade, bin). Disable to require everything upfront.",
)
@click.option(
    "--outdir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Where to write the centromere BEDs. Default: same directory as each --input.",
)
def cmd(
    inputs_raw: tuple[str, ...],
    telo_raw: tuple[str, ...],
    split_haps_regex: str | None,
    db_id: str | None,
    db_root_arg: Path | None,
    centromere_feature_set: str | None,
    coarse_bin_size: int,
    fine_bin_size: int,
    min_scaffold_length: int,
    telo_motif: str | None,
    acrocentrics_raw: tuple[str, ...],
    threads: int,
    bgzip: bool,
    auto: bool,
    output_dir: Path | None,
) -> None:
    """Find per-contig centromere coordinates.

    \b
    This command is for genome assemblies. It cascades through scaffold
    (which itself needs contig names from a FASTA) to produce per-contig
    centromere coordinates; read-level inputs (FASTQ / BAM) have no
    chromosome-scale contig concept and are rejected with a clear error.
    For read-level annotation, use `karyoscope annotate` instead.

    \b
    Examples:
        # Full cascade from FASTA
        karyoscope centromeres -i hap1=HG002.hap1.fa.gz -i hap2=HG002.hap2.fa.gz

        # Coarse-only (skip fine-refinement pass)
        karyoscope centromeres -i CHM13.fa.gz --fine-bin-size 0
    """
    # --- parse named-pair options -----------------------------------
    parsed_inputs: list[tuple[str | None, Path]] = [_parse_named_path(raw) for raw in inputs_raw]
    parsed_telos: dict[str, Path] = {}
    for raw in telo_raw:
        name, path = _parse_named_path(raw)
        if name is None:
            raise click.UsageError(f"--telo requires NAME=PATH form (got {raw!r})")
        if name in parsed_telos:
            raise click.UsageError(f"duplicate --telo for name {name!r}")
        parsed_telos[name] = path

    explicit_names = {name for name, _ in parsed_inputs if name is not None}
    inputs: list[InputSpec] = []
    for name, path in parsed_inputs:
        telo_path = None
        if name is not None and name in parsed_telos:
            telo_path = parsed_telos.pop(name)
        inputs.append(InputSpec(name=name, path=path, telo_path=telo_path))
    if parsed_telos:
        raise click.UsageError(
            f"--telo entries had no matching --input: {sorted(parsed_telos)} "
            f"(known input names: {sorted(explicit_names)})"
        )

    # --- acrocentrics flag ------------------------------------------
    acrocentrics: set[str] | None
    if acrocentrics_raw:
        flattened: list[str] = []
        for entry in acrocentrics_raw:
            flattened.extend(_split_comma(entry))
        if not flattened:
            raise click.UsageError("--acrocentric was given but produced no chromosome names")
        acrocentrics = set(flattened)
    else:
        acrocentrics = None
        logger.warning(
            "no --acrocentric given; falling back to human acrocentrics %s. "
            "Pass --acrocentric explicitly for non-human assemblies.",
            sorted(DEFAULT_HUMAN_ACROCENTRICS),
        )

    # --- bin-size validation ----------------------------------------
    if coarse_bin_size < 1:
        raise click.UsageError("--coarse-bin-size must be a positive integer")
    if fine_bin_size < 0:
        raise click.UsageError("--fine-bin-size must be >= 0 (0 disables refinement)")
    fine: int | None = fine_bin_size if fine_bin_size > 0 else None

    db_root = paths.ensure_db_root(db_root_arg)
    try:
        results = centromeres_run(
            inputs,
            db_root=db_root,
            db_id=db_id,
            centromere_feature_set=centromere_feature_set,
            coarse_bin_size=coarse_bin_size,
            fine_bin_size=fine,
            min_scaffold_length=min_scaffold_length,
            telo_motif=telo_motif,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            bgzip=bgzip,
            auto=auto,
            output_dir=output_dir,
        )
    except (
        CentromereError,
        ScaffoldError,
        DatabaseNotFoundError,
        DatabaseLayoutError,
        ManifestError,
        ToolNotFoundError,
        ExternalToolError,
        KaryoscopeError,
    ) as e:
        raise click.ClickException(str(e)) from e

    click.echo("Wrote:")
    for input_name, result in results.items():
        click.echo(f"  [{result.hap_label or '?'}] {input_name}")
        click.echo(f"    centromeres: {result.centromeres_bed} ({result.num_contigs} contigs)")
