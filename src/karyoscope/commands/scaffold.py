"""``karyoscope scaffold`` -- order, orient, and rename assembly contigs.

Takes one or more FASTA inputs (one per haplotype, typically), reads
the per-feature-set smoothed BEDs that ``annotate`` produces, runs
``seqtk telo`` to flag telomere-bearing contigs, bins the chromosome
and region feature sets at 1 Mb, then classifies each contig to a
chromosome, orients it (reverse-complements if the evidence says it's
reversed), and orders contigs into the canonical chromosome x hap x
category x descending-length sequence the karyotype renderer expects.

Outputs per input (named after the input FASTA's stem):

* ``<stem>.<dbid>.scaffold_map.tsv`` -- the authoritative map from
  scaffolded contig name back to source. Downstream stages
  (centromeres, karyotype) parse this rather than the encoded name.
* ``<stem>.<dbid>.scaffold_stats.tsv`` -- the legacy 2-column
  ``<name>\\t<TPCQT>`` format kept for back-compat.
* ``<stem>.<dbid>.<feature_set>.smoothed.scaffolded.bed[.gz]`` --
  one per feature set, rewritten with the new contig names.

The command auto-derives missing prerequisites: if the annotation
BEDs aren't there it runs ``annotate``; if the telomere file is
missing it runs ``seqtk telo``; if the 1 Mb binned BEDs aren't there
it runs ``bin`` in-process. Pass ``--no-auto`` to turn missing
prerequisites into hard errors instead.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.core.scaffold import DEFAULT_HUMAN_ACROCENTRICS
from karyoscope.core.scaffold_run import InputSpec, scaffold_run
from karyoscope.exceptions import (
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    ManifestError,
    ScaffoldError,
)

logger = logging.getLogger(__name__)


def _parse_named_path(value: str) -> tuple[str | None, Path]:
    """Parse ``NAME=PATH`` or bare ``PATH``.

    Returns ``(None, Path(value))`` when no ``=`` is present.
    """
    if "=" in value:
        name, _, path = value.partition("=")
        name = name.strip()
        if not name:
            raise click.BadParameter(f"empty name in {value!r}; use NAME=PATH or just PATH")
        return name, Path(path)
    return None, Path(value)


def _split_comma(value: str) -> list[str]:
    """Split a comma-separated string into stripped non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


@click.command(
    help="Order, orient, and rename assembly contigs into canonical scaffolds.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "inputs_raw",
    multiple=True,
    required=True,
    help="FASTA input. Repeat per haplotype. Form: 'NAME=PATH' (e.g. 'hap1=hap1.fa.gz') "
    "or bare 'PATH' to auto-infer the label from the filename stem.",
)
@click.option(
    "--telo",
    "telo_raw",
    multiple=True,
    help="Optional precomputed seqtk telo output. Form: 'NAME=PATH'. "
    "Without this, scaffold runs seqtk telo on each input automatically.",
)
@click.option(
    "--split-haps",
    "split_haps_regex",
    type=str,
    default=None,
    help="Optional regex applied per contig name; capture group 1 is the hap label. "
    "Overrides the built-in patterns (hifiasm h[12]tg, hap1/hap2, maternal/paternal).",
)
@click.option(
    "--db",
    "db_id",
    type=str,
    default=None,
    help="Database id to use. Default: the unique installed database if exactly one is installed.",
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
    help="Restrict the scaffolded-BED outputs to this feature set. Repeatable. "
    "Default: every feature set declared in the database's manifest.",
)
@click.option(
    "--bin-size",
    type=int,
    default=1_000_000,
    show_default=True,
    help="Bin size (bp) for the orientation BEDs. The manuscript uses 1 Mb.",
)
@click.option(
    "--min-scaffold-length",
    type=int,
    default=5_000_000,
    show_default=True,
    help="Drop contigs shorter than this that have no telomere.",
)
@click.option(
    "--acrocentric",
    "acrocentrics_raw",
    multiple=True,
    help="Chromosome name to treat as acrocentric in the flip decision. Repeatable; "
    "accepts comma-separated lists. Default: human acrocentrics (chr13/14/15/21/22) "
    "with a warning to set it explicitly for non-human assemblies.",
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
    "--mode",
    type=click.Choice(["fasta", "bed", "both"], case_sensitive=False),
    default="fasta",
    show_default=True,
    help="What to write per input. 'fasta' writes a scaffolded FASTA only; "
    "'bed' writes per-feature-set scaffolded BEDs only (used by `karyoscope karyotype`); "
    "'both' writes both. The map and legacy stats files are always written.",
)
@click.option(
    "--keep-unscaffolded/--drop-unscaffolded",
    default=True,
    show_default=True,
    help="In FASTA mode, append contigs that did not get scaffolded (too short, no leaf "
    "chromosome) at the end of the output under their original names. Disable to drop them entirely.",
)
@click.option(
    "--bgzip/--no-bgzip",
    default=True,
    show_default=True,
    help="bgzip the scaffolded output BEDs and FASTA.",
)
@click.option(
    "--auto/--no-auto",
    default=True,
    show_default=True,
    help="Auto-derive missing inputs (annotate, seqtk telo, bin). Disable to require everything upfront.",
)
@click.option(
    "--outdir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Where to write scaffolded outputs. Default: same directory as each --input.",
)
def cmd(
    inputs_raw: tuple[str, ...],
    telo_raw: tuple[str, ...],
    split_haps_regex: str | None,
    db_id: str | None,
    db_root_arg: Path | None,
    feature_sets_arg: tuple[str, ...],
    bin_size: int,
    min_scaffold_length: int,
    acrocentrics_raw: tuple[str, ...],
    mode: str,
    keep_unscaffolded: bool,
    threads: int,
    bgzip: bool,
    auto: bool,
    output_dir: Path | None,
) -> None:
    """Run the scaffolder.

    \b
    Examples:
        # Pangenome convention: hap1 and hap2 in separate files
        karyoscope scaffold -i hap1=HG002.hap1.fa.gz -i hap2=HG002.hap2.fa.gz

        # HG002 distributed as a single combined FASTA (auto-detected diploid)
        karyoscope scaffold -i HG002.fa.gz

        # Explicit unassigned set
        karyoscope scaffold -i hap1=h1.fa.gz -i hap2=h2.fa.gz -i unassigned=other.fa.gz

        # Bring your own telomere file
        karyoscope scaffold -i hap1=h1.fa.gz --telo hap1=h1.telo
    """
    # --- parse named-pair options -----------------------------------
    parsed_inputs: list[tuple[str | None, Path]] = []
    for raw in inputs_raw:
        parsed_inputs.append(_parse_named_path(raw))

    parsed_telos: dict[str, Path] = {}
    for raw in telo_raw:
        name, path = _parse_named_path(raw)
        if name is None:
            raise click.UsageError(
                f"--telo requires NAME=PATH form (got {raw!r}); the name should match an --input"
            )
        if name in parsed_telos:
            raise click.UsageError(f"duplicate --telo for name {name!r}")
        parsed_telos[name] = path

    # Build InputSpec list. We pass telo paths through alongside the
    # explicit names so scaffold_run can later match telo->input by
    # the same name.
    inputs: list[InputSpec] = []
    explicit_names = {name for name, _ in parsed_inputs if name is not None}
    for name, path in parsed_inputs:
        telo_path = None
        if name is not None and name in parsed_telos:
            telo_path = parsed_telos.pop(name)
        inputs.append(InputSpec(name=name, path=path, telo_path=telo_path))
    if parsed_telos:
        unmatched = sorted(parsed_telos.keys())
        raise click.UsageError(
            f"--telo entries had no matching --input: {unmatched} "
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

    feature_sets = list(feature_sets_arg) if feature_sets_arg else None

    mode_normalised = mode.lower()
    if mode_normalised == "fasta" and feature_sets:
        raise click.UsageError(
            "--feature-set has no effect with --mode fasta (no per-feature-set BEDs are written). "
            "Use --mode bed or --mode both, or drop --feature-set."
        )

    db_root = paths.ensure_db_root(db_root_arg)
    try:
        results = scaffold_run(
            inputs,
            db_root=db_root,
            db_id=db_id,
            feature_sets=feature_sets,
            mode=mode_normalised,
            bin_size=bin_size,
            min_scaffold_length=min_scaffold_length,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            bgzip=bgzip,
            keep_unscaffolded=keep_unscaffolded,
            auto=auto,
            output_dir=output_dir,
        )
    except (
        ScaffoldError,
        DatabaseNotFoundError,
        DatabaseLayoutError,
        ManifestError,
        ToolNotFoundError,
        ExternalToolError,
        KaryoscopeError,
    ) as e:
        raise click.ClickException(str(e)) from e

    # --- summarise on stdout ----------------------------------------
    click.echo("Wrote:")
    for input_name, result in results.items():
        click.echo(f"  [{result.hap_label}] {input_name}")
        click.echo(f"    map:    {result.map_path}")
        click.echo(f"    stats:  {result.stats_path}")
        if result.scaffolded_fasta is not None:
            click.echo(f"    fasta:  {result.scaffolded_fasta}")
        for fs, p in result.scaffolded_beds.items():
            click.echo(f"    {fs}:  {p}")
