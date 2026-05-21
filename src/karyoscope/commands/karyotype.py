"""``karyoscope karyotype`` -- render karyotype SVG visualisations.

Top of the cascade. Given one or more FASTA inputs, this command runs
the whole pipeline through to a karyotype SVG: annotate -> seqtk telo
-> bin -> scaffold (mode='bed') -> (centromere mode only) centromeres
-> render. Existing intermediate files are reused; ``--no-auto``
turns missing inputs into hard errors.

Three render modes:

* ``full`` (default): whole-chromosome view, 1 Mb bins, 10 Mb scale bar.
* ``subtelomere``: zoomed view of p/q-arm telomeric ends, 100 bp
  bins, 10 kb scale bar. Only contigs flagged with at least one
  telomere appear.
* ``centromere``: zoomed view of each contig's centromere, 100 kb
  bins, 1 Mb scale bar. Only contigs with centromere coordinates
  appear.

One SVG is rendered per requested feature set (``--feature-set``,
repeatable; default: every feature set declared in the manifest).
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.commands.scaffold import _parse_named_path, _split_comma
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.core.karyotype_run import karyotype_run
from karyoscope.core.scaffold import DEFAULT_HUMAN_ACROCENTRICS
from karyoscope.core.scaffold_run import InputSpec
from karyoscope.exceptions import (
    CentromereError,
    DatabaseLayoutError,
    DatabaseNotFoundError,
    KaryoscopeError,
    KaryotypeError,
    ManifestError,
    ScaffoldError,
)

logger = logging.getLogger(__name__)


@click.command(
    help="Render karyotype SVG from per-input scaffolded annotations.",
    no_args_is_help=True,
)
@click.option(
    "--input",
    "-i",
    "inputs_raw",
    multiple=True,
    required=True,
    help="FASTA input. Repeat per haplotype. Form: 'NAME=PATH' or bare 'PATH'.",
)
@click.option(
    "--telo",
    "telo_raw",
    multiple=True,
    help="Optional precomputed seqtk telo output. Form: 'NAME=PATH'.",
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
    help="Override the database root directory.",
)
@click.option(
    "--feature-set",
    "feature_sets_arg",
    multiple=True,
    help="Feature set to render. Repeatable. Default: every feature set in the manifest "
    "(one SVG per set).",
)
@click.option(
    "--mode",
    "modes_raw",
    type=click.Choice(["genome", "subtelomere", "centromere"], case_sensitive=False),
    multiple=True,
    help="Which view(s) to render. Repeatable. Default: render every mode "
    "(genome, centromere, subtelomere).",
)
@click.option(
    "--sex",
    type=click.Choice(["male", "female", "reference", "unknown"], case_sensitive=False),
    default="unknown",
    show_default=True,
    help="Sample sex. 'unknown' draws sex-chromosome haps only where data is present.",
)
@click.option(
    "--sex-determination-system",
    type=click.Choice(["XY", "X0", "ZW", "ZO"], case_sensitive=False),
    default="XY",
    show_default=True,
    help="Sex-determination system.",
)
@click.option(
    "--background-color",
    type=click.Choice(["white", "black"], case_sensitive=False),
    default="white",
    show_default=True,
    help="Background colour. 'white' draws sequence outlines; 'black' uses light text.",
)
@click.option(
    "--bin-size",
    type=int,
    default=None,
    help="Bin size (bp) for the SVG. Default depends on --mode: 1Mb for genome, "
    "100Kb for centromere, 100bp for subtelomere. Only valid with exactly one --mode.",
)
@click.option(
    "--subtelomere-boundary",
    type=int,
    default=250_000,
    show_default=True,
    help="Subtelomere window size (bp). Only used in --mode subtelomere.",
)
@click.option(
    "--min-scaffold-length",
    type=int,
    default=5_000_000,
    show_default=True,
    help="Drop contigs shorter than this (no telomere) during the scaffold step.",
)
@click.option(
    "--acrocentric",
    "acrocentrics_raw",
    multiple=True,
    help="Chromosome name to treat as acrocentric during scaffold's flip decision. "
    "Repeatable; accepts comma-separated lists. Default: human acrocentrics with a warning.",
)
@click.option(
    "--no-human-chroms",
    "no_human_chroms",
    is_flag=True,
    default=False,
    help="Don't seed the chromosome list with the standard human set "
    "(chr1..chr22, chrX, chrY). Use for non-human assemblies so the SVG "
    "shows only the chromosomes actually in the data.",
)
@click.option(
    "--format",
    "formats_raw",
    type=click.Choice(["svg", "pdf", "png"], case_sensitive=False),
    multiple=True,
    help="Output format(s). Repeatable. Default: svg only. "
    "PDF and PNG are generated by converting the SVG via cairosvg.",
)
@click.option(
    "--sample-label",
    "sample_label",
    type=str,
    default=None,
    help="Sample label rendered in the SVG title band. Default: joined "
    "stems of the input FASTAs (e.g. 'HG002.maternal + HG002.paternal').",
)
@click.option(
    "--no-title",
    "no_title",
    is_flag=True,
    default=False,
    help="Don't draw the title band at the top of the SVG.",
)
@click.option(
    "--no-legend",
    "no_legend",
    is_flag=True,
    default=False,
    help="Don't draw the color legend in the right margin of the SVG.",
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
    "--auto/--no-auto",
    default=True,
    show_default=True,
    help="Auto-derive missing inputs. Disable to require everything upfront.",
)
@click.option(
    "--bgzip/--no-bgzip",
    default=True,
    show_default=True,
    help="bgzip the intermediate scaffolded BEDs (and the centromeres BED, "
    "when --mode includes centromere) that the cascade produces. The SVG / "
    "PDF / PNG outputs themselves are unaffected -- this only controls "
    "compression of the on-disk BED intermediates.",
)
@click.option(
    "--scaffolding/--no-scaffolding",
    default=True,
    show_default=True,
    help="Write full-resolution scaffolded BEDs to disk during the cascade. "
    "Pass --no-scaffolding to skip the expensive per-feature-set "
    "rewrite_bed step (saves ~5-10 min on whole-genome HG002): the "
    "scaffold_map.tsv is still written, and the map is applied at bin "
    "time so the binned-scaffolded BEDs and the final karyotype SVGs "
    "are equivalent. The full-resolution scaffolded BEDs simply aren't "
    "materialised as a side artifact.",
)
@click.option(
    "--outdir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Where to write the SVGs. Default: same directory as the first --input.",
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Explicit output path base. The mode and feature_set will be appended; "
    "with --output foo.svg you get foo.<dbid>.<mode>.<fs>.karyotype.svg. "
    "Conflicts with --outdir when both are set.",
)
def cmd(
    inputs_raw: tuple[str, ...],
    telo_raw: tuple[str, ...],
    split_haps_regex: str | None,
    db_id: str | None,
    db_root_arg: Path | None,
    feature_sets_arg: tuple[str, ...],
    modes_raw: tuple[str, ...],
    sex: str,
    sex_determination_system: str,
    background_color: str,
    bin_size: int | None,
    subtelomere_boundary: int,
    min_scaffold_length: int,
    acrocentrics_raw: tuple[str, ...],
    no_human_chroms: bool,
    formats_raw: tuple[str, ...],
    sample_label: str | None,
    no_title: bool,
    no_legend: bool,
    threads: int,
    auto: bool,
    bgzip: bool,
    scaffolding: bool,
    output_dir: Path | None,
    output_path: Path | None,
) -> None:
    """Render karyotype SVGs.

    By default, every (mode, feature_set) combination is rendered --
    three modes (genome, centromere, subtelomere) times every feature
    set declared in the database's manifest. Restrict either axis
    with ``--mode`` and/or ``--feature-set`` (both repeatable).

    \b
    Examples:
        # Diploid male, every mode for every feature set (the default)
        karyoscope karyotype -i hap1=HG002.hap1.fa.gz -i hap2=HG002.hap2.fa.gz \\
                              --sex male

        # Just the genome view for the chromosome feature set
        karyoscope karyotype -i hap1=HG002.hap1.fa.gz -i hap2=HG002.hap2.fa.gz \\
                              --sex male --mode genome --feature-set chromosome

        # All modes but only the region feature set
        karyoscope karyotype -i hap1=HG002.hap1.fa.gz -i hap2=HG002.hap2.fa.gz \\
                              --sex male --feature-set region
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
            "no --acrocentric given; falling back to human acrocentrics %s.",
            sorted(DEFAULT_HUMAN_ACROCENTRICS),
        )

    if output_dir is not None and output_path is not None:
        raise click.UsageError("--outdir and --output cannot both be set")

    sex_resolved: str | None = None if sex.lower() == "unknown" else sex.lower()
    feature_sets = list(feature_sets_arg) if feature_sets_arg else None
    modes = [m.lower() for m in modes_raw] if modes_raw else None
    formats = [f.lower() for f in formats_raw] if formats_raw else None

    db_root = paths.ensure_db_root(db_root_arg)
    try:
        results = karyotype_run(
            inputs,
            db_root=db_root,
            db_id=db_id,
            feature_sets=feature_sets,
            modes=modes,
            sex=sex_resolved,
            sex_determination_system=sex_determination_system.upper(),
            background_color=background_color.lower(),
            formats=formats,
            sample_label=sample_label,
            show_title=not no_title,
            show_legend=not no_legend,
            bin_size=bin_size,
            subtelomere_boundary=subtelomere_boundary,
            min_scaffold_length=min_scaffold_length,
            acrocentrics=acrocentrics,
            split_haps_regex=split_haps_regex,
            threads=threads,
            auto=auto,
            bgzip=bgzip,
            scaffolding=scaffolding,
            output_dir=output_dir,
            output_path=output_path,
            seed_human_chromosomes=not no_human_chroms,
        )
    except (
        KaryotypeError,
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
    for r in results:
        click.echo(f"  [{r.mode}/{r.feature_set}]")
        for p in r.output_paths:
            click.echo(f"    {p}")
