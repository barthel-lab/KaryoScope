"""``karyoscope prep-bed`` — turn a source annotation into a feature-set BED.

One subcommand per *source format*, because unrelated formats produce the same
kind of feature set: RepeatMasker output and an EDTA GFF3 both yield a
``repeat`` set but share no parsing. Keying on the format means each subcommand
carries exactly the options that apply to it, with no flags that are silently
ignored depending on the input.

Each subcommand writes the BED (and, where the format implies one, a hierarchy)
and prints the matching ``feature_sets:`` stanza on **stdout**, with progress
and warnings on stderr — so the stanza can be redirected into a build spec
without capturing the commentary.

What these subcommands deliberately do *not* do is gap-fill, flatten overlaps,
or drop sequences: ``karyoscope build`` already owns all three, via
``background:``, ``flatten:`` and ``exclude:``.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import click

from karyoscope.core.prep import asat as asat_prep
from karyoscope.core.prep import cytoband as cytoband_prep
from karyoscope.core.prep import genes as genes_prep
from karyoscope.core.prep import pave as pave_prep
from karyoscope.core.prep import repeats as repeats_prep
from karyoscope.core.prep import structural as structural_prep
from karyoscope.core.prep.common import PrepError, PrepResult, SeqidRewriter, render_stanza

logger = logging.getLogger(__name__)

_IN = dict(exists=True, dir_okay=False, path_type=Path)
_OUT = dict(dir_okay=False, writable=True, path_type=Path)


def _seqid_options(func: Callable[..., Any]) -> Callable[..., Any]:
    """Attach the two seqid-rewriting options every converter shares."""
    func = click.option(
        "--seqid-map",
        type=click.Path(**_IN),
        default=None,
        help="Two-column 'old new' table renaming annotation seqids to match the assembly "
        "(e.g. NC_060925.1 -> chr1). Takes precedence over --rename-prefix.",
    )(func)
    return click.option(
        "--rename-prefix",
        metavar="OLD:NEW",
        default=None,
        help="Rewrite a leading seqid prefix, e.g. 'Col-CEN_chr:Chr'. Names not starting "
        "with OLD are left alone.",
    )(func)


def _rewriter(rename_prefix: str | None, seqid_map: Path | None) -> SeqidRewriter:
    """Build the rewriter, reporting a malformed option as a usage error.

    This runs while assembling ``_run``'s arguments, outside its handler, so it
    has to translate the error itself or the user gets a bare traceback.
    """
    try:
        return SeqidRewriter.build(rename_prefix=rename_prefix, seqid_map=seqid_map)
    except PrepError as e:
        raise click.UsageError(str(e)) from e


def _report(result: PrepResult) -> None:
    """Commentary to stderr, the pasteable stanza to stdout."""
    for note in result.notes:
        click.echo(f"note: {note}", err=True)
    click.echo(result.summary(), err=True)
    click.echo(err=True)
    click.echo(render_stanza(result))


def _guard(paths: list[Path | None], force: bool) -> None:
    existing = [str(p) for p in paths if p is not None and p.exists()]
    if existing and not force:
        raise click.ClickException(
            "refusing to overwrite (pass --force): " + ", ".join(sorted(existing))
        )


def _run(
    fn: Callable[..., PrepResult], outputs: list[Path | None], force: bool, **kwargs: Any
) -> None:
    _guard(outputs, force)
    try:
        result = fn(**kwargs)
    except PrepError as e:
        raise click.ClickException(str(e)) from e
    _report(result)


@click.group(
    "prep-bed",
    help="Convert a source annotation (RepeatMasker, EDTA, GFF3/GTF, UCSC cytoband, .fai) "
    "into a feature-set BED for `karyoscope build`.",
    no_args_is_help=True,
)
def cmd() -> None:
    """Prepare feature-set BEDs from source annotation formats."""


# -- repeat sets ------------------------------------------------------


@cmd.command("repeatmasker")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="RepeatMasker output: the native .out table or the UCSC BED repackaging "
    "(auto-detected; plain or gzipped).",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy.",
)
@click.option(
    "--colors",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a colours file using the reference KaryoScope repeat palette, "
    "with the five RNA leaves sharing one legend row.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file, for a priority-mode build.",
)
@click.option("--name", default="repeat", show_default=True, help="Feature-set name in the stanza.")
@click.option(
    "--background",
    default="nonrepeat",
    show_default=True,
    help="Gap-fill label build should use for the bases no repeat covers.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def repeatmasker(
    input_path: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None,
    priority: Path | None,
    name: str,
    background: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """RepeatMasker output -> a `repeat` feature set.

    \b
    Leaves are the RepeatMasker classes; the hierarchy and palette reproduce the
    shipped HKS_human_CHM13_v2 `repeat` set. Unrecognised classes are labelled
    `Unknown` rather than dropped, so their bases are never mistaken for
    non-repeat by the gap-fill.

    \b
    Example:
        karyoscope prep-bed repeatmasker --input rm.bed \\
            --output repeat.bed --hierarchy repeat.tsv --colors repeat_colors.tsv
    """
    _run(
        repeats_prep.from_repeatmasker,
        [output, hierarchy, colors, priority],
        force,
        input_path=input_path,
        output=output,
        hierarchy=hierarchy,
        colors=colors,
        priority=priority,
        name=name,
        background=background,
        rename=_rewriter(rename_prefix, seqid_map),
    )


@cmd.command("edta")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="EDTA TE GFF3 (plain or gzipped).",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file, for a priority-mode build.",
)
@click.option("--name", default="repeat", show_default=True, help="Feature-set name in the stanza.")
@click.option(
    "--background",
    default="nonrepeat",
    show_default=True,
    help="Gap-fill label build should use for the bases no TE covers.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def edta(
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None,
    name: str,
    background: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """An EDTA TE GFF3 -> a `repeat` feature set.

    \b
    EDTA's `Classification=` vocabulary aliases each superfamily under both
    spelled-out and Wicker three-letter names (DNA/HAT and DNA/DTA are both
    hAT); those are normalised to one leaf per superfamily so the BED labels and
    the hierarchy leaves agree by construction.

    \b
    No colours are written: there is no established palette for the EDTA
    vocabulary, so colours are left to build's automatic assignment.
    """
    _run(
        repeats_prep.from_edta,
        [output, hierarchy, priority],
        force,
        input_path=input_path,
        output=output,
        hierarchy=hierarchy,
        priority=priority,
        name=name,
        background=background,
        rename=_rewriter(rename_prefix, seqid_map),
    )


# -- gene set ---------------------------------------------------------


@cmd.command("gff-gene")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="Gene annotation in GFF3 or GTF (plain or gzipped; dialect auto-detected).",
)
@click.option(
    "--lengths",
    type=click.Path(**_IN),
    required=True,
    help="samtools .fai (or 2-column sizes file) for the target assembly.",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy.",
)
@click.option(
    "--colors",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a colours file using the reference exon/intron/intergenic palette.",
)
@click.option("--name", default="gene", show_default=True, help="Feature-set name in the stanza.")
@click.option(
    "--feature-type",
    default="exon",
    show_default=True,
    help="GFF/GTF column-3 type to read as an exon.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def gff_gene(
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None,
    name: str,
    feature_type: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A GFF3/GTF gene annotation -> an `exon`/`intron`/`intergenic` set.

    \b
    Introns are derived per transcript as the gaps between its exons. Where
    transcripts disagree the more specific label wins (exon > intron >
    intergenic), so alternative splicing never double-labels a base.

    \b
    The result tiles every sequence in --lengths, so the stanza sets
    `background: null`. If the annotation uses accessions rather than
    chromosome names, map them with --seqid-map.

    \b
    Example:
        karyoscope prep-bed gff-gene --input genes.gtf.gz --lengths ref.fa.fai \\
            --output gene.bed --hierarchy gene.tsv
    """
    _run(
        genes_prep.from_gff,
        [output, hierarchy, colors],
        force,
        input_path=input_path,
        lengths=lengths,
        output=output,
        hierarchy=hierarchy,
        colors=colors,
        name=name,
        feature_type=feature_type,
        rename=_rewriter(rename_prefix, seqid_map),
    )


@cmd.command("pave")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="PaVE genome records: a JSON array of /api/genome/{id} responses (plain or gzipped).",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file. The reading frames overlap, "
    "so a priority-mode build needs this.",
)
@click.option(
    "--colors",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a colours file: early genes cool, late genes warm, URR grey, "
    "with the E5 variants sharing one legend row.",
)
@click.option(
    "--fasta",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the genome sequences the records carry, for build's `sequence:`. "
    "bgzip and index the result.",
)
@click.option(
    "--taxonomy",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the ICTV lineage as a hierarchy, for a `type` feature set whose "
    "BED comes from `prep-bed fai`.",
)
@click.option("--name", default="gene", show_default=True, help="Feature-set name in the stanza.")
@click.option(
    "--background",
    default="intergenic",
    show_default=True,
    help="Gap-fill label build should use for the bases no ORF covers.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def pave(
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None,
    colors: Path | None,
    fasta: Path | None,
    taxonomy: Path | None,
    name: str,
    background: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """PaVE genome records -> a papillomavirus `gene` feature set.

    \b
    Leaves are the ORFs — E6, E7, E1, E2, E5, E10, L2, L1 — plus the URR,
    grouped early/late as in the standard genome map. The three spliced
    transcripts (E1^E4, E8^E2, E6*) are dropped because they lie wholly inside
    the ORFs they are spliced from, and the E1BS/E2BS binding motifs because at
    12-20 bp they are shorter than any usable k.

    \b
    A record carries the sequence and the ICTV lineage as well as the
    coordinates, so --fasta and --taxonomy write the other two files a
    papillomavirus build needs. Nothing else has to be downloaded.

    \b
    Example:
        karyoscope prep-bed pave --input pave_human_ref.json \\
            --output hpv_gene.bed --hierarchy hpv_gene.tsv \\
            --priority hpv_gene.priority.txt --fasta HPV.fasta \\
            --taxonomy hpv_type.hierarchy.txt
    """
    _run(
        pave_prep.from_pave,
        [output, hierarchy, priority, colors, fasta, taxonomy],
        force,
        input_path=input_path,
        output=output,
        hierarchy=hierarchy,
        priority=priority,
        colors=colors,
        fasta=fasta,
        taxonomy=taxonomy,
        name=name,
        background=background,
        rename=_rewriter(rename_prefix, seqid_map),
    )


# -- cytoband set -----------------------------------------------------


@cmd.command("cytoband")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="UCSC cytoBand.txt(.gz) (5 columns) or cytoBandMapped BED (6 columns).",
)
@click.option(
    "--lengths",
    type=click.Path(**_IN),
    required=True,
    help="samtools .fai (or 2-column sizes file) for the target assembly.",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy (chromosome -> band group -> band).",
)
@click.option(
    "--colors",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a colours file keyed on Giemsa stain, grouping the legend by "
    "stain so several hundred bands collapse to a handful of rows.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file, for a priority-mode build.",
)
@click.option(
    "--name", default="cytoband", show_default=True, help="Feature-set name in the stanza."
)
@click.option(
    "--primary-pattern",
    default=cytoband_prep.DEFAULT_PRIMARY_PATTERN,
    show_default=True,
    help="Regex of seqids that carry banding. The default drops alt/random/fix "
    "scaffolds, unplaced contigs and the mitochondrion.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def cytoband(
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None,
    priority: Path | None,
    name: str,
    primary_pattern: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A UCSC cytoband table -> a `cytoband` feature set.

    \b
    Band labels keep their chromosome (chr1 + p36.33 -> 1p36.33) so no label is
    ambiguous, and nest three deep: chromosome -> band group -> band.

    \b
    Sequences with no banding are not given a placeholder label. They are
    reported for the spec's top-level `exclude:`, so no feature set claims them.

    \b
    Example:
        karyoscope prep-bed cytoband --input cytoBand.txt.gz --lengths ref.fa.fai \\
            --output cytoband.bed --hierarchy cytoband.tsv --colors cytoband_colors.tsv
    """
    _run(
        cytoband_prep.from_ucsc,
        [output, hierarchy, colors, priority],
        force,
        input_path=input_path,
        lengths=lengths,
        output=output,
        hierarchy=hierarchy,
        colors=colors,
        priority=priority,
        name=name,
        primary_pattern=primary_pattern,
        rename=_rewriter(rename_prefix, seqid_map),
    )


# -- structural sets --------------------------------------------------


@cmd.command("fai")
@click.option(
    "--lengths",
    type=click.Path(**_IN),
    required=True,
    help="samtools .fai (or 2-column sizes file) for the target assembly.",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a flat hierarchy with every sequence under the root, as a "
    "starting point for hand-curated grouping.",
)
@click.option(
    "--name", default="chromosome", show_default=True, help="Feature-set name in the stanza."
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def fai(
    lengths: Path,
    output: Path,
    hierarchy: Path | None,
    name: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A samtools .fai -> a `chromosome` feature set.

    \b
    One whole-length record per sequence, labelled with its own name. No
    grouping hierarchy is derived: how sequences group (autosome vs sex,
    metacentric vs acrocentric, haplotype) is organism-specific curation that a
    .fai cannot supply. Use --hierarchy to get a flat starting point to edit.

    \b
    Keep non-karyotype sequences out with build's `exclude:` rather than by
    filtering here, so every feature set agrees about what exists.
    """
    _run(
        structural_prep.from_fai,
        [output, hierarchy],
        force,
        lengths=lengths,
        output=output,
        hierarchy=hierarchy,
        name=name,
        rename=_rewriter(rename_prefix, seqid_map),
    )


@cmd.command("censat")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="CenSat annotation BED (e.g. chm13v2.0.cenSatv2.1.bed; plain or gzipped).",
)
@click.option(
    "--lengths",
    type=click.Path(**_IN),
    required=True,
    help="samtools .fai (or 2-column sizes file) for the target assembly.",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy (CenSat v2.1 tree).",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file, ranking centromeric over rDNA over arm.",
)
@click.option("--name", default="region", show_default=True, help="Feature-set name in the stanza.")
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def censat(
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None,
    name: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A CenSat annotation -> a fully-tiled `region` feature set.

    \b
    CenSat qualifies every label with the specific arrays it contains
    (`gSat(TAR1)`); the part before the parenthesis is the leaf, so hundreds of
    distinct values collapse to the 14 features the hierarchy names.

    \b
    Every remaining base is labelled by which arm it falls on, split at the
    centromere -- located from the `ct` (centromeric transition) features that
    bracket it, or from the extent of all centromeric features where a sequence
    has no `ct`. Pericentromeric remnants far out on an arm keep their own
    labels; only the gaps around them become arm.

    \b
    Example:
        karyoscope prep-bed censat --input chm13v2.0.cenSatv2.1.bed \\
            --lengths CHM13.fa.gz.fai --output region.bed --hierarchy region.tsv
    """
    _run(
        structural_prep.from_censat,
        [output, hierarchy, priority],
        force,
        input_path=input_path,
        lengths=lengths,
        output=output,
        hierarchy=hierarchy,
        priority=priority,
        name=name,
        rename=_rewriter(rename_prefix, seqid_map),
    )


@cmd.command("satellite")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="Centromeric satellite monomers as GFF/GFF3 (1-based) or BED (0-based); "
    "the coordinate convention follows the file suffix.",
)
@click.option(
    "--lengths",
    type=click.Path(**_IN),
    required=True,
    help="samtools .fai (or 2-column sizes file) for the target assembly.",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file, for a priority-mode build.",
)
@click.option("--name", default="region", show_default=True, help="Feature-set name in the stanza.")
@click.option(
    "--satellite",
    default="satellite",
    show_default=True,
    help="Leaf label for the satellite bands, e.g. 'CEN180' or 'aSat'.",
)
@click.option(
    "--merge-gap",
    type=int,
    default=10,
    show_default=True,
    help="Merge monomers separated by at most this many bases into one band. The "
    "default bridges 1-2 bp monomer-boundary artefacts without swallowing real "
    "interior insertions; 0 keeps strict monomer bands.",
)
@click.option(
    "--cluster-gap",
    type=int,
    default=500_000,
    show_default=True,
    help="Gap for clustering bands when locating the centromere core.",
)
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def satellite(
    input_path: Path,
    lengths: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None,
    name: str,
    satellite: str,
    merge_gap: int,
    cluster_gap: int,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A centromeric-satellite annotation -> a fully-tiled `region` set.

    \b
    Monomers merge into bands, the densest cluster of bands is taken as the
    centromere core, and every remaining base is labelled by which side of the
    core it falls on: p_arm, the satellite itself, cen_gap inside the core, and
    q_arm. This is the one converter that tiles, because only the satellite
    annotation knows where the centromere is — build's gap-fill has just one
    label and could not tell the arms apart.

    \b
    p and q are assigned by coordinate, so this assumes the assembly is oriented
    short-arm-first.

    \b
    Example:
        karyoscope prep-bed satellite --input CEN180.gff --lengths ref.fa.fai \\
            --satellite CEN180 --output region.bed --hierarchy region.tsv
    """
    _run(
        structural_prep.from_satellite,
        [output, hierarchy, priority],
        force,
        input_path=input_path,
        lengths=lengths,
        output=output,
        hierarchy=hierarchy,
        priority=priority,
        name=name,
        satellite=satellite,
        merge_gap=merge_gap,
        cluster_gap=cluster_gap,
        rename=_rewriter(rename_prefix, seqid_map),
    )


@cmd.command("asat")
@click.option(
    "--input",
    "input_path",
    type=click.Path(**_IN),
    required=True,
    help="CenSat annotation BED (e.g. chm13v2.0.cenSatv2.1.bed; plain or gzipped).",
)
@click.option("--output", type=click.Path(**_OUT), required=True, help="Output feature-set BED.")
@click.option(
    "--hierarchy",
    type=click.Path(**_OUT),
    required=True,
    help="Output 'child<TAB>parent' hierarchy: the asat/alpha_hor/dhor/mon scaffold, "
    "with every array flat under its class.",
)
@click.option(
    "--priority",
    type=click.Path(**_OUT),
    default=None,
    help="Also write the tree as a 3-column priority file ranking alpha_hor over dhor over "
    "mon, so a k-mer an array shares with divergent or monomeric alpha-satellite resolves "
    "to the array. Names the background explicitly at the same priority as asat.",
)
@click.option(
    "--colors",
    type=click.Path(**_OUT),
    default=None,
    help="Also write a colours file, one colour per suprachromosomal family with a "
    "legend_group so the legend collapses to nine rows instead of one per array.",
)
@click.option(
    "--background",
    default="background",
    show_default=True,
    help="Gap-fill label named in the stanza and in the priority file.",
)
@click.option(
    "--class",
    "classes",
    type=click.Choice(asat_prep.CLASSES),
    multiple=True,
    default=asat_prep.CLASSES,
    show_default=True,
    help="alpha-satellite class to include. Repeatable; the default takes all four.",
)
@click.option("--name", default="asat", show_default=True, help="Feature-set name in the stanza.")
@_seqid_options
@click.option("--force", is_flag=True, help="Overwrite existing output files.")
def asat(
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None,
    colors: Path | None,
    background: str,
    classes: tuple[str, ...],
    name: str,
    rename_prefix: str | None,
    seqid_map: Path | None,
    force: bool,
) -> None:
    """A CenSat annotation -> a per-array `asat` feature set.

    \b
    Where `prep-bed censat` collapses CenSat to its broad classes (one `aSat`
    label for every array), this keeps the per-array detail in the
    parenthetical: `hor_1_5(S1C1/5/19H1L)` becomes a record labelled
    `S1C1_5_19H1L`.

    \b
    CenSat names more than one array on intervals where two arrays' sequence is
    interleaved. Each is emitted as its own record rather than picking a
    winner, letting HKS resolve the shared k-mers to their common ancestor --
    so do NOT pass `build --flatten` for this set. Bare continuation suffixes
    are expanded first: `hor_1_1(S3C1H2-A,B,C)` is three records for S3C1H2-A,
    S3C1H2-B and S3C1H2-C, not leaves called `A`, `B` and `C`.

    \b
    Every alpha-satellite class is included by default. Dropping one leaves that
    sequence to build's gap-fill, whose leaf sits at the hierarchy root, so every
    k-mer a named array shares with it resolves to the root. On CHM13, excluding
    `mon` and `dhor` puts 36.9% of array bases on the root, against 4.6% with
    them included.

    \b
    Arrays are left flat under `alpha_hor`. Structure among them is a phylogeny
    no annotation file contains; derive it separately (e.g. mashtree over the
    per-array sequences) and replace that star.

    \b
    Example:
        karyoscope prep-bed asat --input chm13v2.0.cenSatv2.1.bed \\
            --output asat.bed --hierarchy asat.tsv
    """
    _run(
        asat_prep.from_censat,
        [output, hierarchy, priority, colors],
        force,
        input_path=input_path,
        output=output,
        hierarchy=hierarchy,
        priority=priority,
        colors=colors,
        background=background,
        classes=tuple(classes),
        name=name,
        rename=_rewriter(rename_prefix, seqid_map),
    )
