"""``karyoscope build`` — build an HKS index database from a genome + BEDs.

Two forms feed the same pipeline (:func:`karyoscope.core.build.build_database`):

* **Simple:** ``--sequence`` plus one or more ``--feature-set NAME=annot.bed``.
  Good for one or two feature sets.
* **Spec file:** ``--spec build.yaml`` describing the genome and every feature
  set (BED / hierarchy / priorities / colours / background). Scales to the
  multi-feature-set databases in the registry.

See ``docs`` / README for the BED contract and the spec schema.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from karyoscope import paths
from karyoscope.commands.scaffold import _parse_named_path
from karyoscope.core.build import BuildResult, build_database
from karyoscope.core.buildspec import BuildSpec
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.exceptions import (
    BuildError,
    DatabaseLayoutError,
    KaryoscopeError,
    ManifestError,
)

logger = logging.getLogger(__name__)


def _named_path_map(values: tuple[str, ...], flag: str) -> dict[str, Path]:
    """Parse repeated ``NAME=PATH`` options into ``{name: Path}``."""
    out: dict[str, Path] = {}
    for raw in values:
        name, path = _parse_named_path(raw)
        if name is None:
            raise click.UsageError(f"{flag} requires NAME=PATH form (got {raw!r})")
        if name in out:
            raise click.UsageError(f"duplicate {flag} for feature set {name!r}")
        out[name] = path
    return out


def _named_str_map(values: tuple[str, ...], flag: str) -> dict[str, str]:
    """Parse repeated ``NAME=VALUE`` options into ``{name: value}`` (value is text)."""
    out: dict[str, str] = {}
    for raw in values:
        if "=" not in raw:
            raise click.UsageError(f"{flag} requires NAME=VALUE form (got {raw!r})")
        name, _, value = raw.partition("=")
        name = name.strip()
        if not name:
            raise click.UsageError(f"empty name in {flag} value {raw!r}")
        if name in out:
            raise click.UsageError(f"duplicate {flag} for feature set {name!r}")
        out[name] = value
    return out


@click.command(
    help="Build an HKS index database from a genome and per-feature-set BED annotations.",
    no_args_is_help=True,
)
@click.option(
    "--id", "db_id", type=str, default=None, help="Database id (also the directory name)."
)
@click.option(
    "--sequence",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Genome FASTA (plain or bgzipped). Required for BED (mode A) feature sets.",
)
@click.option(
    "--feature-set",
    "feature_sets_raw",
    multiple=True,
    metavar="NAME=BED",
    help="A feature set as NAME=annotation.bed (4th column = leaf label). Repeatable.",
)
@click.option(
    "--background",
    "backgrounds_raw",
    multiple=True,
    metavar="NAME=LABEL",
    help="Gap-fill label for a feature set (default 'background'). Repeatable.",
)
@click.option(
    "--hierarchy",
    "hierarchies_raw",
    multiple=True,
    metavar="NAME=PATH",
    help="Edge-list ('child parent') hierarchy for a feature set. Repeatable.",
)
@click.option(
    "--priority",
    "priorities_raw",
    multiple=True,
    metavar="NAME=PATH",
    help="Priority file ('child priority parent' or 'name priority') for a feature set; "
    "enables priority mode. Repeatable.",
)
@click.option(
    "--colors",
    "colors_raw",
    multiple=True,
    metavar="NAME=PATH",
    help="Colours file ('feature<tab>color') for a feature set. Repeatable.",
)
@click.option(
    "--flatten-order",
    "flatten_orders_raw",
    metavar="NAME=PATH",
    multiple=True,
    help="Ranking used to flatten a feature set's overlaps, when that ranking differs from "
    "--priority. Implies --flatten for that set. Repeatable.",
)
@click.option(
    "--flatten",
    is_flag=True,
    help="Pre-flatten overlapping BED regions to one label per base (usually unnecessary "
    "under HKS). Applies to all feature sets in the simple form.",
)
@click.option(
    "--variable-k",
    "variable_k",
    is_flag=True,
    help="Build a variable-k index (queryable at any k <= s from one build, e.g. a k-sweep). "
    "The base is built from the generated per-feature FASTAs. Applies to all feature sets "
    "in the simple form. Not combinable with --priority.",
)
@click.option(
    "--spec",
    "spec_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Build-spec YAML describing the genome and every feature set (alternative to the "
    "--id/--sequence/--feature-set flags).",
)
@click.option(
    "--db-version", type=str, default="1.0.0", show_default=True, help="Database version (semver)."
)
@click.option(
    "-s",
    "--s",
    "s",
    type=int,
    default=31,
    show_default=True,
    help="Maximum query length / k-mer size.",
)
@click.option(
    "-t", "--threads", type=int, default=4, show_default=True, help="Threads for HKS construction."
)
@click.option(
    "--mem-gigas",
    type=int,
    default=8,
    show_default=True,
    help="RAM budget (GB) for base-index construction.",
)
@click.option(
    "--external-memory",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Run base-index construction in external-memory mode using this scratch directory.",
)
@click.option("--forward-only", is_flag=True, help="Do not add reverse-complemented k-mers.")
@click.option(
    "--exclude",
    "exclude_raw",
    multiple=True,
    help="Sequence name to exclude from the whole build (e.g. an organelle 'ChrM'). "
    "Repeatable; accepts comma-separated lists. Excluded sequences are dropped from every "
    "feature BED and the gap-fill index, so they read as 'none' everywhere and never appear "
    "as karyotype chromosomes. Keep non-karyotype sequences out of the chromosome set this way.",
)
@click.option(
    "--db-root",
    "db_root_arg",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override the database root (default: $KARYOSCOPE_DB or ~/.karyoscope/db/).",
)
@click.option("--no-register", is_flag=True, help="Build only; do not record in installed.json.")
@click.option(
    "--force", is_flag=True, help="Overwrite an existing database directory / install record."
)
@click.option(
    "--keep-intermediates", is_flag=True, help="Keep the per-feature FASTA working directory."
)
def cmd(
    db_id: str | None,
    sequence: Path | None,
    feature_sets_raw: tuple[str, ...],
    backgrounds_raw: tuple[str, ...],
    flatten_orders_raw: tuple[str, ...],
    hierarchies_raw: tuple[str, ...],
    priorities_raw: tuple[str, ...],
    colors_raw: tuple[str, ...],
    flatten: bool,
    variable_k: bool,
    spec_path: Path | None,
    db_version: str,
    s: int,
    threads: int,
    mem_gigas: int,
    external_memory: Path | None,
    forward_only: bool,
    exclude_raw: tuple[str, ...],
    db_root_arg: Path | None,
    no_register: bool,
    force: bool,
    keep_intermediates: bool,
) -> None:
    """Build (and by default register) an HKS database.

    \b
    Examples:
        karyoscope build --id HKS_mygenome --sequence genome.fa.gz \\
            --feature-set repeat=repeat.bed --background repeat=nonrepeat
        karyoscope build --spec build.yaml
    """
    db_root = paths.default_db_root(db_root_arg)

    try:
        spec = _build_spec(
            spec_path=spec_path,
            db_id=db_id,
            sequence=sequence,
            feature_sets_raw=feature_sets_raw,
            backgrounds_raw=backgrounds_raw,
            flatten_orders_raw=flatten_orders_raw,
            hierarchies_raw=hierarchies_raw,
            priorities_raw=priorities_raw,
            colors_raw=colors_raw,
            flatten=flatten,
            variable_k=variable_k,
            db_version=db_version,
            s=s,
            threads=threads,
            mem_gigas=mem_gigas,
            external_memory=external_memory,
            forward_only=forward_only,
            exclude_raw=exclude_raw,
        )
    except BuildError as e:
        raise click.ClickException(str(e)) from e

    click.echo(f"Building {spec.id} v{spec.version} ({len(spec.feature_sets)} feature set(s))...")
    try:
        result = build_database(
            spec,
            db_root,
            register=not no_register,
            force=force,
            keep_intermediates=keep_intermediates,
        )
    except (
        BuildError,
        ToolNotFoundError,
        ExternalToolError,
        DatabaseLayoutError,
        ManifestError,
        KaryoscopeError,
    ) as e:
        raise click.ClickException(str(e)) from e

    _report(result)


def _build_spec(
    *,
    spec_path: Path | None,
    db_id: str | None,
    sequence: Path | None,
    feature_sets_raw: tuple[str, ...],
    backgrounds_raw: tuple[str, ...],
    flatten_orders_raw: tuple[str, ...],
    hierarchies_raw: tuple[str, ...],
    priorities_raw: tuple[str, ...],
    colors_raw: tuple[str, ...],
    flatten: bool,
    variable_k: bool,
    db_version: str,
    s: int,
    threads: int,
    mem_gigas: int,
    external_memory: Path | None,
    forward_only: bool,
    exclude_raw: tuple[str, ...],
) -> BuildSpec:
    if spec_path is not None:
        conflicting = (
            feature_sets_raw
            or db_id
            or sequence
            or backgrounds_raw
            or exclude_raw
            or flatten_orders_raw
        )
        if conflicting:
            raise click.UsageError(
                "--spec cannot be combined with --id/--sequence/--feature-set/--background/"
                "--exclude; put everything in the spec file or use the flags."
            )
        return BuildSpec.from_yaml(spec_path)

    if not db_id:
        raise click.UsageError("--id is required (or use --spec)")
    if not feature_sets_raw:
        raise click.UsageError("at least one --feature-set NAME=bed is required (or use --spec)")

    return BuildSpec.from_flags(
        db_id=db_id,
        version=db_version,
        sequence=sequence,
        feature_beds=_named_path_map(feature_sets_raw, "--feature-set"),
        backgrounds=_named_str_map(backgrounds_raw, "--background"),
        flatten_orders=_named_path_map(flatten_orders_raw, "--flatten-order"),
        hierarchies=_named_path_map(hierarchies_raw, "--hierarchy"),
        priorities=_named_path_map(priorities_raw, "--priority"),
        colors=_named_path_map(colors_raw, "--colors"),
        flatten=flatten,
        variable_k=variable_k,
        s=s,
        threads=threads,
        mem_gigas=mem_gigas,
        external_memory=external_memory,
        forward_only=forward_only,
        exclude=[s.strip() for raw in exclude_raw for s in raw.split(",") if s.strip()],
    )


def _report(result: BuildResult) -> None:
    click.echo(f"Built {result.db_id} at {result.db_dir}")
    for fr in result.feature_sets:
        bg = f", background={fr.background}" if fr.background else ""
        click.echo(f"  {fr.name}: {len(fr.leaves)} leaves, {fr.mode} mode{bg}")
    if result.registered:
        click.echo(
            f"Registered {result.db_id}. Run `karyoscope info {result.db_id}` to inspect it."
        )
    else:
        click.echo(f"Not registered. Run `karyoscope register {result.db_id}` to make it usable.")
