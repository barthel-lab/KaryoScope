"""Shared helpers for KaryoScope CLI option handling."""

from __future__ import annotations

import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)


def parse_named_path(value: str) -> tuple[str | None, Path]:
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


def split_comma(value: str) -> list[str]:
    """Split a comma-separated string into stripped non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def resolve_db_root_flag(
    db_root: Path | None,
    legacy_db: Path | None,
    *,
    command: str,
) -> Path | None:
    """Reconcile the canonical ``--db-root`` flag with the deprecated ``--db`` alias.

    ``info`` and ``download`` historically spelled the database-root override
    ``--db``. That collides with the data commands (``annotate``, ``bin``,
    ``scaffold``, ``centromeres``, ``karyotype``), where ``--db`` selects a
    database *id* and ``--db-root`` is the root directory. The root override is
    now ``--db-root`` on every command; ``--db`` remains a hidden, deprecated
    alias on ``info``/``download`` for one release.

    Returns the effective root override (``None`` if neither flag was given).
    Raises :class:`click.UsageError` if both flags are supplied, and emits a
    deprecation warning when the legacy ``--db`` flag is used.
    """
    if legacy_db is None:
        return db_root
    if db_root is not None:
        raise click.UsageError(
            f"`{command}`: pass --db-root, not both --db-root and --db "
            "(--db is a deprecated alias for --db-root)."
        )
    logger.warning(
        "`--db` is deprecated for `%s` and will be removed in a future release; "
        "use `--db-root` instead.",
        command,
    )
    return legacy_db


def resolve_resource_check_flag(
    no_resource_check: bool,
    no_space_check: bool,
    *,
    command: str,
) -> bool:
    """Reconcile ``--no-resource-check`` with the deprecated ``--no-space-check``.

    The escape hatch originally skipped one check, on free disk. ``annotate``
    now also refuses to start when there is not enough *memory* to hold the
    index, and one flag governs both — so "space" had become the wrong word
    for what it does. ``--no-space-check`` remains a hidden alias for one
    release, per the CLI stability promise in the README.

    Returns True if the resource checks should be skipped. Raises
    :class:`click.UsageError` if both flags are supplied, and warns when the
    legacy one is used.
    """
    if not no_space_check:
        return no_resource_check
    if no_resource_check:
        raise click.UsageError(
            f"`{command}`: pass --no-resource-check, not both it and "
            "--no-space-check (the latter is a deprecated alias)."
        )
    logger.warning(
        "`--no-space-check` is deprecated for `%s` and will be removed in a "
        "future release; use `--no-resource-check` instead. It now covers "
        "memory as well as disk.",
        command,
    )
    return True
