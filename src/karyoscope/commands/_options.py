"""Shared helpers for KaryoScope CLI option handling."""

from __future__ import annotations

import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)


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
