"""Shared pytest fixtures for the KaryoScope test suite."""

from __future__ import annotations

import pytest
from click.testing import CliRunner


@pytest.fixture
def cli_runner() -> CliRunner:
    """A click CliRunner for invoking subcommands in tests."""
    return CliRunner()
