"""Integration tests for the ``get_featureIDs`` C++ binary and its Python wrapper.

Marked ``@pytest.mark.integration`` so a plain ``pytest`` run skips them
on systems where the binary hasn't been built. CI builds the binary
explicitly and then runs ``pytest -m integration``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.core.io.kmc import (
    BINARY_NAME,
    ENV_OVERRIDE,
    get_featureids_binary,
    run_get_featureids,
)

pytestmark = pytest.mark.integration


# --- Fixture: locate (or skip on missing) the built binary ----------


def _find_binary_for_tests() -> str | None:
    """Skip-aware version of ``get_featureids_binary``.

    Returns the binary path if found, else ``None`` (the test will skip).
    """
    # Try the production lookup first.
    try:
        return get_featureids_binary()
    except ToolNotFoundError:
        return None


@pytest.fixture(scope="module")
def featureids_binary() -> str:
    path = _find_binary_for_tests()
    if path is None:
        pytest.skip(
            f"{BINARY_NAME} binary not found — build with "
            "`make -C native/get_featureIDs` or set "
            f"{ENV_OVERRIDE} to its location."
        )
    return path


@pytest.fixture
def query_fasta(tmp_path: Path) -> Path:
    """A small FASTA we can query against the dummy db.

    The dummy db's k-mer index is k=21, built from synthetic sequences.
    These query sequences are chosen so the run-length-encoded output
    has a predictable structure (each sequence becomes one or more BED
    records).
    """
    fa = tmp_path / "query.fa"
    fa.write_text(
        ">test_seq_1\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">test_seq_2\n"
        "GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT\n"
    )
    return fa


# --- Tests ----------------------------------------------------------


def test_binary_runs_help_without_arguments(featureids_binary: str) -> None:
    """Smoke test: --help exits 0 and prints recognisable text."""
    import subprocess

    result = subprocess.run(
        [featureids_binary, "--help"], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0
    assert "get_featureIDs" in result.stdout
    assert "--db" in result.stdout
    assert "--input" in result.stdout


def test_wrapper_runs_against_dummy_db(
    featureids_binary: str,
    unpacked_dummy_db: Path,
    query_fasta: Path,
    tmp_path: Path,
) -> None:
    """End-to-end: Python wrapper → binary → BED output."""
    output_dir = tmp_path / "out"
    kmc_db = unpacked_dummy_db / "index" / "features"

    bed_path = run_get_featureids(
        db_path=kmc_db,
        input_path=query_fasta,
        output_dir=output_dir,
        threads=1,
        capture=True,
    )

    assert bed_path.is_file(), f"expected output BED at {bed_path}"
    assert bed_path.parent == output_dir
    # Expected naming: <query_basename>.<kmc_basename>.combined.presmoothed.featureIDs.bed
    assert bed_path.name == "query.features.combined.presmoothed.featureIDs.bed"


def test_bed_output_is_well_formed(
    featureids_binary: str,
    unpacked_dummy_db: Path,
    query_fasta: Path,
    tmp_path: Path,
) -> None:
    """The BED produced should be a valid 4-column tab-separated file."""
    output_dir = tmp_path / "out"
    bed_path = run_get_featureids(
        db_path=unpacked_dummy_db / "index" / "features",
        input_path=query_fasta,
        output_dir=output_dir,
        threads=1,
        capture=True,
    )

    lines = bed_path.read_text().splitlines()
    assert lines, "BED file is empty"

    seen_seq_names: set[str] = set()
    for i, line in enumerate(lines, start=1):
        fields = line.split("\t")
        assert len(fields) == 4, (
            f"line {i}: expected 4 tab-separated fields, got {len(fields)}: {line!r}"
        )
        seq_name, start, end, feature_id = fields

        # start, end, feature_id must all be integers
        s, e, fid = int(start), int(end), int(feature_id)
        assert s >= 0
        assert e > s
        assert fid >= 0
        seen_seq_names.add(seq_name)

    # Both query sequences should appear in the output.
    assert seen_seq_names == {"test_seq_1", "test_seq_2"}


def test_wrapper_with_custom_prefix(
    featureids_binary: str,
    unpacked_dummy_db: Path,
    query_fasta: Path,
    tmp_path: Path,
) -> None:
    """The --prefix flag overrides the auto-derived output basename."""
    output_dir = tmp_path / "out"
    bed_path = run_get_featureids(
        db_path=unpacked_dummy_db / "index" / "features",
        input_path=query_fasta,
        output_dir=output_dir,
        prefix="my_custom_run",
        threads=1,
        capture=True,
    )

    assert bed_path.name == "my_custom_run.combined.presmoothed.featureIDs.bed"
    assert bed_path.is_file()


def test_wrapper_creates_output_dir(
    featureids_binary: str,
    unpacked_dummy_db: Path,
    query_fasta: Path,
    tmp_path: Path,
) -> None:
    """The output_dir is created if it doesn't exist."""
    nested = tmp_path / "deeply" / "nested" / "out"
    assert not nested.exists()
    bed_path = run_get_featureids(
        db_path=unpacked_dummy_db / "index" / "features",
        input_path=query_fasta,
        output_dir=nested,
        threads=1,
        capture=True,
    )
    assert bed_path.is_file()
    assert nested.is_dir()


def test_wrapper_propagates_subprocess_failure(
    featureids_binary: str,
    unpacked_dummy_db: Path,
    tmp_path: Path,
) -> None:
    """An invalid KMC db path causes the binary to fail; we surface ExternalToolError."""
    output_dir = tmp_path / "out"
    nonexistent_input = tmp_path / "does_not_exist.fa"

    with pytest.raises(ExternalToolError) as exc_info:
        run_get_featureids(
            db_path=unpacked_dummy_db / "index" / "features",
            input_path=nonexistent_input,
            output_dir=output_dir,
            threads=1,
            capture=True,
        )
    assert exc_info.value.returncode != 0


# --- Tests for the binary-lookup logic itself ------------------------


def test_env_override_pointing_at_missing_file_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """If $KARYOSCOPE_GET_FEATUREIDS is set but bogus, we tell the user."""
    monkeypatch.setenv(ENV_OVERRIDE, str(tmp_path / "no-such-file"))
    with pytest.raises(ToolNotFoundError, match=ENV_OVERRIDE):
        get_featureids_binary()


def test_env_override_pointing_at_real_file_is_used(
    monkeypatch: pytest.MonkeyPatch,
    featureids_binary: str,
) -> None:
    """A valid override path takes precedence over $PATH lookup."""
    monkeypatch.setenv(ENV_OVERRIDE, featureids_binary)
    found = get_featureids_binary()
    assert found == featureids_binary


def test_lookup_falls_back_when_env_not_set(
    monkeypatch: pytest.MonkeyPatch, featureids_binary: str
) -> None:
    """Without the env var, the lookup still resolves via $PATH or source tree."""
    monkeypatch.delenv(ENV_OVERRIDE, raising=False)
    found = get_featureids_binary()
    # Should at least be the same file as our test fixture (resolved differently).
    assert Path(found).is_file()
