"""Shared pytest fixtures for the KaryoScope test suite."""

from __future__ import annotations

import sys

# Friendly diagnostic for a common environment mistake — especially on
# macOS where multiple Python installations coexist (Python.org, Homebrew,
# conda, etc.). If `pytest` is running from a different Python than the
# one where `pip install -e .[dev]` was run, the user sees this message
# instead of a wall of `ModuleNotFoundError` tracebacks.
#
# This must run BEFORE the `from karyoscope...` imports below.
try:
    import karyoscope  # noqa: F401
except ImportError as _import_err:
    _msg = (
        "\n"
        "    KaryoScope is not importable from this Python interpreter.\n"
        "\n"
        f"    Python being used:  {sys.executable}\n"
        "\n"
        "    This usually means pytest is running from a different Python\n"
        "    environment than the one where 'pip install -e .[dev]' was run.\n"
        "    Common on macOS where multiple Python installations coexist.\n"
        "\n"
        "    To fix, either:\n"
        "      (a) Run pytest through the right Python (recommended):\n"
        "              python -m pytest\n"
        "      (b) Or install karyoscope in this Python:\n"
        f"              {sys.executable} -m pip install -e '.[dev]'\n"
    )
    raise RuntimeError(_msg) from _import_err

import gzip
import shutil
import tarfile
from collections import OrderedDict
from pathlib import Path

import pytest
from click.testing import CliRunner

# Paths to committed test fixtures.
TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
DUMMY_TARBALL = DATA_DIR / "dummy_db.tar.gz"
DUMMY_REGISTRY = DATA_DIR / "dummy_registry.yaml"
DUMMY_SHA256_FILE = DATA_DIR / "dummy_db.sha256"
DUMMY_DB_ID = "KS_dummy_test_v1"


def _extractall_compat(tar: tarfile.TarFile, dest: Path) -> None:
    """Extract a tar's members to ``dest``, using ``filter='data'`` where available.

    The ``filter`` keyword argument to ``TarFile.extractall`` was added in
    Python 3.12 (and backported to 3.10.12 / 3.11.4 as a CVE-2007-4559
    mitigation). Older patch-level releases of 3.10 and 3.11 — including
    the Python versions that some GitHub Actions runners ship — don't have
    it and raise ``TypeError``. Production code (``download.py``) has the
    same try/except; this helper exists so the test fixtures don't need
    to repeat it.

    Test fixtures extract from tarballs we built ourselves, so the
    ``filter='data'`` is for forward-compatibility (Python 3.14 will
    require a filter) rather than for security per se.
    """
    try:
        tar.extractall(dest, filter="data")  # type: ignore[arg-type]
    except TypeError:  # pragma: no cover — exercised only on older Pythons
        tar.extractall(dest)


def read_fasta_records(path: Path) -> OrderedDict[str, str]:
    """Read a FASTA into ``{name: sequence}`` -- a whole-file test helper.

    Used by scaffold/FASTA tests to read a (small) output assembly back and
    assert on sequence content. It deliberately lives here, not in
    production ``karyoscope.core.io.fasta``, because it holds the whole file
    in memory: the pipeline itself streams FASTA (``read_fasta_lengths`` /
    per-contig spill), so no production code should ever reach for a
    whole-file reader. Test fixtures are tiny, so the in-memory read is fine.

    Name = first whitespace token of the ``>`` header; multi-line bodies are
    concatenated; blank lines skipped; CR/LF stripped; insertion order kept.
    Plain and ``.gz`` inputs both supported.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    records: OrderedDict[str, str] = OrderedDict()
    current_name: str | None = None
    current_chunks: list[str] = []
    with opener(path, "rt") as h:
        for raw in h:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_chunks)
                head = line[1:].lstrip()
                current_name = head.split()[0] if head else ""
                current_chunks = []
            elif current_name is not None:
                current_chunks.append(line)
    if current_name is not None:
        records[current_name] = "".join(current_chunks)
    return records


@pytest.fixture
def cli_runner() -> CliRunner:
    """A click CliRunner for invoking subcommands in tests."""
    return CliRunner()


@pytest.fixture
def dummy_db_tarball() -> Path:
    """Path to the committed dummy database tarball."""
    assert DUMMY_TARBALL.is_file(), (
        f"missing {DUMMY_TARBALL}; regenerate it with `python tests/data/build_dummy_db.py`"
    )
    return DUMMY_TARBALL


@pytest.fixture
def dummy_db_sha256() -> str:
    """The SHA-256 of the committed dummy database tarball."""
    text = DUMMY_SHA256_FILE.read_text().strip()
    # Format: "<hex>  <filename>"
    return text.split()[0]


@pytest.fixture
def dummy_db_url(dummy_db_tarball: Path) -> str:
    """A ``file://`` URL pointing at the dummy db tarball."""
    return dummy_db_tarball.absolute().as_uri()


@pytest.fixture
def unpacked_dummy_db(tmp_path: Path, dummy_db_tarball: Path) -> Path:
    """A freshly extracted copy of the dummy database in a temp directory.

    Use this when a test needs to operate on a real on-disk KaryoScope
    database (e.g., to call ``validate_database_layout``). The extraction
    is per-test, so tests can mutate the directory without affecting one
    another.
    """
    target_root = tmp_path / "unpacked"
    target_root.mkdir()
    with tarfile.open(dummy_db_tarball, "r:gz") as tar:
        _extractall_compat(tar, target_root)
    return target_root / DUMMY_DB_ID


@pytest.fixture
def dummy_registry_url(tmp_path: Path, dummy_db_tarball: Path, dummy_db_sha256: str) -> str:
    """A ``file://`` URL pointing at a customized dummy_registry.yaml.

    The committed ``dummy_registry.yaml`` has a placeholder URL. This
    fixture rewrites the placeholder to point at the local dummy db
    tarball, then writes the result to a temp file and returns its
    ``file://`` URL. Use this as the ``--registry-url`` (or
    ``registry_url=`` argument) in tests that need a full registry.
    """
    text = DUMMY_REGISTRY.read_text()
    rewritten = text.replace(
        "PLACEHOLDER_WILL_BE_REWRITTEN_BY_TEST_FIXTURE",
        dummy_db_tarball.absolute().as_uri(),
    )
    # Sanity check: confirm the SHA-256 matches what's in the registry.
    assert dummy_db_sha256 in rewritten, (
        "the SHA-256 in dummy_registry.yaml does not match the actual tarball; "
        "regenerate fixtures with `python tests/data/build_dummy_db.py`"
    )
    out = tmp_path / "dummy_registry.yaml"
    out.write_text(rewritten)
    return out.absolute().as_uri()


@pytest.fixture
def isolated_db_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """An isolated, empty database root for one test.

    Sets ``$KARYOSCOPE_DB`` so any code path that calls
    :func:`karyoscope.paths.ensure_db_root` without an explicit override
    picks up this temporary directory. Cleans up the env var via monkeypatch.
    """
    db_root = tmp_path / "karyoscope_db"
    monkeypatch.setenv("KARYOSCOPE_DB", str(db_root))
    return db_root


@pytest.fixture
def populated_db_root(
    tmp_path: Path,
    dummy_db_tarball: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    """A database root with the dummy db already extracted and recorded."""
    from karyoscope.installed import InstalledRecord, now_iso, record_install

    db_root = tmp_path / "karyoscope_db"
    db_root.mkdir()
    monkeypatch.setenv("KARYOSCOPE_DB", str(db_root))

    # Extract the dummy db directly into db_root.
    with tarfile.open(dummy_db_tarball, "r:gz") as tar:
        _extractall_compat(tar, db_root)
    assert (db_root / DUMMY_DB_ID).is_dir()

    # Pretend it was installed.
    record_install(
        db_root,
        DUMMY_DB_ID,
        InstalledRecord(
            version="1.0.0",
            installed_at=now_iso(),
            source_url="file://test/dummy_db.tar.gz",
            source_sha256="0" * 64,
            directory=DUMMY_DB_ID,
        ),
    )
    return db_root


# Defensive cleanup: if a test mucks up by writing to the real ~/.karyoscope,
# we don't want it to bleed into other tests. The isolated_db_root and
# populated_db_root fixtures already isolate via $KARYOSCOPE_DB; this is a
# belt-and-suspenders check.
@pytest.fixture(autouse=True)
def _block_real_home_db(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Make sure tests don't accidentally touch the real ~/.karyoscope."""
    fake_home = tmp_path / "_fake_home"
    fake_home.mkdir(exist_ok=True)
    monkeypatch.setenv("HOME", str(fake_home))


# Keep shutil import warm so the module loads cleanly even if unused.
_ = shutil
