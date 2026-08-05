"""Tests for the upgraded ``karyoscope info`` command."""

from __future__ import annotations

import tarfile
from pathlib import Path

from click.testing import CliRunner

from karyoscope.cli import main


def _invoke(runner: CliRunner, *args: str) -> object:
    return runner.invoke(main, list(args), catch_exceptions=False)


# --- info (no arguments) --------------------------------------------


def test_info_empty_db_root(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    """With nothing installed, `info` should report cleanly."""
    result = _invoke(cli_runner, "info", "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert str(isolated_db_root) in result.output
    # Either "No installed databases" or the "root does not exist yet" message.
    assert "No installed" in result.output or "does not exist" in result.output


def test_info_lists_installed_databases(cli_runner: CliRunner, populated_db_root: Path) -> None:
    result = _invoke(cli_runner, "info", "--db-root", str(populated_db_root))
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "Version:" in result.output
    assert "Size:" in result.output


# --- info <database_id> ---------------------------------------------


def test_info_detailed_view_of_installed_database(
    cli_runner: CliRunner, populated_db_root: Path
) -> None:
    result = _invoke(cli_runner, "info", "KS_dummy_test_v1", "--db-root", str(populated_db_root))
    assert result.exit_code == 0, result.output
    assert "KS_dummy_test_v1" in result.output
    assert "k-mer:" in result.output
    assert "size=21" in result.output
    assert "Feature sets:" in result.output
    assert "chromosome:" in result.output
    assert "region:" in result.output


def test_info_unknown_database_id(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(cli_runner, "info", "not_a_real_database", "--db-root", str(isolated_db_root))
    assert result.exit_code != 0
    assert "not installed" in result.output


# --- info <path> ----------------------------------------------------


def test_info_existing_database_directory(
    cli_runner: CliRunner,
    unpacked_dummy_db: Path,
    isolated_db_root: Path,
) -> None:
    """Pointing at a database directory on disk works even if not installed."""
    result = _invoke(cli_runner, "info", str(unpacked_dummy_db), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "KaryoScope database directory" in result.output
    assert "KS_dummy_test_v1" in result.output


def test_info_existing_regular_file(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    a_file = tmp_path / "something.bed"
    a_file.write_text("chrom\tstart\tend\tname\n")
    result = _invoke(cli_runner, "info", str(a_file), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Type: file" in result.output


def test_info_missing_path(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    nope = tmp_path / "does-not-exist"
    result = _invoke(cli_runner, "info", str(nope), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0
    assert "Does not exist" in result.output


# --- info <archive> -------------------------------------------------


def _tar_from(src_root: Path, dest: Path) -> Path:
    """Pack every entry under ``src_root`` into a gzipped tarball."""
    with tarfile.open(dest, "w:gz") as tar:
        for item in sorted(src_root.iterdir()):
            tar.add(item, arcname=item.name)
    return dest


def test_info_validates_a_database_archive(
    cli_runner: CliRunner,
    dummy_db_tarball: Path,
    isolated_db_root: Path,
) -> None:
    """The archive gets its layout checked, not just its size reported.

    This is the check the registry's CONTRIBUTING tells contributors to
    run before opening a PR.
    """
    result = _invoke(cli_runner, "info", str(dummy_db_tarball), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Type: database archive" in result.output
    assert "Layout valid: yes" in result.output
    assert "Top-level directory: KS_dummy_test_v1" in result.output
    assert "Database id:  KS_dummy_test_v1" in result.output
    # The manifest and hierarchy summaries are the point of validating.
    assert "Index type:" in result.output
    assert "Feature sets:" in result.output


def test_info_archive_leaves_nothing_behind(
    cli_runner: CliRunner,
    dummy_db_tarball: Path,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """Inspecting an archive must not unpack it into the working tree.

    Staging happens in a TemporaryDirectory that is removed on the way
    out; nothing may be left beside the archive or under the db root.
    """
    result = _invoke(cli_runner, "info", str(dummy_db_tarball), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert not list(tmp_path.rglob("manifest.yaml"))
    assert not list(tmp_path.rglob("KS_dummy_test_v1"))
    assert not list(dummy_db_tarball.parent.glob("KS_dummy_test_v1"))


def test_info_archive_without_manifest(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    src = tmp_path / "src"
    (src / "NotADatabase").mkdir(parents=True)
    (src / "NotADatabase" / "readme.txt").write_text("nothing to see\n")
    archive = _tar_from(src, tmp_path / "not_a_db.tar.gz")

    result = _invoke(cli_runner, "info", str(archive), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Layout valid: NO" in result.output
    assert "no manifest.yaml" in result.output


def test_info_archive_with_multiple_top_level_dirs(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """download installs one top-level directory; two is unusable."""
    src = tmp_path / "src"
    for name in ("db_one", "db_two"):
        (src / name).mkdir(parents=True)
        (src / name / "manifest.yaml").write_text("id: whatever\n")
    archive = _tar_from(src, tmp_path / "two_roots.tar.gz")

    result = _invoke(cli_runner, "info", str(archive), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Layout valid: NO" in result.output
    assert "exactly one top-level directory" in result.output
    assert "db_one, db_two" in result.output


def test_info_archive_reports_missing_index_files(
    cli_runner: CliRunner,
    unpacked_dummy_db: Path,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """A truncated archive is caught rather than reported as valid."""
    for stray in unpacked_dummy_db.glob("index/*"):
        stray.unlink()
    archive = _tar_from(unpacked_dummy_db.parent, tmp_path / "broken.tar.gz")

    result = _invoke(cli_runner, "info", str(archive), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Layout valid: NO" in result.output


def test_info_archive_flags_id_directory_mismatch(
    cli_runner: CliRunner,
    unpacked_dummy_db: Path,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """The manifest id and the directory name have to agree."""
    renamed = unpacked_dummy_db.parent / "some_other_name"
    unpacked_dummy_db.rename(renamed)
    archive = _tar_from(renamed.parent, tmp_path / "mismatch.tar.gz")

    result = _invoke(cli_runner, "info", str(archive), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Layout valid: yes" in result.output
    assert "does not match the top-level directory" in result.output


def test_info_archive_rejects_path_traversal_members(
    cli_runner: CliRunner,
    unpacked_dummy_db: Path,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """A member escaping the root is skipped, not written outside staging."""
    evil = tmp_path / "evil.txt"
    evil.write_text("pwned\n")
    archive = tmp_path / "traversal.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(unpacked_dummy_db, arcname=unpacked_dummy_db.name)
        tar.add(evil, arcname="../../escaped.txt")

    result = _invoke(cli_runner, "info", str(archive), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    # The legitimate database is still reported...
    assert "Layout valid: yes" in result.output
    # ...and the traversal member did not land anywhere near the tree.
    assert not (tmp_path.parent / "escaped.txt").exists()
    assert not (tmp_path / "escaped.txt").exists()


def test_info_non_archive_file_is_unchanged(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """Only archive-suffixed files take the new path."""
    a_file = tmp_path / "reads.fa.gz"
    a_file.write_bytes(b"\x1f\x8b\x08\x00not-a-tar")
    result = _invoke(cli_runner, "info", str(a_file), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Type: file" in result.output


def test_info_unreadable_archive(
    cli_runner: CliRunner,
    tmp_path: Path,
    isolated_db_root: Path,
) -> None:
    """A file named like an archive but not one fails gracefully."""
    bogus = tmp_path / "corrupt.tar.gz"
    bogus.write_bytes(b"this is not a tar archive at all")
    result = _invoke(cli_runner, "info", str(bogus), "--db-root", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert "Layout valid: NO" in result.output
    assert "cannot read archive" in result.output


# --- --db / --db-root flag handling ---------------------------------


def test_info_db_is_deprecated_alias_for_db_root(
    cli_runner: CliRunner, isolated_db_root: Path
) -> None:
    """The legacy --db flag still resolves the database root."""
    result = _invoke(cli_runner, "info", "--db", str(isolated_db_root))
    assert result.exit_code == 0, result.output
    assert str(isolated_db_root) in result.output


def test_info_db_and_db_root_conflict(cli_runner: CliRunner, isolated_db_root: Path) -> None:
    result = _invoke(
        cli_runner, "info", "--db", str(isolated_db_root), "--db-root", str(isolated_db_root)
    )
    assert result.exit_code != 0
    assert "not both" in result.output
