"""Unit tests for :mod:`karyoscope.core.io.kmc` wrapper helpers.

These tests don't need the ``get_featureIDs`` binary -- they exercise
the pure-Python helper logic (OOM-hint augmentation, prefix
inference). The subprocess-driving tests live in ``test_kmc.py``
and are marked integration.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from karyoscope.core.external import ExternalToolError
from karyoscope.core.io.kmc import _augment_with_oom_hint, _infer_prefix

# --- _augment_with_oom_hint ----------------------------------------


class TestAugmentWithOomHint:
    """The OOM hint is appended only for SIGKILL-like exit codes
    (-9, 137). Everything else passes through with no hint, because
    pointing the user at memory when their input is malformed (or
    similar) would be misleading.
    """

    cmd: ClassVar[list[str]] = ["get_featureIDs", "--db", "x"]

    def test_sigkill_python_form_gets_hint(self) -> None:
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=-9, stderr="upstream error\n")
        assert isinstance(e, ExternalToolError)
        assert "OOM-killer" in e.stderr or "OOM" in e.stderr.upper()
        assert "upstream error" in e.stderr  # original preserved
        assert e.returncode == -9

    def test_sigkill_shell_form_gets_hint(self) -> None:
        """``137`` is the shell-style ``128 + 9`` SLURM reports."""
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=137)
        assert "OOM" in e.stderr.upper()
        assert "137" in e.stderr  # the actual code is interpolated

    def test_non_oom_exit_code_does_not_get_hint(self) -> None:
        """Plain failures (malformed input, missing file, etc.) shouldn't
        be misattributed to memory issues."""
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=1, stderr="parse error at line 3")
        assert "OOM" not in e.stderr.upper()
        assert "KaryoScope hint" not in e.stderr
        assert e.stderr == "parse error at line 3"

    def test_sigterm_does_not_get_oom_hint(self) -> None:
        """SIGTERM (-15) is normally a user / scheduler ``timeout`` kill,
        not OOM. Pointing at memory would be misleading."""
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=-15)
        assert "KaryoScope hint" not in e.stderr

    def test_returncode_zero_does_not_get_hint(self) -> None:
        """Defensive: even success-code construction (shouldn't happen
        normally) doesn't get a hint."""
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=0)
        assert e.stderr == ""

    def test_hint_mentions_actionable_fixes(self) -> None:
        """The hint should mention the concrete next steps users care
        about: --mem on SLURM, --threads, moving off a login node."""
        e = _augment_with_oom_hint(cmd=self.cmd, returncode=-9)
        lower = e.stderr.lower()
        assert "--mem" in lower
        assert "--threads" in lower or "threads" in lower
        assert "login node" in lower or "compute node" in lower

    def test_stdout_passthrough(self) -> None:
        e = _augment_with_oom_hint(
            cmd=self.cmd, returncode=-9, stdout="some output", stderr="some err"
        )
        assert e.stdout == "some output"
        assert "some err" in e.stderr


# --- _infer_prefix -------------------------------------------------


class TestInferPrefix:
    """Mirrors the C++ binary's ``get_fasta_prefix`` derivation logic."""

    def test_fasta_gz(self) -> None:
        assert (
            _infer_prefix(Path("my_assembly.fa.gz"), Path("/db/features")) == "my_assembly.features"
        )

    def test_fastq_gz(self) -> None:
        assert _infer_prefix(Path("reads.fastq.gz"), Path("/db/features")) == "reads.features"

    def test_bam(self) -> None:
        assert _infer_prefix(Path("aln.bam"), Path("/db/features")) == "aln.features"

    def test_db_path_with_kmc_pre_suffix_stripped(self) -> None:
        """Passing the actual KMC index file (not the basename) should
        still produce the basename-style prefix."""
        assert _infer_prefix(Path("x.fa"), Path("/db/features.kmc_pre")) == "x.features"
