"""Tests for :mod:`karyoscope.core.external`."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from karyoscope.core.external import (
    ExternalToolError,
    ToolNotFoundError,
    describe_returncode,
    require_tool,
    run_tool,
)

# --- require_tool ----------------------------------------------------


def test_require_tool_finds_existing_binary() -> None:
    # `python` is always available (it's literally running us).
    path = require_tool("python") if sys.platform != "win32" else require_tool("python.exe")
    assert path  # non-empty path returned
    assert Path(path).is_file() or Path(path).is_symlink()


def test_require_tool_missing_raises() -> None:
    with pytest.raises(ToolNotFoundError, match=r"not found"):
        require_tool("totally_definitely_does_not_exist_xyz123")


def test_require_tool_missing_includes_hint() -> None:
    with pytest.raises(ToolNotFoundError, match=r"conda install"):
        require_tool(
            "totally_definitely_does_not_exist_xyz123",
            install_hint="Install with: conda install -c bioconda foo",
        )


# --- run_tool happy paths --------------------------------------------


def test_run_tool_success_capture(tmp_path: Path) -> None:
    result = run_tool(
        [sys.executable, "-c", "print('hello')"],
        capture=True,
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "hello"


def test_run_tool_success_no_check() -> None:
    result = run_tool(
        [sys.executable, "-c", "import sys; sys.exit(0)"],
        check=False,
        capture=True,
    )
    assert result.returncode == 0


def test_run_tool_passes_cwd(tmp_path: Path) -> None:
    result = run_tool(
        [sys.executable, "-c", "import os; print(os.getcwd())"],
        capture=True,
        cwd=tmp_path,
    )
    # On macOS /tmp may be symlinked to /private/tmp; compare resolved paths.
    assert Path(result.stdout.strip()).resolve() == tmp_path.resolve()


def test_run_tool_passes_env() -> None:
    result = run_tool(
        [sys.executable, "-c", "import os; print(os.environ['KS_TEST_VAR'])"],
        capture=True,
        env={**os.environ, "KS_TEST_VAR": "hello_world"},
    )
    assert result.stdout.strip() == "hello_world"


def test_run_tool_pipes_input() -> None:
    result = run_tool(
        [sys.executable, "-c", "import sys; print(sys.stdin.read().upper())"],
        capture=True,
        input_text="quiet",
    )
    assert result.stdout.strip() == "QUIET"


# --- run_tool failure modes ------------------------------------------


def test_run_tool_failure_raises_with_check() -> None:
    with pytest.raises(ExternalToolError) as exc_info:
        run_tool(
            [sys.executable, "-c", "import sys; sys.exit(7)"],
        )
    assert exc_info.value.returncode == 7
    assert exc_info.value.cmd[0] == sys.executable


def test_run_tool_failure_includes_stderr_tail() -> None:
    with pytest.raises(ExternalToolError) as exc_info:
        run_tool(
            [
                sys.executable,
                "-c",
                "import sys; print('oops happened', file=sys.stderr); sys.exit(1)",
            ],
        )
    assert "oops happened" in str(exc_info.value)
    assert exc_info.value.stderr.strip() == "oops happened"


def test_run_tool_failure_no_check_returns_result() -> None:
    """check=False: failure doesn't raise; caller inspects the return code."""
    result = run_tool(
        [sys.executable, "-c", "import sys; sys.exit(3)"],
        check=False,
        capture=True,
    )
    assert result.returncode == 3


def test_external_tool_error_carries_command_for_inspection() -> None:
    """The exception exposes the original argv so callers can branch on it."""
    try:
        run_tool([sys.executable, "-c", "import sys; sys.exit(1)"])
    except ExternalToolError as e:
        assert e.cmd[0] == sys.executable
        assert e.returncode == 1
    else:  # pragma: no cover
        pytest.fail("expected ExternalToolError")


def test_run_tool_message_truncates_long_stderr() -> None:
    """A long stderr is shown as the last ~10 lines in the error message."""
    code = "import sys\nfor i in range(50):\n    print(f'line {i}', file=sys.stderr)\nsys.exit(1)\n"
    with pytest.raises(ExternalToolError) as exc_info:
        run_tool([sys.executable, "-c", code])
    msg = str(exc_info.value)
    # Last 10 lines: line 40 .. line 49
    assert "line 49" in msg
    # Earlier lines should not be in the message preview (but should be on .stderr)
    assert "line 0\n" not in msg
    assert "line 0" in exc_info.value.stderr


# --- signal / OOM reporting -------------------------------------------


def test_negative_returncode_is_reported_as_a_signal() -> None:
    """-9 is not an exit code, and calling it one sends users nowhere.

    This is the exact message a colleague hit: `exit code -9` is
    unsearchable, while SIGKILL points straight at the cause.
    """
    assert describe_returncode(-9) == "killed by SIGKILL (signal 9)"
    assert describe_returncode(-15) == "killed by SIGTERM (signal 15)"


def test_shell_style_signal_codes_are_translated() -> None:
    """SLURM and Docker report 128 + N for a signal death."""
    assert "SIGKILL" in describe_returncode(137)
    assert "137" in describe_returncode(137)


def test_ordinary_exit_codes_are_left_alone() -> None:
    assert describe_returncode(1) == "exit code 1"
    assert describe_returncode(2) == "exit code 2"


def test_oom_hint_is_shown_only_for_oom_like_codes() -> None:
    hint = "Request at least 16 GB.\n"
    for code in (-9, 137):
        e = ExternalToolError(["hks", "lookup"], code, oom_hint=hint)
        assert "Request at least 16 GB." in str(e), code
        assert "KaryoScope hint" in str(e), code
    for code in (1, 2, -15):
        e = ExternalToolError(["hks", "lookup"], code, oom_hint=hint)
        assert "Request at least 16 GB." not in str(e), code


def test_the_hint_never_pollutes_captured_stderr() -> None:
    """stderr must stay exactly what the tool wrote.

    An earlier version spliced the hint into stderr, so a caller inspecting
    it programmatically saw KaryoScope's words mixed into the tool's output.
    """
    e = ExternalToolError(["hks"], -9, stderr="real tool output\n", oom_hint="advice\n")
    assert e.stderr == "real tool output\n"
    assert "advice" in str(e)


def test_a_tool_with_no_hint_still_explains_the_signal() -> None:
    """Every tool benefits from the translation, even without advice."""
    e = ExternalToolError(["bgzip", "-f", "x.bed"], -9)
    assert "SIGKILL" in str(e)
    assert "out of memory" in str(e) or "too much memory" in str(e)


def test_killed_by_signal_predicate() -> None:
    assert ExternalToolError(["x"], -9).killed_by_signal
    assert ExternalToolError(["x"], 137).killed_by_signal
    assert not ExternalToolError(["x"], 1).killed_by_signal
