"""Tests for :mod:`karyoscope.core.external`."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from karyoscope.core.external import (
    ExternalToolError,
    ToolNotFoundError,
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
