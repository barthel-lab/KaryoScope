"""Wrappers for running external command-line tools.

KaryoScope relies on a handful of external binaries — the KMC k-mer
counter most prominently, and others (bgzip, tabix, ...) likely soon.
This module centralizes the patterns for finding and invoking them so
that subprocess error handling, logging, and "tool not found" messages
are consistent everywhere.

Two entry points:

* :func:`require_tool` — locate a tool on ``$PATH`` and return its path.
  Raises :class:`ToolNotFoundError` with an actionable message if it's
  missing.
* :func:`run_tool` — invoke a subprocess. Logs the command at DEBUG.
  Raises :class:`ExternalToolError` on non-zero exit with the command,
  exit code, and tail of stderr in the message.

Neither function is "pluggable" or designed for extension. They're a
hundred lines of utility code that should stay that way.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)


class ToolNotFoundError(KaryoscopeError):
    """A required external tool was not found on ``$PATH``."""


class ExternalToolError(KaryoscopeError):
    """An external command exited with a non-zero status.

    The exception carries the original command, exit code, and captured
    stderr/stdout for inspection by callers that want to handle specific
    failure modes.
    """

    def __init__(
        self,
        cmd: Sequence[str],
        returncode: int,
        stderr: str = "",
        stdout: str = "",
    ) -> None:
        self.cmd = list(cmd)
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = stdout
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        cmd_str = " ".join(str(c) for c in self.cmd)
        msg = f"command failed with exit code {self.returncode}: {cmd_str}"
        if self.stderr:
            # Include the tail of stderr to help diagnose; full stderr is
            # available on the exception object for callers that want it.
            tail = "\n".join(self.stderr.strip().splitlines()[-10:])
            if tail:
                msg += f"\n--- stderr (last 10 lines) ---\n{tail}"
        return msg


def require_tool(
    name: str,
    *,
    install_hint: str | None = None,
) -> str:
    """Locate ``name`` on ``$PATH`` and return its absolute path.

    Parameters
    ----------
    name
        The bare command name (e.g., ``"kmc"``).
    install_hint
        Optional one-line suggestion for how the user could install the
        tool, appended to the error message on failure. For example,
        ``"Install with: conda install -c bioconda kmc"``.

    Raises
    ------
    ToolNotFoundError
        If the tool is not on the user's ``$PATH``.
    """
    path = shutil.which(name)
    if path is None:
        msg = f"required tool {name!r} was not found on $PATH"
        if install_hint:
            msg += f". {install_hint}"
        raise ToolNotFoundError(msg)
    logger.debug("located tool %r at %s", name, path)
    return path


def run_tool(
    cmd: Sequence[str | Path],
    *,
    check: bool = True,
    capture: bool = False,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run an external command and return its :class:`CompletedProcess`.

    Behavior:

    * The exact command is logged at DEBUG.
    * If ``check`` is True (the default) and the command exits non-zero,
      raises :class:`ExternalToolError` (always — even if ``capture`` is
      False, ``capture_output`` is enabled internally so the error
      message can include stderr).
    * If ``capture`` is False (the default) and the command succeeds,
      stdout and stderr go to this process's stdout/stderr — useful for
      tools that emit their own progress.
    * If ``capture`` is True, stdout and stderr are captured and
      available on the returned object as text.

    Parameters
    ----------
    cmd
        The argv-style command to run. Each element is stringified.
    check
        Raise :class:`ExternalToolError` on non-zero exit.
    capture
        Capture stdout/stderr instead of passing them through.
    cwd
        Working directory for the subprocess.
    env
        Environment overrides. If ``None``, the subprocess inherits this
        process's environment. If provided, it *replaces* the
        environment entirely — for incremental additions, the caller
        should merge with :data:`os.environ` first.
    input_text
        Optional string to pipe to the subprocess's stdin.
    """
    cmd_list = [str(c) for c in cmd]
    logger.debug("running: %s", " ".join(cmd_list))

    # When `check=True` we need stderr captured even if `capture=False`,
    # so that ExternalToolError can include it. Run with capture in that
    # case and re-emit stdout/stderr if capture wasn't requested.
    must_capture = capture or check

    result = subprocess.run(
        cmd_list,
        check=False,
        capture_output=must_capture,
        cwd=str(cwd) if cwd is not None else None,
        env=dict(env) if env is not None else None,
        input=input_text,
        text=True,
    )

    if check and result.returncode != 0:
        raise ExternalToolError(
            cmd_list,
            result.returncode,
            stderr=result.stderr or "",
            stdout=result.stdout or "",
        )

    if must_capture and not capture:
        # We captured for error-reporting only; pass through on success.
        import sys

        if result.stdout:
            sys.stdout.write(result.stdout)
        if result.stderr:
            sys.stderr.write(result.stderr)

    return result
