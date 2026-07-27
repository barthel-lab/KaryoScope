"""Milestone reporting to stdout for long-running commands.

KaryoScope keeps two output channels apart on purpose (see
:mod:`karyoscope.cli`): *program output* on stdout via ``click.echo``, and
*diagnostics* on stderr via ``logging``, hidden unless the user passes
``-v``. That split is sound, but it left the two slowest commands silent:
``download`` and ``build`` announced themselves and their results, while
``annotate`` (7-22 minutes) and ``karyotype`` (tens of minutes) printed
nothing at all until their closing ``Wrote:`` block. A user watching a
blank terminal for twenty minutes cannot tell a working run from a hung
one.

This module adds a third, deliberately thin thing: a handful of
*milestone* lines on stdout, in the same shape ``download`` and ``build``
already use. It is not a second logging system. Milestones are the events
a person waiting on the command cares about — "this feature set is done,
five to go" — and nothing else. The detailed per-step timings stay at INFO
where they were, so ``-v`` remains the way to get the full picture and
nothing is printed twice.

Core modules take a :class:`Progress` explicitly rather than reaching for
a global, so a library caller gets silence by default (:data:`SILENT`) and
tests can capture output without touching global state.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TextIO

import click


def format_duration(seconds: float) -> str:
    """Render a duration the way someone reading a progress line wants it.

    Sub-minute durations keep a decimal (``45.3s``) because at that scale
    the difference matters; longer ones don't (``4m05s``, ``1h02m``).
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        minutes, secs = divmod(int(seconds), 60)
        return f"{minutes}m{secs:02d}s"
    hours, remainder = divmod(int(seconds), 3600)
    return f"{hours}h{remainder // 60:02d}m"


#: Spaces per nesting level. One level is enough to read as "this
#: happened inside the thing above it" without marching off the page.
_INDENT = "  "


class Tracker:
    """Counter for one sequence of comparable items, from :meth:`Progress.track`.

    The counter lives here rather than on :class:`Progress` because the
    commands nest: ``karyotype`` tracks its renders and, part-way through,
    calls ``annotate``, which tracks its feature sets. Sharing one counter
    made the inner sequence overwrite the outer one's total and produced
    impossible lines like ``[4/3]``. A tracker per sequence makes nesting
    correct by construction.
    """

    def __init__(self, owner: Progress, labels: Sequence[str]) -> None:
        self._owner = owner
        self._total = len(labels)
        self._index = 0
        self._width = max((len(label) for label in labels), default=0)

    def step(self, label: str, seconds: float) -> None:
        """Report one item finishing, as ``[2/6] region  3m45s``."""
        self._index += 1
        counter = f"[{self._index}/{self._total}]"
        self._owner._emit(f"{counter} {label.ljust(self._width)}  {format_duration(seconds)}")


class Progress:
    """Emits milestone lines to stdout, or nowhere when disabled.

    Parameters
    ----------
    enabled
        When False every method is a no-op. This is what ``--quiet``
        and the library default use.
    stream
        Where to write. Defaults to stdout via ``click.echo``. Tests
        pass a :class:`io.StringIO`.
    depth
        Nesting level, controlling indentation. Use :meth:`child` rather
        than setting this directly.
    """

    def __init__(
        self, *, enabled: bool = True, stream: TextIO | None = None, depth: int = 0
    ) -> None:
        self.enabled = enabled
        self._stream = stream
        self._depth = depth

    def child(self) -> Progress:
        """A reporter one level deeper, for work invoked by this command.

        ``karyotype`` cascades into ``annotate``; indenting the inner
        command's lines keeps its headline from reading as a second
        top-level command that started on its own.
        """
        return Progress(enabled=self.enabled, stream=self._stream, depth=self._depth + 1)

    def _write(self, line: str) -> None:
        if not self.enabled:
            return
        if self._stream is None:
            click.echo(line)
        else:
            self._stream.write(line + "\n")

    def _emit(self, body: str) -> None:
        """Write an indented detail line, one level in from this reporter's headline."""
        self._write(f"{_INDENT * (self._depth + 1)}{body}")

    def start(self, headline: str, *details: str) -> None:
        """Announce what the command is about to do.

        ``details`` are indented continuation lines — the run's shape
        (feature sets, threads, estimated output), which is exactly what
        a user checks before walking away from a long job.
        """
        self._write(f"{_INDENT * self._depth}{headline}")
        for detail in details:
            if detail:
                self._emit(detail)

    def track(self, labels: Sequence[str]) -> Tracker:
        """Open a counter over ``labels`` for ``[i/N]`` reporting.

        Column alignment needs the longest label up front; deriving it per
        line instead would make the output jitter as it goes.
        """
        return Tracker(self, labels)

    def stage(self, label: str, seconds: float) -> None:
        """Report a named phase finishing, with no ``[i/N]`` counter.

        For work that isn't a clean sequence of equivalent items — the KMC
        backend smooths every feature set in a single streaming pass, so
        there is no per-feature-set moment to report.
        """
        self._emit(f"{label}  {format_duration(seconds)}")

    def note(self, text: str) -> None:
        """Report something that isn't timed — a decision or a warning."""
        self._emit(text)


#: The no-op reporter. Used as the default everywhere in ``core`` so that
#: importing KaryoScope as a library never writes to stdout.
SILENT = Progress(enabled=False)


def from_context() -> Progress:
    """Build a reporter honouring the group-level ``-q/--quiet`` flag.

    Reads the flag off the click context rather than adding a parameter to
    every command signature, mirroring how ``-v`` reaches the whole
    program through one call to ``logging`` config in
    :func:`karyoscope.cli.main`.

    Returns an enabled reporter outside a click context (a library caller
    that asked for one explicitly wants output).
    """
    ctx = click.get_current_context(silent=True)
    quiet = bool(ctx is not None and isinstance(ctx.obj, dict) and ctx.obj.get("quiet"))
    return Progress(enabled=not quiet)
