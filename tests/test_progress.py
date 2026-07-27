"""Tests for :mod:`karyoscope.progress`."""

from __future__ import annotations

import io

from karyoscope.progress import SILENT, Progress, format_duration


def _reporter() -> tuple[Progress, io.StringIO]:
    buf = io.StringIO()
    return Progress(stream=buf), buf


# --- format_duration --------------------------------------------------


def test_format_duration_keeps_a_decimal_under_a_minute() -> None:
    assert format_duration(0.4) == "0.4s"
    assert format_duration(24.94) == "24.9s"


def test_format_duration_switches_to_minutes_and_zero_pads_seconds() -> None:
    assert format_duration(60) == "1m00s"
    assert format_duration(245.3) == "4m05s"


def test_format_duration_switches_to_hours() -> None:
    assert format_duration(3600) == "1h00m"
    assert format_duration(3720) == "1h02m"


# --- basic emission ---------------------------------------------------


def test_start_writes_a_headline_and_indented_details() -> None:
    p, buf = _reporter()
    p.start("Annotating x.fa against DB", "6 feature set(s), 16 threads")
    assert buf.getvalue() == "Annotating x.fa against DB\n  6 feature set(s), 16 threads\n"


def test_empty_details_are_skipped() -> None:
    p, buf = _reporter()
    p.start("Headline", "")
    assert buf.getvalue() == "Headline\n"


def test_step_numbers_and_pads_to_the_longest_label() -> None:
    p, buf = _reporter()
    tracker = p.track(["chromosome", "gene"])
    tracker.step("chromosome", 245.3)
    tracker.step("gene", 12.0)
    assert buf.getvalue() == ("  [1/2] chromosome  4m05s\n  [2/2] gene        12.0s\n")


def test_stage_has_no_counter() -> None:
    """The KMC backend has phases, not a countable sequence of equals."""
    p, buf = _reporter()
    p.stage("k-mer query", 50.0)
    assert buf.getvalue() == "  k-mer query  50.0s\n"


def test_note_is_untimed() -> None:
    p, buf = _reporter()
    p.note("reusing the combined BED from a previous run")
    assert buf.getvalue() == "  reusing the combined BED from a previous run\n"


# --- disabled reporter ------------------------------------------------


def test_disabled_reporter_writes_nothing() -> None:
    buf = io.StringIO()
    p = Progress(enabled=False, stream=buf)
    p.start("Headline", "detail")
    p.track(["a"]).step("a", 1.0)
    p.stage("phase", 1.0)
    p.note("note")
    assert buf.getvalue() == ""


def test_silent_is_disabled() -> None:
    """core defaults to SILENT so importing KaryoScope never prints."""
    assert SILENT.enabled is False


def test_child_of_a_disabled_reporter_stays_disabled() -> None:
    """--quiet must survive the cascade, not just the outermost command."""
    assert Progress(enabled=False).child().enabled is False


# --- nesting ----------------------------------------------------------


def test_nested_trackers_keep_separate_counters() -> None:
    """Regression: a nested sequence used to clobber the outer counter.

    ``karyotype`` tracks its renders, then calls ``annotate`` part-way
    through, which tracks its feature sets. With one shared counter the
    outer sequence resumed from the inner one's index against the inner
    one's total and printed impossible lines like ``[4/3]``.
    """
    p, buf = _reporter()
    outer = p.track(["genome/chromosome", "genome/repeat"])
    inner = p.child().track(["chromosome", "repeat", "region"])
    for label in ("chromosome", "repeat", "region"):
        inner.step(label, 1.0)
    outer.step("genome/chromosome", 1.0)
    outer.step("genome/repeat", 1.0)

    lines = buf.getvalue().splitlines()
    assert [line.strip().split()[0] for line in lines] == [
        "[1/3]",
        "[2/3]",
        "[3/3]",
        "[1/2]",
        "[2/2]",
    ]


def test_child_indents_one_level_deeper() -> None:
    p, buf = _reporter()
    p.start("Rendering karyotypes", "2 render(s)")
    child = p.child()
    child.start("Annotating", "3 feature set(s)")
    child.track(["chromosome"]).step("chromosome", 1.0)
    p.track(["genome/chromosome"]).step("genome/chromosome", 1.0)

    assert buf.getvalue() == (
        "Rendering karyotypes\n"
        "  2 render(s)\n"
        "  Annotating\n"
        "    3 feature set(s)\n"
        "    [1/1] chromosome  1.0s\n"
        "  [1/1] genome/chromosome  1.0s\n"
    )


def test_track_of_an_empty_sequence_does_not_crash() -> None:
    p, buf = _reporter()
    p.track([])
    assert buf.getvalue() == ""
