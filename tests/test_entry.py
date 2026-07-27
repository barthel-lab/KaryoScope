"""Tests for :mod:`karyoscope._entry`, the console-script entry point.

Its whole job is the failure path: a dependency missing at *import* time
cannot be caught by anything inside ``karyoscope.cli``, because the
generated console script imports that module directly. Routing through
``_entry`` is what makes the failure reportable at all.
"""

from __future__ import annotations

import builtins
import sys

import pytest

from karyoscope import _entry


@pytest.fixture
def block_import(monkeypatch: pytest.MonkeyPatch):
    """Make a named module unimportable on the next ``karyoscope.cli`` import.

    Every ``karyoscope.*`` module has to leave ``sys.modules``, not just
    ``karyoscope.cli``. ``drawsvg`` is imported by
    ``karyoscope.core.karyotype``, several levels down the chain; if that
    module is still cached -- which it is as soon as any other test in the
    suite has touched it -- re-importing ``karyoscope.cli`` is served from
    cache and never re-executes the ``import drawsvg`` line at all. An
    earlier version of this fixture purged only ``karyoscope.cli`` and so
    passed in isolation but silently tested nothing in a full run.
    """

    def _block(name: str) -> None:
        real = builtins.__import__

        def fake(mod: str, *args: object, **kwargs: object):
            if mod == name or mod.startswith(name + "."):
                raise ModuleNotFoundError(f"No module named {name!r}", name=name)
            return real(mod, *args, **kwargs)

        doomed = [
            m for m in sys.modules if m == "karyoscope" or m.startswith(("karyoscope.", name))
        ]
        for mod in doomed:
            monkeypatch.delitem(sys.modules, mod, raising=False)
        monkeypatch.setattr(builtins, "__import__", fake)

    return _block


def test_missing_dependency_reports_cleanly_and_exits_nonzero(
    block_import, capsys: pytest.CaptureFixture[str]
) -> None:
    """A broken install must not surface as a traceback.

    Before this, every command -- including --help and version, the two a
    user would reach for to diagnose it -- raised a bare
    ModuleNotFoundError out of the console script.
    """
    block_import("drawsvg")
    assert _entry.main() == 1
    err = capsys.readouterr().err
    assert "installation is incomplete" in err
    assert "drawsvg" in err
    assert "Traceback" not in err


def test_message_names_the_interpreter(block_import, capsys: pytest.CaptureFixture[str]) -> None:
    """ "Wrong Python" is the most common cause and the least obvious one."""
    block_import("drawsvg")
    _entry.main()
    assert sys.executable in capsys.readouterr().err


def test_message_suggests_a_fix(block_import, capsys: pytest.CaptureFixture[str]) -> None:
    block_import("drawsvg")
    _entry.main()
    assert "pip install" in capsys.readouterr().err


def test_report_uses_the_exception_name_when_present(
    capsys: pytest.CaptureFixture[str],
) -> None:
    _entry._report_broken_install(ModuleNotFoundError("boom", name="cairosvg"))
    assert "'cairosvg'" in capsys.readouterr().err


def test_report_degrades_gracefully_without_a_name(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A bare ImportError carries no .name; still say something useful."""
    _entry._report_broken_install(ImportError("something went wrong"))
    err = capsys.readouterr().err
    assert "a required package" in err
    assert "something went wrong" in err


def test_healthy_install_delegates_to_the_cli(monkeypatch: pytest.MonkeyPatch) -> None:
    """The happy path must be a plain pass-through.

    Substitutes a stub module rather than patching an attribute on the real
    one: ``_entry`` resolves ``karyoscope.cli`` through ``sys.modules`` at
    call time, so stubbing there is what actually intercepts it regardless
    of what else in the suite has imported.
    """
    import types

    called: list[bool] = []
    stub = types.ModuleType("karyoscope.cli")
    stub.main = lambda: called.append(True) or 0  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "karyoscope.cli", stub)

    assert _entry.main() == 0
    assert called == [True]
