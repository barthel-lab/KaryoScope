"""Tests for :mod:`karyoscope._entry`, the console-script entry point.

Its whole job is the failure path: a dependency missing at *import* time
cannot be caught by anything inside ``karyoscope.cli``, because the
generated console script imports that module directly. Routing through
``_entry`` is what makes the failure reportable at all.

The check used to be inferred from a side-effect -- ``karyoscope.cli``
imported every command module, so a missing dependency blew up during that
import. Subcommands are now lazy, so the dependencies are checked directly
instead, and these tests exercise that: both the up-front check and the
backstop for a dependency that is present but does not import.
"""

from __future__ import annotations

import importlib.util
import sys
import types

import pytest

from karyoscope import _entry


@pytest.fixture
def hide_module(monkeypatch: pytest.MonkeyPatch):
    """Make one module invisible to ``find_spec``, as an uninstalled one is."""

    def _hide(name: str) -> None:
        real = importlib.util.find_spec

        def fake(mod: str, *args: object, **kwargs: object):
            if mod == name:
                return None
            return real(mod, *args, **kwargs)

        monkeypatch.setattr(importlib.util, "find_spec", fake)

    return _hide


class TestMissingDependency:
    def test_reports_cleanly_and_exits_nonzero(
        self, hide_module, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """A broken install must not surface as a traceback.

        Before this, every command -- including --help and version, the two a
        user would reach for to diagnose it -- raised a bare
        ModuleNotFoundError out of the console script.
        """
        hide_module("drawsvg")
        assert _entry.main() == 1
        err = capsys.readouterr().err
        assert "installation is incomplete" in err
        assert "drawsvg" in err
        assert "Traceback" not in err

    def test_message_names_the_interpreter(
        self, hide_module, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """ "Wrong Python" is the most common cause and the least obvious one."""
        hide_module("drawsvg")
        _entry.main()
        assert sys.executable in capsys.readouterr().err

    def test_message_suggests_a_fix(self, hide_module, capsys: pytest.CaptureFixture[str]) -> None:
        hide_module("drawsvg")
        _entry.main()
        assert "pip install" in capsys.readouterr().err

    def test_does_not_reach_the_cli(self, hide_module, monkeypatch) -> None:
        """The point of checking up front is to fail before dispatch."""
        called: list[bool] = []
        stub = types.ModuleType("karyoscope.cli")
        stub.main = lambda: called.append(True) or 0  # type: ignore[attr-defined]
        monkeypatch.setitem(sys.modules, "karyoscope.cli", stub)
        hide_module("cairosvg")
        assert _entry.main() == 1
        assert called == []


class TestRequiredModules:
    """The dependency list is read from metadata so it cannot drift."""

    def test_covers_every_declared_dependency(self) -> None:
        mods = set(_entry._required_modules())
        # Every hard dependency in pyproject.toml, by module name.
        assert {"click", "drawsvg", "cairosvg", "requests", "yaml", "tqdm"} <= mods

    def test_maps_distribution_names_to_module_names(self) -> None:
        # Installed as "pyyaml", imported as "yaml" -- checking the wrong one
        # would report a healthy install as broken.
        mods = _entry._required_modules()
        assert "yaml" in mods
        assert "pyyaml" not in mods

    def test_falls_back_when_metadata_is_unreadable(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A source checkout that was never pip-installed still gets checked."""
        import importlib.metadata

        def boom(_name: str) -> None:
            raise importlib.metadata.PackageNotFoundError("karyoscope")

        monkeypatch.setattr(importlib.metadata, "requires", boom)
        assert set(_entry._required_modules()) == set(_entry._FALLBACK_REQUIREMENTS)

    def test_skips_optional_extras(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A missing dev extra is not a broken install."""
        import importlib.metadata

        monkeypatch.setattr(
            importlib.metadata,
            "requires",
            lambda _n: ["click>=8.1", 'pytest>=8; extra == "dev"'],
        )
        assert _entry._required_modules() == ["click"]


class TestImportBackstop:
    """find_spec proves a module is present, not that it imports."""

    def test_installed_but_unimportable_still_reports(
        self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        # The realistic case: a compiled dependency built against the wrong
        # ABI. find_spec finds it; importing raises.
        stub = types.ModuleType("karyoscope.cli")

        def explode() -> int:
            raise ImportError("libcairo.so.2: cannot open shared object file", name="cairosvg")

        stub.main = explode  # type: ignore[attr-defined]
        monkeypatch.setitem(sys.modules, "karyoscope.cli", stub)

        assert _entry.main() == 1
        err = capsys.readouterr().err
        assert "installation is incomplete" in err
        assert "cairosvg" in err
        assert "Traceback" not in err

    def test_report_degrades_gracefully_without_a_name(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """A bare ImportError carries no .name; still say something useful."""
        _entry._report_broken_install("a required package", "something went wrong")
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
    called: list[bool] = []
    stub = types.ModuleType("karyoscope.cli")
    stub.main = lambda: called.append(True) or 0  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "karyoscope.cli", stub)

    assert _entry.main() == 0
    assert called == [True]
