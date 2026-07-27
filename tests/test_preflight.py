"""Tests for :mod:`karyoscope.preflight`.

The point of this module is that it must agree with the backends about
what "installed" means -- a second, simpler lookup here would reject
working installs. Several tests below pin that agreement.
"""

from __future__ import annotations

import pytest

from karyoscope import preflight
from karyoscope.exceptions import MissingDependencyError

# --- resolve_binary ---------------------------------------------------


def test_resolve_binary_finds_something_on_path() -> None:
    # Not a KaryoScope dependency, but guaranteed present and not in the
    # resolver table, so it exercises the plain $PATH branch.
    assert preflight.resolve_binary("python") is not None


def test_resolve_binary_returns_none_for_a_missing_tool() -> None:
    assert preflight.resolve_binary("karyoscope_definitely_not_a_real_tool") is None


def test_resolve_binary_delegates_to_the_backend_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    """``hks`` and ``get_featureIDs`` must not be resolved with a bare which().

    Both accept an environment override, and ``get_featureIDs`` also
    resolves out of the source tree for editable installs. Re-implementing
    that here would report "not installed" for a perfectly good install.
    """
    assert "hks" in preflight._RESOLVERS
    assert "get_featureIDs" in preflight._RESOLVERS

    called: list[str] = []

    def fake() -> str:
        called.append("hks")
        return "/somewhere/hks"

    monkeypatch.setitem(preflight._RESOLVERS, "hks", fake)
    assert preflight.resolve_binary("hks") == "/somewhere/hks"
    assert called == ["hks"]


def test_resolve_binary_treats_a_raising_resolver_as_absent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def boom() -> str:
        raise RuntimeError("not found")

    monkeypatch.setitem(preflight._RESOLVERS, "hks", boom)
    assert preflight.resolve_binary("hks") is None


def test_get_featureids_resolves_from_the_source_tree() -> None:
    """Regression: an editable install has no get_featureIDs on $PATH.

    It lives at ``native/get_featureIDs/build/get_featureIDs`` instead. A
    preflight that only checked $PATH blocked every source install.
    """
    from karyoscope.core.io import kmc

    try:
        expected = kmc.get_featureids_binary()
    except Exception:
        pytest.skip("get_featureIDs is not built in this checkout")
    assert preflight.resolve_binary("get_featureIDs") == expected


# --- check / require --------------------------------------------------


def test_check_reports_only_the_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        preflight, "resolve_binary", lambda name: None if name == "seqtk" else "/bin/" + name
    )
    missing = preflight.check(["bgzip", "seqtk", "tabix"])
    assert [d.name for d in missing] == ["seqtk"]


def test_check_ignores_unknown_names() -> None:
    """An unrecognised key is a KaryoScope bug, not a user dependency problem."""
    assert preflight.check(["not_a_dependency"]) == []


def test_require_is_silent_when_everything_is_present(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(preflight, "resolve_binary", lambda name: "/bin/" + name)
    preflight.require(["bgzip", "seqtk"], context="testing")


def test_require_reports_every_missing_dependency_at_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One message, not one per rerun.

    Discovering missing tools one at a time means a user pays a full
    pipeline run for each.
    """
    monkeypatch.setattr(preflight, "resolve_binary", lambda name: None)
    with pytest.raises(MissingDependencyError) as excinfo:
        preflight.require(["bgzip", "seqtk", "samtools"], context="annotate against DB_x")
    message = str(excinfo.value)
    assert "annotate against DB_x" in message
    for name in ("bgzip", "seqtk", "samtools"):
        assert name in message
    # Each one carries an actionable install hint.
    assert message.count("conda install") >= 3


def test_require_mentions_karyoscope_version_for_diagnosis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(preflight, "resolve_binary", lambda name: None)
    with pytest.raises(MissingDependencyError) as excinfo:
        preflight.require(["bgzip"], context="testing")
    assert "karyoscope version" in str(excinfo.value)


def test_python_dependencies_are_checked_by_import_not_path() -> None:
    """cairosvg is a Python package; $PATH would never find it."""
    assert preflight.DEPENDENCIES["cairosvg"].kind == "python"
    # click is definitely importable, so the python branch resolves it.
    assert preflight._python_available(
        preflight.Dependency(name="click", kind="python", purpose="x", install_hint="y")
    )
    assert not preflight._python_available(
        preflight.Dependency(
            name="no_such_module_xyz", kind="python", purpose="x", install_hint="y"
        )
    )


def test_every_dependency_declares_a_kind_and_a_hint() -> None:
    for name, dep in preflight.DEPENDENCIES.items():
        assert dep.name == name
        assert dep.kind in ("binary", "python")
        assert dep.purpose and dep.install_hint
