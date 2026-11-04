"""Wrappers for the ``karyozip`` CLI, which reads and writes ``.kst`` tracks.

A ``.kst`` file is a compressed annotation track: the same rows as a
presmoothed BED, 12-44x smaller than the bgzipped BEDs ``annotate`` writes
by default. ``annotate --format kst`` produces one; every reader of
annotate output streams it back through ``karyozip decompress`` (see
:mod:`karyoscope.core.io.track`). The smoothed track is not stored: karyozip
derives it on read with ``--smooth``, using the database's ``hierarchy.tsv``.

Two shapes of track:

* **per-set** -- one feature set's presmoothed BED, labels as names. Read
  with ``karyozip decompress -i F``.
* **merged** -- every feature set of a database in one file, labels stored
  as ``features.tsv`` featureIDs. Written with ``compress -d features.tsv``
  and read back one set at a time with ``-d features.tsv --set S``.

KaryoScope locates ``karyozip`` on ``$PATH``, like ``bgzip``.
"""

from __future__ import annotations

import logging
import re
import subprocess
from collections.abc import Sequence
from pathlib import Path

from karyoscope.core.external import require_tool, run_tool
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)

#: Effort levels ``annotate --kst-effort`` offers. ``high`` exists in karyozip
#: but makes every region query decode the whole file, so annotate does not
#: offer it.
KST_EFFORTS: tuple[str, ...] = ("low", "medium")

INSTALL_HINT = (
    "Build karyozip from its repository and put the binary on PATH, or re-run with --format bed."
)

#: The format magic of every ``.kst`` this module knows how to inspect.
_KST_MAGIC = b"KST00001"

#: Bit of the header's second word that marks labels as featureIDs, i.e. a
#: track written with ``compress -d``. Mirrors ``DBFLAG`` in karyozip's
#: ``rust/src/kst.rs`` (and the C build's equivalent).
_DBFLAG = 1 << 58

_VERSION_RE = re.compile(r"karyozip\s+(\d[\w.]*)")


class KstError(KaryoscopeError):
    """A ``.kst`` track could not be read the way it was asked for."""


def karyozip_binary() -> str:
    """Absolute path of ``karyozip`` on ``$PATH``."""
    return require_tool("karyozip", install_hint=INSTALL_HINT)


def karyozip_version(binary: str | None = None) -> str | None:
    """The first line of ``karyozip --version`` (e.g. ``karyozip 0.2.0, format KST00001``)."""
    try:
        exe = binary or karyozip_binary()
        proc = subprocess.run(
            [exe, "--version"], capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.TimeoutExpired, KaryoscopeError) as exc:
        logger.debug("could not read the karyozip version: %s", exc)
        return None
    text = (proc.stdout or proc.stderr or "").strip()
    return text.splitlines()[0] if text else None


def parse_karyozip_version(line: str | None) -> str | None:
    """The bare version number out of :func:`karyozip_version`'s line."""
    if not line:
        return None
    m = _VERSION_RE.search(line)
    return m.group(1) if m else None


def is_kst(path: Path | str) -> bool:
    """True when ``path`` names a ``.kst`` track (by extension)."""
    return str(path).endswith(".kst")


def kst_holds_feature_ids(path: Path) -> bool:
    """True when ``path`` is a merged track (written with ``compress -d``).

    Read straight from the header rather than by asking karyozip, which would
    have to start decoding a block to say so. A file with a magic this module
    does not know falls back to asking karyozip.
    """
    with Path(path).open("rb") as h:
        head = h.read(24)
    if len(head) < 24:
        raise KstError(f"{path} is too short to be a .kst track")
    if head[:8] != _KST_MAGIC:
        return _probe_feature_ids(path)
    word = int.from_bytes(head[16:24], "little")
    return bool(word & _DBFLAG)


def _probe_feature_ids(path: Path) -> bool:
    """Ask karyozip whether ``path`` holds featureIDs, by its refusal to say otherwise."""
    proc = subprocess.Popen(
        [karyozip_binary(), "decompress", "-i", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        assert proc.stdout is not None
        first = proc.stdout.read(1)
    finally:
        proc.kill()
        _, err = proc.communicate()
    return not first and b"holds featureIDs" in (err or b"")


def _threads_args(threads: int) -> list[str]:
    return ["-t", str(threads)] if threads and threads > 1 else []


def compress(
    *,
    inputs: Sequence[Path],
    output: Path,
    effort: str = "low",
    features: Path | None = None,
    hierarchy: Path | None = None,
    threads: int = 0,
) -> Path:
    """``karyozip compress`` one BED (per-set) or one BED per set (merged).

    With ``features`` the track is merged: ``inputs`` are one BED per column of
    ``features.tsv``, in column order (or KMC's single combined featureID BED).
    karyozip stages ``-o`` and renames it into place, so a killed run never
    leaves a partial ``.kst`` under the final name.
    """
    if effort not in ("low", "medium", "high"):
        raise KaryoscopeError(f"unknown karyozip effort {effort!r}")
    cmd: list[str | Path] = [karyozip_binary(), "compress"]
    if features is not None:
        cmd += ["-d", features]
        if hierarchy is not None:
            cmd += ["--hier", hierarchy]
    for p in inputs:
        cmd += ["-i", p]
    cmd += ["-o", output, "-e", effort, *_threads_args(threads)]
    run_tool(cmd, capture=True)
    return output


def index(path: Path) -> Path:
    """``karyozip index``: write ``<path>.kzi`` so single-read queries stay fast."""
    run_tool([karyozip_binary(), "index", "-i", path], capture=True)
    return Path(str(path) + ".kzi")


def decompress_args(
    path: Path,
    *,
    smoothed: bool = False,
    feature_set: str | None = None,
    features: Path | None = None,
    hierarchy: Path | None = None,
    merged: bool | None = None,
) -> list[str]:
    """The ``karyozip decompress`` argv that prints one feature set's BED.

    ``merged`` is whether ``path`` holds featureIDs; ``None`` reads it from the
    header. A merged track needs ``feature_set`` and ``features``; smoothing
    any track needs ``hierarchy``. The resulting rows are exactly the BED
    ``annotate --format bed`` would have written for that set.
    """
    if merged is None:
        merged = kst_holds_feature_ids(path)
    cmd = [karyozip_binary(), "decompress", "-i", str(path)]
    if merged:
        if feature_set is None:
            raise KstError(
                f"{path.name} holds every feature set of its database; say which one to read"
            )
        if features is None:
            raise KstError(
                f"{path.name} holds every feature set of its database, so reading it needs "
                f"that database's features.tsv (give the database)"
            )
        cmd += ["-d", str(features)]
        if hierarchy is not None:
            cmd += ["--hier", str(hierarchy)]
        cmd += ["--set", feature_set]
    if smoothed:
        if hierarchy is None:
            raise KstError(
                f"the smoothed track of {path.name} is derived on read and needs the "
                f"database's hierarchy.tsv (give the database)"
            )
        if not merged:
            cmd += ["--hier", str(hierarchy)]
            if feature_set is not None:
                cmd += ["--set", feature_set]
        cmd.append("--smooth")
    return cmd


def decompress_to_file(
    path: Path,
    output: Path,
    *,
    smoothed: bool = False,
    feature_set: str | None = None,
    features: Path | None = None,
    hierarchy: Path | None = None,
    merged: bool | None = None,
    threads: int = 0,
) -> Path:
    """Write one feature set of ``path`` to the plain BED ``output``."""
    cmd = decompress_args(
        path,
        smoothed=smoothed,
        feature_set=feature_set,
        features=features,
        hierarchy=hierarchy,
        merged=merged,
    )
    run_tool([*cmd, "-o", str(output), *_threads_args(threads)], capture=True)
    return output
