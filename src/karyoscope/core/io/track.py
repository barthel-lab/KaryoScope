"""One reader for every annotation track ``annotate`` writes.

``annotate`` writes a track as a BED (``.bed`` or bgzipped ``.bed.gz``, the
default) or, with ``--format kst``, as a compressed ``.kst``. Every command
that reads annotate output goes through :func:`open_track_text` (a text
handle of BED lines) or :func:`open_track` (parsed rows), so it accepts all
three without knowing which it got:

* ``.bed`` / ``.bed.gz`` -- opened directly. ``smoothed``, ``feature_set``
  and ``db_dir`` are ignored: a BED is already one feature set, smoothed or
  not, and the caller picked which file.
* ``.kst`` -- streamed from ``karyozip decompress`` through a pipe, so no
  BED touches the disk. ``smoothed=True`` adds ``--smooth`` (the smoothed
  track is derived on read, never stored). A merged track (every set of a
  database in one file) also needs ``feature_set`` and ``db_dir``, for
  ``--db features.tsv --set``; smoothing needs ``db_dir`` for its
  ``hierarchy.tsv``.

A :class:`Track` bundles a path with those options, so a pipeline that
locates annotate output by name (:func:`annotation_track`) can hand the
result to code that only knows how to open a path.
"""

from __future__ import annotations

import gzip
import logging
import os
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from karyoscope.core.external import ExternalToolError
from karyoscope.core.io.karyozip import decompress_args, is_kst
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Track:
    """A track file plus how to read it.

    Path-like (``os.fspath``, ``name``, ``suffix``, ``is_file``), so code that
    used to take the BED's :class:`~pathlib.Path` can take a :class:`Track`
    unchanged as long as it opens it with :func:`open_track_text`.
    """

    path: Path
    smoothed: bool = False
    feature_set: str | None = None
    db_dir: Path | None = None

    def __fspath__(self) -> str:
        return os.fspath(self.path)

    def __str__(self) -> str:
        if not is_kst(self.path):
            return str(self.path)
        what = "smoothed" if self.smoothed else "presmoothed"
        fs = f" {self.feature_set}" if self.feature_set else ""
        return f"{self.path} ({what}{fs})"

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def suffix(self) -> str:
        return self.path.suffix

    @property
    def parent(self) -> Path:
        return self.path.parent

    def is_file(self) -> bool:
        return self.path.is_file()


TrackSource = Path | Track | str


def as_track(
    src: TrackSource,
    *,
    smoothed: bool | None = None,
    feature_set: str | None = None,
    db_dir: Path | None = None,
) -> Track:
    """Normalise ``src`` to a :class:`Track`; keyword arguments override its fields."""
    base = src if isinstance(src, Track) else Track(Path(src))
    return Track(
        path=base.path,
        smoothed=base.smoothed if smoothed is None else smoothed,
        feature_set=feature_set if feature_set is not None else base.feature_set,
        db_dir=db_dir if db_dir is not None else base.db_dir,
    )


def is_plain_bed(src: TrackSource) -> bool:
    """True for an uncompressed, seekable BED on disk (not ``-``, ``.gz`` or ``.kst``)."""
    s = os.fspath(src) if not isinstance(src, str) else src
    return s != "-" and not s.endswith(".gz") and not is_kst(s)


def _db_files(db_dir: Path | None) -> tuple[Path | None, Path | None]:
    """``(features.tsv, hierarchy.tsv)`` of a database directory, from its manifest."""
    if db_dir is None:
        return None, None
    from karyoscope.manifest import parse_manifest

    manifest = parse_manifest(Path(db_dir) / "manifest.yaml")
    features = Path(db_dir) / manifest.features if manifest.features else None
    return features, Path(db_dir) / manifest.hierarchy


def kst_command(track: Track) -> list[str]:
    """The ``karyozip decompress`` argv that prints ``track`` as BED."""
    features, hierarchy = _db_files(track.db_dir)
    return decompress_args(
        track.path,
        smoothed=track.smoothed,
        feature_set=track.feature_set,
        features=features,
        hierarchy=hierarchy,
    )


@contextmanager
def _kst_stream(track: Track) -> Iterator[IO[str]]:
    """Yield the stdout of ``karyozip decompress`` as a text handle.

    karyozip's failure is the caller's failure: a non-zero exit raises
    :class:`ExternalToolError` with its stderr once the handle is done, so a
    track that stops decoding part-way never reads as a short one. A caller
    that stops reading early is not a failure: karyozip is then still writing,
    so it is terminated and its exit status ignored.
    """
    cmd = kst_command(track)
    logger.debug("streaming %s: %s", track, " ".join(cmd))
    with tempfile.TemporaryFile() as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True)
        assert proc.stdout is not None
        finished = False
        try:
            yield proc.stdout
            # Read past what the caller consumed. Nothing left means karyozip
            # closed its output, so its exit status is the verdict on the
            # track; anything left means the caller stopped early and the rest
            # is unwanted.
            finished = proc.stdout.read(1) == ""
        finally:
            if not finished:
                proc.kill()
            proc.stdout.close()
            rc = proc.wait()
        if finished and rc != 0:
            err.seek(0)
            raise ExternalToolError(cmd, rc, stderr=err.read().decode(errors="replace"))


@contextmanager
def open_track_text(
    src: TrackSource,
    *,
    smoothed: bool | None = None,
    feature_set: str | None = None,
    db_dir: Path | None = None,
) -> Iterator[IO[str]]:
    """Open an annotation track as a text handle of 4-column BED lines.

    See the module docstring for what each option means per format.
    """
    track = as_track(src, smoothed=smoothed, feature_set=feature_set, db_dir=db_dir)
    path = track.path
    if is_kst(path):
        if not path.is_file():
            raise KaryoscopeError(f"track not found: {path}")
        with _kst_stream(track) as h:
            yield h
        return
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as h:
            yield h
    else:
        with path.open("r") as h:
            yield h


def iter_bed_rows(
    handle: IO[str], *, source: object = "track"
) -> Iterator[tuple[str, int, int, str]]:
    """Parse ``(name, start, end, label)`` rows from a BED handle.

    Skips blank lines and ``#``/``track``/``browser`` headers; raises
    :class:`KaryoscopeError` naming the line on a malformed row.
    """
    for i, raw in enumerate(handle, start=1):
        line = raw.rstrip("\n")
        if not line or line.startswith(("#", "track", "browser")):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            raise KaryoscopeError(
                f"{source}:{i}: expected 4+ tab-separated columns, got {len(parts)}: {raw!r}"
            )
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError as e:
            raise KaryoscopeError(f"{source}:{i}: non-integer coordinates: {raw!r}") from e
        yield parts[0], start, end, parts[3]


def open_track(
    src: TrackSource,
    *,
    smoothed: bool | None = None,
    feature_set: str | None = None,
    db_dir: Path | None = None,
) -> Iterator[tuple[str, int, int, str]]:
    """Yield ``(name, start, end, label)`` for every row of an annotation track.

    A generator: the file (or the karyozip child) is closed when it is
    exhausted or closed.
    """
    track = as_track(src, smoothed=smoothed, feature_set=feature_set, db_dir=db_dir)
    with open_track_text(track) as h:
        yield from iter_bed_rows(h, source=track)


def annotation_track(
    out_dir: Path,
    stem: str,
    db_id: str,
    fs: str,
    variant: str = "smoothed",
    *,
    db_dir: Path | None = None,
) -> Path | Track:
    """Locate the ``annotate`` output for one feature set, whichever format it is in.

    Looks, in order, for ``<stem>.<db>.<fs>.<variant>.bed.gz``, then ``.bed``,
    then a per-set ``<stem>.<db>.<fs>.kst``, then the merged ``<stem>.<db>.kst``.
    A BED is returned as its :class:`~pathlib.Path` (exactly as before
    ``.kst`` existed); a ``.kst`` as a :class:`Track` that reads ``fs``'s
    ``variant`` from it. With none on disk, returns the ``.bed.gz`` path, so
    the caller's "missing" message names the default.
    """
    gz = out_dir / f"{stem}.{db_id}.{fs}.{variant}.bed.gz"
    if gz.is_file():
        return gz
    plain = out_dir / f"{stem}.{db_id}.{fs}.{variant}.bed"
    if plain.is_file():
        return plain
    smoothed = variant == "smoothed"
    for kst in (out_dir / f"{stem}.{db_id}.{fs}.kst", out_dir / f"{stem}.{db_id}.kst"):
        if kst.is_file():
            return Track(kst, smoothed=smoothed, feature_set=fs, db_dir=db_dir)
    return gz
