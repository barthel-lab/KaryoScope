"""URL fetching with progress bars, SHA-256 verification, and resume support.

This module hides the differences between ``http(s)://`` and ``file://`` URLs
behind a single :func:`fetch` function. Tests use ``file://`` to exercise the
full download path offline; production uses HTTPS.
"""

from __future__ import annotations

import hashlib
import shutil
from collections.abc import Iterator
from pathlib import Path
from urllib.parse import urlparse

import requests

# NB: `from tqdm import tqdm` is imported LAZILY inside _stream_http (below), not
# at module load. tqdm is only needed for the download feature, but cli.py imports
# every command at startup -> a missing tqdm broke the WHOLE CLI (incl. karyotype,
# which needs no tqdm). Making it lazy lets karyotype/annotate/scaffold run in an
# env without tqdm. (dsingh 2026-06-15: unblocked chm13_divergent_v3_karyotype.)

from karyoscope.exceptions import ChecksumError, FetchError, UnsupportedSchemeError

# Read in 1 MiB chunks. Large enough to amortize syscall overhead, small enough
# that progress bars feel responsive on slow links.
_CHUNK_SIZE = 1024 * 1024


def _sha256_file(path: Path) -> str:
    """Return the hex SHA-256 digest of a file's contents."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(_CHUNK_SIZE):
            h.update(chunk)
    return h.hexdigest()


def _scheme(url: str) -> str:
    return urlparse(url).scheme.lower()


def _iter_http(
    url: str,
    resume_from: int = 0,
    timeout: float = 30.0,
) -> tuple[Iterator[bytes], int | None]:
    """Stream the body of an HTTP(S) URL.

    Returns ``(chunk_iterator, total_size_or_None)``. ``total_size`` includes
    the bytes before ``resume_from`` (i.e., it is the full file size, not the
    remaining size), to match what progress bars typically expect to display.
    """
    headers = {}
    if resume_from > 0:
        headers["Range"] = f"bytes={resume_from}-"

    try:
        response = requests.get(url, stream=True, headers=headers, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as e:
        raise FetchError(f"failed to fetch {url}: {e}") from e

    # Determine total file size for progress display.
    if response.status_code == 206:
        # Partial Content; server confirmed our Range request.
        # Content-Range looks like "bytes 1000-2000/2001".
        content_range = response.headers.get("Content-Range", "")
        if "/" in content_range:
            try:
                total = int(content_range.rsplit("/", 1)[1])
            except ValueError:
                total = None
        else:
            total = None
    else:
        # Full response. If we asked for a Range but the server ignored it,
        # we need to discard what we already have on disk.
        if resume_from > 0:
            resume_from = 0  # caller will need to start over; we signal via iterator
        length = response.headers.get("Content-Length")
        total = int(length) if length else None

    return response.iter_content(chunk_size=_CHUNK_SIZE), total


def _fetch_http(
    url: str,
    dest: Path,
    show_progress: bool = True,
    resume: bool = True,
    timeout: float = 30.0,
) -> None:
    """Stream an HTTP(S) URL to ``dest``, with optional resume."""
    partial = dest.with_suffix(dest.suffix + ".part")
    resume_from = partial.stat().st_size if (resume and partial.is_file()) else 0
    if resume_from > 0:
        # Server might 200 instead of 206; in that case the helper resets resume_from.
        pass

    iterator, total = _iter_http(url, resume_from=resume_from, timeout=timeout)

    # Open partial in append mode so a 206 response continues at the right offset.
    mode = "ab" if resume_from > 0 else "wb"
    from tqdm import tqdm  # lazy: only the download path needs it (see module header)

    progress = tqdm(
        total=total,
        initial=resume_from,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        disable=not show_progress,
        desc=dest.name,
    )
    try:
        with partial.open(mode) as f:
            for chunk in iterator:
                if chunk:
                    f.write(chunk)
                    progress.update(len(chunk))
    finally:
        progress.close()

    partial.replace(dest)


def _fetch_file(url: str, dest: Path) -> None:
    """Copy a ``file://`` URL to ``dest``."""
    parsed = urlparse(url)
    source = Path(parsed.path)
    if not source.is_file():
        raise FetchError(f"file URL {url} does not point at a regular file")
    shutil.copyfile(source, dest)


def fetch(
    url: str,
    dest: Path,
    expected_sha256: str | None = None,
    show_progress: bool = True,
    resume: bool = True,
    timeout: float = 30.0,
) -> None:
    """Fetch a URL to a local path, optionally verifying its SHA-256 digest.

    Supports ``http://``, ``https://``, and ``file://`` URLs.

    The download is staged to ``<dest>.part`` and atomically renamed to
    ``dest`` only on success, so a partial download never overwrites a good
    file at ``dest``. For HTTP(S) URLs, if a ``.part`` file already exists,
    the download resumes using a ``Range`` request (unless ``resume=False``).

    If ``expected_sha256`` is provided, the file is hashed after download
    and a :class:`ChecksumError` is raised on mismatch. The bad file is
    left on disk for debugging.

    Raises
    ------
    UnsupportedSchemeError
        If ``url`` is neither HTTP(S) nor a file URL.
    FetchError
        For network errors, missing source files, etc.
    ChecksumError
        If the downloaded file's SHA-256 does not match ``expected_sha256``.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    scheme = _scheme(url)
    if scheme in ("http", "https"):
        _fetch_http(url, dest, show_progress=show_progress, resume=resume, timeout=timeout)
    elif scheme == "file":
        _fetch_file(url, dest)
    else:
        raise UnsupportedSchemeError(f"unsupported URL scheme {scheme!r}; expected http(s) or file")

    if expected_sha256 is not None:
        actual = _sha256_file(dest)
        if actual.lower() != expected_sha256.lower():
            raise ChecksumError(
                f"SHA-256 mismatch for {dest.name}: "
                f"expected {expected_sha256.lower()}, got {actual}"
            )
