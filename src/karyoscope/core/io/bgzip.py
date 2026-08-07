"""In-place ``bgzip`` compression of output files.

Shared by every pipeline stage that compresses its BED/FASTA outputs
(``annotate``, ``scaffold``, ``centromeres``, ``karyotype``).
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

from karyoscope.core.external import require_tool, run_tool
from karyoscope.diskspace import format_bytes

logger = logging.getLogger(__name__)


def bgzip_file(path: Path, threads: int = 1) -> Path:
    """Compress ``path`` in-place with ``bgzip``, returning the new path.

    ``bgzip`` removes the source file by default (matches gzip's behaviour).
    Returns ``Path(str(path) + ".gz")``. Logs per-file start + completion
    at INFO so a long bgzip pass (12 files for a 6-feature-set human
    database) doesn't look like the pipeline has hung.

    ``threads`` is forwarded as ``bgzip -@``; the htslib bgzip compresses
    a single file in parallel when given more than one thread. We
    process files sequentially within the bgzip pass, so passing the
    user's full ``--threads`` here is the right call (no contention
    with concurrent file compressions). ``threads=1`` (the default)
    omits ``-@`` entirely for cleanest subprocess invocation.
    """
    bgzip = require_tool(
        "bgzip",
        install_hint="Install htslib (`conda install -c bioconda htslib`), "
        "or rerun with --no-bgzip to skip compression.",
    )
    orig_size = path.stat().st_size
    logger.info("bgzipping %s (%s, threads=%d)", path.name, format_bytes(orig_size), threads)
    t0 = time.perf_counter()
    cmd = [bgzip, "-f"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(path))
    run_tool(cmd)
    out_path = Path(str(path) + ".gz")
    out_size = out_path.stat().st_size if out_path.is_file() else 0
    dt = time.perf_counter() - t0
    logger.info(
        "bgzipped %s (%s -> %s) in %.1fs",
        out_path.name,
        format_bytes(orig_size),
        format_bytes(out_size),
        dt,
    )
    return out_path
