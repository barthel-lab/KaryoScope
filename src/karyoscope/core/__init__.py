"""Core algorithmic logic and I/O shared across KaryoScope commands.

This package houses code used by multiple commands (and command-specific
orchestration that isn't pure presentation). Stage 4 introduces two
pieces:

* :mod:`karyoscope.core.external` — a small wrapper around
  :mod:`subprocess` for invoking external tools like ``kmc`` with
  consistent error handling and logging.
* :mod:`karyoscope.core.io` — parsers and writers for the file formats
  KaryoScope reads and produces.

Future stages will add command-specific modules
(``annotate.py``, ``bin.py``, ``scaffold.py``, etc.) and a couple of
cross-cutting algorithm modules (``smooth.py``). See
``docs/architecture.md`` (or the closest equivalent) for the design
rationale.
"""
