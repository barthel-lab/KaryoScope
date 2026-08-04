"""Source-format converters behind ``karyoscope prep-bed``.

Each converter turns one *source annotation format* into the two files
``karyoscope build`` consumes for a feature set: a 4-column BED whose 4th
column is a leaf label, and (where the format implies one) a ``child parent``
hierarchy. Converters are keyed by input format rather than by feature-set
name because unrelated formats can produce the same kind of set — RepeatMasker
output and an EDTA GFF3 both yield a ``repeat`` set but share no parsing at all.

Converters do **not** tile, gap-fill, flatten overlaps, or drop sequences:
``build`` already does all of that (``background:``, ``flatten:``, ``exclude:``).
The one exception is semantic tiling that only the source format can do — the
p-arm/q-arm split around a satellite core in :mod:`.structural`.
"""

from karyoscope.core.prep.common import PrepResult, render_stanza

__all__ = ["PrepResult", "render_stanza"]
