"""Repeat annotations -> a ``repeat`` feature set.

Two unrelated sources produce the same kind of set, which is why ``prep-bed``
keys its subcommands on the input format rather than on the output name:

* :func:`from_repeatmasker` — RepeatMasker, in either the native ``.out`` table
  or the UCSC BED repackaging of it. Leaves are the RepeatMasker *classes*, and
  the hierarchy and palette reproduce the shipped ``HKS_human_CHM13_v2``
  ``repeat`` set exactly.
* :func:`from_edta` — an EDTA TE GFF3, whose ``Classification=`` vocabulary
  aliases each superfamily under both spelled-out and Wicker three-letter names.

Neither gap-fills: ``build`` adds the background leaf (conventionally
``nonrepeat``) itself.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from karyoscope.core.prep.common import (
    ColorRow,
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    open_text,
    write_bed,
    write_colors,
    write_hierarchy,
    write_priority,
)

# -- RepeatMasker -----------------------------------------------------

#: RepeatMasker class -> leaf label. Keys are the class part of ``class/family``
#: with any trailing ``?`` (RepeatMasker's "uncertain" marker) already stripped.
#: The values are the 15 non-background leaves of the shipped CHM13 v2 set.
RM_CLASS_TO_LEAF = {
    "DNA": "DNA",
    "LINE": "LINE",
    "Low_complexity": "Low_complexity",
    "LTR": "LTR",
    "RC": "RC",
    "Retroposon": "Retroposon",
    "rRNA": "rRNA",
    "Satellite": "Satellite",
    "scRNA": "scRNA",
    "Simple_repeat": "Simple_repeat",
    "SINE": "SINE",
    "snRNA": "snRNA",
    "srpRNA": "srpRNA",
    "tRNA": "tRNA",
    "Unknown": "Unknown",
}

#: Classes RepeatMasker emits that carry no repeat signal worth a leaf of their
#: own. Dropped rather than labelled, so ``build``'s gap-fill covers them.
RM_DROP_CLASSES = frozenset({"ARTEFACT", "Artefact"})

#: Leaf for a class absent from :data:`RM_CLASS_TO_LEAF`. RepeatMasker's own
#: catch-all is ``Unknown``, so unrecognised classes join it rather than being
#: skipped — skipping would leave those bases to the ``nonrepeat`` gap-fill,
#: which asserts the opposite of what the annotation says.
RM_FALLBACK_LEAF = "Unknown"

#: ``child parent`` edges of the shipped CHM13 v2 ``repeat`` tree, root last.
RM_HIERARCHY: list[Edge] = [
    ("repeat", "categorized"),
    ("Interspersed_Repeat", "repeat"),
    ("RNA", "Interspersed_Repeat"),
    ("rRNA", "RNA"),
    ("scRNA", "RNA"),
    ("snRNA", "RNA"),
    ("srpRNA", "RNA"),
    ("tRNA", "RNA"),
    ("Transposable_Element", "Interspersed_Repeat"),
    ("Class_I_Retrotransposition", "Transposable_Element"),
    ("LINE", "Class_I_Retrotransposition"),
    ("LINE-dependent_Retroposon", "Class_I_Retrotransposition"),
    ("Retroposon", "LINE-dependent_Retroposon"),
    ("SINE", "LINE-dependent_Retroposon"),
    ("LTR", "Class_I_Retrotransposition"),
    ("Class_II_DNA_Transposition", "Transposable_Element"),
    ("DNA", "Class_II_DNA_Transposition"),
    ("RC", "Class_II_DNA_Transposition"),
    ("Unknown", "Interspersed_Repeat"),
    ("Satellite", "repeat"),
    ("Noninterspersed", "repeat"),
    ("Low_complexity", "Noninterspersed"),
    ("Simple_repeat", "Noninterspersed"),
]

#: ``node -> (colour, legend_group)`` for the shipped palette. The five RNA
#: leaves share one colour, so grouping them collapses five legend rows into
#: one without breaking the "same group ⇒ same colour" invariant ``build``
#: enforces. Every other leaf has a colour of its own and stays ungrouped.
_INTERIOR = "#B0C4DE"
RM_COLORS: dict[str, tuple[str, str]] = {
    "repeat": (_INTERIOR, ""),
    "Interspersed_Repeat": (_INTERIOR, ""),
    "RNA": ("#90D5FF", "RNA"),
    "rRNA": ("#90D5FF", "RNA"),
    "scRNA": ("#90D5FF", "RNA"),
    "snRNA": ("#90D5FF", "RNA"),
    "srpRNA": ("#90D5FF", "RNA"),
    "tRNA": ("#90D5FF", "RNA"),
    "Transposable_Element": (_INTERIOR, ""),
    "Class_I_Retrotransposition": (_INTERIOR, ""),
    "LINE": ("#00FF00", ""),
    "LINE-dependent_Retroposon": (_INTERIOR, ""),
    "Retroposon": ("#800080", ""),
    "SINE": ("#FF0000", ""),
    "LTR": ("#0000FF", ""),
    "Class_II_DNA_Transposition": (_INTERIOR, ""),
    "DNA": ("#FFC0CB", ""),
    "RC": ("#FFEE8C", ""),
    "Unknown": ("#FF5C00", ""),
    "Satellite": ("#008B8B", ""),
    "Noninterspersed": (_INTERIOR, ""),
    "Low_complexity": ("#6B8E23", ""),
    "Simple_repeat": ("#DAA520", ""),
}

#: Colour for the ``nonrepeat`` gap-fill leaf ``build`` adds.
RM_BACKGROUND_COLOR = "#808080"


def _rm_leaf(class_family: str) -> str | None:
    """Map a RepeatMasker ``class/family`` to a leaf, or ``None`` to drop the row."""
    klass = class_family.split("/", 1)[0].rstrip("?")
    if klass in RM_DROP_CLASSES:
        return None
    return RM_CLASS_TO_LEAF.get(klass, RM_FALLBACK_LEAF)


def _looks_like_bed(fields: list[str]) -> bool:
    """True if a tab-split line is the UCSC BED repackaging rather than ``.out``."""
    return len(fields) >= 4 and fields[1].isdigit() and fields[2].isdigit() and "#" in fields[3]


def _parse_repeatmasker(
    path: Path, rename: SeqidRewriter
) -> tuple[list[Record], Counter[str], str]:
    """Parse either RepeatMasker dialect into records, keeping a class tally.

    Returns ``(records, class_counts, dialect)``. Sniffing here is safe in a way
    it would not be in ``build``: both dialects are the *same tool's* output, and
    the choice is reported rather than silently assumed.
    """
    records: list[Record] = []
    seen_classes: Counter[str] = Counter()
    dialect = ""
    with open_text(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith(("#", "track", "browser")):
                continue
            tab_fields = line.split("\t")
            if not dialect:
                dialect = "ucsc-bed" if _looks_like_bed(tab_fields) else "out"
            if dialect == "ucsc-bed":
                if not _looks_like_bed(tab_fields):
                    raise PrepError(f"{path}:{lineno}: not a UCSC RepeatMasker BED line: {line!r}")
                seqid, start, end = tab_fields[0], int(tab_fields[1]), int(tab_fields[2])
                # Column 4 is ``name#class/family``, e.g. ``TAR1#Satellite/subtelo``.
                class_family = (
                    tab_fields[3].split("#", 1)[1] if "#" in tab_fields[3] else tab_fields[3]
                )
            else:
                cols = line.split()
                # The 3-line ``.out`` banner has no leading integer SW score.
                if len(cols) < 11 or not cols[0].lstrip("-").isdigit():
                    continue
                seqid = cols[4]
                # ``.out`` positions are 1-based inclusive.
                start, end = int(cols[5]) - 1, int(cols[6])
                class_family = cols[10]
            seen_classes[class_family.split("/", 1)[0].rstrip("?")] += 1
            leaf = _rm_leaf(class_family)
            if leaf is None:
                continue
            records.append((rename(seqid), start, end, leaf))
    if not dialect:
        raise PrepError(f"{path}: no RepeatMasker records found")
    return records, seen_classes, dialect


def from_repeatmasker(
    *,
    input_path: Path,
    output: Path,
    hierarchy: Path,
    colors: Path | None = None,
    priority: Path | None = None,
    name: str = "repeat",
    background: str = "nonrepeat",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Convert RepeatMasker output into a ``repeat`` feature set."""
    rename = rename or SeqidRewriter()
    records, seen_classes, dialect = _parse_repeatmasker(input_path, rename)

    n_records = write_bed(output, records)
    n_edges = write_hierarchy(hierarchy, RM_HIERARCHY)
    n_colors = 0
    if colors is not None:
        rows: list[ColorRow] = [(node, *RM_COLORS[node]) for node, _p in RM_HIERARCHY]
        rows.append((background, RM_BACKGROUND_COLOR, ""))
        n_colors = write_colors(colors, name, rows)
    if priority is not None:
        write_priority(priority, RM_HIERARCHY)

    notes = [f"read {input_path} as the RepeatMasker '{dialect}' dialect"]
    unmapped = sorted(
        c for c in seen_classes if c not in RM_CLASS_TO_LEAF and c not in RM_DROP_CLASSES
    )
    if unmapped:
        total = sum(seen_classes[c] for c in unmapped)
        notes.append(
            f"{len(unmapped)} unrecognised class(es) ({total:,} rows) labelled "
            f"{RM_FALLBACK_LEAF}: {', '.join(unmapped)}"
        )
    dropped = sorted(c for c in seen_classes if c in RM_DROP_CLASSES)
    if dropped:
        total = sum(seen_classes[c] for c in dropped)
        notes.append(f"dropped {total:,} row(s) of non-repeat class: {', '.join(dropped)}")

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=background,
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        colors=colors,
        n_colors=n_colors,
        notes=notes,
    )


# -- EDTA -------------------------------------------------------------

#: EDTA ``Classification=`` value -> normalised leaf. EDTA aliases the same
#: superfamily under both spelled-out and Wicker three-letter names (``DNA/HAT``
#: and ``DNA/DTA`` are both hAT), so several keys collapse onto one leaf.
EDTA_ALIAS = {
    "LTR/Gypsy": "Gypsy",
    "LTR/Copia": "Copia",
    "LTR/unknown": "LTR_other",
    "LTR?": "LTR_other",
    "LINE/L1": "L1",
    "LINE/unknown": "LINE_other",
    "LINE?": "LINE_other",
    "LINE_element": "LINE_other",
    "SINE": "SINE_other",
    "SINE_element": "SINE_other",
    "RathE1_cons": "RathE",
    "RathE2_cons": "RathE",
    "RathE3_cons": "RathE",
    "DNA/HAT": "hAT",
    "DNA/DTA": "hAT",
    "DNA/En-Spm": "CACTA",
    "DNA/DTC": "CACTA",
    "DNA/MuDR": "Mutator",
    "DNA/DTM": "Mutator",
    "DNA/Harbinger": "PIF_Harbinger",
    "DNA/DTH": "PIF_Harbinger",
    "DNA/Mariner": "Tc1_Mariner",
    "DNA/Tc1": "Tc1_Mariner",
    "DNA/DTT": "Tc1_Mariner",
    "DNA/Pogo": "Tc1_Mariner",
    "MITE/DTA": "MITE",
    "MITE/DTM": "MITE",
    "MITE/DTC": "MITE",
    "MITE/DTT": "MITE",
    "MITE/DTH": "MITE",
    "RC/Helitron": "Helitron",
    "DNA/Helitron": "Helitron",
    "helitron": "Helitron",
    "DNA": "DNA_other",
    "Unknown": "TE_unclassified",
    "Unassigned": "TE_unclassified",
}

EDTA_FALLBACK_LEAF = "TE_unclassified"

#: ``child parent`` edges covering every :data:`EDTA_ALIAS` value.
EDTA_HIERARCHY: list[Edge] = [
    ("Transposable_Element", "categorized"),
    ("Class_I_Retro", "Transposable_Element"),
    ("LTR", "Class_I_Retro"),
    ("Gypsy", "LTR"),
    ("Copia", "LTR"),
    ("LTR_other", "LTR"),
    ("LINE", "Class_I_Retro"),
    ("L1", "LINE"),
    ("LINE_other", "LINE"),
    ("SINE", "Class_I_Retro"),
    ("SINE_other", "SINE"),
    ("RathE", "SINE"),
    ("Class_II_DNA", "Transposable_Element"),
    ("TIR", "Class_II_DNA"),
    ("hAT", "TIR"),
    ("CACTA", "TIR"),
    ("Mutator", "TIR"),
    ("PIF_Harbinger", "TIR"),
    ("Tc1_Mariner", "TIR"),
    ("MITE", "TIR"),
    ("Helitron", "Class_II_DNA"),
    ("DNA_other", "Class_II_DNA"),
    ("TE_unclassified", "Transposable_Element"),
]

_CLASSIFICATION = re.compile(r"Classification=([^;]+)")


def from_edta(
    *,
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None = None,
    name: str = "repeat",
    background: str = "nonrepeat",
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Convert an EDTA TE GFF3 into a ``repeat`` feature set.

    No colours are emitted: unlike RepeatMasker there is no established
    KaryoScope palette for the EDTA vocabulary, so ``build``'s automatic
    assignment is a better default than an invented one.
    """
    rename = rename or SeqidRewriter()
    records: list[Record] = []
    unmapped: Counter[str] = Counter()

    with open_text(input_path) as gff:
        for raw in gff:
            if raw.startswith("#"):
                if raw.startswith("##FASTA"):
                    break
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            match = _CLASSIFICATION.search(parts[8])
            classification = match.group(1) if match else ""
            leaf = EDTA_ALIAS.get(classification)
            if leaf is None:
                leaf = EDTA_FALLBACK_LEAF
                if classification:
                    unmapped[classification] += 1
            # GFF3 is 1-based inclusive; BED is 0-based half-open.
            records.append((rename(parts[0]), int(parts[3]) - 1, int(parts[4]), leaf))

    n_records = write_bed(output, records)
    n_edges = write_hierarchy(hierarchy, EDTA_HIERARCHY)
    if priority is not None:
        write_priority(priority, EDTA_HIERARCHY)

    notes: list[str] = []
    if unmapped:
        total = sum(unmapped.values())
        notes.append(
            f"{len(unmapped)} classification(s) ({total:,} rows) fell to "
            f"{EDTA_FALLBACK_LEAF}: {', '.join(sorted(unmapped))}"
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=background,
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        notes=notes,
    )
