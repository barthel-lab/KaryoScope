"""alpha-satellite feature set from a CenSat annotation.

Where :func:`karyoscope.core.prep.structural.from_censat` collapses CenSat to
its broad classes (one ``aSat`` label for every array), this converter keeps the
per-array detail CenSat carries in the parenthetical — ``hor_1_5(S1C1/5/19H1L)``
becomes a record labelled ``S1C1_5_19H1L`` — so each alpha-satellite HOR array is
its own feature.

Two properties of the source drive everything here.

**CenSat names more than one array per interval.** ``hor_16_3(S1C16H1L,S1C1/5/19H1L)``
is one interval that genuinely contains two arrays' sequence interleaved. This
emits *one record per (interval, label) pair* rather than picking a winner:
``build`` writes each record into its own leaf's FASTA, and HKS then resolves
the shared k-mers to ``LCA(A, B)`` on its own. Flattening to a single label per
base would assert an ownership the annotation does not claim — and would be
undone by ``build --flatten`` anyway, which is why that flag must stay off for
this set.

**CenSat ships in two dialects.** The CenSatData release names the class
outright — ``active_hor(...)``, ``hor(...)``, ``dhor(...)``, ``mixedAlpha(...)``,
``mon`` — while the older UCSC track carries an interval index instead
(``hor_1_5(...)``, ``mon_1_2``) and marks a live array only by the trailing
``L`` of its name. Both are read. The ``L`` rule reproduces the release's
``active_hor`` set exactly on CHM13 v2 (20 arrays, no disagreement), so the two
dialects yield the same classification.

**Continuation names are bare suffixes.** ``hor_1_1(S3C1H2-A,B,C)`` means
S3C1H2-A, S3C1H2-**B** and S3C1H2-**C**, not three arrays called ``A``, ``B``
and ``C``. A plain comma split silently produces single-letter leaves;
:func:`expand_labels` reattaches the stem.

The hierarchy mirrors the shipped ``region`` set, so the two describe the same
biology with the same names::

    asat -> categorized
      alpha_hor -> asat
        active_hor  the live centromeric arrays
        hor         inactive higher-order repeats
        dhor        divergent ones
      mon -> asat

Arrays are **flat within their class**. Real array-to-array structure is a
phylogeny, which no annotation file contains; deriving it (mashtree over the
per-array sequences, say) is a separate step that replaces those stars.

The scaffold around them is what matters for correctness: ``mon`` and ``dhor``
sit *inside* ``asat``, so a k-mer from the conserved monomer core shared between
a named array and unnamed monomeric alpha-satellite resolves to ``asat`` rather
than escaping to the root. Leaving that sequence to ``build``'s gap-fill instead
is the one mistake this converter exists to prevent: ``background`` is a leaf at
the hierarchy root, so every k-mer a feature shares with it is labelled with the
root and paints nothing.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from karyoscope.core.prep.common import (
    Edge,
    PrepError,
    PrepResult,
    Record,
    SeqidRewriter,
    open_text,
    write_bed,
    write_colors,
    write_hierarchy,
)

#: The alpha-satellite classes a leaf can belong to, in hierarchy order.
#: ``active_hor`` is the live centromeric array, ``hor`` an inactive higher-order
#: repeat, ``dhor`` a divergent one, ``mon`` monomeric alpha-satellite.
CLASSES = ("active_hor", "hor", "dhor", "mon")

#: Root of the alpha-satellite subtree, and the HOR node beneath it. The shape
#: mirrors the shipped ``region`` set's ``aSat -> alpha_hor -> {active_hor, dhor,
#: hor, mixedAlpha}`` / ``mon -> aSat``, so the two feature sets describe the
#: same biology with the same names.
ASAT = "asat"
ALPHA_HOR = "alpha_hor"

#: Where each class's leaves hang.
CLASS_PARENT = {
    "active_hor": "active_hor",
    "hor": "hor",
    "dhor": "dhor",
    "mon": "mon",
}

#: ``child parent`` edges for the classes themselves.
CLASS_EDGES = {
    "active_hor": (("active_hor", ALPHA_HOR),),
    "hor": (("hor", ALPHA_HOR),),
    "dhor": (("dhor", ALPHA_HOR),),
    "mon": (("mon", ASAT),),
}

#: A CenSat class that names arrays without asserting where they belong.
#: ``mixedAlpha`` marks a stretch where two arrays' sequence is interleaved; the
#: arrays it names are classified by the records that do assert a class.
NAME_ONLY_CLASSES = ("mixedAlpha",)

#: Which class decides a leaf's home when it is named in more than one.
_PARENT_PRECEDENCE = ("active_hor", "dhor", "hor", "mon")

#: Priority per node; lower wins. A live array is the most specific call
#: available, so it takes a k-mer it shares with an inactive or divergent array;
#: HOR of any kind takes one it shares with monomeric alpha-satellite. Mirrors
#: the region set's ``centromeric`` < ``rDNA`` < ``arm``.
CLASS_PRIORITY = {
    ALPHA_HOR: 1,
    "mon": 2,
    "active_hor": 1,
    "hor": 2,
    "dhor": 3,
}

#: Priority given to ``asat`` and to the gap-fill leaf. They are deliberately
#: **equal**, so a k-mer shared between a feature and the rest of the genome
#: ties and falls back to the root rather than being claimed by either side.
#:
#: The gap-fill has to be named even though its value changes nothing: ``build``
#: requires a priority for every node, and an unlisted one would take the *best*
#: priority and outrank every feature.
ROOT_PRIORITY = 1

#: Priority for every leaf. Siblings must be all-equal or all-distinct, and
#: there is no basis for ranking one array above another: a tie resolves to
#: their common ancestor, which is the honest answer.
LEAF_PRIORITY = 1

#: Suprachromosomal family for each leaf, used as the legend group. An array
#: name opens with its family — ``S1C10H1L`` is SF1, ``S01C6H1L`` is SF01 — and
#: three arrays are chimeric, naming two families (``S01/1C3H1L``,
#: ``S4/6C13/14/21H1``, ``S4/6C13/14/21/22H8``). Those take the **first**, which
#: is the primary-family reading; the alternative inclusive reading would put one
#: array in two legend rows, and a legend row is one swatch.
_SF_RE = re.compile(r"^S(0?[0-9]+(?:_[0-9]+)*)C")

#: One colour per family, plus monomeric. A legend group is a single swatch, so
#: every array in a family paints the same colour; arrays are still separate
#: leaves in the annotation, and are told apart on a karyotype by which
#: centromere they sit in. Ordered light-to-dark within the warm/cool split so
#: neighbouring families stay distinguishable.
SF_COLORS = {
    "SF01": "#8DD3C7",
    "SF02": "#BC80BD",
    "SF1": "#E41A1C",
    "SF2": "#377EB8",
    "SF3": "#4DAF4A",
    "SF4": "#FF7F00",
    "SF5": "#984EA3",
    "SF6": "#A65628",
}

#: Monomeric alpha-satellite is not a family, so it gets its own row.
MON_COLOR = "#B3B3B3"

#: Interior nodes and the gap-fill. Interior nodes are never drawn (the binner
#: competes among leaves), so their colour is structural.
INTERIOR_COLOR = "#B0C4DE"
BACKGROUND_COLOR = "#808080"


def sf_of(label: str) -> str | None:
    """``S1C10H1L`` -> ``SF1``; ``S01_1C3H1L`` -> ``SF01`` (primary family)."""
    m = _SF_RE.match(label)
    return f"SF{m.group(1).split('_')[0]}" if m else None


#: A live array's name ends in ``L``. That is the only marker the older dialect
#: gives, and it agrees exactly with the newer dialect's ``active_hor`` class.
_LIVE_SUFFIX = "L"

#: The class token at the front of a CenSat label, in either dialect:
#: ``active_hor(...)`` / ``mon`` (release) and ``hor_1_2(...)`` / ``mon_X_3``
#: (UCSC track, which adds a sequence infix and an interval index). Alternatives
#: are ordered longest-first so ``dhor`` and ``active_hor`` are not truncated to
#: ``hor``.
_CLASS_RE = re.compile(r"^(active_hor|mixedAlpha|dhor|hor|mon)(?:_[A-Za-z0-9]+)*(?=\(|$)")

#: A bare continuation suffix: the ``B`` of ``S3C1H2-A,B``.
_BARE_SUFFIX_RE = re.compile(r"[A-Z]")

#: The variant token at the end of an array name — ``-A``, or the ``L`` that
#: marks a live array. Stripping it yields the stem a continuation attaches to.
_VARIANT_TAIL_RE = re.compile(r"(-[A-Z]|L)$")


def censat_class(raw: str) -> str | None:
    """The CenSat class token, or ``None`` for a non-alpha-satellite record.

    ``hor_1_2(S3C1H2-A,B)`` and ``hor(S3C1H2-A,B)`` both give ``hor``;
    ``bsat_13_1`` and ``HSat3_7`` give ``None``.
    """
    m = _CLASS_RE.match(raw)
    return m.group(1) if m else None


def home_class(label: str, record_class: str) -> str:
    """Which class a named array belongs to.

    A name ending in ``L`` is a live array wherever it is seen, which is the
    only signal the older dialect gives and agrees exactly with the newer
    dialect's ``active_hor``. Otherwise the array takes its record's class —
    except under ``mixedAlpha``, which names arrays without classifying them and
    so leaves them to be placed by the records that do.
    """
    if label.endswith(_LIVE_SUFFIX):
        return "active_hor"
    if record_class in NAME_ONLY_CLASSES:
        return "hor"
    return record_class


def safe_label(name: str) -> str:
    """``S1C1/5/19H1L`` -> ``S1C1_5_19H1L``.

    ``build`` writes one FASTA per leaf, named by the label, so a ``/`` would be
    read as a path separator. CenSat uses it for the multi-chromosome arrays.
    """
    return name.replace("/", "_")


def expand_labels(names: list[str]) -> list[str]:
    """Reattach the stem to bare continuation suffixes.

    ``["S3C1H2-A", "B", "C"]`` -> ``["S3C1H2-A", "S3C1H2-B", "S3C1H2-C"]``.
    A single capital letter continues the previous full name with its variant
    token replaced; anything else is a full name and becomes the new stem.
    """
    out: list[str] = []
    stem: str | None = None
    for name in names:
        if stem is not None and _BARE_SUFFIX_RE.fullmatch(name):
            out.append(f"{stem}-{name}")
        else:
            out.append(name)
            stem = _VARIANT_TAIL_RE.sub("", name)
    return out


def parse_names(raw: str) -> list[str]:
    """The array names a CenSat record lists, before expansion.

    ``hor_1_1(S3C1H2-A,B,C)`` -> ``["S3C1H2-A", "B", "C"]``. Empty for a record
    with no parenthetical (every ``mon`` record).
    """
    inner = raw.partition("(")[2].rpartition(")")[0]
    return [part.strip() for part in inner.split(",") if part.strip()]


def parse_labels(raw: str, cls: str) -> list[str]:
    """Leaf labels for one CenSat record.

    The parenthetical lists the arrays; an interval without one is labelled
    with its class.
    """
    names = parse_names(raw)
    if not names:
        return [cls]
    return [safe_label(name) for name in expand_labels(names)]


def from_censat(
    *,
    input_path: Path,
    output: Path,
    hierarchy: Path,
    priority: Path | None = None,
    colors: Path | None = None,
    background: str = "background",
    name: str = "asat",
    classes: tuple[str, ...] = CLASSES,
    rename: SeqidRewriter | None = None,
) -> PrepResult:
    """Turn a CenSat annotation into a per-array alpha-satellite feature set.

    ``classes`` selects which alpha-satellite classes to include. The default takes
    all three; dropping ``mon`` or ``dhor`` leaves that sequence to ``build``'s
    gap-fill, which costs the named arrays every k-mer they share with it.

    ``priority`` additionally writes a 3-column priority file ranking the
    classes (see :data:`CLASS_PRIORITY`), which resolves a k-mer shared between
    a named array and monomeric or divergent alpha-satellite to the array
    rather than to ``asat``.
    """
    rename = rename or SeqidRewriter()

    unknown_classes = [c for c in classes if c not in CLASSES]
    if unknown_classes:
        raise PrepError(
            f"unknown alpha-satellite class(es) {unknown_classes}; known: {list(CLASSES)}"
        )
    if not classes:
        raise PrepError("at least one alpha-satellite class must be selected")

    records: list[Record] = []
    # label -> the classes it was seen in, to decide its parent
    label_classes: dict[str, set[str]] = defaultdict(set)
    n_intervals = 0
    mixed_only: set[str] = set()
    n_multi = 0
    multi_bp = 0
    expansions: set[str] = set()

    with open_text(input_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            cls = censat_class(parts[3])
            if cls is None:
                continue
            if cls in NAME_ONLY_CLASSES:
                # Its arrays are only wanted if some HOR class was selected.
                if not any(c in classes for c in ("active_hor", "hor", "dhor")):
                    continue
            elif cls not in classes:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            seqid = rename(parts[0])

            names = parse_names(parts[3])
            labels = parse_labels(parts[3], cls)
            if len(labels) > 1:
                n_multi += 1
                multi_bp += end - start
            # `names` and `labels` align one-for-one unless the record had no
            # parenthetical, in which case `names` is empty and nothing expanded.
            for original, label in zip(names, labels, strict=len(names) == len(labels)):
                if _BARE_SUFFIX_RE.fullmatch(original):
                    expansions.add(label)
            n_intervals += 1
            for label in labels:
                if cls not in NAME_ONLY_CLASSES:
                    label_classes[label].add(home_class(label, cls))
                else:
                    # Record the sequence, but let a classifying record decide
                    # where the array lives; fall back only if none ever does.
                    label_classes.setdefault(label, set())
                    mixed_only.add(label)
                records.append((seqid, start, end, label))

    if not records:
        raise PrepError(
            f"{input_path}: no alpha-satellite records found for class(es) "
            f"{list(classes)} — is this a CenSat annotation?"
        )

    for label in mixed_only:
        if not label_classes[label]:
            label_classes[label].add(home_class(label, "mixedAlpha"))

    records.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    n_records = write_bed(output, records)
    edges = _edges(label_classes, classes)
    n_edges = write_hierarchy(hierarchy, edges)
    n_priority = 0
    if priority is not None:
        n_priority = _write_priority(priority, edges, background)
    n_colors = 0
    if colors is not None:
        n_colors = _write_colors(colors, edges, name, background)

    covered = _union_bp(records)
    notes = [
        f"{n_intervals:,} CenSat interval(s) over class(es) {', '.join(classes)} "
        f"-> {n_records:,} (interval, label) record(s), {covered:,} bp",
        f"{len(label_classes):,} leaf label(s)",
    ]
    if n_multi:
        notes.append(
            f"{n_multi} interval(s) name more than one array ({multi_bp:,} bp) and are "
            "emitted once per label — do NOT pass build --flatten for this set, or the "
            "shared sequence is forced onto a single leaf"
        )
    if expansions:
        shown = sorted(expansions)
        notes.append(
            f"expanded {len(shown)} bare continuation suffix(es) to full array names: "
            + ", ".join(shown[:6])
            + ("..." if len(shown) > 6 else "")
        )
    if priority is not None:
        # Report the sibling groups separately: priorities are only ever
        # compared within a group, so one flat list would imply an ordering
        # (`mon` against `hor`) that `plca` never evaluates.
        used = {c for v in label_classes.values() for c in v}
        groups = []
        under_asat = [ALPHA_HOR] if used & {"active_hor", "hor", "dhor"} else []
        under_asat += ["mon"] if "mon" in used else []
        if len(under_asat) > 1:
            groups.append(f"under {ASAT}: " + " < ".join(under_asat))
        under_hor = [c for c in ("active_hor", "hor", "dhor") if c in used]
        if len(under_hor) > 1:
            groups.append(f"under {ALPHA_HOR}: " + " < ".join(under_hor))
        notes.append(
            "priority (lower wins) — "
            + "; ".join(groups)
            + f". {background!r} is named explicitly at the same priority as {ASAT!r}, so the "
            "two tie and a k-mer shared with the rest of the genome falls back to the root"
        )
    missing = [c for c in CLASSES if c not in classes]
    if missing:
        notes.append(
            f"class(es) {', '.join(missing)} excluded — that alpha-satellite falls to build's "
            "background leaf, which sits at the hierarchy root, so every k-mer the named "
            "arrays share with it will resolve to the root and paint nothing"
        )

    return PrepResult(
        name=name,
        bed=output,
        n_records=n_records,
        background=background,
        hierarchy=hierarchy,
        n_edges=n_edges,
        priority=priority,
        n_priority=n_priority,
        colors=colors,
        n_colors=n_colors,
        notes=notes,
    )


def _write_priority(path: Path, edges: list[Edge], background: str) -> int:
    """Write the tree as 3-column ``child<TAB>priority<TAB>parent``.

    The gap-fill leaf is emitted even though ``build`` would attach it anyway:
    an unlisted node defaults to priority 0, the best value, which would make
    the background outrank every feature.
    """
    n = 0
    with path.open("w") as out:
        for child, parent in edges:
            prio = ROOT_PRIORITY if child == ASAT else CLASS_PRIORITY.get(child, LEAF_PRIORITY)
            out.write(f"{child}\t{prio}\t{parent}\n")
            n += 1
        out.write(f"{background}\t{ROOT_PRIORITY}\tcategorized\n")
    return n + 1


def _write_colors(path: Path, edges: list[Edge], set_name: str, background: str) -> int:
    """Colour every node, grouping the leaves by suprachromosomal family.

    Eighty-odd arrays cannot each have a legend row, and they cannot each have a
    distinct colour *and* be grouped — a legend group is one swatch, so
    ``build``/``karyotype`` reject a group spanning two colours. Grouping by
    family collapses the legend to nine rows while keeping the distinction that
    survives at karyotype scale: arrays are chromosome-specific, so which
    centromere a band sits in already says which array it is.

    An array whose family cannot be read from its name keeps its own row rather
    than being folded into a wrong one.
    """
    parents = dict(edges)
    interior = set(parents.values()) | {ASAT, ALPHA_HOR, *CLASSES}
    rows: list[tuple[str, str, str]] = []
    for child, _parent in edges:
        if child in interior and child != "mon":
            rows.append((child, INTERIOR_COLOR, ""))
            continue
        if child == "mon":
            rows.append((child, MON_COLOR, "monomeric"))
            continue
        sf = sf_of(child)
        if sf is not None and sf in SF_COLORS:
            rows.append((child, SF_COLORS[sf], sf))
        else:
            rows.append((child, INTERIOR_COLOR, ""))
    rows.append((background, BACKGROUND_COLOR, ""))
    return write_colors(path, set_name, rows)


def _edges(label_classes: dict[str, set[str]], classes: tuple[str, ...]) -> list[Edge]:
    """Scaffold edges plus one edge per leaf.

    A label seen in more than one class is placed by :data:`_PARENT_PRECEDENCE`.
    ``alpha_hor`` is emitted only when some HOR class produced a leaf, so a
    ``--class mon`` run does not carry an empty tier.
    """
    used = {c for v in label_classes.values() for c in v}
    edges: list[Edge] = [(ASAT, "categorized")]
    if used & {"active_hor", "hor", "dhor"}:
        edges.append((ALPHA_HOR, ASAT))
    for cls in CLASSES:
        if cls in used:
            edges.extend(CLASS_EDGES[cls])

    for label in sorted(label_classes):
        seen = label_classes[label]
        cls = next(c for c in _PARENT_PRECEDENCE if c in seen)
        parent = CLASS_PARENT[cls]
        # A class whose own name is the leaf (``mon``) is already an edge.
        if label != parent:
            edges.append((label, parent))
    return edges


def _union_bp(records: list[Record]) -> int:
    """Bases covered, merging the duplicate records multi-label intervals make."""
    by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for seqid, start, end, _label in records:
        by_seqid[seqid].append((start, end))
    total = 0
    for spans in by_seqid.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start > cur_end:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
            else:
                cur_end = max(cur_end, end)
        total += cur_end - cur_start
    return total
