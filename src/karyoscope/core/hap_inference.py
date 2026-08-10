"""Haplotype-label inference for ``karyoscope scaffold``.

Real-world assemblies arrive in many shapes:

* one file per haplotype (the pangenome convention) — easy.
* one combined file with hap-tagged contig names (HG002 distributed as
  one FASTA, with hifiasm ``h1tg...`` / ``h2tg...``, ``haplotype1...`` /
  ``haplotype2...``, or verkko ``MATERNAL`` / ``PATERNAL`` markings) —
  needs splitting.
* one file for a haploid assembly (CHM13) — single label.
* mixed: a hap1.fa.gz + hap2.fa.gz + unassigned.fa.gz triple — needs
  per-file labels.

Scaffold's encoded contig name is always ``<chrom>_<hap>_<contig>[_rc]``
(decided in the Stage 5d-1 design discussion), so every contig has to
end up with a hap label one way or the other. This module implements
the inference rules without ever asking the user beyond ``-i`` and
``--split-haps``.

Rule summary (per :class:`InputSpec` in :mod:`karyoscope.core.scaffold`):

1. **Explicit name on -i** (e.g. ``-i hap1=foo.fa``) → every contig in
   that file carries that label, period. No regex applied.
2. **Multiple inputs, names omitted** → for each input, try matching
   the filename stem against the built-in pattern library; if a clean
   label is detected use it, otherwise fall back to positional
   ``input1`` / ``input2`` / ... labels.
3. **Single input, name omitted** → try the built-in patterns
   against each contig name. The set of matched labels determines what
   happens:

   * No matches → all contigs become ``hap1`` (warn once).
   * Exactly one distinct match → all contigs (including the
     non-matching ones) become that label.
   * Two or more distinct matches → matching contigs keep their
     matched label; non-matching contigs take the lexically first
     matched label, with a warning.

4. **``--split-haps REGEX``** → applied per contig, taking precedence
   over the built-in patterns. The regex's first capture group is the
   label. Contigs that don't match the regex fall through to the
   file-level label.

``unassigned`` is reserved and is only ever the result of an explicit
``-i unassigned=PATH``. The auto-inference rules never produce it.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from pathlib import Path

from karyoscope.core.io.fasta import read_fasta_contig_names

__all__ = [
    "assign_per_input_labels",
    "classify_contigs",
    "display_hap_labels",
    "infer_hap_from_contig",
    "infer_hap_from_filename",
    "read_fasta_contig_names",
    "short_hap_label",
]

logger = logging.getLogger(__name__)


# --- pattern library ------------------------------------------------


def _compile(pat: str) -> re.Pattern[str]:
    return re.compile(pat, re.IGNORECASE)


#: Built-in (regex, label-template) pairs applied against contig names
#: and filename stems. Order matters — first match wins.
#:
#: The label template supports ``\1`` back-references; for example the
#: hifiasm pattern captures the hap digit and substitutes it into
#: ``hap\1`` so ``h1tg000001l`` becomes ``hap1``.
_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    # hifiasm dual/trio: h1tg..., h2tg...
    (_compile(r"^h([12])tg"), r"hap\1"),
    # full word "haplotype1" / "haplotype2" (e.g. "haplotype1-0000001",
    # as emitted by some long-read assemblers). Listed before the
    # generic hap1/hap2 rule below because that rule's `hap([12])`
    # cannot match here ("hap" is followed by "l"), so the full word
    # needs its own pattern. Bounded by the trailing class so
    # "haplotype10" is not read as "hap1".
    (_compile(r"(?:^|[._\-/])haplotype([12])(?:[._\-/]|$)"), r"hap\1"),
    # explicit hap1/hap2 anywhere, bounded by word edges or .-_/
    (_compile(r"(?:^|[._\-/])hap([12])(?:[._\-/]|$)"), r"hap\1"),
    # h1/h2 alone (rare — accepts H1, H2 too)
    (_compile(r"(?:^|[._\-/])h([12])(?:[._\-/]|$)"), r"hap\1"),
    # verkko / common pedigree-phased: maternal / paternal (and short forms)
    (_compile(r"(?:^|[._\-/])(maternal|mat)(?:[._\-/]|$)"), "maternal"),
    (_compile(r"(?:^|[._\-/])(paternal|pat)(?:[._\-/]|$)"), "paternal"),
)


def _match_label(name: str, extra_patterns: Iterable[re.Pattern[str]] = ()) -> str | None:
    """Apply the built-in patterns (and any extras) to ``name``.

    Returns the first matching label, or ``None`` if no pattern hits.
    ``extra_patterns`` accept the same back-reference semantics as the
    built-in patterns; callers use this for ``--split-haps``.
    """
    for pat, template in _PATTERNS:
        m = pat.search(name)
        if m is not None:
            return m.expand(template)
    for pat in extra_patterns:
        m = pat.search(name)
        if m is not None:
            if m.lastindex:
                return m.group(1)
            return m.group(0)
    return None


# --- public entry points --------------------------------------------


def infer_hap_from_contig(name: str) -> str | None:
    """Return the inferred hap label for a contig name, or ``None``."""
    return _match_label(name)


def infer_hap_from_filename(stem: str) -> str | None:
    """Return the inferred hap label for a filename stem, or ``None``."""
    return _match_label(stem)


def assign_per_input_labels(
    paths_with_explicit_names: list[tuple[str | None, Path]],
) -> list[str]:
    """Resolve a per-input hap label for each entry, in input order.

    Used when a user supplied ``-i path1 -i path2 ...`` without
    explicit ``NAME=`` prefixes. For each entry, an explicit name (if
    present) wins; otherwise the filename stem is pattern-matched;
    otherwise the positional fallback ``input{n}`` is used.

    If the same label would be assigned to two inputs, ``input{n}``
    fallbacks are applied to disambiguate (since the encoded contig
    names embed the hap label and must be globally unique).
    """
    out: list[str | None] = [None] * len(paths_with_explicit_names)
    seen: dict[str, int] = {}

    def _accept(idx: int, label: str) -> None:
        if label in seen:
            # Collision; fall back to positional later.
            return
        out[idx] = label
        seen[label] = idx

    for i, (name, p) in enumerate(paths_with_explicit_names):
        if name is not None:
            _accept(i, name)
            continue
        stem = _strip_fasta_ext(p.name)
        guessed = infer_hap_from_filename(stem)
        if guessed is not None:
            _accept(i, guessed)

    # Positional fallback for anything still unset (no explicit name,
    # no pattern match, or pattern collision).
    pos = 1
    for i, label in enumerate(out):
        if label is None:
            while f"input{pos}" in seen:
                pos += 1
            fb = f"input{pos}"
            out[i] = fb
            seen[fb] = i
            pos += 1

    return [label for label in out if label is not None]


def short_hap_label(hap: str) -> str:
    """Compact display form for one hap label.

    ``hap<digits>`` becomes ``h<digits>``; anything else is reduced to its first
    character (so ``maternal``/``paternal`` render as ``m``/``p``).

    This form is LOSSY and may collide across labels -- callers that render more
    than one haplotype must go through :func:`display_hap_labels`, which keeps
    the set distinguishable.
    """
    if hap.startswith("hap") and hap[3:].isdigit():
        return f"h{hap[3:]}"
    return hap[:1]


def display_hap_labels(haps: Iterable[str]) -> dict[str, str]:
    """Map every hap label to a display label, preserving distinguishability.

    Hap labels are globally unique by construction (see
    :func:`assign_per_input_labels`, which falls back to ``input{n}`` precisely
    to keep them so). The compact form from :func:`short_hap_label` throws that
    away: ``HG00097_hap1`` and ``HG00097_hap2`` both shorten to ``H``, which
    rendered every haplotype column of a diploid karyotype with the same letter.

    So: use the compact form only when it stays unique across the whole set,
    otherwise fall back to the full labels for *every* hap. All-short or
    all-full keeps the columns readable as a group, rather than mixing widths.
    """
    ordered = list(dict.fromkeys(haps))
    compact = {hap: short_hap_label(hap) for hap in ordered}
    if len(set(compact.values())) == len(ordered):
        return compact
    logger.debug(
        "compact hap labels collide (%s); falling back to full labels",
        sorted(set(compact.values())),
    )
    return {hap: hap for hap in ordered}


def classify_contigs(
    contig_names: list[str],
    *,
    file_level_label: str,
    split_haps_regex: str | None = None,
    is_only_input: bool = False,
    explicit_name_given: bool = False,
) -> dict[str, str]:
    """Resolve each contig to a hap label.

    Parameters
    ----------
    contig_names
        Every contig name found in the input FASTA.
    file_level_label
        The hap label for this input, as resolved by
        :func:`assign_per_input_labels`.
    split_haps_regex
        Optional user-supplied regex applied per contig (overrides the
        built-in patterns). The regex's first capture group is the
        label. Empty matches and missing groups fall through to the
        file-level label.
    is_only_input
        ``True`` when the scaffold invocation has exactly one input.
        Enables the single-input contig-name inference rule (built-in
        patterns scan every contig name and may split the file into
        multiple haps).
    explicit_name_given
        ``True`` when the user passed ``-i NAME=PATH`` rather than
        ``-i PATH``. Disables per-contig auto-inference (the explicit
        name wins for everything in the file).
    """
    extras: tuple[re.Pattern[str], ...] = ()
    if split_haps_regex is not None:
        extras = (_compile(split_haps_regex),)

    if explicit_name_given:
        # User insisted on a label; honour it for every contig.
        return {name: file_level_label for name in contig_names}

    # When the user did not provide a name and there is only one input
    # file, try the single-input split inference.
    if is_only_input and split_haps_regex is None:
        return _single_input_inference(contig_names, file_level_label)

    # General path: per-contig patterns (built-in + optional
    # --split-haps), falling back to the file-level label when no
    # pattern matches.
    out: dict[str, str] = {}
    unmatched: list[str] = []
    for name in contig_names:
        # --split-haps takes precedence when provided; otherwise the
        # built-in patterns apply.
        label = _match_label(name, extras) if extras else _match_label(name)
        if label is None:
            out[name] = file_level_label
            unmatched.append(name)
        else:
            out[name] = label
    if extras and unmatched:
        logger.info(
            "%d contig(s) did not match --split-haps; assigned to file label %r",
            len(unmatched),
            file_level_label,
        )
    return out


def _single_input_inference(contig_names: list[str], fallback: str) -> dict[str, str]:
    """Single-input inference: scan contigs for built-in patterns.

    No matches → all contigs become ``hap1`` (warn once);
    one distinct match → all contigs become that label;
    two+ distinct matches → matched contigs keep their label,
    non-matching contigs take the lexically-first matched label with a
    warning.
    """
    per_contig: dict[str, str | None] = {}
    matched_labels: set[str] = set()
    for name in contig_names:
        label = _match_label(name)
        per_contig[name] = label
        if label is not None:
            matched_labels.add(label)

    if not matched_labels:
        if fallback == "hap1":
            logger.warning(
                "no haplotype patterns matched any contig in this input; "
                "all contigs will be labelled %r",
                "hap1",
            )
        # Fall back to file-level label (typically "hap1").
        return {name: fallback for name in contig_names}

    if len(matched_labels) == 1:
        only = next(iter(matched_labels))
        # Even contigs that didn't pattern-match in this case probably
        # belong to the one haplotype detected; assign them all.
        return {name: only for name in contig_names}

    # Multiple haps detected. Matched contigs keep their label; the
    # rest go to the lexically-first label, with a warning.
    default = sorted(matched_labels)[0]
    out: dict[str, str] = {}
    unmatched_count = 0
    for name, label in per_contig.items():
        if label is None:
            out[name] = default
            unmatched_count += 1
        else:
            out[name] = label
    if unmatched_count:
        logger.warning(
            "detected %d haplotypes %s; %d contig(s) matched no pattern and "
            "were assigned %r as the fallback",
            len(matched_labels),
            sorted(matched_labels),
            unmatched_count,
            default,
        )
    return out


# --- helpers --------------------------------------------------------


_FASTA_EXTS: tuple[str, ...] = (
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
)


def _strip_fasta_ext(name: str) -> str:
    """Strip a recognised FASTA extension from a filename."""
    name_lower = name.lower()
    for ext in _FASTA_EXTS:
        if name_lower.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


# ``read_fasta_contig_names`` lives in :mod:`karyoscope.core.io.fasta`;
# re-exported at the top of this module for back-compat with callers
# imported under the old name.
