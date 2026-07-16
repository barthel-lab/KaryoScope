"""The parsed specification driving ``karyoscope build``.

A :class:`BuildSpec` is the single input to
:func:`karyoscope.core.build.build_database`. It can be assembled two ways —
:meth:`BuildSpec.from_flags` (the simple CLI form) or :meth:`BuildSpec.from_yaml`
(a build-spec file for multi-feature-set databases) — and both funnel through
the same :meth:`BuildSpec.validate`.

Each feature set is one of two input modes:

* **mode A (BED):** ``bed`` set; the shared ``sequence`` genome is sliced into
  per-leaf FASTAs. Gaps are filled with a named ``background`` leaf by default.
* **mode B (FASTAs):** ``fastas`` or ``per_seq_file`` set; the feature FASTAs are
  used directly, with no coordinate system (so ``background`` is not allowed).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from karyoscope.exceptions import BuildError

#: Default gap-fill label when a mode-A feature set doesn't name one.
DEFAULT_BACKGROUND = "background"

#: Default maximum query length (``-s``) / k-mer size.
DEFAULT_S = 31

#: Sentinel distinguishing "background not mentioned" (→ default) from an
#: explicit ``null`` (→ gap-fill disabled).
_UNSET = object()


@dataclass
class FeatureSetSpec:
    """One feature set within a :class:`BuildSpec`."""

    name: str
    # Mode A:
    bed: Path | None = None
    # Mode B (exactly one of these):
    fastas: list[Path] | None = None
    per_seq_file: Path | None = None
    # Optional per-set inputs:
    background: str | None = None  # None = gap-fill disabled (or mode B)
    hierarchy: Path | None = None
    priority: Path | None = None
    colors: Path | None = None
    flatten: bool = False
    # Enable HKS variable-k queries: one index answers any k <= s (e.g. a k-sweep
    # from a single build). HKS variable-k needs a dummy node for each feature
    # sequence's start, so when *any* set requests it the base index is built
    # from the per-feature FASTAs (for mode-A sets, the genome slices we generate;
    # for mode-B, the input FASTAs) rather than the whole genome — see
    # :func:`karyoscope.core.build._build_base_index`. Mutually exclusive with
    # ``priority`` (an HKS constraint).
    variable_k: bool = False

    @property
    def mode(self) -> str:
        return "bed" if self.bed is not None else "fasta"


@dataclass
class BuildSpec:
    """A fully-resolved build request."""

    id: str
    version: str
    feature_sets: list[FeatureSetSpec]
    sequence: Path | None = None
    s: int = DEFAULT_S
    threads: int = 4
    mem_gigas: int = 8
    external_memory: Path | None = None
    forward_only: bool = False
    roles: dict[str, str] = field(default_factory=dict)
    smoothing: dict[str, object] = field(default_factory=dict)

    # -- constructors --------------------------------------------------

    @classmethod
    def from_flags(
        cls,
        *,
        db_id: str,
        version: str,
        sequence: Path | None,
        feature_beds: dict[str, Path],
        backgrounds: dict[str, str] | None = None,
        hierarchies: dict[str, Path] | None = None,
        priorities: dict[str, Path] | None = None,
        colors: dict[str, Path] | None = None,
        flatten: bool = False,
        variable_k: bool = False,
        s: int = DEFAULT_S,
        threads: int = 4,
        mem_gigas: int = 8,
        external_memory: Path | None = None,
        forward_only: bool = False,
    ) -> BuildSpec:
        """Build a spec from the simple ``--feature-set NAME=bed`` CLI form.

        Every feature set is mode A (BED). Each mode-A set gap-fills by default,
        using ``backgrounds[name]`` if given else :data:`DEFAULT_BACKGROUND`.
        ``flatten`` and ``variable_k`` apply to every set.
        """
        backgrounds = backgrounds or {}
        hierarchies = hierarchies or {}
        priorities = priorities or {}
        colors = colors or {}

        sets: list[FeatureSetSpec] = []
        for name, bed in feature_beds.items():
            sets.append(
                FeatureSetSpec(
                    name=name,
                    bed=bed,
                    background=backgrounds.get(name, DEFAULT_BACKGROUND),
                    hierarchy=hierarchies.get(name),
                    priority=priorities.get(name),
                    colors=colors.get(name),
                    flatten=flatten,
                    variable_k=variable_k,
                )
            )
        spec = cls(
            id=db_id,
            version=version,
            sequence=sequence,
            feature_sets=sets,
            s=s,
            threads=threads,
            mem_gigas=mem_gigas,
            external_memory=external_memory,
            forward_only=forward_only,
        )
        spec.validate()
        return spec

    @classmethod
    def from_yaml(cls, path: Path) -> BuildSpec:
        """Parse and validate a build-spec YAML file."""
        if not path.is_file():
            raise BuildError(f"build spec not found: {path}")
        try:
            data = yaml.safe_load(path.read_text())
        except yaml.YAMLError as e:
            raise BuildError(f"could not parse build spec {path}: {e}") from e
        if not isinstance(data, dict):
            raise BuildError(f"build spec {path} must be a YAML mapping")

        base = path.parent

        def _resolve(value: object) -> Path:
            p = Path(str(value))
            return p if p.is_absolute() else (base / p)

        db_id = data.get("id")
        version = data.get("version")
        if not isinstance(db_id, str) or not db_id:
            raise BuildError(f"{path}: 'id' must be a non-empty string")
        if version is None:
            raise BuildError(f"{path}: 'version' is required")
        version = str(version)

        sequence = _resolve(data["sequence"]) if data.get("sequence") is not None else None

        kmer = data.get("kmer") or {}
        s = int(kmer.get("s", DEFAULT_S)) if isinstance(kmer, dict) else DEFAULT_S

        build = data.get("build") or {}
        if not isinstance(build, dict):
            raise BuildError(f"{path}: 'build' must be a mapping")
        threads = int(build.get("threads", 4))
        mem_gigas = int(build.get("mem_gigas", 8))
        ext = build.get("external_memory")
        external_memory = _resolve(ext) if ext is not None else None
        forward_only = bool(build.get("forward_only", False))

        raw_sets = data.get("feature_sets")
        if not isinstance(raw_sets, list) or not raw_sets:
            raise BuildError(f"{path}: 'feature_sets' must be a non-empty list")

        sets: list[FeatureSetSpec] = []
        for i, entry in enumerate(raw_sets):
            if not isinstance(entry, dict):
                raise BuildError(f"{path}: feature_sets[{i}] must be a mapping")
            name = entry.get("name")
            if not isinstance(name, str) or not name:
                raise BuildError(f"{path}: feature_sets[{i}] missing a non-empty 'name'")

            bed = _resolve(entry["bed"]) if entry.get("bed") is not None else None
            fastas = (
                [_resolve(p) for p in entry["fastas"]] if entry.get("fastas") is not None else None
            )
            per_seq_file = (
                _resolve(entry["per_seq_file"]) if entry.get("per_seq_file") is not None else None
            )

            bg_raw = entry.get("background", _UNSET)
            if bg_raw is _UNSET:
                # Not mentioned: default on for mode A, N/A for mode B.
                background = DEFAULT_BACKGROUND if bed is not None else None
            elif bg_raw is None:
                background = None  # explicitly disabled
            else:
                background = str(bg_raw)

            sets.append(
                FeatureSetSpec(
                    name=name,
                    bed=bed,
                    fastas=fastas,
                    per_seq_file=per_seq_file,
                    background=background,
                    hierarchy=_resolve(entry["hierarchy"]) if entry.get("hierarchy") else None,
                    priority=_resolve(entry["priority"]) if entry.get("priority") else None,
                    colors=_resolve(entry["colors"]) if entry.get("colors") else None,
                    flatten=bool(entry.get("flatten", False)),
                    variable_k=bool(entry.get("variable_k", False)),
                )
            )

        roles = data.get("roles") or {}
        smoothing = data.get("smoothing") or {}
        if not isinstance(roles, dict):
            raise BuildError(f"{path}: 'roles' must be a mapping")
        if not isinstance(smoothing, dict):
            raise BuildError(f"{path}: 'smoothing' must be a mapping")

        spec = cls(
            id=db_id,
            version=version,
            sequence=sequence,
            feature_sets=sets,
            s=s,
            threads=threads,
            mem_gigas=mem_gigas,
            external_memory=external_memory,
            forward_only=forward_only,
            roles=dict(roles),
            smoothing=dict(smoothing),
        )
        spec.validate()
        return spec

    # -- validation ----------------------------------------------------

    def validate(self) -> None:
        """Check the spec is internally consistent and inputs exist.

        Raises :class:`BuildError` on the first problem found.
        """
        if not self.id:
            raise BuildError("build spec: 'id' must be non-empty")
        if not self.version:
            raise BuildError("build spec: 'version' must be non-empty")
        if self.s < 1 or self.s > 256:
            raise BuildError(f"build spec: s must be in [1, 256], got {self.s}")
        if not self.feature_sets:
            raise BuildError("build spec: at least one feature set is required")

        seen: set[str] = set()
        needs_sequence = False
        for fs in self.feature_sets:
            if not fs.name or any(c in fs.name for c in "/\\ \t"):
                raise BuildError(
                    f"feature-set name {fs.name!r} is empty or contains a path/whitespace "
                    "character (it becomes a filename component)"
                )
            if fs.name in seen:
                raise BuildError(f"duplicate feature-set name {fs.name!r}")
            seen.add(fs.name)

            n_inputs = sum(x is not None for x in (fs.bed, fs.fastas, fs.per_seq_file))
            if n_inputs == 0:
                raise BuildError(
                    f"feature set {fs.name!r}: give one of bed / fastas / per_seq_file"
                )
            if n_inputs > 1:
                raise BuildError(
                    f"feature set {fs.name!r}: bed, fastas and per_seq_file are mutually exclusive"
                )

            if fs.bed is not None:
                needs_sequence = True
                _require_file(fs.bed, f"feature set {fs.name!r} bed")
            else:
                if fs.background is not None:
                    raise BuildError(
                        f"feature set {fs.name!r}: 'background' only applies to BED "
                        "(mode A) inputs; there is no coordinate system to gap-fill in "
                        "FASTA (mode B) inputs"
                    )
                if fs.flatten:
                    raise BuildError(
                        f"feature set {fs.name!r}: 'flatten' only applies to BED inputs"
                    )
                if fs.fastas is not None:
                    if not fs.fastas:
                        raise BuildError(f"feature set {fs.name!r}: 'fastas' list is empty")
                    for p in fs.fastas:
                        _require_file(p, f"feature set {fs.name!r} fasta")
                else:
                    _require_file(fs.per_seq_file, f"feature set {fs.name!r} per_seq_file")

            if fs.variable_k and fs.priority is not None:
                raise BuildError(
                    f"feature set {fs.name!r}: variable_k and priority are mutually "
                    "exclusive in HKS"
                )
            if fs.background == "none":
                raise BuildError(
                    f"feature set {fs.name!r}: background label 'none' is reserved by HKS"
                )
            for label, p in (
                ("hierarchy", fs.hierarchy),
                ("priority", fs.priority),
                ("colors", fs.colors),
            ):
                if p is not None:
                    _require_file(p, f"feature set {fs.name!r} {label}")

        if needs_sequence:
            if self.sequence is None:
                raise BuildError(
                    "a BED (mode A) feature set was given but no 'sequence' genome; "
                    "provide --sequence / sequence:"
                )
            _require_file(self.sequence, "sequence genome")


def _require_file(path: Path | None, what: str) -> None:
    if path is None or not Path(path).is_file():
        raise BuildError(f"{what}: file not found: {path}")
