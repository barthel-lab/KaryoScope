"""Build the dummy test database fixtures used by the test suite.

This script is not run by CI. It is run manually when the test fixtures
need to be regenerated — e.g., if we change the on-disk database layout,
or update KMC's index format.

Requirements:
    - kmc binary on $PATH (or set KMC env var)

What it produces:
    tests/data/dummy_db/
        manifest.yaml
        hierarchy.tsv
        features.tsv
        colors.txt
        index/
            features.kmc_pre
            features.kmc_suf
    tests/data/dummy_db.tar.gz   (gzipped tarball of the above)
    tests/data/dummy_db.sha256   (text file containing the tarball's SHA-256)

The build is deterministic: running it twice should produce byte-identical
output (modulo KMC's internal nondeterminism, which is minimal). Commit the
output files so CI can use them without needing kmc.

Usage:
    cd <repo root>
    python tests/data/build_dummy_db.py
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "dummy_db"
TARBALL = HERE / "dummy_db.tar.gz"
SHA256_FILE = HERE / "dummy_db.sha256"

DUMMY_DB_ID = "KS_dummy_test_v1"

# Fake DNA sequences. The point is to give KMC something to index that
# yields a real, valid .kmc_pre / .kmc_suf pair. The biology is meaningless.
#: A 23 bp seed sequence chosen so that:
#:
#: * its three overlapping 21-mers (call them m1, m2, m3) are pairwise
#:   distinct, and
#: * none of m1/m2/m3 is the reverse complement of any other.
#:
#: The second condition matters because KMC stores canonical k-mers
#: (min of forward and revcomp), so if any pair were revcomps their
#: counts would merge into a single bucket and break the
#: count-equals-featureID property. If you change this seed, regenerate
#: the fixture and verify the dump still shows exactly three k-mers
#: with counts 1, 2, 3.
#:
#: Together with the three-sequence layout below this yields a KMC index
#: whose counters are exactly 1, 2, and 3 — which means we can use the
#: counter values directly as the featureIDs in ``FEATURES_TSV``.
DUMMY_SEED = "ACGTGCTAGCTAGGCTATCGTAC"  # 23 bp

#: The dummy FASTA. Three sequences that follow the
#: "subsequence A once, B twice, C three times" pattern (with
#: A = m1, B = m2, C = m3, each a single 21-mer):
#:
#: * ``seq_a`` is the full 23 bp seed and contains all three k-mers
#:   (m1 at pos 0, m2 at pos 1, m3 at pos 2)
#: * ``seq_b`` is the seed minus its first base and contains m2, m3
#: * ``seq_c`` is the seed minus its first two bases and contains only m3
#:
#: Total occurrences across the corpus: m1=1, m2=2, m3=3. With canonical
#: counting (KMC's default), these map to KMC counters of 1, 2, 3
#: respectively, which then become featureIDs 1, 2, 3 in
#: :data:`FEATURES_TSV`.
DUMMY_FASTA = f"""\
>seq_a
{DUMMY_SEED}
>seq_b
{DUMMY_SEED[1:]}
>seq_c
{DUMMY_SEED[2:]}
"""

MANIFEST_YAML = """\
# Manifest for the dummy KaryoScope test database.
# This database is NOT biologically meaningful. It exists only to exercise
# the install and validation code paths in the test suite.
id: KS_dummy_test_v1
version: "1.0.0"
karyoscope_min_version: "0.1.0"
description: >
  Dummy database used by the KaryoScope test suite. Not for real use.

index:
  type: kmc
  basename: index/features

hierarchy: hierarchy.tsv
features: features.tsv
colors: colors.txt

kmer:
  size: 21
  type: fixed
  max_size: 21

feature_sets:
  - chromosome
  - region

roles:
  chromosome_assignment: chromosome

smoothing:
  recommended_window_bp: 1000
"""

#: Dummy hierarchy in the real (feature_set, child, parent) schema.
#:
#: Both feature sets have a three-level structure rooted at
#: ``categorized``, which is enough depth to exercise non-trivial LCA
#: computation in smoothing tests:
#:
#: chromosome
#:   chr1, chr2 → autosome → categorized
#:
#: region
#:   rA, rB → aSat → centromeric → categorized
#:   rC     → HSat → centromeric → categorized
#:
#: With this structure, LCA(rA, rB) = aSat (one level up) and
#: LCA(rA, rC) = centromeric (two levels up), so tests can verify both
#: the "promote to immediate ancestor" and "promote farther up" cases.
HIERARCHY_TSV = """\
feature_set\tchild\tparent
chromosome\tautosome\tcategorized
chromosome\tchr1\tautosome
chromosome\tchr2\tautosome
region\tcentromeric\tcategorized
region\taSat\tcentromeric
region\tHSat\tcentromeric
region\trA\taSat
region\trB\taSat
region\trC\tHSat
"""

FEATURES_TSV = """\
featureID\tchromosome\tregion
1\tchr1\trA
2\tchr1\trB
3\tchr2\trC
"""

COLORS_TXT = """\
feature_set\tfeature\tcolor
chromosome\tchr1\t#1f77b4
chromosome\tchr2\t#ff7f0e
region\trA\t#2ca02c
region\trB\t#d62728
region\trC\t#9467bd
"""


def _find_kmc() -> str:
    """Locate the kmc binary, preferring $KMC, then $PATH, then /tmp/bin."""
    candidate = os.environ.get("KMC")
    if candidate:
        if not Path(candidate).is_file():
            raise SystemExit(f"$KMC is set to {candidate} but file does not exist")
        return candidate
    kmc = shutil.which("kmc")
    if kmc:
        return kmc
    fallback = Path("/tmp/bin/kmc")
    if fallback.is_file():
        return str(fallback)
    raise SystemExit("could not locate the 'kmc' binary; install it or set the KMC env var")


def _run_kmc(kmc: str, fasta: Path, out_basename: Path, work_dir: Path) -> None:
    """Run KMC to produce ``out_basename.kmc_pre`` and ``.kmc_suf``.

    Flag rationale:

    * ``-k21`` matches the manifest's k-mer size.
    * ``-ci1`` includes singletons (our smallest count is 1).
    * ``-cs1000`` raises the counter ceiling well above 3 (the default is
      255, which is plenty for us, but we set this explicitly so any
      future fixture changes don't quietly saturate).
    * ``-fa`` treats input as FASTA.
    * ``-t1`` for reproducibility (no parallel-merge ordering).
    * ``-m2`` is KMC's minimum allowed memory budget (2 GB).

    KMC's default canonical-k-mer behaviour (forward and revcomp share a
    counter) is left on. Real KaryoScope databases are built canonical,
    and this fixture matches that convention. The :data:`DUMMY_SEED`
    docstring explains how the seed avoids revcomp collisions among its
    three k-mers, so canonical counting still yields counts 1, 2, 3.
    """
    cmd = [
        kmc,
        "-k21",
        "-ci1",
        "-cs1000",
        "-fa",
        "-t1",
        "-m2",
        str(fasta),
        str(out_basename),
        str(work_dir),
    ]
    subprocess.run(cmd, check=True)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def _make_tarball(src_dir: Path, archive: Path, top_level: str) -> None:
    """Create a deterministic tarball of ``src_dir`` rooted at ``top_level``.

    The tarball entries are sorted by name; mtimes are clamped to 0; uids
    and gids are zeroed. This makes the output byte-stable across runs and
    machines, which keeps the committed SHA-256 stable.
    """
    if archive.exists():
        archive.unlink()
    entries = sorted(src_dir.rglob("*"))
    # Include the top-level directory itself first.
    with tarfile.open(archive, "w:gz", compresslevel=9) as tar:
        # Add the root directory entry explicitly so empty parents don't
        # surprise consumers (though we have no empty dirs here).
        root_info = tarfile.TarInfo(name=top_level)
        root_info.type = tarfile.DIRTYPE
        root_info.mode = 0o755
        root_info.mtime = 0
        root_info.uid = 0
        root_info.gid = 0
        root_info.uname = ""
        root_info.gname = ""
        tar.addfile(root_info)
        for path in entries:
            rel = path.relative_to(src_dir)
            arcname = f"{top_level}/{rel.as_posix()}"
            info = tar.gettarinfo(str(path), arcname=arcname)
            info.mtime = 0
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            if info.isreg():
                with path.open("rb") as f:
                    tar.addfile(info, f)
            else:
                tar.addfile(info)


def main() -> int:
    kmc = _find_kmc()
    print(f"using kmc at: {kmc}")

    # Clean previous output.
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True)

    # Write text fixtures.
    _write_text(OUT_DIR / "manifest.yaml", MANIFEST_YAML)
    _write_text(OUT_DIR / "hierarchy.tsv", HIERARCHY_TSV)
    _write_text(OUT_DIR / "features.tsv", FEATURES_TSV)
    _write_text(OUT_DIR / "colors.txt", COLORS_TXT)
    (OUT_DIR / "index").mkdir()

    # Build the KMC index from a tiny FASTA in a tempdir, then move the
    # resulting .kmc_pre / .kmc_suf into place.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        fasta = tmp_path / "dummy.fa"
        fasta.write_text(DUMMY_FASTA)
        work_dir = tmp_path / "kmc_work"
        work_dir.mkdir()
        kmc_out = tmp_path / "features"
        _run_kmc(kmc, fasta, kmc_out, work_dir)
        for ext in (".kmc_pre", ".kmc_suf"):
            shutil.copy(kmc_out.with_suffix(ext), OUT_DIR / "index" / f"features{ext}")

    print(f"wrote dummy database to {OUT_DIR}")
    sizes = []
    for p in sorted(OUT_DIR.rglob("*")):
        if p.is_file():
            sizes.append((p.relative_to(OUT_DIR), p.stat().st_size))
    total = sum(s for _, s in sizes)
    for rel, size in sizes:
        print(f"  {rel}: {size} bytes")
    print(f"  total: {total} bytes")

    # Tar it up.
    _make_tarball(OUT_DIR, TARBALL, top_level=DUMMY_DB_ID)
    digest = _sha256(TARBALL)
    SHA256_FILE.write_text(f"{digest}  {TARBALL.name}\n")
    print(f"\nwrote tarball: {TARBALL} ({TARBALL.stat().st_size} bytes)")
    print(f"sha256: {digest}")
    print(f"wrote checksum: {SHA256_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
