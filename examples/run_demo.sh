#!/usr/bin/env bash
#
# KaryoScope quick demo.
#
# Runs the full annotation step end to end in a few seconds against the tiny
# synthetic database that ships with the repository (tests/data/dummy_db.tar.gz).
# It downloads nothing and needs no special hardware, so it is the fastest way
# to confirm that an installation works -- both the Python package and the
# compiled `get_featureIDs` helper.
#
# This is a smoke test on constructed inputs, not a biological example. See the
# "Quick start" section of the README for the real-data workflow.
#
# Usage:
#     bash examples/run_demo.sh
#
# Prerequisites: `karyoscope` installed and on PATH (see the README), and the
# C++ helper built (native/get_featureIDs/build/get_featureIDs). `kmc` is NOT
# required -- it is only needed to *build* a database, not to query one.

set -euo pipefail

# Resolve the repo root from this script's location so the demo can be run
# from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

TARBALL="$REPO_ROOT/tests/data/dummy_db.tar.gz"
DEMO_FASTA="$SCRIPT_DIR/demo.fa"
DB_ID="KS_dummy_test_v1"

if [[ ! -f "$TARBALL" ]]; then
    echo "error: cannot find the synthetic database at $TARBALL" >&2
    echo "       run this from a KaryoScope checkout that includes tests/data/." >&2
    exit 1
fi

if ! command -v karyoscope >/dev/null 2>&1; then
    echo "error: 'karyoscope' is not on PATH. Install it first (see the README)." >&2
    exit 1
fi

# Work in a throwaway directory so the demo never touches a real database root
# or leaves files behind.
WORK="$(mktemp -d)"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

export KARYOSCOPE_DB="$WORK/db"
mkdir -p "$KARYOSCOPE_DB"
OUTDIR="$WORK/out"

echo "==> Installing the bundled synthetic database into a throwaway root"
tar -xzf "$TARBALL" -C "$KARYOSCOPE_DB"
karyoscope register "$DB_ID"

echo
echo "==> Annotating examples/demo.fa"
karyoscope annotate -i "$DEMO_FASTA" -o "$OUTDIR" --no-bgzip

echo
echo "==> Smoothed chromosome track:"
cat "$OUTDIR/demo.${DB_ID}.chromosome.smoothed.bed"
echo
echo "==> Smoothed region track:"
cat "$OUTDIR/demo.${DB_ID}.region.smoothed.bed"

echo
echo "Demo complete. Expected: seq_with_features annotated as chr1/chr2 (rA/rB/rC),"
echo "and seq_novel labelled 'novel' (its k-mers are not in the database)."
