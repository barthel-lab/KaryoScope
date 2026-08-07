"""Fetch PaVE's human reference papillomavirus genomes into one JSON array.

PaVE serves one record per genome from ``/api/genome/{id}``, and that record
carries the sequence, the ORF coordinates and the ICTV lineage. There is no
bulk endpoint that returns all three, so this walks the reference list and
writes ``pave_human_ref.json`` — the single input to
``karyoscope prep-bed pave``.

Records are sorted by locus id and serialised with fixed formatting, so two
runs against an unchanged database produce byte-identical output.

Usage:
    python fetch_pave_human_ref.py [OUTPUT]
"""

import json
import sys
import urllib.request

API = "https://pave.niaid.nih.gov/api"


def get(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.load(response)


def main(output="pave_human_ref.json"):
    listing = get(f"{API}/genome?limit=5000&includeNonRef=false")["data"]
    ids = sorted(g["locus_id"] for g in listing if "Human" in g["virus_tags"] and g["is_ref"])
    print(f"{len(ids)} human reference genomes", file=sys.stderr)

    records = []
    for n, locus in enumerate(ids, 1):
        records.append(get(f"{API}/genome/{locus}"))
        if n % 25 == 0 or n == len(ids):
            print(f"  {n}/{len(ids)}", file=sys.stderr)

    with open(output, "w") as out:
        json.dump(records, out, indent=1, sort_keys=True)
    print(f"wrote {output}", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
