# Examples

## [`karyotypes/`](karyotypes/) — reference plots to compare against

Karyotype plots for six assemblies — CHM13, HG002, the HG008 tumour/normal
pair, an HPRC population sample, and Arabidopsis — with notes on what each one
shows and what to look for in your own output. Start there if you have run
`karyoscope karyotype` and want to know whether the result looks right.

## `run_demo.sh` — quick demo / install smoke test

Runs the full KaryoScope annotation step end to end in a few seconds, using the
tiny synthetic database that ships with the repository
(`tests/data/dummy_db.tar.gz`). It downloads nothing and needs no special
hardware.

```bash
bash examples/run_demo.sh
```

The script installs the synthetic database into a throwaway database root,
registers it, annotates [`demo.fa`](demo.fa), prints the resulting BED tracks,
and cleans up after itself. It is a smoke test on constructed inputs — it
confirms that both the Python package and the compiled `get_featureIDs` helper
work, but it is **not** a biological example. For the real-data workflow, see
the **Quick start** section of the top-level [README](../README.md).

### What you should see

`demo.fa` contains two sequences: one whose 21-mers are in the synthetic
database, and one (poly-T) whose k-mers are not. The smoothed chromosome track
should be:

```
seq_with_features	0	2	chr1
seq_with_features	2	3	chr2
seq_novel	0	10	novel
```

and the smoothed region track:

```
seq_with_features	0	1	rA
seq_with_features	1	2	rB
seq_with_features	2	3	rC
seq_novel	0	10	novel
```

### Notes

- `kmc` is **not** required to run the demo. `kmc` is only needed to *build* a
  database; querying an existing one uses the bundled `get_featureIDs` helper.
- The synthetic database is a real, structurally complete KaryoScope database
  (a genuine KMC index plus the manifest, hierarchy, features, and colors
  files). It is the same fixture the integration test suite runs against; see
  [`tests/data/README.md`](../tests/data/README.md) for how it is built.
