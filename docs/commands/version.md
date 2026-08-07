# karyoscope version

Print the KaryoScope version and environment information, useful for bug reports.

## Synopsis

```
karyoscope version
```

## Description

Prints the KaryoScope version, the Python interpreter and platform, the default database root and how many databases are installed, the versions of the Python dependencies, and whether the external tools (KMC, `get_featureIDs`, HKS, bgzip, tabix, seqtk) are found on PATH. If you only need the version string, `karyoscope --version` / `-V` prints just that. Please include this output in bug reports.

## Example output

```
KaryoScope 1.0.0
  Python: 3.12.12 (/path/to/python)
  Platform: macOS-26.5.1-arm64-arm-64bit
  Default database root: /home/you/.karyoscope/db
  Installed databases: 0

Python dependencies:
  click: 8.4.0
  drawsvg: 2.4.1
  cairosvg: 2.9.0
  requests: 2.34.2
  pyyaml: 6.0.3
  tqdm: 4.67.3
  jsonschema: 4.26.0

External tools:
  KMC: not found on PATH
  get_featureIDs: K-mer feature analysis tool (at /path/to/get_featureIDs)
  HKS: Running hks version 0.4.0 (at /path/to/hks)
  bgzip: bgzip (htslib) 1.22.1 (at /path/to/bgzip)
  tabix: tabix (htslib) 1.22.1 (at /path/to/tabix)
  seqtk: ... (at /path/to/seqtk)
```

## See also

- [`karyoscope info`](info.md) — inspect installed databases in more detail
