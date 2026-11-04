# Future work

Notes on behaviour we know we want to change in KaryoScope itself but have not yet.
Each item records the observation that motivated it, so the fix can be tested against
the case that exposed it.

## Do not smooth across non-ACGT k-mers (scaffold N gaps)

**Observed 2026-09-11** on HPRC release-2 sample HG00128, hap2 chr15 (GenBank scaffold
`HG00128#2#CM090102.1`). The record carries a 2.39 Mb run of `N` at 95.36–97.75 Mb
inside 15q26.3. The smoothed BEDs label the whole run with the flanking feature in
every set (`15q26.3`, `q_arm`, `LINE`, `intron`), so 15q26.3 measures 6.3 Mb against a
4.0 Mb reference and downstream copy-number logic calls a 1.6x duplication. A 782 kb
`N` run on the same contig's p-arm stayed `novel`, so the fill is not deterministic.

**Root cause.** Two different things are currently labelled `novel`:

1. ACGT k-mers that are absent from the database (genuinely novel / divergent
   sequence). Smoothing across short runs of these is desired.
2. Non-standard k-mers (containing `N` or other non-ACGT symbols). These are almost
   always assembly scaffold gaps, and there is no sequence there to impute.

The smoother cannot tell them apart, so it treats a multi-megabase gap like a short
novel run and projects the flanking label across it. (`max_gap` does not prevent this:
it bounds the BED-coordinate distance between *related* intervals, and a contiguous
run of novel intervals is never more than `max_gap` from its neighbour.)

**Fix to make.** Carry the distinction through the pipeline: have `annotate` emit
non-ACGT positions with a separate sentinel (or a flag on the interval) and have
`smooth` refuse to bridge a run of them longer than some small threshold (a few
hundred bp), leaving the span `novel`. Hosted pre-computed annotations would need to be
regenerated to benefit.

**Interim workaround (KaryoScope-ISCN, merged in PR #7 and on by default; `--min-n-gap`
sets the threshold).** After
`remap-bed` and before binning, ISCN runs `seqtk gap` on the input FASTA, remaps the
gap BED through the scaffold map, and relabels the overlapping spans of every scaffolded
BED as `novel` — "unsmoothing" the gaps. That corrects ISCN's copy-number and dup/del
calls without touching the hosted BEDs, but every other consumer of the smoothed BEDs
still sees the filled gaps.

Prevalence in HPRC release 2 (16 haplotypes surveyed): 3–17 N runs ≥ 10 kb per
haplotype; 8 of 16 haplotypes have at least one run ≥ 1 Mb.
