# Test data

Small direct RNA sequencing (DRS) test set derived from the upstream tool repository:
<https://github.com/Dunning-Hotopp-Lab/Direct_RNA_Sequencing/tree/main/test-data>

It is a subset of *E. coli* E2348/69 DRS data from the tp.py manuscript
(Mattick et al. 2024, mBio, <https://doi.org/10.1128/mbio.02359-24>; SRA: SRR18061003,
SRR18061004), restricted to 11 ranges on two plasmids of assembly GCF_014117345.2.

Note: chromosome IDs are shortened (`NZ_CP059841` for `NZ_CP059841.1` and `NZ_CP059842`
for `NZ_CP059842.2`) so that test parameter values do not collide with Galaxy test-data
file lookups. The tool treats chromosome IDs as opaque strings, so this does not affect
the prediction results in any way.

## Inputs

* `E2348_69_plasmids.reads.bed` — BED6 of reads that lie entirely within a range and on
  the same strand (`bamToBed -i merged_all.bam > merged_all.readtable.bed`).
* `E2348_69_plasmids.ranges.bed` — distinct, non-overlapping per-strand ranges
  (chr, start, end, strand, read-count) from
  `bedtools merge -i merged_all.readtable.bed -s -c 6,4 -o distinct,count | awk '$5>20 {print $0}'`.

## Expected outputs (`expected/`)

All expected files are generated with the bioconda package `drs-transcript-prediction`
1.0.4 (verified to produce byte-identical output to the expected files committed upstream,
which were produced with v1.0.5):

```
drs-transcript-prediction -r E2348_69_plasmids.reads.bed -d E2348_69_plasmids.ranges.bed \
    -o expected/E2348_69_plasmids

drs-transcript-prediction -r E2348_69_plasmids.reads.bed -d E2348_69_plasmids.ranges.bed \
    --max_depth 200 -o expected/E2348_69_plasmids.max_depth_200

drs-transcript-prediction -r E2348_69_plasmids.reads.bed -d E2348_69_plasmids.ranges.bed \
    -c NZ_CP059842 --candidates_only -o expected/E2348_69_plasmids.NZ_CP059842

drs-transcript-prediction -r E2348_69_plasmids.reads.bed -d E2348_69_plasmids.ranges.bed \
    --min_region_positional_depth 1 --candidate_offset_len 50 --depth_delta_threshold 2 \
    -o expected/E2348_69_plasmids.tuned

drs-transcript-prediction -r E2348_69_plasmids.reads.bed -d E2348_69_plasmids.ranges.bed \
    --assigned_read_ratio_threshold 0.9 -o expected/E2348_69_plasmids.ratio09
```

Regenerate (and re-commit) these files whenever `@TOOL_VERSION@` is bumped.
