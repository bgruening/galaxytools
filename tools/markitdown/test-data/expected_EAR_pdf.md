ERGA Assembly Report
v24.08.26
Tags: ERGA-BGE
|     | TxID  |     |     | 2934182    |     |     |
| --- | ----- | --- | --- | ---------- | --- | --- |
|     | ToLID |     |     | xgPhyFlav1 |     |     |
Species Phyllidia flava
|                   | Class          |                       |               | Gastropoda |               |     |
| ----------------- | -------------- | --------------------- | ------------- | ---------- | ------------- | --- |
|                   | Order          |                       | Nudibranchia  |            |               |     |
|                   | Genome Traits  |                       | Expected      |            | Observed      |     |
| Haploid size (bp) |                |                       | 3,011,504,214 |            | 3,343,225,297 |     |
|                   | Haploid Number | 12 (source: ancestor) |               |            |               | 28  |
|                   | Ploidy         | 2 (source: ancestor)  |               |            |               | 2   |
|                   | Sample Sex     |                       |               | XX         |               | XX  |
EBP metrics summary and curation notes
Obtained EBP quality metric for hap1: 7.8.Q6
The following metrics were automatically flagged as below EBP recommended standards
or different from expected:
. Observed Haploid Number is different from Expected
. QV value is less than 40 for hap1
. Kmer completeness value is less than 90 for hap1
. Assembly length loss > 3% for hap1
Curator notes
. Interventions/Gb: 2
. Contamination notes: "No presence of contaminants. Mitochondrial genome was removed
from the assembly"
. Other observations: "Large collapsed repeat in chr5, haplotypic inversion in chr12"

Quality metrics table
|              |            | Pre-curation  |             |               | Curated     |
| ------------ | ---------- | ------------- | ----------- | ------------- | ----------- |
|              | Metrics    |               |  hap1       |               |  hap1       |
|              | Total bp   | 3,660,021,265 |             | 3,343,225,297 |             |
|              | GC %       |               | 41.62       |               | 41.09       |
|              | Gaps/Gbp   |               | 34.97       |               | 22.13       |
| Total gap bp |            |               | 1,874,523   |               | 138,686     |
|              | Scaffolds  |               | 973         |               | 104         |
| Scaffold N50 |            |               | 111,683,922 |               | 125,437,426 |
| Scaffold L50 |            |               | 12          |               | 10          |
| Scaffold L90 |            |               | 35          |               | 24          |
|              | Contigs    |               | 1,101       |               | 178         |
|              | Contig N50 |               | 73,275,821  |               | 77,716,104  |
|              | Contig L50 |               | 16          |               | 14          |
|              | Contig L90 |               | 85          |               | 43          |
|              | QV         |               | 68.0694     |               | 6           |
| Kmer compl.  |            |               | 95.5615     |               | 9           |
| BUSCO sing.  |            |               | 95.3%       |               | 95.2%       |
| BUSCO dupl.  |            |               | 0.7%        |               | 0.7%        |
| BUSCO frag.  |            |               | 1.1%        |               | 1.1%        |
| BUSCO miss.  |            |               | 2.9%        |               | 3.0%        |
BUSCO 5.4.7 Lineage: mammalia_odb10 (genomes:24, BUSCOs:9226)

HiC contact map of curated assembly
hap1 [LINK]

K-mer spectra of curated assembly
dataset_a95c0e30-e3ee-4059-830c-e735acee dataset_d0f465df-7b4b-4011-a131-e707c50b
cf0b.dat 6695.dat
dataset_56f39779-5c88-4801-8bc1-b21787fc
aa45.dat

Post-curation contamination screening
hap1. Bubble plot circles are scaled by sequence length, positioned by coverage and
GC proportion, and coloured by taxonomy. Histograms show total assembly length
distribution on each axis.

Data profile
Data HiFi Bionano OmniC
Coverage 40x 10x 90x
Assembly pipeline
- Hifiasm
|_ ver: 0.19.4
|_ key param: HiC
|_ key param: l0
- purge_dups
|_ ver: 1.2.6
|_ key param: NA
- Bionano_solve
|_ ver: Galaxy_3.7.0
|_ key param: NA
- YaHS
|_ ver: 1.1
|_ key param: NA
Curation pipeline
- GRIT_Rapid
|_ ver: 2.0
|_ key param: NA
- HiGlass
|_ ver: 1.0
|_ key param: NA
Submitter: John Doe
Affiliation: Galaxy EU
Date and time: 2024-08-29 10:46:27 CEST