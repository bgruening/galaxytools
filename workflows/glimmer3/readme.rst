This package is a Galaxy workflow for gene prediction using Glimmer3.

It uses the Glimmer3 tool (Delcher et al. 2007) trained on a known set of
genes to generate gene predictions on a new genome, and then calls EMBOSS
(Rice et al. 2000) to translate the predictions into a FASTA file of
predicted protein sequences. The workflow requires two input files:

* Nucleotide FASTA file of know gene sequences (training set)
* Nucleotide FASTA file of genome sequence or assembled contigs

First an interpolated context model (ICM) is built from the set of known
genes, preferably from the closest relative organism(s) available. Next this
ICM model is used to predict genes on the genomic FASTA file. This produces
a FASTA file of the predicted gene nucleotide sequences, which is translated
into protein sequences using the EMBOSS tool transeq.

Glimmer is intended for finding genes in microbial DNA, especially bacteria,
archaea, and viruses.

See http://www.galaxyproject.org for information about the Galaxy Project.


Sample Data
===========

As an example, we will use the first public assembly of the 2011 Shiga-toxin
producing *Escherichia coli* O104:H4 outbreak in Germany. This was part of the
open-source crowd-sourcing analysis described in Rohde et al. (2011) and here:
https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/wiki

You can upload this assembly directly into Galaxy using the "Upload File" tool
with either of these URLs - Galaxy should recognise this is a FASTA file with
3,057 sequences:

* http://static.xbase.ac.uk/files/results/nick/TY2482/TY2482.fasta.txt
* https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/blob/master/strains/TY2482/seqProject/BGI/assemblies/NickLoman/TY2482.fasta.txt

This FASTA file ``TY2482.fasta.txt`` was the initial TY-2482 strain assembled
by Nick Loman from 5 runs of Ion Torrent data released by the BGI, using the
MIRA 3.2 assembler. It was initially released via his blog,
http://pathogenomics.bham.ac.uk/blog/2011/06/ehec-genome-assembly/

We will also need a training set of known *E. coli* genes, for example the
model strain *Escherichia coli* str. K-12 substr. MG1655 which is well
annotated. You can upload the NCBI FASTA file ``NC_000913.ffn`` of the
gene nucleotide sequences directly into Galaxy via this URL, which Galaxy
should recognise as a FASTA file with 4,321 sequences:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.ffn

Then run the workflow, which should produce 2,333 predicted genes for the
TY2482 assembly (two FASTA files, nucleotide and protein sequences).


Citation
========

If you use this workflow directly, or a derivative of it, or the associated
Glimmer wrappers for Galaxy, in work leading to a scientific publication,
please cite:

Cock, P.J.A., Grüning, B., Paszkiewicz, K. and Pritchard, L. (2013)
Galaxy tools and workflows for sequence analysis with applications in
molecular plant pathology. (Submitted).

For Glimmer3 please cite:

Delcher, A.L., Bratke, K.A., Powers, E.C., and Salzberg, S.L. (2007)
Identifying bacterial genes and endosymbiont DNA with Glimmer.
Bioinformatics 23(6), 673-679.
https://doi.org/10.1093/bioinformatics/btm009

For EMBOSS please cite:

Rice, P., Longden, I. and Bleasby, A. (2000)
EMBOSS: The European Molecular Biology Open Software Suite
Trends in Genetics 16(6), 276-277.
https://doi.org/10.1016/S0168-9525(00)02024-2


Additional References
=====================

Rohde, H., Qin, J., Cui, Y., Li, D., Loman, N.J., et al. (2011)
Open-source genomic analysis of shiga-toxin-producing E. coli O104:H4.
New England Journal of Medicine 365, 718-724.
https://doi.org/10.1056/NEJMoa1107643


Availability
============

This workflow is available on the main Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/bgruening/glimmer_gene_calling_workflow

Development is being done on github:

https://github.com/bgruening/galaxytools/workflows/glimmer3/


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/bgruening/glimmer3
* http://toolshed.g2.bx.psu.edu/view/devteam/emboss_5
