This is package is a Galaxy workflow for gene prediction using Glimmer3.

It uses the Glimmer3 tool (Delcher et al. 2007) trained on a known set of
genes to generate gene predictions on a new genome, and then calls EMBOSS to
translate the predictions into a FASTA file of predicted protein sequences.
The workflow requires two input files:

* Nucleotide FASTA file of know gene sequences (training set)
* Nucleotide FASTA file of genome sequence or assembled contigs

First an interpolated context model (ICM) is built from the set of known
genes, preferably from the closest relative organism(s) available. Next this
ICM model is used to predict genes on the genomic FASTA file. This produces
a FASTA file of the predicted gene nucleotide sequences, which is translated
into protein sequences using EMBOSS (Rice et al. 2000).

Glimmer is intended for finding genes in microbial DNA, especially bacteria,
archaea, and viruses.

See http://www.galaxyproject.org for information about the Galaxy Project.


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
http://dx.doi.org/10.1093/bioinformatics/btm009

For EMBOSS please cite:

Rice, P., Longden, I. and Bleasby, A. (2000)
EMBOSS: The European Molecular Biology Open Software Suite
Trends in Genetics 16(6), 276-277.
http://dx.doi.org/10.1016/S0168-9525(00)02024-2


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
