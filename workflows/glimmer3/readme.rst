==============================
Glimmer3 gene calling workflow
==============================

This Tool Shed Repository contains a workflow for the gene prediction of from a given nucleotide FASTA file.

At first an interpolated context model (ICM) is build from a know set of genes, preferable from the closest relative available organism(s). In a following step this ICM model is used to predict genes on the second input. The output is a FASTA file with nucleotide sequences that is further converted to proteins sequences.

To run that worflow glimmer_ und the EMBOSS_ suite is required. Both can be installed from the Tool Shed.

.. _glimmer: http://www.cbcb.umd.edu/software/glimmer/
.. _EMBOSS: http://emboss.sourceforge.net/

| A. L. Delcher, K.A. Bratke, E.C. Powers, and S.L. Salzberg. Identifying bacterial genes and endosymbiont DNA with Glimmer. Bioinformatics (Advance online version) (2007).

EMBOSS: The European Molecular Biology Open Software Suite (2000) 
Rice,P. Longden,I. and Bleasby,A. 
Trends in Genetics 16, (6) pp276--277

************
Availability
************

This workflow is available on the main Galaxy Tool Shed:
http://toolshed.g2.bx.psu.edu/view/bgruening/glimmer_gene_calling_workflow

Development is being done here:
https://github.com/bgruening/galaxytools/workflows/glimmer3/
