Galaxy wrapper for RepeatMasker
=====================================

This wrapper is copyright 2012 by Björn Grüning.

This is a wrapper for the command line tool of GlimmerHMM.
http://www.cbcb.umd.edu/software/GlimmerHMM/

GlimmerHMM is a new gene finder based on a Generalized Hidden Markov Model (GHMM). Although the gene finder conforms to the overall mathematical framework of a GHMM,
additionally it incorporates splice site models adapted from the GeneSplicer program and a decision tree adapted from GlimmerM. It also utilizes
Interpolated Markov Models for the coding and noncoding models.
Currently, GlimmerHMM's GHMM structure includes introns of each phase, intergenic regions, and four types of exons (initial, internal, final, and single).

Majoros, W.H., Pertea, M., and Salzberg, S.L. TigrScan and GlimmerHMM: two open-source ab initio eukaryotic gene-finders Bioinformatics 20 2878-2879.
Pertea, M. and S. L. Salzberg (2002). "Computational gene finding in plants." Plant Molecular Biology 48(1-2): 39-48.
The Arabidopsis Genome Initiative, (2000) "Analysis of the genome sequence of the flowering plant Arabidopsis thaliana", Nature. Dec 14; 408(6814):796-815.
Pertea, M., S. L. Salzberg, et al. (2000). "Finding genes in Plasmodium falciparum." Nature 404(6773): 34; discussion 34-5.
Salzberg, S. L., M. Pertea, et al. (1999). "Interpolated Markov models for eukaryotic gene finding." Genomics 59(1): 24-31. 


Installation
============

To install Glimmer3, please download GlimmerHMM from 

ftp://ftp.cbcb.umd.edu/pub/software/glimmerhmm/

and follow the installation instructions.
To extract the glimmerHMM predicted genes, the GFF Parser from Brad Chapman (ttp://github.com/chapmanb/bcbb/tree/master/gff) was used and is included.

To install the wrapper copy the glimmerHMM folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example:

<tool file="gene_prediction/tools/glimmerHMM/glimmerhmm_predict.xml" />
<tool file="gene_prediction/tools/glimmerHMM/glimmerhmm_to_sequence.xml" />


History
=======

v0.1 - Initial public release


Wrapper Licence (MIT/BSD style)
===============================

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

