=======================================
Galaxy wrapper for Glimmer gene calling
=======================================

This wrapper is copyright 2012-2013 by Björn Grüning.

This is a wrapper for the command line tool of Glimmer3_.

.. _Glimmer3: http://www.cbcb.umd.edu/software/glimmer/

Glimmer is a system for finding genes in microbial DNA, 
especially the genomes of bacteria, archaea, and viruses. 
Glimmer (Gene Locator and Interpolated Markov ModelER) uses interpolated 
Markov models (IMMs) to identify the coding regions and distinguish them from noncoding DNA.

| S. Salzberg, A. Delcher, S. Kasif, and O. White. Microbial gene identification using interpolated Markov models, Nucleic Acids Research 26:2 (1998), 544-548.

| A.L. Delcher, D. Harmon, S. Kasif, O. White, and S.L. Salzberg. Improved microbial gene identification with GLIMMER, Nucleic Acids Research 27:23 (1999), 4636-4641.

| A.L. Delcher, K.A. Bratke, E.C. Powers, and S.L. Salzberg. Identifying bacterial genes and endosymbiont DNA with Glimmer. Bioinformatics (Advance online version) (2007). 


============
Installation
============

Since version 0.2 the recommended installation procedure is via the `Galaxy Tool Shed`.

.. _`Galaxy Tool Shed`: http://toolshed.g2.bx.psu.edu/view/bjoern-gruening/glimmer3

To install Glimmer3 manually, please download Glimmer3 from:: 

    http://www.cbcb.umd.edu/software/glimmer/glimmer302.tar.gz

and follow the installation instructions. You can also use packages from your distribution like http://packages.debian.org/stable/science/tigr-glimmer

To install the wrapper copy the glimmer3 folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example::

	<tool file="gene_prediction/tools/glimmer3/glimmer_w_icm.xml" />
	<tool file="gene_prediction/tools/glimmer3/glimmer_wo_icm.xml" />
	<tool file="gene_prediction/tools/glimmer3/glimmer_build-icm.xml" />

=======
History
=======

- v0.1: Initial public release
- v0.2: Add tool shed integration

===============================
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

