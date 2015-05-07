Galaxy wrapper for t-RNA prediction tools
=========================================

This wrapper is copyright 2012-2013 by Björn Grüning.

This prepository contains wrapper for the command line tools of tRNAscan-SE_ and Arogorn_.

.. _tRNAscan-SE: http://lowelab.ucsc.edu/tRNAscan-SE/
.. _Arogorn: http://mbio-serv2.mbioekol.lu.se/ARAGORN/

Dean Laslett and Bjorn Canback
ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences Nucl. Acids Res. (2004) 32(1): 11-16
doi:10.1093/nar/gkh152

Todd M. Lowe and Sean R. Eddy
tRNAscan-SE: A Program for Improved Detection of Transfer RNA Genes in Genomic Sequence Nucl. Acids Res. (1997) 25(5): 0955-964
doi:10.1093/nar/25.5.0955 


============
Installation
============

The t-RNA prediction wrappers are available through the toolshed_ and can be automatically installed.

.. _toolshed: http://toolshed.g2.bx.psu.edu/view/bjoern-gruening/trna_prediction

For manuel installation, please download tRNAscan-SE from the following URL and follow the install instructions::

	http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz

Arogorn can be download from::

	http://mbio-serv2.mbioekol.lu.se/ARAGORN/aragorn1.2.33.c

With a recent GNU-Compiler (gcc) you can compile it with the following command::

	gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.33.c

Please include aragorn and tRNAscan-SE into your PATH::

	export PATH=$PATH:/home/user/bin/aragorn/bin/


To install the wrappers copy the files aragorn.xml and tRNAscan.xml in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example add the following lines::

	<tool file="trna_prediction/aragorn.xml" />
	<tool file="trna_prediction/tRNAscan.xml" />


=======
History
=======

tRNAscan:

    - v0.1: Initial public release
    - v0.2: add fasta output
    - v0.2.1: added tool-dependency
    - v0.2.2: patch from Nicola Soranzo added
    - v0.3: add unit tests, documentation improvements, bug fixes

aragorn:

    - v0.1: Initial public release
    - v0.2: added options, upgrade to 1.2.36, tool-dependency
    - v0.3: add unit tests, documentation improvements
    - v0.4: added gff3 parser
    - v0.5: improved gff3 parser, now handles introns and full gene model



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

