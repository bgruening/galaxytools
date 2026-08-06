Conversion Suite powered by biopython
=====================================

This suite is copyright 2012 by Björn Grüning, Brad Chapman and Peter Cock

It's aim is to provide a place to host different small and useful scripts, that do not justify an own repository.
The first wrappers are all powered by biopython [1]. 
So if you like it please cite:

P.J.A. Cock, T. Antao, J.T. Chang, B.A. Chapman, C.J. Cox, A. Dalke, I. Friedberg, T. Hamelryck, F. Kauff, B. Wilczynski and M.J.L. de Hoon (2009)
Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, doi:10.1093/bioinformatics/btp163

[1] http://biopython.org/


Installation
============

Biopython is included in every major unix distribution and you can install with your favorite package manager.
If you need to install it manually please refer to:
http://biopython.org/DIST/docs/install/Installation.html


The GFF manipulation tools also need an additional library developed by Brad Chapman. You can find it under:
https://github.com/chapmanb/bcbb/tree/master/gff/BCBio

Make sure to put the BCBio folder somewhere in your PYTHONPATH.
For example with:
export PYTHONPATH=$PYTHONPATH:/home/user/path-to-BCBio/


To install the wrapper copy the converters folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example:

<section name="Converters" id="additional-converters">
    <tool file="converters/translate.xml" />
    <tool file="converters/reverse_comlement.xml" />
    <tool file="converters/gff_to_sequence.xml" />
    <tool file="converters/gff_to_gb_or_embl.xml" />
    <tool file="converters/gb_to_sequence.xml" />
    <tool file="converters/gb_embl_to_gff.xml" />
    <tool file="converters/sequence_converter.xml" />
</section>

It is probably a good idea to include the tools in one of the existing sections like 'Convert Formats'.

Contact
=======

If you encounter any problems, have patches or suggestions for improvments you can send an email to:

bjoern-at-gruenings-dot-eu


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

