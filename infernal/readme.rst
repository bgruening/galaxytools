================================================
Galaxy wrapper for Infernal prediction tools
================================================

Infernal ("INFERence of RNA ALignment") is for searching DNA sequence databases 
for RNA structure and sequence similarities. It is an implementation of a special 
case of profile stochastic context-free grammars called covariance models (CMs). 
A CM is like a sequence profile, but it scores a combination of sequence consensus 
and RNA secondary structure consensus, so in many cases, it is more capable of 
identifying RNA homologs that conserve their secondary structure more than their 
primary sequence. 

This wrapper is copyright 2013 by:
 * Bjoern Gruening


This prepository contains wrapper for the Infernal_ command line tool.

.. _Infernal: http://infernal.janelia.org/


E. P. Nawrocki, D. L. Kolbe, and S. R. Eddy, Infernal 1.0: Inference of RNA alignments , Bioinformatics 25:1335-1337 (2009), . 


============
Installation
============

Please download install Infernal and the tool wrappers with the Galaxy Tool Shed:

=============
Miscellaneous
=============

Included in that repository is a CM datatype for INFERNAL 1.1. If you need that datatype in an additionl package,
I can source it out as separate package. Please contact me in that case.


=======
History
=======

interproscan:

 - v1.1.0: Initial public release




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

