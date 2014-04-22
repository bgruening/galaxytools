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


The recommended citation for using Infernal 1.1 is E. P. Nawrocki and S. R. Eddy, Infernal 1.1: 100-fold faster RNA homology searches , Bioinformatics 29:2933-2935 (2013).


============
Installation
============

Please install Infernal and the tool wrappers with the Galaxy Tool Shed from

http://toolshed.g2.bx.psu.edu/view/bgruening/infernal 


=======
History
=======

 - v1.1.0: Initial public release


Bug Reports
===========

You can file an issue here https://github.com/bgruening/galaxytools/issues or ask
us on the Galaxy development list http://lists.bx.psu.edu/listinfo/galaxy-dev


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
