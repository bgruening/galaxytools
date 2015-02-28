==================================================
Galaxy wrapper for InterProScan 4 prediction tools
==================================================

**Note**:

This wrapper is for InterProScan 4.x if you want to use InterProScan 5 please have a look at:
http://toolshed.g2.bx.psu.edu/view/bgruening/interproscan5 

-----

InterProScan is a tool that combines different protein signature recognition methods native to the InterPro 
member databases into one resource with look up of corresponding InterPro and GO annotation.

This wrapper is copyright 2012-2013 by:
 *  Bjoern Gruening, Pharmaceutical Bioinformatics, University of Freiburg
 *  Konrad Paszkiewicz, Exeter Sequencing Service, University of Exeter


This prepository contains wrapper for the InterProScan_ command line tool.

.. _InterProScan: http://www.ebi.ac.uk/interpro/


Zdobnov E.M. and Apweiler R. "InterProScan - an integration platform for the signature-recognition methods in InterPro" Bioinformatics, 2001, 17(9): p. 847-8.


============
Installation
============

Please download install InterProScan according to:

ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/Installing_InterProScan.txt

Please see also:
ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/index.html

And rebuild the indizes if necessary

.. code:: 

	index_data.pl -f interpro.xml -inx -v -bin -bforce
	index_data.pl -f match_complete.xml -inx -v -bin -bforce
	index_data.pl -f Pfam-A.seed -inx -v -bin -bforce
	index_data.pl -f Pfam-C -inx -v -bin -bforce
	index_data.pl -f prints.pval -inx -v -bin -bforce
	index_data.pl -f sf.seq -inx -v -bin -bforce
	index_data.pl -f sf_hmm -inx -v -bin -bforce
	index_data.pl -f smart.HMMs -inx -v -bin -bforce
	index_data.pl -f superfamily.hmm -inx -v -bin -bforce
	index_data.pl -f TIGRFAMs_HMM.LIB -inx -v -bin -bforce


Add the tool definition to your tool_conf.xml file under Galaxy root:
.. code::

	<tool file="iprscan/interproscan.xml" />

=============
Input formats
=============

The standard interproscan input is either genomic or protein sequences. In the case of genomic sequences Interproscan will run an ORF 
prediction tool. However this tends to lose the ORF information (e.g. start/end co-ordinates) from the header. As such the requirement here is to input ORF 
sequences (e.g. from EMBOSS getorf) and to then replace any spaces in the FASTA header with underscores. This workaround generally preserves the relevant 
positional information. 


=======
History
=======

interproscan:

 - v1.1: Initial public release
 - v1.2: Merge with Konrad Paszkiewicz repository


=============
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

