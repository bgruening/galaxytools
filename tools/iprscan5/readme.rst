==================================================
Galaxy wrapper for InterProScan 5 prediction tools
==================================================

InterProScan is a tool that combines different protein signature recognition methods native to the InterPro 
member databases into one resource with look up of corresponding InterPro and GO annotation.

This wrapper is copyright 2013 by:
 * Bjoern Gruening
 * Konrad Paszkiewicz


This prepository contains a wrapper for the InterProScan_ command line tool.

.. _InterProScan: http://www.ebi.ac.uk/interpro/interproscan.html


Quevillon E., Silventoinen V., Pillai S., Harte N., Mulder N., Apweiler R., Lopez R. (2005). InterProScan: protein domains identifier. Nucleic Acids Res. 33 (Web Server issue): W116-W120


============
Installation
============

Please download install InterProScan according to:

https://code.google.com/p/interproscan/wiki/HowToDownload


========
Citation
========

If you use this Galaxy tool in work leading to a scientific
publication, in addition to citing the invididual underlying tools, please cite:

Peter Cock, Bjoern Gruening, Konrad Paszkiewicz and Leighton Pritchard (2013).
Galaxy tools and workflows for sequence analysis with applications
in molecular plant pathology. PeerJ 1:e167
https://doi.org/10.7717/peerj.167

Full reference information is included in the help text.


=============
Input formats
=============

The standard interproscan input is either genomic or protein sequences. 
In the case of genomic sequences Interproscan will run an ORF prediction tool.


=======
History
=======

interproscan:

 - v5.0: Initial public release of version 5.0


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

