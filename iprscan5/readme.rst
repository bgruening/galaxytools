================================================
Galaxy wrapper for InterProScan prediction tools
================================================

InterProScan is a tool that combines different protein signature recognition methods native to the InterPro 
member databases into one resource with look up of corresponding InterPro and GO annotation.

This wrapper is copyright 2013 by:
 * Bjoern Gruening
 * Konrad Paszkiewicz


This prepository contains wrapper for the InterProScan_ command line tool.

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
http://dx.doi.org/10.7717/peerj.167

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

