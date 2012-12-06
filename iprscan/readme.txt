Galaxy wrapper for InterProScan prediction tools
=========================================

InterProScan is a tool that combines different protein signature recognition methods native to the InterPro member databases into one resource with look up of corresponding InterPro and GO annotation.

This wrapper is copyright 2012 by
*  Bjoern Gruening, Pharmaceutical Bioinformatics, University of Freiburg
*  Konrad Paszkiewicz, Exeter Sequencing Service, University of Exeter


This prepository contains wrapper for the InterProScan command line tool.
http://www.ebi.ac.uk/interpro/


Zdobnov E.M. and Apweiler R. "InterProScan - an integration platform for the signature-recognition methods in InterPro" Bioinformatics, 2001, 17(9): p. 847-8.



Installation
============

Please download install InterProScan according to:

ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/Installing_InterProScan.txt

Please see also:
ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/index.html

And rebuild the indizes if necessary:

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

<tool file="iprscan/interproscan.xml" />



History
=======

interproscan:
v0.1 - initial commit



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

