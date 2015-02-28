===============================
Galaxy wrapper for RepeatMasker
===============================

This wrapper is copyright 2013 by Björn Grüning.

This is a wrapper for the command line tool of RepeatMasker from the Institute for Systems Biology.
http://www.repeatmasker.org/


Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-3.0.
1996-2010 <http://www.repeatmasker.org>. 


Additional Information:
Using RepeatMasker to identify repetitive elements in genomic sequences.
http://www.ncbi.nlm.nih.gov/pubmed/19274634

============
Installation
============

To install RepeatMasker, please use the following instructions:

http://www.repeatmasker.org/RMDownload.html

To install the wrapper copy the file RepeatMasker.xml in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
Add a line like the following:

Add the tool definition to your tool_conf.xml file under Galaxy root.
	<tool file="RepeatMasker/RepeatMasker.xml" />

=======
History
=======

- v1.1: Initial public release
- v0.1.1: patch from Simon Guest, to create empty files if no repeat is found
- v0.1.2: remove trailing semicolon, redirect all output to stdout

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

