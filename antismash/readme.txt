Galaxy wrapper for AntiSmash
=====================================

This wrapper is copyright 2012 by Björn Grüning.

This is a wrapper for the command line tool of antiSMASH.

antiSMASH allows the rapid genome-wide identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genomes.
It integrates and cross-links with a large number of in silico secondary metabolite analysis tools.

http://antismash.secondarymetabolites.org/

Marnix H. Medema, Kai Blin, Peter Cimermancic, Victor de Jager, Piotr Zakrzewski, Michael A. Fischbach, Tilmann Weber, Rainer Breitling & Eriko Takano (2011).
antiSMASH: Rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters. Nucleic Acids Research 39: W339-W346.


Installation
============

Currently these wrapper is tested with version 1.1 and the modified version of antismash.py included in that repository.

Install or downlaod antiSMASH from:

http://antismash.secondarymetabolites.org/download.html

... and follow the instructions.
Please replace the antismash.py file with the one inlcuded in that repository.
Edit the following lines in multi_antiSMASH_wrapper.py and antiSMASH_wrapper.py and adopt it to your installation.

blastdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
pfamdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
antismash_path = '/home/galaxy/bin/antismash-1.1.0/antismash.py'


To install the wrapper copy the antiSMASH folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example:

<section name="PKS-NRPS prediction" id="pks-nrps_prediction">
    <tool file="pks-nrps/tools/antiSMASH/antiSMASH.xml" />
    <tool file="pks-nrps/tools/antiSMASH/multi_antiSMASH.xml" />
</section>


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

