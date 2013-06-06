Galaxy wrapper for Augustus
=====================================

This wrapper is copyright 2012 by Björn Grüning.

This is a wrapper for the command line tool of augustus.
http://bioinf.uni-greifswald.de/augustus/

AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences.

Oliver Keller, Martin Kollmar, Mario Stanke, Stephan Waack (2011)
A novel hybrid gene prediction method employing protein multiple sequence alignments
Bioinformatics, doi: 10.1093/bioinformatics/btr010

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008)
Using native and syntenically mapped cDNA alignments to improve de novo gene finding
Bioinformatics, doi: 10.1093/bioinformatics/btn013

Mario Stanke and Stephan Waack (2003)
Gene Prediction with a Hidden-Markov Model and a new Intron Submodel. 
Bioinformatics, Vol. 19, Suppl. 2, pages ii215-ii225




Installation
============

Install or downlaod augustus from:

http://bioinf.uni-greifswald.de/augustus/binaries/

and follow the installation instructions or copy the binaries into your $PATH

To install the wrapper copy the augustus folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example:

<section name="Gene Prediction" id="gene_prediction">
    <tool file="gene_prediction/tools/augustus/augustus.xml" />
</section>


Set the AUGUSTUS_CONFIG_PATH to /path_to_augustus/augustus/config with
    export AUGUSTUS_CONFIG_PATH=/path_to_augustus/augustus/config
or modify the wrapper and use the following additional commandline switch:
    --AUGUSTUS_CONFIG_PATH=/path_to_augustus/augustus/config




History
=======

v0.1 - Initial public release
v0.2 - Added tool_dependencies.xml file and update the augustus version (thanks to James Johnson)


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

