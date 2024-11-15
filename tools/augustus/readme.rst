Galaxy wrapper for Augustus
===========================

This wrapper is copyright 2012-2013 by Björn Grüning.

This is a wrapper for the command line tool of Augustus_.

.. _augustus: http://bioinf.uni-greifswald.de/augustus/

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

The recommended installation is by means of the toolshed_.
If you need to install it manually here is a short introduction.

.. _toolshed:  http://toolshed.g2.bx.psu.edu/view/bgruening/augustus


Install or downlaod augustus from::

    http://bioinf.uni-greifswald.de/augustus/binaries/

and follow the installation instructions or copy the binaries into your $PATH. To install the wrapper copy the augustus folder in the galaxy tools folder and modify the tools_conf.xml file to make the tool available to Galaxy.

For example::

  <section name="Gene Prediction" id="gene_prediction">
    <tool file="gene_prediction/tools/augustus/augustus.xml" />
  </section>


Set the *AUGUSTUS_CONFIG_PATH* to /path_to_augustus/augustus/config with::

  export AUGUSTUS_CONFIG_PATH=/path_to_augustus/augustus/config

or modify the wrapper and use the following additional commandline switch::

  --AUGUSTUS_CONFIG_PATH=/path_to_augustus/augustus/config


History
=======

- v0.1: Initial public release
- v0.2: Added tool_dependencies.xml file and update the augustus version (thanks to James Johnson)
- v0.3: upgrade to augustus 2.7, added new organisms and new parameters, output additional sequence files
- v0.3.1: added parallelism and changed the output parameters from boolean to a select box

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

