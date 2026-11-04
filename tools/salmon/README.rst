Galaxy wrappers for salmon
=====================================

These wrappers are copyright 2014, updated 2018 by Björn Grüning, Mhd Ramez Alrawas and additional contributors.
All rights reserved. See the licence text below.

Currently tested with salmon 2.11.0.


Automated Installation
======================

Galaxy should be able to automatically install the salmon for you if you have Conda dependency management enabled in Galaxy.
For more information please read the documentation at: https://github.com/galaxyproject/galaxy/blob/dev/doc/source/admin/conda_faq.rst


Upgrading from salmon 1.x to salmon 2.x
=======================================

Salmon 2.0 is a from-scratch rewrite in Rust. Two changes are important for
existing Galaxy users:

* **A new index format.** Indices built with the C++ salmon (up to and
  including 1.12.x) cannot be read by salmon 2.x and have to be rebuilt with
  ``salmon index``. The version tag used in the ``salmon_indexes_versioned``
  data table was therefore bumped from ``q7`` to ``r3`` (the value reported in
  the ``index_version`` field of the new ``info.json``). Existing local data
  tables should be rebuilt and the version column updated accordingly.
* **``salmon alevin`` was removed.** Single-cell quantification has moved to
  the alevin-fry ecosystem (https://github.com/COMBINE-lab/alevin-fry), so the
  alevin wrapper was dropped from this repository.

A few command-line flags that the wrappers used to pass were also removed or
renamed upstream (``--decoy`` became ``--decoys``, the long-read ``--ont``
mode and the alignment sampling options ``--sampleOut``/``--sampleUnaligned``
are gone, and the various ``--mimicBT2``-style scoring presets no longer
exist). New options exposed by the wrapper include ``--sketch``
(pseudoalignment mode) and ``--ignoreTxVersion`` (needed when combining an
Ensembl cDNA FASTA with an Ensembl GTF).

Salmon 2.x also requires read-name-grouped alignments for alignment-based
quantification; the wrapper now runs ``samtools collate`` on the input BAM.
Finally, the Gibbs-sampling and bootstrap options became mutually exclusive,
so they are only passed when set to a non-zero value.


History
=======

======== ======================================================================
Version  Changes
-------- ----------------------------------------------------------------------
v0.6.3.0  - First release
-------- ----------------------------------------------------------------------
v0.7.2    - Second release
-------- ----------------------------------------------------------------------
v0.8.2    - Third release
-------- ----------------------------------------------------------------------
v0.9.1    - Fourth release
-------- ----------------------------------------------------------------------
v0.11.2   - Fifth release
-------- ----------------------------------------------------------------------
v2.11.0   - Sixth release. Update to the salmon 2.x Rust rewrite: new index
           format (rebuild your indices), alevin removed (use alevin-fry),
           long-read quantification removed (use oarfish), renamed decoy flag,
           new sketch and transcript-version options.
-------- ---------------------------------------------------------------------

======== ======================================================================


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
