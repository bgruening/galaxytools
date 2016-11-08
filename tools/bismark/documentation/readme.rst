===============
Bismark Wrapper
===============

Bismark_ uses Bowtie or Bowtie2 to map bisulfite converted sequence reads to a reference genome and determine cytosine methylation states.

Publication: http://www.ncbi.nlm.nih.gov/pubmed/21493656

User Guide: http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.7.12.pdf

.. _bismark: http://www.bioinformatics.babraham.ac.uk/projects/bismark/

Preparation
===========

Create your reference index with *bismark_genome_preparation* in your normal Galaxy Bowtie2/Botwie index directory. It will create a Bisulfite_Genome folder directly in your Bowtie2/Bowtie index directory.
If you follow that approach you do not need to specify or modify an extra .loc file.
That wrapper will extract the path to the Bisulfite_Genome folder from ./tool-data/bowtie2_indices.loc or ./tool-data/bowtie_indices.loc.

=======
History
=======

- v0.7: Initial public release
- v0.7.8: update and add Tool Shed Integration
- v0.7.11.1 change default output to BAM, from now on samtools are required
- v0.7.11.2 added multi-threading to samtools (samtools > 0.1.19 is required)
- v0.7.12 upgrade to bismark 0.7.12 and fix a major slowdown
- v0.7.12.1 define a dependency to samtools 0.1.19


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

