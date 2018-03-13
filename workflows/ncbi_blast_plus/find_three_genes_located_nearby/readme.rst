Galaxy workflow for the identification of candidate genes clusters
------------------------------------------------------------------

This approach screens three proteins against a given genome sequence, leading to a genome position
were all three genes are located nearby. As usual in Galaxy workflows every
parameter, including the proximity distance, can be changed and additional steps
can be easily added. For example additional filtering to refine the initial BLAST
hits, or inclusion of a third query sequence.

.. image:: https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_three_genes_located_nearby/find_three_genes_located_nearby.png


Sample Data
===========

As an example, we will use three protein sequences from *Pan troglodytes* (Chimpanzee)
which are part of the β-globin cluster.

You can upload all sequences directly into Galaxy using the "Upload tool"
with either of these URLs - Galaxy should recognise this is FASTA files.

Query sequences:

* `P61920.fasta <https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_three_genes_located_nearby/P61920.fasta>`_
* `P61921.fasta <https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_three_genes_located_nearby/P61921.fasta>`_
* `Q6LDH1.fasta <https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_three_genes_located_nearby/Q6LDH1.fasta>`_

Genome sequence:

* http://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz


In addition you can find the query sequences at the UniProt server:
 * http://www.uniprot.org/uniprot/P61920 (Hemoglobin subunit gamma-1)
   ::

     >sp|P61920|HBG1_PANTR Hemoglobin subunit gamma-1 OS=Pan troglodytes GN=HBG1 PE=1 SV=2
     MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPK
     VKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFG
     KEFTPEVQASWQKMVTAVASALSSRYH


 * http://www.uniprot.org/uniprot/P61921 (Hemoglobin subunit gamma-2)
   ::

     >sp|P61921|HBG2_PANTR Hemoglobin subunit gamma-2 OS=Pan troglodytes GN=HBG2 PE=1 SV=2
     MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPK
     VKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFG
     KEFTPEVQASWQKMVTGVASALSSRYH


 * http://www.uniprot.org/uniprot/Q6LDH1 (Hemoglobin subunit epsilon)
   ::

     >sp|Q6LDH1|HBE_PANTR Hemoglobin subunit epsilon OS=Pan troglodytes GN=HBE1 PE=2 SV=3
     MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPK
     VKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFG
     KEFTPEVQAAWQKLVSAVAIALAHKYH


Running this workflow with the given inputs gave the following results.

===== ========= ==========
chrom     start        end       
----- --------- ----------
 chr1 169005092 169005184
 chr1 169004698 169004988
 chr1 169003907 169004035
 chr1 169001669 169001761
 chr1 169001086 169001340
 chr1 169000136 169000264
 chr1 168994127 168994219
 chr1 168993733 168993999
 chr1 168992844 168992975
 chr1 168988741 168988833
 chr1 168988508 168988624
 chr1 168988407 168988511
 chr1 168987452 168987583
 chr1 168977112 168977237
 chr1 168972535 168972633
 chr1 168972177 168972452
 chr1 168971401 168971529
 chr1 168965272 168965430
 chr1 168964417 168964698
 chr1 168958546 168958704
 chr1 168957691 168957972
 chr1 168952738 168952893
 chr1 168945810 168946055
 chr1 168945619 168945717
===== ========= ==========


Citation
========

If you use this workflow directly, or a derivative of it, or the associated
NCBI BLAST wrappers for Galaxy, in work leading to a scientific publication,
please cite:

Peter J. A. Cock, John M. Chilton, Björn Grüning, James E. Johnson, Nicola Soranzo
NCBI BLAST+ integrated into Galaxy

* http://biorxiv.org/content/early/2015/01/21/014043
* https://doi.org/10.1101/014043


Availability
============

This workflow is available on the main Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/bgruening/find_three_genes_located_nearby_workflow

Development is being done on github:

https://github.com/bgruening/galaxytools/tree/master/workflows/ncbi_blast_plus/find_three_genes_located_nearby


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus


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


