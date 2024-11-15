Galaxy workflow for the identification of candidate genes clusters
------------------------------------------------------------------

This approach screens two proteins against all nucleotide sequence from the
NCBI nt database within hours on our cluster, leading to all organisms with an inter-
esting gene structure for further investigation. As usual in Galaxy workflows every
parameter, including the proximity distance, can be changed and additional steps
can be easily added. For example additional filtering to refine the initial BLAST
hits, or inclusion of a third query sequence.

.. image:: https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_genes_located_nearby/find_genes_located_nearby.png


Sample Data
===========

As an example, we will use two protein sequences from *Streptomyces aurantiacus*
that are part of a gene cluster, responsible for metabolite producion.

You can upload both sequences directly into Galaxy using the "Upload File" tool
with either of these URLs - Galaxy should recognise this is FASTA files.

* `WP_037658548.fasta <https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_genes_located_nearby/WP_037658548.fasta>`_
* `WP_037658557.fasta <https://raw.githubusercontent.com/bgruening/galaxytools/master/workflows/ncbi_blast_plus/find_genes_located_nearby/WP_037658557.fasta>`_

In addition you can find both sequences at the NCBI server:
 * http://www.ncbi.nlm.nih.gov/protein/739806622 (cytochrome P450)
   ::
   
     >gi|739806622|ref|WP_037658557.1| cytochrome P450 [Streptomyces aurantiacus]
     MQRTCPFSVPPVYTKFREESPITQVVLPDGGKAWLVTKYDDVRAVMANPKLSSDRRAPDFPVVVPGQNAA
     LAKHAPFMIILDGAEHAAARRPVISEFSVRRVAAMKPRIQEIVDGFIDDMLKMPKPVDLNQVFSLPVPSL
     VVSEILGMPYEGHEYFMELAEILLRRTTDEQGRIAVSVELRKYMDKLVEEKIENPGDDLLSRQIELQRQQ
     GGIDRPQLASLCLLVLLAGHETTANMINLGVFSMLTKPELLAEIKADPSKTPKAVDELLRFYTIPDFGAH
     RLALDDVEIGGVLIRKGEAVIASTFAANRDPAVFDDPEELDFGRDARHHVAFGYGPHQCLGQNLGRLELQ
     VVFDTLFRRLPELRLAVPEEELSFKSDALVYGLYELPVTW


 * http://www.ncbi.nlm.nih.gov/protein/739806613 (beta-ACP synthase)
   ::
  
     >gi|739806613|ref|WP_037658548.1| beta-ACP synthase [Streptomyces aurantiacus]
     MSGRRVVVTGMEVLAPGGVGTDNFWSLLSEGRTATRGITFFDPAQFRSRVAAEIDFDPYAHGLTPQEVRR
     MDRAAQFAVVAARGAVADSGLDTDTLDPYRIGVTIGSAVGATMSLDEDYRVVSDAGRLDLVDHTYADPFF
     YNYFVPSSFATEVARLVGAQGPSSVVSAGCTSGLDSVGYAVELIREGTADVMVAGATDAPISPITMACFD
     AIKATTPRHDDPEHASRPFDDTRNGFVLGEGTAVFVLEELESARRRGARIYAEIAGYATRSNAYHMTGLR
     PDGAEMAEAITVALDEARMNPTAIDYINAHGSGTKQNDRHETAAFKRSLGEHAYRTPVSSIKSMVGHSLG
     AIGSIEIAASILAIQHDVVPPTANLHTPDPQCDLDYVPLNAREQIVDAVLTVGSGFGGFQSAMVLAQPER
     NAA


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

http://toolshed.g2.bx.psu.edu/view/bgruening/find_genes_located_nearby_workflow

Development is being done on github:

https://github.com/bgruening/galaxytools/tree/master/workflows/ncbi_blast_plus/find_genes_located_nearby


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

