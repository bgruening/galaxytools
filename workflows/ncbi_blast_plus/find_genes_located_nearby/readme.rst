Galaxy workflow for the identification of candidate genes clusters
------------------------------------------------------------------

Sample Data
===========

As an example, we will use two protein sequences from *Streptomyces aurantiacus*
that are part of a gene cluster, responsible for metabolite producion.

You can upload both sequences directly into Galaxy using the "Upload File" tool
with either of these URLs - Galaxy should recognise this is FASTA files.


http://www.ncbi.nlm.nih.gov/protein/739806622 (cytochrome P450)

```
>gi|739806622|ref|WP_037658557.1| cytochrome P450 [Streptomyces aurantiacus]
MQRTCPFSVPPVYTKFREESPITQVVLPDGGKAWLVTKYDDVRAVMANPKLSSDRRAPDFPVVVPGQNAA
LAKHAPFMIILDGAEHAAARRPVISEFSVRRVAAMKPRIQEIVDGFIDDMLKMPKPVDLNQVFSLPVPSL
VVSEILGMPYEGHEYFMELAEILLRRTTDEQGRIAVSVELRKYMDKLVEEKIENPGDDLLSRQIELQRQQ
GGIDRPQLASLCLLVLLAGHETTANMINLGVFSMLTKPELLAEIKADPSKTPKAVDELLRFYTIPDFGAH
RLALDDVEIGGVLIRKGEAVIASTFAANRDPAVFDDPEELDFGRDARHHVAFGYGPHQCLGQNLGRLELQ
VVFDTLFRRLPELRLAVPEEELSFKSDALVYGLYELPVTW
```

http://www.ncbi.nlm.nih.gov/protein/739806613 (beta-ACP synthase)

```
>gi|739806613|ref|WP_037658548.1| beta-ACP synthase [Streptomyces aurantiacus]
MSGRRVVVTGMEVLAPGGVGTDNFWSLLSEGRTATRGITFFDPAQFRSRVAAEIDFDPYAHGLTPQEVRR
MDRAAQFAVVAARGAVADSGLDTDTLDPYRIGVTIGSAVGATMSLDEDYRVVSDAGRLDLVDHTYADPFF
YNYFVPSSFATEVARLVGAQGPSSVVSAGCTSGLDSVGYAVELIREGTADVMVAGATDAPISPITMACFD
AIKATTPRHDDPEHASRPFDDTRNGFVLGEGTAVFVLEELESARRRGARIYAEIAGYATRSNAYHMTGLR
PDGAEMAEAITVALDEARMNPTAIDYINAHGSGTKQNDRHETAAFKRSLGEHAYRTPVSSIKSMVGHSLG
AIGSIEIAASILAIQHDVVPPTANLHTPDPQCDLDYVPLNAREQIVDAVLTVGSGFGGFQSAMVLAQPER
NAA
```

Citation
========

If you use this workflow directly, or a derivative of it, or the associated
NCBI BLAST wrappers for Galaxy, in work leading to a scientific publication,
please cite:

Peter J. A. Cock, John M. Chilton, Björn Grüning, James E. Johnson, Nicola Soranzo
NCBI BLAST+ integrated into Galaxy

http://biorxiv.org/content/early/2015/01/21/014043
http://dx.doi.org/10.1101/014043


Availability
============

This workflow is available on the main Galaxy Tool Shed:

http://toolshed.g2.bx.psu.edu/view/bgruening/find_genes_located_nearby_workflow

Development is being done on github:

https://github.com/bgruening/galaxytools/workflows/ncbi_blast_plus/


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus
