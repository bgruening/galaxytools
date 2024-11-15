Galaxy wrapper for GlimmerHMM
=====================================

This wrapper has been written in 2012 by Björn Grüning, and updated by Rémi Marenco in 2016.

This is a wrapper for the command line tool of GlimmerHMM.
https://ccb.jhu.edu/software/glimmerhmm/

GlimmerHMM is a gene finder based on a Generalized Hidden Markov Model (GHMM). Although the gene finder conforms to the overall mathematical framework of a GHMM,
additionally it incorporates splice site models adapted from the GeneSplicer program and a decision tree adapted from GlimmerM. It also utilizes
Interpolated Markov Models for the coding and noncoding models.
Currently, GlimmerHMM's GHMM structure includes introns of each phase, intergenic regions, and four types of exons (initial, internal, final, and single).

Majoros, W.H., Pertea, M., and Salzberg, S.L. TigrScan and GlimmerHMM: two open-source ab initio eukaryotic gene-finders Bioinformatics 20 2878-2879.
Pertea, M. and S. L. Salzberg (2002). "Computational gene finding in plants." Plant Molecular Biology 48(1-2): 39-48.
The Arabidopsis Genome Initiative, (2000) "Analysis of the genome sequence of the flowering plant Arabidopsis thaliana", Nature. Dec 14; 408(6814):796-815.
Pertea, M., S. L. Salzberg, et al. (2000). "Finding genes in Plasmodium falciparum." Nature 404(6773): 34; discussion 34-5.
Salzberg, S. L., M. Pertea, et al. (1999). "Interpolated Markov models for eukaryotic gene finding." Genomics 59(1): 24-31. 


Installation
============

To install Glimmer3, please download GlimmerHMM from 

ftp://ccb.jhu.edu/pub/software/glimmerhmm

and follow the installation instructions.
To extract the glimmerHMM predicted genes, the GFF Parser from Brad Chapman (ttp://github.com/chapmanb/bcbb/tree/master/gff) was used and is included.

To install the wrapper copy the glimmerHMM folder in the galaxy tools
folder and modify the $GALAXY_ROOT/config/tool_conf.xml file to make the tool available to Galaxy.
For example:

```xml
<tool file="gene_prediction/tools/glimmerHMM/glimmerhmm_predict.xml" />
<tool file="gene_prediction/tools/glimmerHMM/glimmerhmm_to_sequence.xml" />
```

You also need to use a trained organism by adding them as reference data in Galaxy:

1. Add the *glimmer_hmm_trained_dir* data table to `tool_data_table_conf.xml` in `$GALAXY_ROOT/config/`:
        
    ```xml
    <!-- glimmer_hmm trained_dir -->
    <table name="glimmer_hmm_trained_dir" comment_char="#">
        <columns>value, name, path</columns>
        <file path="tool-data/glimmer_hmm.loc" />
    </table>
    ```
    
2. Add the `glimmer_hmm.loc` file referencing your trained organism, in `tool-data`.
    You have a sample [`glimmer_hmm.loc.sample`] available in the repository to help you configuring it properly
3. Add your data in the chosen folder at step 2. You can get them from the GlimmerHMM tar, `$GLIMMERHMM/trained_dir`

History
=======

- v2.0 - Update by Rémi Marenco to make it work without having to modify the wrapper + add ability to select the species
- v0.1 - Initial public release


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
