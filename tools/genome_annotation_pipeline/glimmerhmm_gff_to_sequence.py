#!/usr/bin/env python
"""Convert GlimmerHMM GFF3 gene predictions into protein sequences.

This works with the GlimmerHMM GFF3 output format:

##gff-version 3
##sequence-region Contig5.15 1 47390
Contig5.15      GlimmerHMM      mRNA    323     325     .       +       .       ID=Contig5.15.path1.gene1;Name=Contig5.15.path1.gene1
Contig5.15      GlimmerHMM      CDS     323     325     .       +       0       ID=Contig5.15.cds1.1;Parent=Contig5.15.path1.gene1;Name=Contig5.15.path1.gene1;Note=final-exon

http://www.cbcb.umd.edu/software/GlimmerHMM/

Modified version of the converter from Brad Chapman: https://github.com/chapmanb/bcbb/blob/master/biopython/glimmergff_to_proteins.py

Usage:
    glimmer_to_proteins.py <glimmer output> <ref fasta> <output file> <convert to protein ... False|True>
"""
import sys
import os
import operator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from BCBio import GFF


def main(glimmer_file, ref_file, out_file, to_protein, ncbi_traslation_table=1):
    if to_protein in ["False", "false", "0"]:
        to_protein = False
    with open(ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    base, ext = os.path.splitext(glimmer_file)

    with open(out_file, "w") as out_handle:
        SeqIO.write(
            protein_recs(glimmer_file, ref_recs, to_protein, ncbi_traslation_table),
            out_handle,
            "fasta",
        )


def protein_recs(glimmer_file, ref_recs, to_protein, ncbi_traslation_table=1):
    """
    Generate protein records from GlimmerHMM gene predictions.
    """
    with open(glimmer_file) as in_handle:
        for rec in glimmer_predictions(in_handle, ref_recs):
            for feature in rec.features:
                seq_exons = []
                for cds in feature.sub_features:
                    if cds.type == "transcript":
                        """
                        Augustus seems to add a special features called transcript.
                        CDS is a sub_feature of that transcript, so go one step deeper
                        """
                        for cds in cds.sub_features:
                            if cds.type == "CDS":
                                seq_exons.append(
                                    rec.seq[
                                        cds.location.nofuzzy_start : cds.location.nofuzzy_end
                                    ]
                                )
                    else:
                        if cds.type == "CDS":
                            seq_exons.append(
                                rec.seq[
                                    cds.location.nofuzzy_start : cds.location.nofuzzy_end
                                ]
                            )
                gene_seq = reduce(operator.add, seq_exons)
                if feature.strand == -1:
                    gene_seq = gene_seq.reverse_complement()

                if to_protein:
                    yield SeqRecord(
                        gene_seq.translate(table=ncbi_traslation_table),
                        feature.qualifiers["ID"][0],
                        "",
                        "",
                    )
                else:
                    yield SeqRecord(gene_seq, feature.qualifiers["ID"][0], "", "")


def glimmer_predictions(in_handle, ref_recs):
    """
    Parse Glimmer output, generating SeqRecord and SeqFeatures for predictions
    """
    for rec in GFF.parse(in_handle, target_lines=1000, base_dict=ref_recs):
        yield rec


if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print __doc__
        sys.exit()
    main(*sys.argv[1:])
