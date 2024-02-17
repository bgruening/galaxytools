#!/usr/bin/env python
"""
Extract sequences from a GFF3 file.

This works for examlpe with the GlimmerHMM GFF3 output format:

#gff-version 3
#sequence-region Glarea 1 38755942
Glarea	GlimmerHMM	mRNA	1017	1644	.	-	.	ID=Glarea.path1.gene1;Name=Glarea.path1.gene1
Glarea	GlimmerHMM	CDS	1017	1126	.	-	2	ID=Glarea.cds1.1;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=final-exon
Glarea	GlimmerHMM	CDS	1258	1317	.	-	2	ID=Glarea.cds1.2;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=internal-exon
Glarea	GlimmerHMM	CDS	1638	1644	.	-	0	ID=Glarea.cds1.3;Parent=Glarea.path1.gene1;Name=Glarea.path1.gene1;Note=initial-exon
Glarea	GlimmerHMM	mRNA	2755	5399	.	-	.	ID=Glarea.path1.gene2;Name=Glarea.path1.gene2
Glarea	GlimmerHMM	CDS	2755	4826	.	-	2	ID=Glarea.cds2.1;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=final-exon
Glarea	GlimmerHMM	CDS	5071	5145	.	-	2	ID=Glarea.cds2.2;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=internal-exon
Glarea	GlimmerHMM	CDS	5202	5214	.	-	0	ID=Glarea.cds2.3;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=internal-exon
Glarea	GlimmerHMM	CDS	5337	5399	.	-	0	ID=Glarea.cds2.4;Parent=Glarea.path1.gene2;Name=Glarea.path1.gene2;Note=initial-exon
Glarea	GlimmerHMM	mRNA	6580	7583	.	+	.	ID=Glarea.path1.gene3;Name=Glarea.path1.gene3

http://www.cbcb.umd.edu/software/GlimmerHMM/

Its based on code from from Brad Chapman: https://github.com/chapmanb/bcbb/blob/master/biopython/glimmergff_to_proteins.py

"""
import os

from BCBio import GFF
from Bio import Seq, SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from utils import check_gff


def main(gff_file, ref_file, ofile, seq_type="CDS"):
    with open(ref_file) as in_handle:
        fasta_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    base, ext = os.path.splitext(gff_file)

    gff_iter = GFF.parse(gff_file, fasta_recs)
    recs = protein_recs(check_gff(gff_iter), fasta_recs, seq_type)
    SeqIO.write(recs, ofile, "fasta")


def protein_recs(gff_iter, ref_recs, seq_type, to_protein=False):
    """
    Extract Introns, Exons and mRNA sequences from the glimmer prediction.
    """
    seq_counter = 0
    for rec in gff_iter:
        print(rec)
        for feature in rec.features:
            if feature.type.lower() in seq_type.lower():
                seq_counter += 1
                iid = feature.qualifiers.get(
                    "ID", feature.qualifiers.get("Name", [str(seq_counter)])
                )[0]
                desc = "%s_%s" % (rec.name, iid)

                # Augustus special cases
                if "gene_id" in feature.qualifiers.keys():
                    iid = feature.qualifiers["gene_id"][0]
                    print(iid)
                    desc = "%s_%s %s" % (rec.name, iid, seq_type)
                    print(desc)
                else:
                    for qualifier, value in feature.qualifiers.items():
                        if value and qualifier.startswith("g"):
                            iid = qualifier
                            desc = "%s_%s %s" % (rec.name, iid, seq_type)
                if not to_protein:
                    seq = SeqRecord(
                        feature.extract(rec.seq), id=iid, name=iid, description=desc
                    )
                else:
                    seq = SeqRecord(
                        Seq.Seq(feature.extract(rec.seq), generic_dna).translate(
                            table=codon_table, to_stop=True, cds=is_complete_cds  # noqa F821  # SBCHECK
                        ),
                        id=iid,
                        name=iid,
                        description=desc,
                    )
                yield seq


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract Introns, Exons, mRNA and CDS from sequences from the glimmer prediction."
    )
    parser.add_argument("-o", "--output", dest="ofile", help="Output file path.")
    parser.add_argument(
        "-t",
        "--seq-type",
        dest="type",
        help="GFF sequence type you want to extract. For example CDS, gene ...",
    )
    parser.add_argument(
        "-s",
        "--sequence",
        dest="sequence_path",
        required=True,
        help="FASTA sequence file. The file which was used for the glimmer prediction.",
    )
    parser.add_argument(
        "-g",
        "--gff",
        dest="gff_path",
        required=True,
        help="GFF File from the glimmer prediction.",
    )
    options = parser.parse_args()

    main(options.gff_path, options.sequence_path, options.ofile, options.type)
