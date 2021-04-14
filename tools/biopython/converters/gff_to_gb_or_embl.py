#!/usr/bin/env python
"""
    Convert a GFF and the associated sequence (FASTA) into GenBank, embl or FASTA format.

    https://raw.github.com/chapmanb/bcbb/master/gff/Scripts/gff/gff_to_genbank.py

    Usage:
    gff_to_gb_or_embl.py <GFF annotation file> <FASTA sequence file> <Output file> <Output format genbank or embl>
"""
import argparse

from BCBio import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from utils import check_gff, fix_ncbi_id


def main(gff_file, fasta_file, outfile, oformat):

    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    gff_iter = add_translation(gff_iter)
    if oformat in ["genbank", "gb"]:
        SeqIO.write(check_gff(fix_ncbi_id(gff_iter)), outfile, oformat)
    else:
        SeqIO.write(check_gff(gff_iter), outfile, oformat)


def add_translation(recs):
    """"""
    for rec in recs:
        for f in rec.features:
            f.qualifiers.update({"translation": f.extract(rec.seq).translate(table=4)})
        yield rec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="GFF to GenBank converter.")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    parser.add_argument("--fasta", "-f", dest="fasta", help="FASTA sequence file")
    parser.add_argument("--gff", "-g", dest="gff", help="GFF file")
    parser.add_argument(
        "--oformat",
        "-k",
        dest="oformat",
        help="output format, embl or genbank",
        choices=["embl", "genbank"],
    )

    parser.add_argument("--outfile", "-o", dest="outfile", help="output file")
    options = parser.parse_args()
    main(options.gff, options.fasta, options.outfile, options.oformat)
