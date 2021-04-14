#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#!/usr/bin/env python

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = ""
__license__ = "GLP3+"

import argparse
import os
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def remove_vector_contamination(query, univec_db, ofile="noncontamination.fasta"):

    tempf = tempfile.NamedTemporaryFile(delete=False)

    cmd = (
        "blastn -query %s -db %s -out %s -outfmt 6 -max_target_seqs 1 -evalue 0.001"
        % (query, univec_db, tempf.name)
    )
    subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

    seq_ids = dict()
    with open(tempf.name) as handle:
        for line in handle:
            if line.strip():
                columns = line.split("\t")
                start = columns[6]
                end = columns[7]
                seq_ids[columns[0]] = (int(start), int(end))

    os.remove(tempf.name)
    new_sequence_array = []
    with open(query, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in seq_ids.keys():
                sequence = str(record.seq)
                sequence = sequence.replace(
                    sequence[seq_ids[record.id][0] - 1 : seq_ids[record.id][1]], ""
                )
                record.seq = sequence
                d = Seq(sequence, IUPAC.IUPACUnambiguousDNA)
                temp = SeqRecord(
                    d, id=record.id, name=record.name, description=record.description
                )
                if len(sequence) > 200:
                    new_sequence_array.append(temp)
            else:
                new_sequence_array.append(record)

    output_file = open(ofile, "w+")
    SeqIO.write(new_sequence_array, output_file, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Removes contamination from FASTA sequences, utilizing the UniVec Database."
    )
    parser.add_argument(
        "-u",
        "--univec",
        dest="univec_path",
        required="True",
        help="Path to the UniVec BLAST Database. It is used for contamination checking.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_sequence",
        required=True,
        help="FASTA input file, can contain gaps in the sequence",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_sequence",
        required=True,
        help="FASTA output file, cleaned from contaminations.",
    )

    options = parser.parse_args()

    remove_vector_contamination(
        options.input_sequence, options.univec_path, ofile=options.output_sequence
    )
