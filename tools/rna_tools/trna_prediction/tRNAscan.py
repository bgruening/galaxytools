#!/usr/bin/env python

"""
    Converts tRNAScan output back to fasta-sequences.
"""
import subprocess
import sys

from Bio import SeqIO


def main(args):
    """
    Call from galaxy:
    tRNAscan.py $organism $mode $showPrimSecondOpt $disablePseudo $showCodons $tabular_output $inputfile $fasta_output

        tRNAscan-SE $organism $mode $showPrimSecondOpt $disablePseudo $showCodons -Q -y -q -b -o $tabular_output $inputfile;
    """
    cmd = f"tRNAscan-SE -Q -y -q --brief {' '.join(args[:-1])}"
    child = subprocess.Popen(
        cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    stdout, stderr = child.communicate()
    return_code = child.returncode
    if return_code:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write(f"Return error code {return_code} from command:\n")
        sys.stderr.write(f"{cmd}\n")
    else:
        sys.stdout.write(stdout)
        sys.stdout.write(stderr)

    outfile = args[-1]
    sequence_file = args[-2]
    tRNAScan_file = args[-3]

    with open(sequence_file) as sequences:
        sequence_recs = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))

    tRNAs = []
    with open(tRNAScan_file) as tRNA_handle:
        for line in tRNA_handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            iid = cols[0].strip()
            start = int(cols[2])
            end = int(cols[3])
            aa = cols[4]
            codon = cols[5]
            rec = sequence_recs.get(iid)
            if start > end:
                new_rec = rec[end:start]
                new_rec.seq = new_rec.seq.reverse_complement()
                tRNAs.append(new_rec)
            else:
                new_rec = rec[start:end]

            new_rec.id = rec.id
            new_rec.name = rec.name
            new_rec.description = f"{rec.description} {aa} {codon} {start} {end}"
            tRNAs.append(new_rec)

    SeqIO.write(tRNAs, open(outfile, "w+"), "fasta")


if __name__ == "__main__":
    main(sys.argv[1:])
