#!/usr/bin/env python

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = ""
__license__ = "GLP3+"

"""
    Housten we have a problem. NCBI does not accept sequences with gaps in it.
    Scaffold on the other hands can have large amounts of N-runs in it.
    These script, remove all N-runs and creates an additional sequence out of
    it. After that these script will give proper names for each entry with a
    unique ID.
    Usage:
    scipt.py old_scaffold_file.fasta new_contig_file.fasta
"""

import argparse
import os
import re

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run(
    ifile,
    cleand_sequence_file,
    agp_file="scaffold.agp",
    locus_tag="XYX",
    sequence_description="",
):
    # read input scaffold file
    """
    Reviewer suggestion:
    Seq1 [organism=Glarea lozoyensis 74030] [strain=74030]

    """

    # sequence_description = '[organism=Glarea lozoyensis 74030] [strain=74030]'

    if agp_file != False:
        agp_file = open(agp_file, "w+")

    new_sequence_array = []
    scaffold_counter = 1
    scaffold_name = locus_tag + "_scaffold%(scaffold_counter)s"
    agp_contig_string = (
        scaffold_name
        + "\t%(start)s\t%(end)s\t%(linecount)s\tW\t%(sequence_iid)s\t1\t%(end_part_sequence)s\t+\n"
    )
    agp_gap_string = (
        scaffold_name
        + "\t%(start)s\t%(end)s\t%(linecount)s\tN\t%(gap_len)s\tfragment\tyes\t\n"
    )
    sequence_counter = 0
    for record in SeqIO.parse(open(ifile), "fasta"):
        sequence_start = 0
        sequence_end = 0
        linecount = 0

        # if the whole sequence has a length under 200, we can skip it
        if len(record.seq) < 200:
            continue

        """
            Filter out all contigs from a scaffold that are shorter than 200
            nucleotides. Two steps, at first replace all short contigs in the
            sequence and afterwards a special case for the beginning of the
            line.
        """

        record_seq = str(record.seq)

        # Check if an arbitrary contig is shorter than 200 nucleotides
        # replace them with gaps
        def match_func(match):
            return "N" * len(match.group())

        record_seq = re.sub(r"(?=([nN][ACGTactg]{1,199}[nN]))", match_func, record_seq)

        # Check if the sequence starts with a contig shorter than 200
        # replace it with N's and remove it
        record_seq = re.sub(r"^[ACGTactg]{1,199}[nN]", match_func, record_seq)

        def match_func(match):
            if len(match.group()) < 200:
                return "N" * len(match.group())
            return match.group()

        record_seq = re.sub(r"[ACGTactg]+", match_func, record_seq)

        # remove first occurences of N's
        record_seq = re.sub("^[nN]+", "", record_seq)
        # remove last occurences of N's
        record_seq = re.sub("[nN]+$", "", record_seq)
        if len(record_seq) < 200:
            continue
        # find all N-runs with at least two N's
        # a single N letter will not be touched
        for match in re.finditer("[nN]{2,}", record_seq):
            seq_before_match = match.string[sequence_start : match.start()]

            linecount += 1
            sequence_counter += 1
            gap_len = match.end() - match.start()
            sequence_end += len(seq_before_match)

            data = {
                "start": sequence_start + 1,  # start relative to the scaffold
                "end": sequence_end,  # end relative to the scaffold
                "linecount": linecount,  # linecount in the scaffold
                "sequence_iid": "Seq%s" % sequence_counter,  # unique contig iid
                "scaffold_counter": scaffold_counter,  # counting the scaffold
                "end_part_sequence": len(
                    seq_before_match
                ),  # end position and len of the contig
                "gap_len": gap_len,  # Gap length
            }
            if agp_file != False:
                agp_file.write(agp_contig_string % data)
            linecount += 1
            sequence_start += len(seq_before_match)
            sequence_end += gap_len
            data.update(
                {
                    "linecount": linecount,
                    "start": sequence_start + 1,
                    "end": sequence_end,
                }
            )

            if len(record_seq) < sequence_end:
                # GAP is at the end of the scaffold, do not track that, in theory that should never happen
                continue
            # write to AGP file
            if agp_file != False:
                agp_file.write(agp_gap_string % data)
            sequence_start = sequence_end

            seq_obj = Seq(seq_before_match, IUPAC.IUPACUnambiguousDNA)
            temp = SeqRecord(
                seq_obj, id="Seq%s" % sequence_counter, description=sequence_description
            )
            assert len(temp.seq) > 199
            new_sequence_array.append(temp)

        """
            sequence after the last gap
        """
        if len(record_seq) > sequence_start:
            seq_after_last_match = record_seq[sequence_start : len(record_seq)]
            if len(seq_after_last_match) >= 200:
                linecount += 1
                sequence_end = len(record_seq)
                sequence_counter += 1
                data = {
                    "linecount": linecount,
                    "start": sequence_start + 1,
                    "end": sequence_end,
                    "sequence_iid": "Seq%s" % sequence_counter,
                    "scaffold_counter": scaffold_counter,  # counting the scaffold
                    "end_part_sequence": sequence_end - sequence_start,
                }
                if agp_file != False:
                    agp_file.write(agp_contig_string % data)
                seq_obj = Seq(seq_after_last_match, IUPAC.IUPACUnambiguousDNA)
                temp = SeqRecord(
                    seq_obj,
                    id="Seq%s" % sequence_counter,
                    name=record.name,
                    description=sequence_description,
                )
                assert len(temp.seq) > 199
                new_sequence_array.append(temp)

        scaffold_counter += 1

    output_file = open(cleand_sequence_file, "w+")
    SeqIO.write(new_sequence_array, output_file, "fasta")
    output_file.close()

    # Test if all contigs bigger than 199
    for record in SeqIO.parse(open(cleand_sequence_file), "fasta"):
        assert len(record.seq) > 199


def run_test(testdir):
    import filecmp

    for filename in os.listdir(testdir):
        path = os.path.join(testdir, filename)
        if os.path.isfile(path):
            run(
                path,
                "/tmp/scaffold2fasta_test_output.fsa",
                "/tmp/scaffold2fasta_test_output.agp",
                "AFG",
                "BBBB",
            )
            assert (
                filecmp.cmp(
                    os.path.join(testdir, "results", filename),
                    "/tmp/scaffold2fasta_test_output.fsa",
                )
                == True
            )
            assert (
                filecmp.cmp(
                    os.path.join(testdir, "results", filename.replace(".fsa", ".agp")),
                    "/tmp/scaffold2fasta_test_output.agp",
                )
                == True
            )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Splitting a scaffold file in separate contigs and create a agp file, preserving the right arragments of the contigs."
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
    parser.add_argument(
        "--agp-file",
        dest="agp_file",
        default=False,
        help="""The output path where the AGP file
                            should be stored""",
    )
    parser.add_argument(
        "--sequence-description",
        dest="seq_description",
        default=False,
        required=True,
        help="""The sequence description will be inserted in each FASTA header and will be
                      included in the genbank file. The NCBI reviever suggested something like [organism=Glarea lozoyensis 74030] [strain=74030]""",
    )

    """
        If you want to choose a different locus_tag prefix, it must be at least 3 characters.
        See http://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf.
        You can check to see if a prefix is available here:  http://www.ncbi.nlm.nih.gov/genomes/lltp.cgi
    """
    parser.add_argument(
        "--locus-tag",
        dest="locus_tag",
        default=False,
        required=True,
        help="""The locus tag each scaffold and gene gets. It should be at least 3 characters.""",
    )

    options = parser.parse_args()

    # run_test( testdir = '/usr/local/galaxy/galaxy-dist/tools/test/genome_annotation_pipeline/scaffold_tests' )

    run(
        options.input_sequence,
        options.output_sequence,
        options.agp_file,
        options.locus_tag,
        options.seq_description,
    )
