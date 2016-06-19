__author__ = 'gpratt'
"""

barcode_collapse.py  read in a .bam file where the
first 9 nt of the read name
are the barcode and merge reads mapped to the same position that have the same barcode

From: https://github.com/YeoLab/gscripts

"""

from collections import Counter
import itertools
from optparse import OptionParser
import sys
import pysam


def stranded_read_start(read):
    if read.is_reverse:
        return read.positions[-1]
    else:
        return read.pos


def output_metrics(metrics_file, total_count, removed_count):
    with open(metrics_file, 'w') as metrics:
        metrics.write("\t".join(["randomer", "total_count", "removed_count"]) + "\n")
        for barcode in total_count.keys():
            metrics.write("\t".join(map(str, [barcode, total_count[barcode], removed_count[barcode]])) + "\n")


def barcode_collapse(in_bam, out_bam):
    number_of_unmapped_mate_pairs = 0
    different_chroms = 0
    removed_count = Counter()
    total_count = Counter()
    result_dict = {}

    with pysam.Samfile(in_bam, 'r') as samfile1:
        with pysam.Samfile(in_bam, 'r') as samfile2:
            samfile_read1 = itertools.islice(samfile1, 0, None, 2)
            samfile_read2 = itertools.islice(samfile2, 1, None, 2)
            for read1, read2 in itertools.izip(samfile_read1, samfile_read2):
                if not read1.qname == read2.qname:
                    print read1.qname, read2.qname
                    raise Exception("Read Names don't match")
                if read1.is_unmapped and read1.is_unmapped:
                    #Both reads don't map, don't even both saving them.
                    continue
                if (not read1.is_unmapped and read2.is_unmapped) or (read1.is_unmapped and read2.is_unmapped):
                    number_of_unmapped_mate_pairs += 1
                    continue
                if read1.rname != read2.rname:
                    different_chroms += 1
                    continue
                randomer = read1.qname.split(":")[0]

                start = stranded_read_start(read1)
                stop = stranded_read_start(read2)
                strand = "-" if read1.is_reverse else "+"
                unique_location = (read1.rname, start, stop, strand, randomer)
                total_count[randomer] += 1
                if unique_location in result_dict:
                    removed_count[randomer] += 1
                    continue

                result_dict[(read1.rname, start, stop, strand, randomer)] = (read1, read2)

        with pysam.Samfile(out_bam, 'wb', template=samfile1) as out_bam:
            for key, (read1, read2) in result_dict.items():
                out_bam.write(read1)
                out_bam.write(read2)

    return total_count, removed_count

if __name__ == "__main__":
    description=""""Paired End randomer aware duplciate removal algorithm."""
    usage="""Assumes paired end reads are adjectent to each other in output file (ie only provide unsorted bams)
             Also assumes no multimappers in the bam file, if there are multimappers behavior is undefined"""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-m", "--metrics_file", dest="metrics_file")
    (options, args) = parser.parse_args()

    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % options.bam)

    total_count, removed_count = barcode_collapse(options.bam, options.out_file)
    output_metrics(options.metrics_file, total_count, removed_count)

    sys.exit(0)
