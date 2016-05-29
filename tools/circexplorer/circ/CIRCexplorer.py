#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CIRCexplorer.py 1.1.7 -- circular RNA analysis toolkit.

Usage: CIRCexplorer.py [options]

Options:
    -h --help                      Show this screen.
    --version                      Show version.
    -f FUSION --fusion=FUSION      TopHat-Fusion fusion BAM file. (used in \
TopHat-Fusion mapping)
    -j JUNC --junc=JUNC            STAR Chimeric junction file. (used in STAR \
mapping)
    -g GENOME --genome=GENOME      Genome FASTA file.
    -r REF --ref=REF               Gene annotation.
    -o PREFIX --output=PREFIX      Output prefix [default: CIRCexplorer].
    --tmp                          Keep temporary files.
    --no-fix                       No-fix mode (useful for species \
with poor gene annotations)
"""

__author__ = 'Xiao-Ou Zhang (zhangxiaoou@picb.ac.cn)'
__version__ = '1.1.7'

from docopt import docopt
import sys
import pysam
from collections import defaultdict
from genomic_interval import Interval
import tempfile
import os
from pkg_resources import parse_version as norm_vr


def convert_fusion(fusion_bam, output_f):
    """
    Extract fusion junction reads from the BAM file
    """
    print('Start to convert fustion reads...')
    fusions = defaultdict(int)
    for i, read in enumerate(parse_bam(fusion_bam)):
        chrom, strand, start, end = read
        segments = [start, end]
        if (i + 1) % 2 == 1:  # first fragment of the fusion junction read
            interval = [start, end]
        else:  # second fragment of the fusion junction read
            sta1, end1 = interval
            sta2, end2 = segments
            if end1 < sta2 or end2 < sta1:  # no overlap between fragments
                sta = sta1 if sta1 < sta2 else sta2
                end = end1 if end1 > end2 else end2
                fusions['%s\t%d\t%d' % (chrom, sta, end)] += 1
    total = 0
    with open(output_f, 'w') as outf:
        for i, pos in enumerate(fusions):
            outf.write('%s\tFUSIONJUNC_%d/%d\t0\t+\n' % (pos, i, fusions[pos]))
            total += fusions[pos]
    print('Converted %d fusion reads!' % total)


def annotate_fusion(ref_f, input_f, output_f):
    """
    Align fusion juncrions to gene annotations
    """
    print('Start to annotate fusion junctions...')
    genes, gene_info = parse_ref1(ref_f)  # gene annotations
    fusions, fusion_index = parse_bed(input_f)  # fusion junctions
    total = set()
    with open(output_f, 'w') as outf:
        for chrom in genes:
            # overlap gene annotations with fusion juncrions
            result = Interval.overlapwith(genes[chrom].interval,
                                          fusions[chrom])
            for itl in result:
                # extract gene annotations
                iso = list(filter(lambda x: x.startswith('iso'), itl[2:]))
                # for each overlapped fusion junction
                for fus in itl[(2 + len(iso)):]:
                    reads = fus.split()[1]
                    fus_start, fus_end = fusion_index[fus]
                    edge_annotations = []  # first or last exon flag
                    for iso_id in iso:
                        g, i, c, s = iso_id.split()[1:]
                        start = gene_info[iso_id][0][0]
                        end = gene_info[iso_id][-1][-1]
                        # fusion junction excesses boundaries of gene
                        # annotation
                        if fus_start < start - 10 or fus_end > end + 10:
                            continue
                        (fusion_info,
                         index,
                         edge) = map_fusion_to_iso(fus_start,
                                                   fus_end, s,
                                                   gene_info[iso_id])
                        if fusion_info:
                            fus_start_str = str(fus_start)
                            fus_end_str = str(fus_end)
                            bed_info = '\t'.join([chrom, fus_start_str,
                                                  fus_end_str,
                                                  'FUSIONJUNC/%s' % reads,
                                                  '0', s, fus_start_str,
                                                  fus_start_str, '0,0,0'])
                            bed = '\t'.join([bed_info, fusion_info, g, i,
                                             index])
                            if not edge:  # not first or last exon
                                outf.write(bed + '\n')
                                total.add(fus)
                            else:  # first or last exon
                                edge_annotations.append(bed)
                    if edge_annotations:  # first or last exon
                        for bed in edge_annotations:
                            outf.write(bed + '\n')
                        total.add(fus)
    print('Annotated %d fusion junctions!' % len(total))


def fix_fusion(ref_f, genome_fa, input_f, output_f, no_fix):
    """
    Realign fusion juncrions
    """
    print('Start to fix fusion junctions...')
    fa = genome_fa
    ref = parse_ref2(ref_f)
    fusions, fusion_names, fixed_flag = fix_bed(input_f, ref, fa, no_fix)
    total = 0
    with open(output_f, 'w') as outf:
        for fus in fusion_names:
            reads = str(fusions[fus])
            fixed = fixed_flag[fus]
            if fixed > 0:
                total += 1
            fixed = str(fixed)
            name = 'circular_RNA/' + reads
            gene, iso, chrom, strand, index = fus.split()
            starts, ends = ref['\t'.join([gene, iso, chrom, strand])]
            exon_num = len(starts)
            intron_num = exon_num - 1
            if ',' in index:  # back spliced exons
                s, e = [int(x) for x in index.split(',')]
                if strand == '+':
                    index_info = ','.join(str(x + 1)
                                          for x in xrange(s, e + 1))
                else:
                    index_info = ','.join(str(exon_num - x)
                                          for x in xrange(s, e + 1))
                start = str(starts[s])
                end = str(ends[e])
                length = str(e - s + 1)
                sizes, offsets = generate_bed(int(start), starts[s:(e + 1)],
                                              ends[s:(e + 1)])
                if s == 0:
                    left_intron = 'None'
                else:
                    left_intron = '%s:%d-%d' % (chrom, ends[s - 1], starts[s])
                if e == len(ends) - 1:
                    right_intron = 'None'
                else:
                    right_intron = '%s:%d-%d' % (chrom, ends[e], starts[e + 1])
                intron = '|'.join([left_intron, right_intron])
                bed = '\t'.join([chrom, start, end, name, fixed, strand, start,
                                 start, '0,0,0', length, sizes, offsets,
                                 reads, 'circRNA', gene, iso, index_info,
                                 intron])
            else:  # ciRNAs
                index, start, end = index.split('|')
                index = int(index)
                if strand == '+':
                    index_info = str(index + 1)
                else:
                    index_info = str(intron_num - index)
                size = str(int(end) - int(start))
                intron = '%s:%d-%d' % (chrom, ends[index], starts[index + 1])
                bed = '\t'.join([chrom, start, end, name, fixed, strand, start,
                                 start, '0,0,0', '1', size, '0',
                                 reads, 'ciRNA', gene, iso, index_info,
                                 intron])
            outf.write(bed + '\n')
    print('Fixed %d fusion junctions!' % total)


def parse_bam(bam):
    fusions = {}
    for read in bam:
        if read.is_secondary:  # not the primary alignment
            continue
        if not read.has_tag('XF'):  # not fusion junctions
            continue
        chr1, chr2 = read.get_tag('XF').split()[1].split('-')
        if chr1 != chr2:  # not on the same chromosome
            continue
        strand = '+' if not read.is_reverse else '-'
        if read.query_name not in fusions:  # first fragment
            fusions[read.query_name] = [chr1, strand, read.reference_start,
                                        read.reference_end]
        else:  # second fragment
            if chr1 == fusions[read.query_name][0] \
               and strand == fusions[read.query_name][1]:
                yield [chr1, strand, read.reference_start, read.reference_end]
                yield fusions[read.query_name]


def parse_ref1(ref_file):
    genes = defaultdict(list)
    gene_info = {}
    with open(ref_file, 'r') as f:
        for line in f:
            gene_id, iso_id, chrom, strand = line.split()[:4]
            total_id = '\t'.join(['iso', gene_id, iso_id, chrom, strand])
            starts = [int(x) for x in line.split()[9].split(',')[:-1]]
            ends = [int(x) for x in line.split()[10].split(',')[:-1]]
            start = starts[0]
            end = ends[-1]
            genes[chrom].append([start, end, total_id])
            gene_info[total_id] = [starts, ends]
    for chrom in genes:
        genes[chrom] = Interval(genes[chrom])
    return (genes, gene_info)


def parse_ref2(ref_file):
    genes = {}
    with open(ref_file, 'r') as f:
        for line in f:
            gene_id, iso_id, chrom, strand = line.split()[:4]
            starts = [int(x) for x in line.split()[9].split(',')[:-1]]
            ends = [int(x) for x in line.split()[10].split(',')[:-1]]
            genes['\t'.join([gene_id, iso_id, chrom, strand])] = [starts, ends]
    return genes


def parse_bed(fus):
    fusions = defaultdict(list)
    fusion_index = {}
    with open(fus, 'r') as f:
        for line in f:
            chrom, start, end, name = line.split()[:4]
            start = int(start)
            end = int(end)
            reads = name.split('/')[1]
            fusion_id = '%s\t%s' % (name, reads)
            fusions[chrom].append([start, end, fusion_id])
            fusion_index[fusion_id] = [start, end]
    return (fusions, fusion_index)


def map_fusion_to_iso(start, end, strand, iso_info):
    starts = iso_info[0]
    ends = iso_info[1]
    # check sequnence within +/-10bp
    start_points = list(range(start - 10, start + 11))
    end_points = list(range(end - 10, end + 11))
    start_index, end_index = None, None
    start_intron_flag, end_intron_flag = False, False
    # check starts
    for i, s in enumerate(starts):
        if s in start_points:
            start_index = i
            break
    else:
        for j, e in enumerate(ends):
            if e in start_points:
                if j != len(ends) - 1:
                    start_index = j
                    start_intron_flag = True
                    break
    # check ends
    for j, e in enumerate(ends):
        if e in end_points:
            end_index = j
            break
    else:
        for i, s in enumerate(starts):
            if s in end_points:
                if i != 0:
                    end_index = i - 1
                    end_intron_flag = True
                    break
    # ciRNAs
    if start_intron_flag and strand == '+' and end < starts[start_index + 1]:
        return ('\t'.join(['1', str(end - start), '0', 'ciRNA']),
                str(start_index), False)
    elif end_intron_flag and strand == '-' and start > ends[end_index]:
        return ('\t'.join(['1', str(end - start), '0', 'ciRNA']),
                str(end_index), False)
    # back spliced exons
    elif (start_index is not None and end_index is not None and
          not start_intron_flag and not end_intron_flag):
        if start_index != 0 and end_index != len(ends) - 1:
            edge_flag = False
        else:
            edge_flag = True
        return convert_to_bed(start, end,
                              starts[start_index:(end_index + 1)],
                              ends[start_index:(end_index + 1)],
                              str(start_index), str(end_index), edge_flag)
    else:
        return(None, None, False)


def convert_to_bed(start, end, starts, ends, start_index, end_index, edge):
    new_starts = [start] + starts[1:]
    new_ends = ends[:-1] + [end]
    block_starts, block_sizes = [], []
    for s, e in zip(new_starts, new_ends):
        block_starts.append(str(s - start))
        block_sizes.append(str(e - s))
    length = len(block_sizes)
    block_starts = ','.join(block_starts)
    block_sizes = ','.join(block_sizes)
    index = ','.join([start_index, end_index])
    return ('\t'.join([str(length), block_sizes, block_starts, 'circRNA']),
            index, edge)


def fix_bed(fusion_file, ref, fa, no_fix):
    fusions = defaultdict(int)
    # make sure order of fusion names according to fusion_file
    fusion_names = []
    fusion_set = set()
    fixed_flag = defaultdict(int)  # flag to indicate realignment
    junctions = set()
    with open(fusion_file, 'r') as f:
        for line in f:
            chrom = line.split()[0]
            strand = line.split()[5]
            start, end = [int(x) for x in line.split()[1:3]]
            junction_info = '%s\t%d\t%d' % (chrom, start, end)
            if junction_info in junctions:
                continue
            reads = int(line.split()[3].split('/')[1])
            flag, gene, iso, index = line.split()[-4:]
            flag = True if flag == 'ciRNA' else False
            name = '\t'.join([gene, iso, chrom, strand, index])
            iso_starts, iso_ends = ref['\t'.join([gene, iso, chrom, strand])]
            if not flag:  # back spliced exons
                s, e = [int(x) for x in index.split(',')]
                # not realign
                if start == iso_starts[s] and end == iso_ends[e]:
                    fusions[name] += reads
                    if name not in fusion_set:
                        fusion_set.add(name)
                        fusion_names.append(name)
                    junctions.add(junction_info)
                # no fix mode
                elif no_fix:
                    fusions[name] += reads
                    if name not in fusion_set:
                        fusion_set.add(name)
                        fusion_names.append(name)
                    fixed_flag[name] += 1
                    junctions.add(junction_info)
                # realign
                elif check_seq(chrom, [start, iso_starts[s], end, iso_ends[e]],
                               fa):
                    fusions[name] += reads
                    if name not in fusion_set:
                        fusion_set.add(name)
                        fusion_names.append(name)
                    fixed_flag[name] += 1
                    junctions.add(junction_info)
            else:  # ciRNAs
                index = int(index)
                if strand == '+':
                    # not realign
                    if start == iso_ends[index]:
                        name += '|'.join(['', str(start), str(end)])
                        fusions[name] += reads
                        if name not in fusion_set:
                            fusion_set.add(name)
                            fusion_names.append(name)
                        junctions.add(junction_info)
                    # realign
                    elif check_seq(chrom, [start, iso_ends[index], end], fa,
                                   intron_flag=True):
                        fixed_start = iso_ends[index]
                        fixed_end = end + fixed_start - start
                        name += '|'.join(['', str(fixed_start),
                                          str(fixed_end)])
                        fusions[name] += reads
                        if name not in fusion_set:
                            fusion_set.add(name)
                            fusion_names.append(name)
                        fixed_flag[name] += 1
                        junctions.add(junction_info)
                else:
                    if end == iso_starts[index + 1]:
                        # not realign
                        name += '|'.join(['', str(start), str(end)])
                        fusions[name] += reads
                        if name not in fusion_set:
                            fusion_set.add(name)
                            fusion_names.append(name)
                        junctions.add(junction_info)
                        # realign
                    elif check_seq(chrom, [end, iso_starts[index + 1], start],
                                   fa, intron_flag=True):
                        fixed_end = iso_starts[index + 1]
                        fixed_start = start + fixed_end - end
                        name += '|'.join(['', str(fixed_start),
                                          str(fixed_end)])
                        fusions[name] += reads
                        if name not in fusion_set:
                            fusion_set.add(name)
                            fusion_names.append(name)
                        fixed_flag[name] += 1
                        junctions.add(junction_info)
    return (fusions, fusion_names, fixed_flag)


def check_seq(chrom, pos, fa, intron_flag=False):
    if not intron_flag:  # back spliced exons
        if pos[0] - pos[1] != pos[2] - pos[3]:
            return False
        if pos[0] < pos[1]:
            seq1 = fa.fetch(chrom, pos[0], pos[1])
            seq2 = fa.fetch(chrom, pos[2], pos[3])
        else:
            seq1 = fa.fetch(chrom, pos[1], pos[0])
            seq2 = fa.fetch(chrom, pos[3], pos[2])
    else:  # ciRNAs
        if abs(pos[0] - pos[1]) <= 5:  # permit mismatches within 5bp
            return True
        elif pos[0] < pos[1]:
            seq1 = fa.fetch(chrom, pos[0], pos[1])
            seq2 = fa.fetch(chrom, pos[2], pos[2] + pos[1] - pos[0])
        else:
            seq1 = fa.fetch(chrom, pos[1], pos[0])
            seq2 = fa.fetch(chrom, pos[2] - pos[0] + pos[1], pos[2])
    if seq1 == seq2:
        return True
    else:
        return False


def generate_bed(start, starts, ends):
    sizes, offsets = [], []
    for s, e in zip(starts, ends):
        sizes.append(str(e - s))
        offsets.append(str(s - start))
    sizes = ','.join(sizes)
    offsets = ','.join(offsets)
    return (sizes, offsets)


def create_temp(tmp_flag, prefix, flag=1):
    if tmp_flag:
        temp_dir = os.getcwd()
        if flag:
            temp1 = temp_dir + '/%s_fusion_junction_info.txt' % prefix
        temp2 = temp_dir + '/%s_annotated_junction_info.txt' % prefix
    else:
        temp_dir = tempfile.mkdtemp()
        if flag:
            temp1 = temp_dir + '/tmp1'
        temp2 = temp_dir + '/tmp2'
    if flag:
        return (temp_dir, temp1, temp2)
    else:
        return(temp_dir, temp2)


def delete_temp(temp_dir, temp1, temp2, flag=1):
    if flag:
        os.remove(temp1)
    os.remove(temp2)
    os.rmdir(temp_dir)


def main():
    pysam_vr = norm_vr(pysam.__version__)
    if pysam_vr < norm_vr('0.8.2'):
        sys.exit('The version of pysam is too low. It should be >= 0.8.2.')
    if len(sys.argv) == 1:
        sys.exit(__doc__)
    options = docopt(__doc__, version=__version__)
    parameters = ('--genome', '--ref')
    miss_parameters = []
    for arg in parameters:
        if not options[arg]:
            miss_parameters.append(arg)
    if miss_parameters:
        sys.exit('Miss required option: ' + ' '.join(miss_parameters))
    if options['--fusion'] and not options['--junc']:
        try:
            fusion_bam = pysam.AlignmentFile(options['--fusion'], 'rb')
        except:
            sys.exit('Please make sure %s is BAM file!' % options['--fusion'])
    elif not options['--junc'] and not options['--fusion']:
        sys.exit('--fusion or --junc should be used!')
    elif options['--junc'] and options['--fusion']:
        sys.exit('Could not use --fusion and --junc simultaneously!')

    if not os.path.isfile(options['--genome'] + '.fai'):
        pysam.faidx(options['--genome'])
    genome_fa = pysam.FastaFile(options['--genome'])
    ref_f = options['--ref']
    output_prefix = options['--output']
    output = output_prefix + '_circ.txt'
    print('Start CIRCexplorer %s' % __version__)
    if options['--junc']:
        temp1 = options['--junc']
        temp_dir, temp2 = create_temp(options['--tmp'], output_prefix, 0)
    else:
        temp_dir, temp1, temp2 = create_temp(options['--tmp'], output_prefix)
        convert_fusion(fusion_bam, temp1)
    annotate_fusion(ref_f, temp1, temp2)
    fix_fusion(ref_f, genome_fa, temp2, output, options['--no-fix'])
    if not options['--tmp']:
        if options['--junc']:
            delete_temp(temp_dir, temp1, temp2, 0)
        else:
            delete_temp(temp_dir, temp1, temp2)


if __name__ == '__main__':
    main()
