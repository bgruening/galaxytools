#!/usr/bin/env python
import sys

full_gene_model = False
if '--full' in sys.argv:
    full_gene_model = True

genome_id = None
stdin_data = []
KEY_ORDER = ('parent', 'source', 'type', 'start', 'end', 'score', 'strand',
             '8', 'quals')

# Table of amino acids
aa_table = {
    'Ala' : 'A',
    'Arg' : 'R',
    'Asn' : 'N',
    'Asp' : 'D',
    'Cys' : 'C',
    'Gln' : 'Q',
    'Glu' : 'E',
    'Gly' : 'G',
    'His' : 'H',
    'Ile' : 'I',
    'Leu' : 'L',
    'Lys' : 'K',
    'Met' : 'M',
    'Phe' : 'F',
    'Pro' : 'P',
    'Ser' : 'S',
    'Thr' : 'T',
    'Trp' : 'W',
    'Tyr' : 'Y',
    'Val' : 'V',
    'Pyl' : 'O',
    'seC' : 'U',
    '???' : 'X' }

def output_line(gff3):
    print '\t'.join(str(gff3[x]) for x in KEY_ORDER)

print '##gff-version 3'
for line in sys.stdin:
    if line.startswith('>'):
        genome_id = line[1:].strip()
        if ' ' in genome_id:
            genome_id = genome_id[0:genome_id.index(' ')]
    else:
        data = line.split()
        if len(data) == 5:
            # Parse data
            strand = '-' if data[2].startswith('c') else '+'
            start, end = data[2][data[2].index('[') + 1:-1].split(',')

            gff3 = {
                'parent': genome_id,
                'source': 'aragorn',
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'score': '.',
                '8': '.',
            }

            aa_long = data[1][5:]
            aa_short = aa_table[aa_long]
            anticodon = data[4][1:data[4].index(")")].upper().replace("T", "U")
            name = 'trn{}-{}'.format(aa_short, anticodon)

            if not full_gene_model:
                gff3.update({
                    'type': 'tRNA',
                    'quals': 'ID=tRNA{0}.{1};Name={name};product={2}'.format(genome_id, *data, name = name),
                })
                output_line(gff3)
            else:
                gff3.update({
                    'type': 'gene',
                    'quals': 'ID=gene{0}.{1};Name={name};product={2}'.format(genome_id, *data, name = name),
                })
                output_line(gff3)
                gff3.update({
                    'type': 'tRNA',
                    'quals': 'ID=tRNA{0}.{1};Parent=gene{0}.{1};Name={name};product={2}'.format(genome_id, *data, name = name),
                })
                output_line(gff3)

                # If no introns
                if ')i(' not in data[4]:
                    gff3['type'] = 'exon'
                    gff3['quals'] = 'Parent=tRNA{0}.{1}'.format(genome_id, *data)
                    output_line(gff3)
                else:
                    intron_location = data[4][data[4].rindex('(') + 1:-1].split(',')
                    intron_start, intron_length = map(int, intron_location)
                    if strand == '+':
                        original_end = gff3['end']
                    else:
                        original_end = gff3['start']

                    # EXON
                    gff3.update({
                        'type': 'exon',
                        'quals': 'Parent=tRNA{0}.{1}'.format(genome_id, *data),
                    })
                    if strand == '+':
                        gff3['end'] = gff3['start'] + intron_start - 2
                    else:
                        gff3['start'] = gff3['end'] - intron_start + 2

                    output_line(gff3)

                    # INTRON
                    gff3.update({
                        'type': 'intron',
                        'quals': 'Parent=tRNA{0}.{1}'.format(genome_id, *data),
                    })
                    if strand == '+':
                        gff3['start'] = gff3['end'] + 1
                        gff3['end'] = gff3['start'] + intron_length + 2
                    else:
                        gff3['end'] = gff3['start'] - 1
                        gff3['start'] = gff3['end'] - intron_length + 1

                    output_line(gff3)

                    # EXON
                    gff3.update({
                        'type': 'exon',
                        'quals': 'Parent=tRNA{0}.{1}'.format(genome_id, *data),
                    })
                    if strand == '+':
                        gff3['start'] = gff3['end'] + 1
                        gff3['end'] = original_end
                    else:
                        gff3['end'] = gff3['start'] - 1
                        gff3['start'] = original_end

                    output_line(gff3)
