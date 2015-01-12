#!/usr/bin/env python
import re
import sys

def start_pattern(string):
    return re.match(r'^[0-9]+\.$', string) \
        or string.startswith('Number of possible') \
        or string.startswith('Searching for')

def blank_line(string):
    return re.match(r'^\s*$', string)

def blocks(iterable):
    accumulator = []
    run_of_blanklines = 0
    for line in iterable:
        # Count blank lines
        if blank_line(line):
            run_of_blanklines += 1
        else:
            run_of_blanklines = 0

        if start_pattern(line) or run_of_blanklines > 3 or 'Mean G+C' in line:
            if accumulator:
                yield accumulator
                accumulator = [line]
        else:
            accumulator.append(line)
    if accumulator:
        yield accumulator

IMPORTANT_INFO = {
    'trna': re.compile(r'tRNA-(?P<codon>[A-Za-z]{3})\((?P<anticodon>[A-Za-z]{3})\)'),
    'trna-pseudo': re.compile(r'tRNA-\?\((?P<codon>[^\)]+)\)\((?P<anticodon>[A-Za-z]{3,})\)'),
    'bases': re.compile(r'(?P<bases>[0-9]+) bases, %GC = (?P<gc>[0-9.]+)'),
    'sequence': re.compile(r'Sequence (?P<complement>[c]{0,1})\[(?P<start>\d+),(?P<end>\d+)\]'),
    'possible_pseudogene': re.compile(r'(?P<pseudo>Possible Pseudogene)'),
}
INFO_GROUPS = ('codon', 'anticodon', 'bases', 'gc', 'complement', 'start', 'end', 'pseudo')

def important_info(block):
    info = {}
    for line in block:
        for matcher in IMPORTANT_INFO:
            matches = IMPORTANT_INFO[matcher].search(line)
            if matches:
                for group in INFO_GROUPS:
                    try:
                        info[group] = matches.group(group)
                    except:
                        pass
    return info


input_raw = open(sys.argv[1],'r').readlines()
possible_blocks = [line for line in blocks(input_raw)]

seqid = None
print '##gff-version-3'
for block in possible_blocks:
    data = None
    fasta_defline = None
    if start_pattern(block[0]):
        if block[0].startswith('Searching for') or 'nucleotides in sequence' in block[-1]:
            try:
                fasta_defline = block[-2].strip()
            except:
                pass
        else:
            data = important_info(block)
    else:
        if 'nucleotides in sequence' in block[-1]:
            try:
                fasta_defline = block[-2].strip()
            except:
                pass
        pass

    if fasta_defline is not None:
        seqid = fasta_defline[0:fasta_defline.index(' ')]

    if data is not None:
        notes = {
            'Codon': data['codon'],
            'Anticodon': data['anticodon'],
        }
        if 'pseudo' in data:
            notes['Note'] = 'Possible pseudogene'

        notestr = ';'.join(['%s=%s' % (k,v) for k,v in notes.iteritems()])

        print '\t'.join([
            seqid,
            'aragorn',
            'tRNA',  #TODO: tmRNA?
            data['start'],
            data['end'],
            '.',
            '.',
            '.',
            notestr
        ])
