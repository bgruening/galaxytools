#!/usr/bin/env python
import re

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

        if start_pattern(line) or run_of_blanklines > 2 or 'Mean G+C' in line:
            if accumulator:
                yield accumulator
                accumulator = [line]
        else:
            accumulator.append(line)
    if accumulator:
        yield accumulator

IMPORTANT_INFO = {
    'trna': re.compile(r'tRNA-(?P<codon>[A-Za-z]{3})\((?P<anticodon>[A-Za-z]{3})\)'),
    'trna-alt': re.compile(r'tRNA-\?\((?P<codon>[^\)]+)\)\((?P<anticodon>[A-Za-z]{3,})\)'),
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

IMPORTANT_INFO_TMRNA = {
    'tag_peptide': re.compile(r'Tag peptide:\s+(?P<pep>[A-Z*]*)'),
    'location': re.compile(r'Location \[(?P<start>\d+),(?P<end>\d+)\]'),
}
INFO_GROUPS_TMRNA = ('start', 'end', 'pep')

def important_info_tmrna(block):
    info = {}
    for line in block:
        for matcher in IMPORTANT_INFO_TMRNA:
            matches = IMPORTANT_INFO_TMRNA[matcher].search(line)
            if matches:
                for group in INFO_GROUPS_TMRNA:
                    try:
                        info[group] = matches.group(group)
                    except:
                        pass
    return info

import fileinput
stdin_data = []
for line in fileinput.input():
    stdin_data.append(line)

possible_blocks = [line for line in blocks(stdin_data)]

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

            # I am not proud of any of this.
            if len(data.keys()) == 0:
                data = important_info_tmrna(block)
                if len(data.keys()) == 0:
                    data = {}
                else:
                    data['type'] = 'tmrna'
            else:
                data['type'] = 'trna'
    else:
        if 'nucleotides in sequence' in block[-1]:
            try:
                fasta_defline = block[-2].strip()
            except:
                pass
        pass

    if fasta_defline is not None:
        try:
            seqid = fasta_defline[0:fasta_defline.index(' ')]
        except:
            seqid = fasta_defline

    if data is not None and len(data.keys()) > 0:

        if data['type'] == 'trna':
            notes = {
                'Codon': data['codon'],
                'Anticodon': data['anticodon'],
            }
            if 'pseudo' in data:
                notes['Note'] = 'Possible pseudogene'
        else:
            notes = {
                'Note': '"Tag peptide: ' + data['pep'] + '"'
            }

        notestr = ';'.join(['%s=%s' % (k,v) for k,v in notes.iteritems()])

        print '\t'.join([
            seqid,
            'aragorn',
            data['type'],
            data['start'],
            data['end'],
            '.',
            '.',
            '.',
            notestr
        ])
