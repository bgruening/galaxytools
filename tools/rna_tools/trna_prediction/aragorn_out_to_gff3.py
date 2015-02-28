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
    'trna-alt': re.compile(r'tRNA-\?\((?P<codon>[^\)]+)\)\((?P<anticodon>[A-Za-z]{2,})\)'),
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
    'location': re.compile(r'Location (?P<complement>[c]{0,1})\[(?P<start>\d+),(?P<end>\d+)\]'),
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
# We're off to a GREAT start, if I'm accessing by index you just know that I'm going to do terrible
# awful things
for block_idx in range(len(possible_blocks)):
    block = possible_blocks[block_idx]
    data = None
    fasta_defline = None

    if block[0].startswith('Searching for') or 'nucleotides in sequence' in block[-1]:
        # Try and get a sequence ID out of it
        try:
            fasta_defline = block[-2].strip()
        except:
            # Failing that, ignore it.
            pass
    else:
        # They DUPLICATE results in multiple places, including a fasta version
        # in the 'full report'.
        possible_ugliness = [x for x in block if x.startswith('>t')]
        if len(possible_ugliness) > 0:
            continue

        # However, if it didn't have one of those all important pieces of
        # information, then it's either a different important piece of
        # information, or complete junk
        data = important_info(block)

        # I am not proud of any of this. We essentially say "if that block
        # didn't come up with useful info, then try making it a tmrna"
        if len(data.keys()) == 0:
            data = important_info_tmrna(block)
            # And if that fails, just none it.
            if len(data.keys()) == 0:
                data = None
            else:
                # But if it didn't, confirm that we're a tmRNA
                data['type'] = 'tmRNA'
        else:
            # If we did have keys, and didn't pass through any of the tmRNA
            # checks, we're tRNA
            data['type'] = 'tRNA'

    # If we got a sequence ID in this block, set the defline
    if 'nucleotides in sequence' in block[-1]:
        try:
            fasta_defline = block[-2].strip()
        except:
            pass

    # if a defline is available, try and extract the fasta header ID
    if fasta_defline is not None:
        try:
            seqid = fasta_defline[0:fasta_defline.index(' ')]
        except:
            seqid = fasta_defline

    # If there's data
    if data is not None and len(data.keys()) > 1:

        # Deal with our flags/notes.
        if data['type'] == 'tRNA':
            # Are these acceptable GFF3 tags?
            notes = {
                'Codon': data['codon'],
                'Anticodon': data['anticodon'],
            }
            if 'pseudo' in data:
                notes['Note'] = 'Possible pseudogene'
        else:
            notes = {
                'Note': 'Tag peptide: ' + data['pep'] + ''
            }

        notestr = ';'.join(['%s="%s"' % (k,v) for k,v in notes.iteritems()])

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
