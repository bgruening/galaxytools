#!/usr/bin/env python

from src.sequence import Sequence

class FastaReader:

    def __init__(self):
        self.seqs = []

    def read(self, io_buffer):
        header = ''
        bases = ''
        for line in io_buffer:
            if line[0] == '>':
                if len(header) > 0:
                    # Save the data
                    self.seqs.append(Sequence(header, bases))
                header = line[1:].strip().split()[0] # Get the next header
                bases = ''
            else:
                bases += line.strip()
        # Add the last sequence
        self.seqs.append(Sequence(header, bases))
        return self.seqs

