#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import re

def add_read_pair_num(input_filepath, output_filepath, read_pair_num):
    to_add = '.' + str(read_pair_num)
    with open(input_filepath,'r') as input_file:
        with open(output_filepath,'w') as output_file:
            for line in input_file:
                if line.startswith('>'):
                    split_line = line.split()
                    seq_id = split_line[0]
                    if seq_id.rfind(to_add) != (len(seq_id)-len(to_add)):
                        split_line[0] = seq_id + to_add
                    output_file.write(' '.join(split_line) + '\n')
                else:
                    output_file.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1_input', required=True)
    parser.add_argument('--r2_input', required=True)
    parser.add_argument('--r1_output', required=True)
    parser.add_argument('--r2_output', required=True)
    args = parser.parse_args()

    add_read_pair_num(args.r1_input, args.r1_output, 1)
    add_read_pair_num(args.r2_input, args.r2_output, 2)