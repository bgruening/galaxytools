### overcoming the problem of SHAPE data working with a single line.
### creating multiple multiple files containg SHAPE data for a single sequence and running RNAfold for every
### single sequence.

import os
import sys
from os import system
from Bio import SeqIO
import re
from subprocess import Popen, PIPE

params_list = sys.argv[1:]
param_list_no_shape = [s for s in params_list if not "--shape=" in s ]
shape_file = [s for s in params_list if "--shape=" in s ]
assert (len(shape_file) == 1)

shape_file = shape_file[0]
shape_file = shape_file.replace('--shape=', '')

params_no_shape  =  " ".join(str(x) for x in param_list_no_shape)

pattern = re.compile("^>.*$")
id_line = ""
with open(shape_file, 'r') as f:
    content = f.read()
    lines = content.split('\n')
    for line in lines:
        if pattern.match(line):
            id_line = line.split()[0]
            id_line = id_line[1:]
            continue
        else:
            with open(id_line +'.tmp', "a") as clFile:
                clFile.write(line + "\n")
                
input_file = sys.stdin

for record in SeqIO.parse(input_file, "fasta"):
    seq = ">{}\n{}".format(record.id,record.seq)
    cmd =  " RNAfold  --shape=" + record.id + '.tmp ' + params_no_shape
    p = Popen(cmd , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate(seq.encode())
    if err:
        raise RuntimeError("Error in calling RNAfold\n{}\n{}\n".format(out, err))
    print (out.decode('utf-8').strip())
