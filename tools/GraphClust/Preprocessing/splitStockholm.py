#!/usr/bin/env python

########
# This script reads multiple alignments merged in single stockholm file
# and splits the alignment blocks according to data.names table
# The first sequence of each alignment file assumed to match to names table entries
# Author: M. Miladi
########
import os
import re
import sys

from Bio import AlignIO, SeqIO
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

stk_file = sys.argv[1]
print ("Parsing and splitting stk file:{}".format(stk_file))
target_f = "alignment_data_split.stk"
pattern = re.compile("^>.*$")
toWriteID = ""

count_for_id = 1
seq_counter = 0
new_id = ""

seq_id = []
seq_string = []
orig_id = []
name_file = "FASTA/data.names"
array_all_chunks = []
with open(name_file, 'r') as f:
    for line in f:
        if len(line.strip()) == 0:
            continue
        seq_id.append(int(line.split()[0]))
        seq_string.append(line.split()[1])
        orig_id_srt = line.split()[3]
        orig_id_srt = orig_id_srt.rsplit('_',1)[0]
        orig_id.append(orig_id_srt)



with open(stk_file) as stk_in:
    alignments = AlignIO.parse(stk_in, "stockholm")#, alphabet=IUPAC.ambiguous_rna)  
    alignments_dic = {(a[0].id):a for a in alignments}


regx_gaps = '[-.~_]'  # valid gap symbols all be converted to "-"
str_gaps = '-.~_'  # valid gap symbols all be converted to "-"


chunks = []
with open(target_f, 'w') as out_stk_handle:
    for i in range(len(orig_id)):
        
        #----------------------
        # We need to map ungapped positions of the chunks to gapped positions of first sequence 
        gap_count = 0
        ungap_ind = 0
        dic_gap_counts = dict()
        cur_alignment = alignments_dic[orig_id[i]]
        for c in cur_alignment[0].seq:
            #print ungap_ind
            if c in str_gaps:
                gap_count += 1
            else:
                dic_gap_counts[ungap_ind] = gap_count
                ungap_ind += 1
        ID =  str(seq_id[i]) + " " + seq_string[i] 
        chunks = re.findall(r'\d+', seq_string[i])
        print (ID,chunks)

        index_start, index_end =int(chunks[1])-1, int(chunks[2])-1
        subalign = cur_alignment[:, index_start + dic_gap_counts[index_start]:
                           index_end+dic_gap_counts[index_end]+1]
        
        #----------------------
        # BioPython does not handel the GF ID entry for alignment
        # So we add entry in the second line manually
        siotmp = StringIO()
        AlignIO.write(subalign, siotmp, format="stockholm")
        stk_lines = siotmp.getvalue().split('\n')
        out_stk_handle.writelines('{}\n'.format(stk_lines[:1]))
        out_stk_handle.write('#=GF ID {}\n'.format(ID))
        out_stk_handle.writelines('\n'.join(stk_lines[1:]))
        #print out_stk_handle.getvalue()

        
