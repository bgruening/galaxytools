#! /usr/bin/python

import sys, os
from collections import defaultdict

def __main__():
    input_file = open(sys.argv[1], 'r')
    top_hits = int( sys.argv[2] )
    output_file = open(sys.argv[3], 'w')
    
    remember_ids = defaultdict(int) # defaultdict with 0
    for line in input_file:
        """
            Write comments on the fly out
        """
        if line.startswith('#'):
            output_file.write(line)
        else:
            """
            a typical gff_line from the blast2gff converter (bo_search2gff.pl):
            7180000016251	TBLASTN	match_part	26257	26601	491	.	0	Target=Sequence:tr|A5JMD6|A5JMD6_9ACTO 2 116;score=491
            We need the QueryID, since we only whant TOP X hits.
            """
            query_id = line.split('Target=Sequence:')[-1].split()[0]
            remember_ids[query_id] += 1
            if remember_ids[query_id] <= top_hits:
                output_file.write(line)

    input_file.close()
    output_file.close()
if __name__ == "__main__" : __main__()
