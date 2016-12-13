import sys
import os
wild_fa = sys.argv[1]
snp_tags = sys.argv[2]
outfile = sys.argv[3]
# Write RNA sequence and SNP tags to a outfile
out = open(outfile,'w')
with open (wild_fa) as in_fa_handle:
    for line in in_fa_handle:
        out.write(line)
with open (snp_tags) as in_snp_handle:
    out.write('*'+in_snp_handle.readline())
out.close()
