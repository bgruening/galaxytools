import re
import glob
import sys
from os import system

def sh(script):
    system("bash -c '%s'" % script)

cdhitcluster = sys.argv[1]

cluster_seqs_stats_path = "RESULTS/*.cluster.all"
cluster_seqs_stats_files = glob.glob(cluster_seqs_stats_path)

repSeqRedSeqdict = {}
repLine = ""
count = 0
first = False
add_FullId = ""
k = 0

with open(cdhitcluster, 'r+') as f:
    content = f.read()
    reps = re.compile("^.*\*$", re.MULTILINE).findall(content)
    lines = content.split('\n')

    for i in range(0, len(lines)):
        line = lines[i]
        if ">Cluster" in line:
            first = True
            count = 0
            repLine = reps[k]
            k = k+1
            continue
        elif not first:
            count += 1
            first = False
        else:
            first = False
            lineArr = []
        if count > 0:
            repLine = repLine.strip()
            rep_FullId = repLine.split()[2]
            rep_FullId = rep_FullId.replace(">","")
            rep_FullId = rep_FullId.replace("...","")
            if "*" in line or not line.strip():
                continue
            line = line.strip()
            add_FullId = line.split()[2]
            add_FullId = add_FullId.replace(">","")
            add_FullId = add_FullId.replace("...","")
            lineArr.append(add_FullId)
            repSeqRedSeqdict[rep_FullId] = lineArr

toWrite = ""
for singleFile in sorted(cluster_seqs_stats_files):
    toWrite = ""
    with open(singleFile, "r+") as clFile:
        file_lines = clFile.readlines()
        for line in file_lines:
            line = '\t'.join(line.split())
            toWrite += line + '\n'
        clFile.seek(0)
        clFile.write(toWrite)
        clFile.truncate()
        first_line = file_lines[0]
        toWrite = ""
        cols = first_line.split()
        file_content =  '\n'.join(file_lines)
        for key, val in repSeqRedSeqdict.items():
            if key in file_content:

                for i in val:
                    cols[3] = "---"
                    cols[4] = "CD-Hit"
                    cols[7] = str(i)
                    if len(first_line.split()) > 9:
                        cols[9] = str(i.rsplit("_",1)[0])
                    toWrite += '\t'.join(cols)
                    toWrite +="\n"
        clFile.write(toWrite)
