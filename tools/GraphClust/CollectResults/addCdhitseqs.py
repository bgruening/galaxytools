import re
import glob
import sys

cdhitcluster = sys.argv[1]
#clusters = sys.argv[2]

cluster_seqs_stats_path = "RESULTS/*.cluster.all"
cluster_seqs_stats_files = glob.glob(cluster_seqs_stats_path)

#clusterFiles = clusters.split(',')
repSeqRedSeqdict = {}
repLine = ""
count = 0
first = False

with open(cdhitcluster, 'r+') as f:
    lines = f.readlines()
    for i in range(0, len(lines)):
        line = lines[i]
        if ">Cluster" in line:
            first = True
            count = 0
            if i+1 < len(lines):
                repLine = lines[i+1]
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
            rep_FullId = rep_FullId.replace(">", "")
            #rep_short_id = re.findall("_".join(["[^_]+"] * 2), rep_FullId)[0]
            rep_FullId = rep_FullId.replace("...", "")
            line = line.strip()
            add_FullId = line.split()[2]
            add_FullId = add_FullId.replace(">", "")
            add_FullId = add_FullId.replace("...", "")
            #add_short_id = re.findall("_".join(["[^_]+"] * 2), add_FullId)[0]
            lineArr.append(add_FullId)
            repSeqRedSeqdict[rep_FullId] = lineArr
            #lineArr.append(add_short_id)
            #repSeqRedSeqdict[rep_short_id] = lineArr

toWrite = ""

for singleFile in sorted(cluster_seqs_stats_files):
    with open(singleFile, "a+") as clFile:
        file_content = clFile.read()
        first_line = file_content.split('\n')[0]
        for key, val in repSeqRedSeqdict.items():
            if key in file_content:
                for i in val:
                    toWrite += first_line.split()[0] + "  " + first_line.split()[1] + "  " + first_line.split()[2] + "  " + " - " + "   " + "CD-Hit" + "    " + first_line.split()[5] + "  " + "ORIGID" + "  "  + str(i) + "\n"
        clFile.write(toWrite)
