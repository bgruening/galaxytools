#!/usr/bin/env python
import fnmatch
import os
import re
import sys
from os import system
from shutil import make_archive

from sklearn import metrics


def sh(script):
    system("bash -c '%s'" % script)


fasta_dir = sys.argv[1]
results_dir = sys.argv[2]
dataNames = os.path.join(fasta_dir, "data.names")

listOfClusters = []
listOfHeaders = []
headersNames = set()
idsNames = set()


names = os.listdir(results_dir)
cluster_seqs_stats_files = fnmatch.filter(names, "*.cluster.all")
with open(dataNames, "r") as names:
    for line2 in names:
        splits2 = line2.split()
        fullHeader = ""
        if len(splits2) >= 6:
            fullHeader = splits2[5]
            headersNames.add(fullHeader)
            fullID = splits2[3]
            idsNames.add(fullID)

blackList = []
numberOfClusters = 0
for singleFile in sorted(cluster_seqs_stats_files):
    singleFile = os.path.join(results_dir, singleFile)
    numberOfClusters += 1
    with open(singleFile, "r") as f:
        for line in f:
            splits = line.split()
            header = ""
            idd = ""
            if len(splits) >= 11:
                header = splits[10]
                idd = splits[8]
            clustNum = splits[2]
            listOfHeaders.append(header)
            listOfClusters.append(clustNum)
            if idd in idsNames:  # header in headersNames:
                blackList.append(idd)

numberOfClusters += 1  # 1 cluster for all unassigned seqs
ignoreBlackList = False
with open(dataNames, "r") as names:
    for line in names:
        splits = line.split()
        fullUniqeId = splits[3]
        fullHeader = ""
        fullID = ""
        if len(splits) >= 6:
            fullHeader = line.split()[5]
            fullID = line.split()[3]
        if ignoreBlackList or (
            fullID not in blackList  # fullHeader not in blackList
            or len(fullHeader) == 0
        ):
            listOfHeaders.append(fullHeader)
            listOfClusters.append(str(numberOfClusters))
            numberOfClusters += 1  # separate cluster for all unassigned seqs
        # else:
        #     print ("Skip header", fullHeader)

toWrite = ""
for i in range(len(listOfClusters)):
    toWrite += "%s\t%s\n" % (listOfHeaders[i], listOfClusters[i])

with open(os.path.join(results_dir, "fullTab.tabular"), "w") as full:
    full.write(toWrite)


pattern = re.compile("^RF.*$")

if len(listOfHeaders) > 1:  # and  pattern.match(str(listOfHeaders[0])):

    completeness_score = metrics.completeness_score(listOfHeaders, listOfClusters)
    homogeneity_score = metrics.homogeneity_score(listOfHeaders, listOfClusters)
    adjusted_rand_score = metrics.adjusted_rand_score(listOfHeaders, listOfClusters)
    adjusted_mutual_info_score = metrics.adjusted_mutual_info_score(
        listOfHeaders, listOfClusters
    )
    v_measure_score = metrics.v_measure_score(listOfHeaders, listOfClusters)

    toWrite = "completeness_score : {}\n".format(completeness_score)
    toWrite += "homogeneity_score : {}\n".format(homogeneity_score)
    toWrite += "adjusted_rand_score : {}\n".format(adjusted_rand_score)
    toWrite += "adjusted_mutual_info_score : {}\n".format(adjusted_mutual_info_score)
    toWrite += "v_measure_score : {}\n".format(v_measure_score)


else:
    toWrite = "completeness_score : NA \nhomogeneity_score : NA \nadjusted_rand_score : NA \nadjusted_mutual_info_score : NA \nv_measure_score : NA"

with open(os.path.join(results_dir, "evaluation.txt"), "w") as fOut:
    fOut.write(toWrite)


make_archive("RESULTS", "zip", root_dir=results_dir)
