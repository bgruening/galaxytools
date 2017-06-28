#!/usr/bin/env python2
import glob
from os import system
import re
from sklearn import metrics
from shutil import make_archive

def sh(script):
    system("bash -c '%s'" % script)

dataNames = "FASTA/data.names"

listOfClusters = []
listOfHeaders = []
headersNames = set()
cluster_seqs_stats_path = "RESULTS/*.cluster.all"
cluster_seqs_stats_files = glob.glob(cluster_seqs_stats_path)

with open(dataNames, "r") as names:
    for line2 in names:
        splits2 = line2.split()
        fullHeader = ''
        if len(splits2) >= 6:
            fullHeader = splits2[5]
            headersNames.add(fullHeader)

blackList = []
numberOfClusters = 0
for singleFile in sorted(cluster_seqs_stats_files):
    numberOfClusters += 1
    with open(singleFile, "r") as f:
        for line in f:
            splits = line.split()
            header = ''
            if len(splits) >= 11:
                header = splits[10]
            clustNum = splits[2]
            listOfHeaders.append(header)
            listOfClusters.append(clustNum)
            if header in headersNames:
                blackList.append(header)

numberOfClusters += 1  # 1 cluster for all unassigned seqs
with open(dataNames, "r") as names:
    for line in names.readlines():
        splits = line.split() 
        fullUniqeId = splits[3]
        fullHeader = ''
        if len(splits) >= 6:
            fullHeader = line.split()[5]
        if fullHeader not in blackList or len(fullHeader) == 0:
            listOfHeaders.append(fullHeader)
            listOfClusters.append(str(numberOfClusters))
            numberOfClusters += 1  # separate cluster for all unassigned seqs

toWrite = ""
for i in range(len(listOfClusters)):
    toWrite += listOfHeaders[i] + "\t" + listOfClusters[i] + '\n'
with open("RESULTS/fullTab.tabular", "w") as full:
    full.write(toWrite)


pattern = re.compile("^RF.*$")

if len(listOfHeaders) > 1: # and  pattern.match(str(listOfHeaders[0])):

    completeness_score = metrics.completeness_score(listOfHeaders, listOfClusters)
    homogeneity_score = metrics.homogeneity_score(listOfHeaders, listOfClusters)
    adjusted_rand_score = metrics.adjusted_rand_score(listOfHeaders, listOfClusters)
    adjusted_mutual_info_score = metrics.adjusted_mutual_info_score(listOfHeaders, listOfClusters)
    v_measure_score = metrics.v_measure_score(listOfHeaders, listOfClusters)

    toWrite = "completeness_score : " + str(completeness_score) + "\n" + "homogeneity_score : " + str(homogeneity_score) + "\n" + "adjusted_rand_score : " +str(adjusted_rand_score)  + "\n" + "adjusted_mutual_info_score : " + str(adjusted_mutual_info_score)+ "\n" + "v_measure_score : " + str(v_measure_score)

else:
    toWrite = "completeness_score : NA \nhomogeneity_score : NA \nadjusted_rand_score : NA \nadjusted_mutual_info_score : NA \nv_measure_score : NA"

with open("RESULTS/evaluation.txt", "w") as fOut:
    fOut.write(toWrite)


make_archive('RESULTS', 'zip', root_dir='RESULTS')
