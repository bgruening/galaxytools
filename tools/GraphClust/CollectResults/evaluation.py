import glob
from os import system
import re


def sh(script):
    system("bash -c '%s'" % script)


dataNames = "FASTA/data.names"

listOfClusters = []
listOfClasses = []
cluster_seqs_stats_path = "RESULTS/*.cluster.all"
cluster_seqs_stats_files = glob.glob(cluster_seqs_stats_path)

blackList = []
numberOfClusters = 0
for singleFile in sorted(cluster_seqs_stats_files):
    numberOfClusters += 1
    with open(singleFile, "r") as f:
        for line in f.readlines():
            uniqueId = line.split()[7]
            clustNum = line.split()[1]
            rnaClass, sep, tail = uniqueId.partition("_")
            listOfClasses.append(rnaClass)
            listOfClusters.append(clustNum)
            with open(dataNames, "r") as names:
                for line in names.readlines():
                    fullUniqeId = line.split()[3]
                    rnaClass, sep, tail = fullUniqeId.partition("_")
                    if fullUniqeId == uniqueId:
                        blackList.append(uniqueId)

numberOfClusters += 1  # 1 cluster for all unassigned seqs
with open(dataNames, "r") as names:
    for line in names.readlines():
        fullUniqeId = line.split()[3]
        rnaClass, sep, tail = fullUniqeId.partition("_")
        rnaClass, sep, tail = fullUniqeId.partition("_")
        if fullUniqeId not in blackList:
            listOfClasses.append(rnaClass)
            listOfClusters.append(str(numberOfClusters))
            numberOfClusters += 1  # separate cluster for all unassigned seqs

toWrite = ""
for i in range(len(listOfClusters)):
    toWrite += listOfClasses[i] + "\t" + listOfClusters[i] + '\n'
with open("RESULTS/fullTab.tabular", "w") as full:
    full.write(toWrite)
