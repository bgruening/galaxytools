import glob
from os import system
import re
from sklearn import metrics

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


pattern = re.compile("^RF.*$")


if len(listOfClasses) > 0 and  pattern.match(str(listOfClasses[0])):

    completeness_score = metrics.completeness_score(listOfClasses, listOfClusters)
    homogeneity_score = metrics.homogeneity_score(listOfClasses, listOfClusters)
    adjusted_rand_score = metrics.adjusted_rand_score(listOfClasses, listOfClusters)
    adjusted_mutual_info_score = metrics.adjusted_mutual_info_score(listOfClasses, listOfClusters)
    v_measure_score = metrics.v_measure_score(listOfClasses, listOfClusters)

    toWrite = "completeness_score : " + str(completeness_score) + "\n" + "homogeneity_score : " + str(homogeneity_score) + "\n" + "adjusted_rand_score : " +str(adjusted_rand_score)  + "\n" + "adjusted_mutual_info_score : " + str(adjusted_mutual_info_score)+ "\n" + "v_measure_score : " + str(v_measure_score)

else:
    toWrite = "completeness_score : NA \nhomogeneity_score : NA \nadjusted_rand_score : NA \nadjusted_mutual_info_score : NA \nv_measure_score : NA"

with open("RESULTS/evaluation.txt", "w") as fOut:
    fOut.write(toWrite)
