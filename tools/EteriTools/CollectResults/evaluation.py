import glob, os, sys
import shutil
from shutil import copyfile
from os import system
import re
#import sklearn
#from sklearn import metrics
#from tabulate import tabulate



def sh(script):
    system("bash -c '%s'" % script)

#dataNames = "test-data/FASTA/data.names"
dataNames = "FASTA/data.names"
#path = sys.argv[2]


listOfClusters = []
listOfClasses = []

#cluster_seqs_stats_path = "test-data/RESULTS/*.cluster.all"

cluster_seqs_stats_path = "RESULTS/*.cluster.all"
cluster_seqs_stats_files = glob.glob(cluster_seqs_stats_path)

blackList =[]
numberOfClusters = 0
for singleFile in sorted(cluster_seqs_stats_files):
    numberOfClusters += 1
    with open(singleFile, "r") as f:
        for line in f.readlines():
            uniqueId = line.split()[6] #was 8
            print "unique id = ", uniqueId
            clustNum = line.split()[1]
            rnaClass, sep, tail = uniqueId.partition("_")
            #classToClusterDict[rnaClass] = clustNum
            listOfClasses.append(rnaClass)
            listOfClusters.append(clustNum)
            print "rnaClass = ", rnaClass
            print "clustNum = ", clustNum


            with open(dataNames, "r") as names:
                for nameLine in names.readlines():
                    fullUniqeId = nameLine.split()[3] #waas 6
                    print "full uniqueId = ", fullUniqeId
                    rnaClass, sep, tail = fullUniqeId.partition("_")

                    n = 2
                    short_unique =  re.findall("_".join(["[^_]+"] * n), fullUniqeId)[0]

                    if short_unique == uniqueId:
                        blackList.append(uniqueId)

print "blackList = ", blackList
print "numberOfClusters = " ,numberOfClusters
numberOfClusters += 1 ### 1 cluster for all unassigned seqs
with open(dataNames, "r") as names:
    for nameLine in names.readlines():
        fullUniqeId = nameLine.split()[3] # was 6
        rnaClass, sep, tail = fullUniqeId.partition("_")
        print "fullUniqeId = ", fullUniqeId
        n = 2
        short_unique =  re.findall("_".join(["[^_]+"] * n), fullUniqeId)[0]
        rnaClass, sep, tail = fullUniqeId.partition("_")
        if short_unique not in blackList:

            listOfClasses.append(rnaClass)
            #listOfClasses.append(fullUniqeId)
            listOfClusters.append(str(numberOfClusters))
            numberOfClusters += 1 ### separate cluster for all unassigned seqs

toWrite = ""
#table=[]
for i in range(len(listOfClusters)):
    toWrite += listOfClasses[i] + "\t" + listOfClusters[i] + '\n'
    #table += [[listOfClasses[i],listOfClusters[i]]]

#f = open('table.txt', 'w')
#f.write(tabulate(table))
#f.close()

#with open("fullTab.tabular", "w") as full:
#    full.write(toWrite)

with open("RESULTS/fullTab.tabular", "w") as full:
    full.write(toWrite)

# with open("listOfClusters.txt", "w") as fclt:
#     fclt.write("" + str(listOfClusters))
#
# with open("listOfClasses.txt", "w") as fcls:
#     fcls.write("" + str(listOfClasses))

#completeness_score = metrics.completeness_score(listOfClasses, listOfClusters)
# homogeneity_score = metrics.homogeneity_score(listOfClasses, listOfClusters)
# adjusted_rand_score = metrics.adjusted_rand_score(listOfClasses, listOfClusters)
# adjusted_mutual_info_score = metrics.adjusted_mutual_info_score(listOfClasses, listOfClusters)
# v_measure_score = metrics.v_measure_score(listOfClasses, listOfClusters)
#
# toWrite = "completeness_score : " + str(completeness_score) + "\n" + "homogeneity_score : " + str(homogeneity_score) + "\n" + "adjusted_rand_score : " +str(adjusted_rand_score)  + "\n" + "adjusted_mutual_info_score : " + str(adjusted_mutual_info_score)+ "\n" + "v_measure_score : " + str(v_measure_score)
#
# with open("RESULTS/evaluation.txt", "w") as fOut:
#     fOut.write(toWrite)
#
#
#
#print "compl_score = ", completeness_score, "\n"
# print "homog_score = ", homogeneity_score , "\n"
# print "v_measure_score = ", v_measure_score , "\n"
# print "adjusted_rand_score = ", adjusted_rand_score , "\n"
# print "adjusted_mutual_info_score = ", adjusted_mutual_info_score , "\n"
