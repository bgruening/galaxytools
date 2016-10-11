import glob, os, sys
import shutil
from shutil import copyfile
from os import system
import re

def sh(script):
    system("bash -c '%s'" % script)


#zipFile = sys.argv[1]
tabularFiles = sys.argv[1]
dataNames = sys.argv[2]

#sh("unzip " + zipFile + " -d ./")

#sh("rm -r RESULTS")
sh("mkdir -p RESULTS")

#sh("rm  RESULTS/*/*.stats")

#print "tab = ", tabularFiles
#print "names = ", dataNames

tabfiles = tabularFiles.split(',')
#print "1 = " , tabfiles[0]
#print "2 = " , tabfiles[1]



#tabularPath = "CLUSTER/*/CMSEARCH/model.tree.cm.tabresult"
#tabularFiles = glob.glob(tabularPath)


for i in range(len(tabfiles)):
    seqNum = 0
    stats = ""
    #clustNum=re.findall(r'\.(.+?)\.',tab)
    #clustNum = clustNum[0]
    clustNum = i+1
    clustNum = str(clustNum)
    uniqueIds = []
#    print "################CLUSTER ", clustNum, " ################CLUSTER "

    writeToClusterSeqStats = ""

    #sh("mkdir -p RESULTS/"+clustNum)


    with open(tabfiles[i], "r") as f:
        for l in f.readlines():
            #print l
            li=l.strip()
            if not li.startswith("#"):
                seqNum = seqNum +1
                origId = l.split()[18]
                #print "unique id = ", origId
                if origId not in uniqueIds:
                    uniqueIds.append(origId)
                with open(dataNames, "r") as names:
                    for line in names.readlines():
                        seq, sep, tail = line.partition("#")
                        seq =  seq.split()[1]
                        if l.split()[0] == seq:
                            #print line.split()[1]
                            #print l.split()[0]
                            writeToClusterSeqStats += "CLUSTER  " + clustNum + "  " + str(line.split()[1]) + " " + str(line.split()[2]) + " " + str(line.split()[3]) + " " +  str(line.split()[4]) + " " +  str(line.split()[5]) +"\n"
        with open("RESULTS/"+clustNum + ".cluster.seqs.stats.txt", "w") as fOut:
            fOut.write(writeToClusterSeqStats)
        writeToClusterStats = ""
        writeToClusterStats = "CLUSTER  " + str(clustNum) + "  " + "SEQS  " + str(seqNum)
        #writeToClusterSeqStats
        writeToClusterStats += "  IDS_UNIQUE  " + str(len(uniqueIds)) + "\n"
        with open("RESULTS/final.clusters.stats.txt", "a") as fOut:
    		fOut.write(writeToClusterStats)


sh("python evaluation.py " + dataNames)





# for tab in tabularFiles:
#     seqNum = 0
#     stats = ""
#     clustNum=re.findall(r'\.(.+?)\.',tab)
#     clustNum = clustNum[0]
#     uniqueIds = []
#     print "################CLUSTER ", clustNum, " ################CLUSTER "
#
#     writeToClusterSeqStats = ""
#
#     sh("mkdir -p RESULTS/"+clustNum)
#
#
#     with open(tab, "r") as f:
#         for l in f.readlines():
#             li=l.strip()
#             if not li.startswith("#"):
#                 seqNum = seqNum +1
#                 origId = l.split()[18]
#                 if origId not in uniqueIds:
#                     uniqueIds.append(origId)
#                 with open(dataNames, "r") as names:
#                     for line in names.readlines():
#                         seq, sep, tail = line.partition("#")
#                         seq =  seq.split()[1]
#                         if l.split()[0] == seq:
#                             #print line.split()[1]
#                             #print l.split()[0]
#                             writeToClusterSeqStats += "CLUSTER  " + clustNum + "  " + str(line.split()[1]) + " " + str(line.split()[2]) + " " + str(line.split()[3]) + " " +  str(line.split()[4]) + " " +  str(line.split()[5]) +"\n"
#         with open("RESULTS/"+clustNum + "/cluster.seqs.stats", "w") as fOut:
#             fOut.write(writeToClusterSeqStats)
#
#         writeToClusterStats = "CLUSTER  " + str(clustNum) + "  " + "SEQS  " + str(seqNum)
#         #writeToClusterSeqStats
#         writeToClusterStats += "  IDS_UNIQUE  " + str(len(uniqueIds)) + "\n"
#         with open("RESULTS/"+clustNum + "/cluster.stats", "w") as fOut:
#     		fOut.write(writeToClusterStats)
#
#
# sh("python evaluation.py " + dataNames)
