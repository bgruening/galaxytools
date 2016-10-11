import glob, os, sys
import shutil
from shutil import copyfile
from os import system
import re

def sh(script):
    system("bash -c '%s'" % script)

seqStatPath = "RESULTS/*.seqs.stats.txt"
seqStatPathFiles = glob.glob(seqStatPath)
listofSets = []

for i in range(len(seqStatPathFiles)):
    writeToClusterSeqStats = ""
    #print seqStatPathFiles[i]
    clustNum =  re.findall(r'\d+', seqStatPathFiles[i])
    #print clustNum
    for j in range(len(seqStatPathFiles)):
     #   print j
        j = j + 1
        if j < len(seqStatPathFiles) and i<j :
            idSet1 = set([])
            idSet2 = set([])
            if i != j :
                print "i = ", seqStatPathFiles [i]
                print "j = ", seqStatPathFiles [j]
                with open(seqStatPathFiles[i], "r") as f1:
                    for l1 in f1.readlines():
                        #print l
                        li1=l1.strip()
                        #seqNum = seqNum +1
                        origId1 = l1.split()[4]
                        #print origId
                        idSet1.add(origId1)
                with open(seqStatPathFiles[j], "r") as f2:
                    for l2 in f2.readlines():
                        #print l
                        li2=l2.strip()
                        #seqNum = seqNum +1
                        origId2 = l2.split()[4]
                        #print origId
                        idSet2.add(origId2)
                lenOfSet = len(idSet1)
                lenOfOverlap = len(idSet1.intersection(idSet2))
                percentOfOverlap = (lenOfOverlap*100)/lenOfSet
                if percentOfOverlap < 50:
                    #print "smaller than 50 = ", percentOfOverlap
                    copyfile(seqStatPathFiles[i], "RESULTS/"+str(clustNum[0]) + ".cluster.seqs.stats.SOFT.txt")
                else :
                    print "bigger  than 50  = ", percentOfOverlap
                    print "len of first =  ", len(idSet1)
                    print "len of second =  ", len(idSet2)
                    print "len of diff =  ", len(idSet1.difference(idSet2) )
                #    print "diff = " , idSet1.difference(idSet2)
                    diff = idSet1.difference(idSet2)
                    diff = list(diff)
                #    print diff
                    with open(seqStatPathFiles[j], "r") as f1:
                        for l1 in f1.readlines():
                            writeToClusterSeqStats += l1

                    with open(seqStatPathFiles[i], "r") as f2:
                        for l2 in f2.readlines():
                            for index in range(len(diff)):
                                #print l2
                                #print diff[index]
                                #print l2.split()[4]
                                if diff[index] == l2.split()[4]:
                                    writeToClusterSeqStats += l2 


                    with open("RESULTS/"+str(clustNum[0]) + ".cluster.seqs.stats.SOFT.txt", "w") as fOut:
                        fOut.write(writeToClusterSeqStats)


        #listofSets.append(idSet)

# for i in range(len(listofSets)):
#     for j in range(len(listofSets)):
#         if i != j :
#             lenOfSet =  len(listofSets[i])
#             lenOfOverlap =  len(listofSets[i].intersection(listofSets[j]))
#             percentOfOverlap = (lenOfOverlap*100)/lenOfSet
#             print "lenOfSet " , i , " = ", lenOfSet
#             print "lenOfOverlap  with ", j , " = ", lenOfOverlap
#             print "percentOfOverlap = ", percentOfOverlap
