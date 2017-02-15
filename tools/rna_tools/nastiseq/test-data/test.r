library(NASTIseq)

### generation of test set
# data(WholeRoot)
#
# write.table(WholeRoot$genepos, file = "input_TAIR10_annotation.tsv",row.names = FALSE ,  sep = "\t", quote = FALSE)
#
# write.table(WholeRoot$smat, file = "input_read_count_smt.tsv",  col.names = FALSE,  sep = "\t", quote = FALSE)
# write.table(WholeRoot$asmat, file = "input_read_count_asmt.tsv", col.names = FALSE, sep = "\t", quote = FALSE)
#
# write.table(WholeRoot$pospairs, file = "input_positive_pair.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

mydata = list()

mydata$genepos = read.table("input_TAIR10_annotation.tsv", header = TRUE, sep = "\t")
mydata$genepos$attributes = as.character(mydata$genepos$attributes)
rownames(mydata$genepos) = mydata$genepos$attributes

mydata$pospairs = read.table("input_positive_pair.tsv", sep = "\t", as.is = TRUE)

mydata$smat = as.matrix(read.table("input_read_count_smt.tsv",  sep = "\t",  row.names = 1))
colnames(mydata$smat) = NULL

mydata$asmat = as.matrix(read.table("input_read_count_asmt.tsv",  sep = "\t",  row.names = 1))
colnames(mydata$asmat) = NULL

WRscore<-getNASTIscore(mydata$smat, mydata$asmat)

negpairs<-getnegativepairs(mydata$genepos)

WRpred<-NASTIpredict(mydata$smat,mydata$asmat, mydata$pospairs, negpairs)

WRpred_rocr<-prediction(WRpred$predictions,WRpred$labels)

thr<-defineFDR(WRpred_rocr,0.05)

WR_names<-FindNATs(WRscore, thr, mydata$pospairs, mydata$genepos)

write.table(WR_names$newpairs, file = "output_newpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$neworphan, file = "output_neworphan.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$knownpairs, file = "output_knownpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

