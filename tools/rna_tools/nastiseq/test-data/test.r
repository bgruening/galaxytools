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

genepos = read.table("input_TAIR10_annotation.tsv", header = TRUE, sep = "\t")
genepos$attributes = as.character(genepos$attributes)
rownames(genepos) = genepos$attributes

pospairs = read.table("input_positive_pair.tsv", sep = "\t", as.is = TRUE)

smat = as.matrix(read.table("input_read_count_smt.tsv",  sep = "\t",  row.names = 1))
colnames(smat) = NULL

asmat = as.matrix(read.table("input_read_count_asmt.tsv",  sep = "\t",  row.names = 1))
colnames(asmat) = NULL

WRscore = getNASTIscore(smat, asmat)

negpairs = getnegativepairs(genepos)

WRpred = NASTIpredict(smat,asmat, pospairs, negpairs)

WRpred_rocr = prediction(WRpred$predictions,WRpred$labels)

thr = defineFDR(WRpred_rocr,0.05)

WR_names = FindNATs(WRscore, thr, pospairs, genepos)

write.table(WR_names$newpairs, file = "output_newpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$neworphan, file = "output_neworphan.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$knownpairs, file = "output_knownpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

