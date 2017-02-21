library(NASTIseq)

## generation of test set
# data(WholeRoot)
# WholeRoot$genepos$feature <- 'gene'
#
# set_attri <- function(attri){
#   attri = paste('gene_id ', '"', attri, '"', ';', sep = '')
#   return(attri)
# }
#
# WholeRoot$genepos$attributes = as.character(lapply(as.character(WholeRoot$genepos$attributes), set_attri))
#
# write.table(WholeRoot$genepos, file = "input_TAIR10_annotation.gtf", row.names = FALSE,  col.names = FALSE,  sep = "\t",  quote = FALSE)
# write.table(WholeRoot$smat, file = "input_read_count_smt.tsv",  col.names = FALSE,  sep = "\t", quote = FALSE)
# write.table(WholeRoot$asmat, file = "input_read_count_asmt.tsv", col.names = FALSE, sep = "\t", quote = FALSE)
#
# write.table(WholeRoot$pospairs, file = "input_positive_pair.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

genepos = read.delim("input_TAIR10_annotation.gtf", header=FALSE, comment.char="#")
colnames(genepos) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
genepos = subset(genepos, feature=="gene")

get_id = function(attri){
  gene_info = strsplit(attri, ";")[[1]][1]
  gene_id = strsplit(gene_info, " ")[[1]][2]
  gene_id = gsub("\"", "", gene_id)
  return(gene_id)
}

genepos$attributes = as.character(lapply(as.character(genepos$attributes), get_id))

pospairs = read.table("input_positive_pair.tsv", sep = "\t", as.is = TRUE)

smat = as.matrix(read.table("input_read_count_smt.tsv",  sep = "\t",  row.names = 1))

asmat = as.matrix(read.table("input_read_count_asmt.tsv",  sep = "\t",  row.names = 1))

WRscore = getNASTIscore(smat, asmat)

negpairs = getnegativepairs(genepos)

WRpred = NASTIpredict(smat,asmat, pospairs, negpairs)

WRpred_rocr = prediction(WRpred$predictions,WRpred$labels)

thr = defineFDR(WRpred_rocr,0.05)

WR_names = FindNATs(WRscore, thr, pospairs, genepos)

write.table(WR_names$newpairs, file = "output_newpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$neworphan, file = "output_neworphan.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
