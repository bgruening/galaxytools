library(NASTIseq)
library(data.table)

## generation of test set
data(WholeRoot)
WholeRoot$genepos$feature <- 'gene'

set_attri <- function(attri){
  attri = paste('gene_id ', '"', attri, '"', ';', sep = '')
  return(attri)
}

WholeRoot$genepos$attributes = as.character(lapply(as.character(WholeRoot$genepos$attributes), set_attri))

#
# write.table(WholeRoot$genepos, file = "input_TAIR10_annotation.tsv",row.names = FALSE ,  sep = "\t", quote = FALSE)
#
# write.table(WholeRoot$smat, file = "input_read_count_smt.tsv",  col.names = FALSE,  sep = "\t", quote = FALSE)
# write.table(WholeRoot$asmat, file = "input_read_count_asmt.tsv", col.names = FALSE, sep = "\t", quote = FALSE)
#
# write.table(WholeRoot$pospairs, file = "input_positive_pair.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# genepos = read.table("Arabidopsis_thaliana.TAIR10.34.gtf", sep = "\t")
# genepos$attributes = as.character(genepos$attributes)

genepos = read.delim("Arabidopsis_thaliana.TAIR10.34.gtf", header=FALSE, comment.char="#")

# genepos = fread('Arabidopsis_thaliana.TAIR10.34.gtf', stringsAsFactors=TRUE)
colnames(genepos) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

genepos = subset(genepos, feature=="gene")

get_id <- function(attri){
  gene_info = strsplit(attri, ";")[[1]][1]
  gene_id = strsplit(gene_info, " ")[[1]][2]
  gene_id = gsub("\"", "", gene_id)
  return(gene_id)
}

genepos$attributes = as.character(lapply(as.character(genepos$attributes), get_id))
# rownames(genepos) = genepos$attributes

pospairs = read.table("input_positive_pair.tsv", sep = "\t", as.is = TRUE)

smat = as.matrix(read.table("input_read_count_smt.tsv",  sep = "\t",  row.names = 1))
colnames(smat) = NULL

asmat = as.matrix(read.table("input_read_count_asmt.tsv",  sep = "\t",  row.names = 1))
colnames(asmat) = NULL

WRscore = getNASTIscore(smat, asmat)

negpairs = getnegativepairs(genepos)

WRpred = NASTIpredict(sm,asmat, pospairs, negpairs)

WRpred_rocr = prediction(WRpred$predictions,WRpred$labels)

thr = defineFDR(WRpred_rocr,0.05)

WR_names = FindNATs(WRscore, thr, pospairs, genepos)

write.table(WR_names$newpairs, file = "output_newpairs.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(WR_names$neworphan, file = "output_neworphan.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

