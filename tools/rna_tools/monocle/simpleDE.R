library(monocle)
library(argparser)
#library(MASS)

p <- arg.parser("Simple differential expression analysis tool for a single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
#p <- add.argument(p, "orderedflag", help="ordered flag");
p <- add.argument(p, "genes", help="Genes short names");
p <- add.argument(p, "formulaString", help="Model Formula String");
p <- add.argument(p, "qvalue", help="q-value");
p <- add.argument(p, "output", help="Output File");


#get input files
argv <- parse.args(p)

SCD <- readRDS(argv[["dataset"]])

marker_genes <- row.names(subset(fData(SCD), gene_short_name %in% unlist(strsplit(argv[["genes"]]," "))))

# Test for differential expression
diff_test_res <- differentialGeneTest(SCD[marker_genes,], fullModelFormulaStr=argv[["formulaString"]])

# Select those genes that are significant at an FDR < qvalue
sig_genes <- subset(diff_test_res, qval < as.numeric(argv[["qvalue"]]))
# Add more information for those genes
sig_genes <- merge(fData(SCD), sig_genes, by="row.names")

# push to output
write.table(sig_genes,file=argv[["output"]],sep="\t",row.names=FALSE)
