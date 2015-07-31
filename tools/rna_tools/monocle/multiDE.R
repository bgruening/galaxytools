library(monocle)
library(argparser)

p <- arg.parser("Multifactorial differential expression analysis tool for a single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
#p <- add.argument(p, "orderedflag", help="ordered flag");
p <- add.argument(p, "genes", help="Genes short names");
p <- add.argument(p, "fullFormulaString", help="Full Model Formula String");
p <- add.argument(p, "reducedFormulaString", help="Reduced Model Formula String");
p <- add.argument(p, "output", help="Output File");

#parse arguments
argv <- parse.args(p)

#get input files
SCD <- readRDS(argv[["dataset"]])

marker_genes <- row.names(subset(fData(SCD), gene_short_name %in% unlist(strsplit(argv[["genes"]]," "))))

# Test for differential expression
diff_test_res <- differentialGeneTest(SCD[marker_genes,], fullModelFormulaStr=argv[["formulaString"]], reducedModelFormulaStr=argv[["reducedFormulaString"]])

# push to output
write.table(diff_test_res, file=argv[["output"]],sep="\t",row.names=FALSE)
