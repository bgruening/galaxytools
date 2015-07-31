library(monocle)
library(argparser)

p <- arg.parser("Plot genes' expression within a group for a single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
p <- add.argument(p, "genes", help="Genes short names");
p <- add.argument(p, "colorby", help="Color by this attribute");
p <- add.argument(p, "output", help="Output File");


#get input files
argv <- parse.args(p)

SCD <- readRDS(argv[["dataset"]])

# get data only for the specified genes
SSCD <- SCD[row.names(subset(fData(SCD), gene_short_name %in% unlist(strsplit(argv[["genes"]]," ")))),]


pdf(argv[["output"]])

# print in pdf
plot_genes_in_pseudotime(SSCD, color_by=argv[["colorby"]])
