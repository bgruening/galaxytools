library(monocle)
library(argparser)

p <- arg.parser("Cluster genes' by expression pattern.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
p <- add.argument(p, "formulaString", help="Model Formula String");
p <- add.argument(p, "samples", help="Number of samples to randomly choose");
p <- add.argument(p, "clusters", help="Number of clusters the samples should be attributed to");
p <- add.argument(p, "output", help="Output File");


argv <- parse.args(p)

#get input files
SCD <- readRDS(argv[["dataset"]])

# Fit model to the expression of randomly chosen data
samples <- sample(nrow(fData(SCD)),as.numeric(argv[["samples"]]))
full_model_fits <- fitModel(SCD[samples, ], modelFormulaStr=argv[["formulaString"]])

# Calculate expression curve matrix
expression_curve_matrix <- responseMatrix(full_model_fits)

clusters <- clusterGenes(expression_curve_matrix, k=as.numeric(argv[["clusters"]]))

pdf(argv[["output"]])

# print in pdf
plot_clusters(SCD[samples,], clusters)
