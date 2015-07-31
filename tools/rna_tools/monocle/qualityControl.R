library(monocle)
library(argparser)
library(reshape2)

p <- arg.parser("Quality control for single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
p <- add.argument(p, "num_cells", help="Number of Cells each gene must be expressed in");
p <- add.argument(p, "output", help="Output File");


#get input files
#CoolValues <- parse.args(p, argv = "-input_sample_sheet -input_gene_annotation -input_expression_matrix -output")
argv <- parse.args(p)
pdf(argv[["output"]])

sessionInfo()

#SCD <- readRDS(argv[[1]])
SCD <- readRDS(argv[["dataset"]])
print(argv[["dataset"]])
print(argv[[1]])


#Filter data only taking expressed genes into account
expressed_genes <- row.names(subset(fData(SCD), num_cells_expressed >= strtoi(argv[["num_cells"]])))

#Filter more
#valid_cells <- row.names(subset(pData(SCD), Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1e+06))
#SCD <- SCD[, valid_cells]

L <- log(exprs(SCD[expressed_genes, ]))

melted_dens_df <- melt(t(scale(t(L))))

qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = "red") + xlab("Standardized log(FPKM)") + ylab("Density")

