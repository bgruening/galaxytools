library(monocle)
library(argparser)
#library(MASS)

p <- arg.parser("Order cells by certain genes' expression for a single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
p <- add.argument(p, "genes", help="Genes short names");
p <- add.argument(p, "formulaString", help="Model Formula String");
p <- add.argument(p, "qvalue", help="Threshold q-value for significantly differentially expressed genes.");
p <- add.argument(p, "components", help="Number of independent components to reduce to");
p <- add.argument(p, "num_paths", help="Number of paths");
p <- add.argument(p, "reverse_path", help="Reverse the path (True/False)");
p <- add.argument(p, "ordered_dataset", help="Ordered Single Cell Dataset");

# get input parameters
argv <- parse.args(p)

# get input file
SCD <- readRDS(argv[["dataset"]])

# filter for genes that are expressed in at least so many cells
#expressed_genes <- row.names(subset(fData(SCD), num_cells_expressed >= argv[["num_cells"]]))

# Filter for a specific gene list
marker_genes <- row.names(subset(fData(SCD), gene_short_name %in% unlist(strsplit(argv[["genes"]]," "))))

# test for differential expression
diff_test_res <- differentialGeneTest(SCD[marker_genes,], fullModelFormulaStr=argv[["formulaString"]])

# filter by q-value
ordering_genes <- row.names(subset(diff_test_res, as.numeric(qval) < as.numeric(argv[["qvalue"]])))

# use both filters
ordering_genes <- intersect(ordering_genes,marker_genes)
SCD <- setOrderingFilter(SCD, ordering_genes)


# reduce dimensions. IRLBA package is not being used (could be set as option..)
# Number of components can be set. 
SCD <- reduceDimension(SCD, max_components = as.numeric(argv[["components"]]), use_irlba = F)

# order cells and save the result
SCD <- orderCells(SCD, num_paths = as.numeric(argv[["num_paths"]]), reverse = argv[["reverse_path"]])
saveRDS(SCD,file=argv[["ordered_dataset"]])



