library(monocle)
library(argparser)

p <- arg.parser("Create Cell Data Set object from text files.");
# add input/output arguments
p <- add.argument(p, "input_expression_matrix", help="Expression Matrix Input");
p <- add.argument(p, "input_sample_sheet", help="Sample Sheet Input");
p <- add.argument(p, "input_gene_annotation", help="Gene Annotation Input");
p <- add.argument(p, "output", help="Output File");

#get input files
#CoolValues <- parse.args(p, argv = "-input_sample_sheet -input_gene_annotation -input_expression_matrix -output")
argv <- parse.args(p)

#argv["output"] 

# create CellDataSet

expr_matrix <- read.table(argv[["input_expression_matrix"]])
sample_sheet <- read.table(argv[["input_sample_sheet"]])
gene_annotation <- read.table(argv[["input_gene_annotation"]])

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
HSMM <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd)

# push to output
saveRDS(HSMM,file=argv[["output"]])

