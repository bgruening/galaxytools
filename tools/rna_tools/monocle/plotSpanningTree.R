library(monocle)
library(argparser)

p <- arg.parser("Plot spanning tree for an ordered single cell dataset.");
# add input/output arguments
p <- add.argument(p, "ordered_dataset", help="Ordered Single Cell Dataset");
p <- add.argument(p, "xComponent", help="X axis component");
p <- add.argument(p, "yComponent", help="Y axis component");
p <- add.argument(p, "color_by", help="Color by this attribute (Pseudotime/State)");
p <- add.argument(p, "output", help="Output File");

# get input parameters
argv <- parse.args(p)

# get input file
OSCD <- readRDS(argv[["ordered_dataset"]])

# print in pdf
pdf(argv[["output"]])
plot_spanning_tree(OSCD, x = as.numeric(argv[["xComponent"]]), y = as.numeric(argv[["yComponent"]]), color_by = argv[["color_by"]])

