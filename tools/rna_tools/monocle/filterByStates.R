library(monocle)
library(argparser)

p <- arg.parser("Filter dataset by state for a single cell dataset.");
# add input/output arguments
p <- add.argument(p, "dataset", help="Single Cell Dataset");
p <- add.argument(p, "states", help="States that should be filtered out");
p <- add.argument(p, "filtered_dataset", help="Filtered Single Cell Dataset");

# get input parameters
argv <- parse.args(p)

# get input file
SCD <- readRDS(argv[["dataset"]])

# Filter for a specific gene list
SCD_filtered <- SCD[, !(pData(SCD)$State %in% as.numeric(unlist(strsplit(argv[["states"]]," "))))]

saveRDS(SCD_filtered,file=argv[["filtered_dataset"]])



