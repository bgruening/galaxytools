suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))

## Assay Data
## F rows of features and S columns of samples
exprs_file <- "$exprs_file"   ## tabular file
exprs <- as.matrix(read.table(exprs_file, header = T, sep = "\t",
                  row.names = 1, as.is = T))

## Phenotype Data
## S rows of samples, and V columns of covariates (e.g. sex, age, etc.)
pdata_file <- "$pdata_file"
pdata <- read.table(pdata_file, row.names = 1, header = T, sep = "\t")

stop(all(rownames(pdata) != colnames(exprs)),
     "Number of Samples between phenotypes and assays are not the same")

metadata <- data.frame(
    labelDescription = c(
        "Patient gender",
        "Case/Control",
        "Tumor progress"),
    row.names = c("gender", "type", "score"))

pheno_data <- new("AnnotatedDataFrame", data = pdata, varMetadata = metadata)

## Annotation and Feature Data, or just a string for type of chip used
annotation <- "hgu95av2"

## Experiment Description
experiment_data <- new(
    "MIAME",
    name = "Pierre Fermat",
    lab = "Francis Galton Lab",
    contact = "pfermat@lab.not.exist",
    title = "Smoking-Cancer Experiment",
    abstract = "An example ExpressionSet",
    url = "www.lab.not.exist",
    other = list(
        notes = "Created from text files"
    ))

e_set <- ExpressionSet(assayData = exprs,
                      phenoData = pheno_data,
                      experimentData = experiment_data,
                      annotation = annotation)

## print to file
capture.output(print(e_set), file = outfile_tab)  # nolint
