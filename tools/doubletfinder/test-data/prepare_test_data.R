#!/usr/bin/env Rscript

# Regenerates test-data/pbmc_small.rda for the DoubletFinder Galaxy wrapper.
#
# The test data is the pbmc_small SeuratObject (230 features x 80 cells,
# already pre-processed with NormalizeData, FindVariableFeatures, ScaleData
# and RunPCA) shipped with the DoubletFinder package itself:
#   https://github.com/chris-mcginnis-ucsf/DoubletFinder/blob/master/data/pbmc_small.rda
# which in turn originates from the SeuratObject package:
#   https://github.com/satijalab/seurat-object/blob/main/data/pbmc_small.rda

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
})

data("pbmc_small", package = "SeuratObject")

stopifnot(inherits(pbmc_small, "Seurat"))
stopifnot(ncol(pbmc_small) == 80)

## Convert to the v5 assay structure: DoubletFinder 2.0.6 calls the defunct
## GetAssayData(slot = ...) on objects reporting a version < 5.0, which fails
## with SeuratObject >= 5.5.
pbmc_small <- UpdateSeuratObject(pbmc_small)
stopifnot(as.numeric(SeuratObject::Version(pbmc_small)[1, 1]) >= 5)

saveRDS(pbmc_small, file = "test-data/pbmc_small.rds")
