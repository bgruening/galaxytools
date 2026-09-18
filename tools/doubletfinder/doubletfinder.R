#!/usr/bin/env Rscript

# Galaxy wrapper for DoubletFinder
# Implements the standard 4-step workflow:
#   paramSweep -> summarizeSweep -> find.pK -> doubletFinder
# with optional homotypic adjustment (modelHomotypic) and pANN reuse.

suppressPackageStartupMessages({
    library(DoubletFinder)
    library(Seurat)
    library(argparse)
})

main <- function() {
    parser <- ArgumentParser(description = "DoubletFinder: doublet detection in single-cell RNA-seq data using artificial nearest neighbors")
    parser$add_argument("--input-object", type = "character", required = TRUE,
        help = "Path to pre-processed Seurat object saved as RDS or RData")
    parser$add_argument("--output-object", type = "character", required = TRUE,
        help = "Path to write annotated Seurat object (RDS)")
    parser$add_argument("--output-table", type = "character", required = TRUE,
        help = "Path to write per-cell classification table (tabular)")
    parser$add_argument("--output-sweep", type = "character", required = TRUE,
        help = "Path to write pK sweep statistics table (tabular)")
    parser$add_argument("--output-plot", type = "character", required = FALSE, default = NULL,
        help = "Path to write BCmvn vs pK plot (PNG)")
    parser$add_argument("--pcs", type = "integer", required = TRUE,
        help = "Number of statistically significant PCs to use")
    parser$add_argument("--pn", type = "double", required = TRUE, default = 0.25,
        help = "Proportion of artificial doublets (pN)")
    parser$add_argument("--pk-mode", type = "character", required = TRUE,
        choices = c("auto", "manual"),
        help = "Whether to select pK automatically via paramSweep/find.pK, or use a user-supplied value")
    parser$add_argument("--pk", type = "double", required = FALSE, default = NULL,
        help = "User-supplied pK value (manual mode)")
    parser$add_argument("--reuse-pann", type = "character", required = FALSE, default = NULL,
        help = "Name of an existing pANN metadata column to reuse (manual mode) instead of recomputing")
    parser$add_argument("--nexp-mode", type = "character", required = TRUE,
        choices = c("expected_rate", "expected_count"),
        help = "Whether the expected number of doublets is given as a rate or a count")
    parser$add_argument("--doublet-rate", type = "double", required = FALSE, default = 0.075,
        help = "Expected doublet rate (expected_rate mode)")
    parser$add_argument("--homotypic-adjust", action = "store_true",
        help = "Adjust the expected doublet count for homotypic doublets using modelHomotypic")
    parser$add_argument("--annotation-col", type = "character", required = FALSE, default = NULL,
        help = "Seurat meta.data column with cell-type annotations used for homotypic adjustment")
    parser$add_argument("--nexp", type = "integer", required = FALSE, default = NULL,
        help = "Expected number of doublets (expected_count mode)")
    parser$add_argument("--sct", action = "store_true",
        help = "The input Seurat object was normalized with SCTransform")
    parser$add_argument("--save-sweep-table", action = "store_true",
        help = "Save the pK sweep statistics table")
    parser$add_argument("--save-plot", action = "store_true",
        help = "Save the BCmvn vs pK plot")
    parser$add_argument("--seed", type = "integer", required = FALSE, default = NULL,
        help = "Random seed for reproducibility")

    args <- parser$parse_args()

    if (!is.null(args$seed)) {
        set.seed(args$seed)
    }

    ## ---------------------------------------------------------------
    ## Load Seurat object (RDS, or RData containing a single Seurat object)
    ## ---------------------------------------------------------------
    read_seurat_object <- function(path) {
        if (!file.exists(path)) {
            stop("Input object file not found: ", path)
        }
        obj <- tryCatch(readRDS(path), error = function(e) NULL)
        if (is.null(obj)) {
            env <- new.env(parent = emptyenv())
            load(path, envir = env)
            candidates <- ls(env)
            seurat_objs <- candidates[vapply(candidates, function(nm) {
                inherits(get(nm, envir = env), "Seurat")
            }, logical(1))]
            if (length(seurat_objs) == 0) {
                stop("No Seurat object found in input file: ", path)
            }
            obj <- get(seurat_objs[1], envir = env)
        }
        if (!inherits(obj, "Seurat")) {
            stop("Input file does not contain a Seurat object.")
        }
        obj
    }

    message("Reading Seurat object...")
    seu <- read_seurat_object(args$input_object)

    n_cells <- ncol(seu)
    message("Loaded Seurat object with ", n_cells, " cells")

    ## ---------------------------------------------------------------
    ## Determine pK
    ## ---------------------------------------------------------------
    sweep_ran <- FALSE
    if (args$pk_mode == "auto") {
        message("Running paramSweep...")
        sweep_res <- DoubletFinder::paramSweep(seu, PCs = 1:args$pcs, sct = args$sct)

        message("Summarizing sweep results...")
        sweep_stats <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)

        message("Finding optimal pK...")
        bcmvn <- DoubletFinder::find.pK(sweep_stats)
        sweep_ran <- TRUE

        if (all(is.na(bcmvn$BCmetric) | is.nan(bcmvn$BCmetric))) {
            stop("Could not determine an optimal pK: all BCmvn values are NA. ",
                 "Consider using manual pK mode instead.")
        }
        opt_row <- which.max(bcmvn$BCmetric)
        pK <- bcmvn$pK[opt_row]
        message("Optimal pK (max BCmvn): ", pK)

        if (args$save_sweep_table) {
            utils::write.table(bcmvn, file = args$output_sweep, sep = "\t",
                quote = FALSE, row.names = FALSE)
        }

        if (args$save_plot && !is.null(args$output_plot)) {
            grDevices::png(args$output_plot, width = 800, height = 600)
            graphics::plot(bcmvn$ParamID, bcmvn$BCmetric,
                type = "b", pch = 16, col = "#41b6c4",
                xlab = "ParamID", ylab = "BCmvn",
                main = paste0("Optimal pK = ", pK))
            graphics::abline(v = bcmvn$ParamID[opt_row], lty = 2, col = "grey50")
            grDevices::dev.off()
        }
    } else {
        if (is.null(args$pk)) {
            stop("--pk is required when --pk-mode is 'manual'")
        }
        pK <- args$pk
        message("Using user-supplied pK: ", pK)
    }

    ## ---------------------------------------------------------------
    ## Determine nExp
    ## ---------------------------------------------------------------
    if (args$nexp_mode == "expected_rate") {
        nExp <- round(args$doublet_rate * n_cells)
        message("Expected number of doublets (rate-based): ", nExp)

        if (args$homotypic_adjust) {
            if (is.null(args$annotation_col)) {
                stop("--annotation-col is required for homotypic adjustment")
            }
            if (!args$annotation_col %in% colnames(seu@meta.data)) {
                stop("Annotation column '", args$annotation_col,
                     "' not found in Seurat meta.data. Available columns: ",
                     paste(colnames(seu@meta.data), collapse = ", "))
            }
            annotations_vec <- as.character(seu@meta.data[[args$annotation_col]])
            if (any(is.na(annotations_vec))) {
                stop("Annotation column '", args$annotation_col, "' contains NA values")
            }
            homotypic_prop <- DoubletFinder::modelHomotypic(annotations_vec)
            nExp <- round(nExp * (1 - homotypic_prop))
            message("Homotypic doublet proportion: ", homotypic_prop)
            message("Adjusted expected number of doublets: ", nExp)
        }
    } else {
        if (is.null(args$nexp)) {
            stop("--nexp is required when --nexp-mode is 'expected_count'")
        }
        nExp <- as.integer(args$nexp)
        message("Expected number of doublets (count-based): ", nExp)
    }

    if (nExp <= 0) {
        stop("Expected number of doublets must be positive, got: ", nExp)
    }
    if (nExp >= n_cells) {
        stop("Expected number of doublets (", nExp,
             ") must be smaller than the number of cells (", n_cells, ")")
    }

    ## ---------------------------------------------------------------
    ## Run doubletFinder
    ## ---------------------------------------------------------------
    if (!is.null(args$reuse_pann) && !args$reuse_pann %in% colnames(seu@meta.data)) {
        stop("pANN column '", args$reuse_pann,
             "' not found in Seurat meta.data. Available columns: ",
             paste(colnames(seu@meta.data), collapse = ", "))
    }

    message("Running doubletFinder with pN = ", args$pn, ", pK = ", pK, ", nExp = ", nExp)
    seu <- DoubletFinder::doubletFinder(
        seu,
        PCs = 1:args$pcs,
        pN = args$pn,
        pK = pK,
        nExp = nExp,
        reuse.pANN = args$reuse_pann,
        sct = args$sct
    )

    ## ---------------------------------------------------------------
    ## Write outputs
    ## ---------------------------------------------------------------
    message("Writing outputs...")

    # Identify the DF result columns just created
    class_col <- grep("^DF\\.classifications", colnames(seu@meta.data), value = TRUE)
    pann_col <- grep("^pANN", colnames(seu@meta.data), value = TRUE)

    if (length(class_col) == 0) {
        stop("DoubletFinder did not add classification columns to meta.data")
    }

    # Use the most recently added classification columns
    class_col <- class_col[length(class_col)]
    pann_col <- if (length(pann_col) > 0) pann_col[length(pann_col)] else NA_character_

    cell_table <- data.frame(
        cell = rownames(seu@meta.data),
        stringsAsFactors = FALSE
    )
    if (!is.na(pann_col)) {
        cell_table$pANN <- seu@meta.data[[pann_col]]
    }
    cell_table$classification <- seu@meta.data[[class_col]]

    utils::write.table(cell_table, file = args$output_table, sep = "\t",
        quote = FALSE, row.names = FALSE)

    # Create an (empty) sweep file if the sweep was not run, so the
    # declared output path always exists
    if (!sweep_ran && args$save_sweep_table) {
        file.create(args$output_sweep)
    }

    saveRDS(seu, file = args$output_object)

    message("Done. Classified ", sum(seu@meta.data[[class_col]] == "Doublet"),
            " doublets out of ", n_cells, " cells.")
}

# Route messages and warnings to stdout (instead of stderr) so that only
# genuine errors produce stderr output; errors still abort with exit code 1.
withCallingHandlers(
    main(),
    message = function(m) {
        cat(conditionMessage(m), "\n")
        invokeRestart("muffleMessage")
    },
    warning = function(w) {
        cat("Warning:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
    }
)
