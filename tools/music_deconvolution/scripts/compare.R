suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to estimate the effectiveness of proportion methods

## Load Conf
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

method_key <- list("MuSiC" = "est_music",
                   "NNLS" = "est_nnls")[[est_method]]


scale_yaxes <- function(gplot, value) {
    if (is.na(value)) {
        gplot
    } else {
        gplot + scale_y_continuous(lim = c(0, value))
    }
}


setFactorData <- function(bulk_data, factor_name = NULL){
    if (is.null(factor_name)){
        factor_name = "None" ## change to something plottable
    }
    pdat <- pData(bulk_data)
    sam_fact <- NULL
    if (factor_name %in% colnames(pdat)){
        sam_fact = cbind(rownames(pdat),
                         as.character(pdat[[factor_name]]))
        cat(paste0("     - factor: ", factor_name, " found in phenotypes\n"))
    } else {
        ## We assign this as the factor for the entire dataset
        sam_fact = cbind(rownames(pdat),
                         factor_name)
        cat(paste0("     - factor: assigning \"", factor_name, "\" to whole dataset\n"))
    }
    colnames(sam_fact) = c("Samples", "Factors")
    return(as.data.frame(sam_fact))
}

## Due to limiting sizes, we need to load and unload
## possibly very large datasets.
processPair <- function(sc_data, bulk_data,
                        ctypes_label, samples_label, ctypes,
                        factor_group){
    ## - Generate
    est_prop <- music_prop(
        bulk.eset = bulk_data, sc.eset = sc_data,
        clusters = ctypes_label,
        samples = samples_label, select.ct = ctypes, verbose = T)
    ## -
    estimated_music_props <- est_prop$Est.prop.weighted
    estimated_nnls_props <- est_prop$Est.prop.allgene
    ## -
    fact_data <- setFactorData(bulk_data, factor_group)
    ## -
    return(list(est_music = estimated_music_props,
                est_nnls = estimated_nnls_props,
                bulk_sample_totals = colSums(exprs(bulk_data)),
                plot_groups = fact_data))
}

musicOnAll <- function (files){
    results <- list()
    for (sc_name in names(files)){
        cat(paste0("sc-group:", sc_name, "\n"))
        scgroup = files[[sc_name]]
        ## - sc Data
        sc_est = readRDS(scgroup$dataset)
        ## - params
        celltypes_label = scgroup$label_cell
        samples_label = scgroup$label_sample
        celltypes = scgroup$celltype

        results[[sc_name]] = list()
        for (bulk_name in names(scgroup$bulk)){
            cat(paste0(" - bulk-group:", bulk_name, "\n"))
            bulkgroup <- scgroup$bulk[[bulk_name]]
            ## - bulk Data
            bulk_est <- readRDS(bulkgroup$dataset)
            ## - bulk params
            pheno_facts <- bulkgroup$pheno_facts
            pheno_excl <- bulkgroup$pheno_excl
            ##
            results[[sc_name]][[bulk_name]] <- processPair(
                sc_est, bulk_est,
                celltypes_label, samples_label,
                celltypes, bulkgroup$factor_group)
            ##
            rm(bulk_est) ## unload
        }
        rm(sc_est) ## unload
    }
    return(results)
}

plotAllIndividualHeatmaps <- function(results){
    pdf(out_heatmulti_pdf, width=8, height=8)
    for (sc_name in names(results)){
        for (bk_name in names(results[[sc_name]])){
            res <- results[[sc_name]][[bk_name]]
            plot_hmap <- Prop_heat_Est(
                data.matrix(res[[method_key]]),
                method.name = est_method) +
                ggtitle(paste0("[", est_method, "Cell type ",
                               "proportions in ",
                               bk_name, " (Bulk) based on ",
                               sc_name, " (scRNA)")) +
                xlab("Cell Types (scRNA)") +
                ylab("Samples (Bulk)") +
                theme(axis.text.x = element_text(angle = -90),
                      axis.text.y = element_text(size = 6))
            print(plot_hmap)
        }
    }
    dev.off()
}

mergeFactorsSpread <- function(grudat_spread, factor_groups){
    ## Generated
    mergeIt <- function(matr, plot_groups, valname){
        ren <- melt(lapply(matr, function(mat){mat["ct"]=rownames(mat); return(mat)}))
        ## - Grab factors and merge into list
        ren_new <- merge(ren, plot_groups, by.x="variable", by.y="Samples")
        colnames(ren_new) <- c("Sample", "Cell", valname, "Bulk", "Factors")
        return(ren_new)
    }
    tab <- merge(mergeIt(grudat$spread$prop, factor_groups, "value.prop"),
                 mergeIt(grudat$spread$scale, factor_groups, "value.scale"),
                 by=c("Sample", "Cell", "Bulk", "Factors"))
    return(tab)
}


plotGroupedHeatmaps <- function(results){
    pdf(out_heatmulti_pdf, width=8, height=8)
    for (sc_name in names(results)){
        named_list <- sapply(
            names(results[[sc_name]]),
            function (n){
                ## We transpose the data here, because
                ## the plotting function omits by default
                ## the Y-axis which are the samples.
                ##   Since the celltypes are the common factor
                ## these should be the Y-axis instead.
                t(data.matrix(results[[sc_name]][[n]][[method_key]]))
            }, simplify=F, USE.NAMES=T)
        named_methods <- names(results[[sc_name]])
        ##
        plot_hmap <- Prop_heat_Est(
            named_list,
            method.name = named_methods) +
            ggtitle(paste0("[", est_method, "] Cell type ",
                           "proportions of ",
                           "Bulk Datasets based on ",
                           sc_name, " (scRNA)")) +
            xlab("Samples (Bulk)") +
            ylab("Cell Types (scRNA)") +
            theme(axis.text.x = element_text(angle = -90),
                  axis.text.y = element_text(size = 6))
        print(plot_hmap)
    }
    dev.off()
}

## Desired plots
## 1. Pie chart:
##    - Per Bulk dataset (using just normalised proportions)
##    - Per Bulk dataset (multiplying proportions by nreads)

unlistNames <- function(results, method, prepend_bkname=FALSE){
    unique(sort(
        unlist(lapply(names(results), function (scname) {
            lapply(names(results[[scname]]), function (bkname) {
                res <- get(method)(results[[scname]][[bkname]][[method_key]])
                if (prepend_bkname){
                    ## We *do not* assume unique bulk sample names
                    ## across different bulk datasets.
                    res <- paste0(bkname, "::", res)
                }
                return(res)
            })
        }))
    ))
}

## convertProportionsToCounts <- function(prop_matrix,

summarizedMatrix <- function(results){
    ## We assume that cell types MUST be unique, but that sample
    ## names do not need to be. For this reason, we must prepend
    ## the bulk dataset name to the individual sample names.
    all_celltypes <- unlistNames(results, "colnames")
    all_samples <- unlistNames(results, "rownames", prepend_bkname=TRUE)

    ## Iterate through all possible samples and populate a table.
    ddff <- data.frame()
    ddff_scale <- data.frame()
    for (cell in all_celltypes){
        for (sample in all_samples){
            group_sname <- unlist(strsplit(sample, split="::"))
            bulk <- group_sname[1]
            id_sample <- group_sname[2]
            for (scgroup in names(results)){
                if (bulk %in% names(results[[scgroup]])){
                    mat_prop <- results[[scgroup]][[bulk]][[method_key]]
                    vec_counts <- results[[scgroup]][[bulk]]$bulk_sample_totals
                    ## - We use sample instead of id_sample because we need to extract
                    ##   bulk sets from the complete matrix later. It's messy, yes.
                    ddff[cell, sample] <- mat_prop[id_sample,cell]
                    ddff_scale[cell, sample] <- mat_prop[id_sample,cell] * vec_counts[[id_sample]]
                }
            }
        }
    }
    return(list(prop=ddff, scaled=ddff_scale))
}

flattenFactorList <- function(results){
    ## Get a 2d DF of all factors across all bulk samples.
    res <- c()
    for (scgroup in names(results)){
        for (bulkgroup in names(results[[scgroup]])){
            dat = results[[scgroup]][[bulkgroup]]$plot_groups
            dat$Samples = paste0(bulkgroup, "::", dat$Samples)
            res <- rbind(res, dat)
        }
    }
    return(res)
}

groupByDataset <- function(summat){
    bulk_names = unlist(
        lapply(names(files),
               function(x) names(files[[x]]$bulk)))
    mat_names = colnames(summat$prop)
    bd <- list()
    bd_scale <- list()
    bd_spread_scale <- list()
    bd_spread_prop <- list()
    for (bname in bulk_names){
        subs <- mat_names[startsWith(mat_names, paste0(bname, "::"))]
        ##print(bname)
        ## -
        bd[[bname]] = rowSums(summat$prop[,subs])
        bd_scale[[bname]] = rowSums(summat$scaled[,subs])
        bd_spread_scale[[bname]] = summat$scaled[,subs]
        bd_spread_prop[[bname]] = summat$prop[,subs]
    }
    return(list(prop=as.data.frame(bd),
                scaled=as.data.frame(bd_scale),
                spread=list(scale=bd_spread_scale,
                            prop=bd_spread_prop)))
}

summarizeHeatmaps <- function(grudat_spread_melt, do.factors){
    # -
    doSingle <- function(grudat_melted, yaxis, xaxis, fillval, title,
                         ylabs=element_blank(), xlabs=element_blank(), USE.LOG=TRUE){
        ## Convert from matrix to long format
        melted = grudat_melted  ## copy?
        if (USE.LOG){
            melted[[fillval]] = log10(melted[[fillval]] + 1)
        }
        return(ggplot(melted) +
               geom_tile(aes_string(y=yaxis, x=xaxis, fill=fillval), colour="white") +
               scale_fill_gradient2(low="steelblue", high="red", mid="white", name=element_blank()) +
               theme(axis.text.x = element_text(angle=-50, hjust=0)) +
               ggtitle(label=title) + xlab(xlabs) + ylab(ylabs))
    }

    doGridPlot <- function(title, xvar, plot="both", ncol=2){
        do.logged = (plot %in% c("log", "both"))
        do.normal = (plot %in% c("normal", "both"))
        plist = list()
        if (do.logged){
            plist[["1"]] = doSingle(grudat_spread_melt, "Cell", xvar, "value.scale", "Reads (log10+1)")
            plist[["2"]] = doSingle(grudat_spread_melt, "Cell", xvar, "value.prop", "Sample (log10+1)")
        }
        if (do.normal){
            plist[["A"]] = doSingle(grudat_spread_melt, "Cell", xvar, "value.scale", "Reads", USE.LOG=F)
            plist[["B"]] = doSingle(grudat_spread_melt, "Cell", xvar, "value.prop", "Sample", USE.LOG=F)
        }
        return(plot_grid(ggdraw() + draw_label(title, fontface="bold"),
                         plot_grid(plotlist=plist, ncol=ncol),
                         ncol=1, rel_heights=c(0.05,0.95)))

    }
    p1 <- doGridPlot("Cell Types against Bulk Datasets", "Bulk", "both")
    p2a <- doGridPlot("Cell Types against Samples", "Sample", "normal", 1)
    p2b <- doGridPlot("Cell Types against Samples (log10+1)", "Sample", "log", 1)
    p3 <- ggplot + theme_void()
    if (do.factors){
        p3 <- doGridPlot("Cell Types against Factors", "Factors", "both")
    }
    return(list(bulk=p1, samples=list(log=p2b, normal=p2a), factors=p3))
}

summarizeBoxPlots <- function(grudat_spread, do.factors){
    common1 <- ggplot(grudat_spread, aes(x=value.prop)) + ggtitle("Sample") +
        xlab(element_blank()) + ylab(element_blank())
    common2 <- ggplot(grudat_spread, aes(x=value.scale)) + ggtitle("Reads") +
        xlab(element_blank()) + ylab(element_blank())

    A=B=list()
    ## Cell type by sample
    A$p1 = common2 + geom_boxplot(aes(y=Cell, color=Bulk))
    A$p2 = common1 + geom_boxplot(aes(y=Cell, color=Bulk))
    ## Sample by Cell type
    B$p1 = common1 + geom_boxplot(aes(y=Bulk, color=Cell))  + ylab("Bulk Dataset")
    B$p2 = common2 + geom_boxplot(aes(y=Bulk, color=Cell))  + ylab("Bulk Dataset")
    ## -- Factor plots are optional
    A$p3 = B$p3 = A$p4 = B$p4 = ggplot() + theme_void()

    if (do.factors){
        A$p3 = common1 + geom_boxplot(aes(y=Cell, color=Factors))
        A$p4 = common2 + geom_boxplot(aes(y=Cell, color=Factors))
        B$p3 = common2 + geom_boxplot(aes(y=Bulk, color=Factors)) + ylab("Bulk Dataset")
        B$p4 = common1 + geom_boxplot(aes(y=Bulk, color=Factors)) + ylab("Bulk Dataset")
    }

    titleA = "Cell Types against Bulk"
    titleB = "Bulk Datasets against Cells"
    if (do.factors){
        titleA = paste0(titleA, " and Factors")
        titleb = paste0(titleB, " and Factors")
    }

    A_all = plot_grid(ggdraw() + draw_label(titleA, fontface="bold"),
                      plot_grid(plotlist=A, ncol=2),
                      ncol=1, rel_heights=c(0.05,0.95))
    B_all = plot_grid(ggdraw() + draw_label(titleB, fontface="bold"),
                      plot_grid(plotlist=B, ncol=2),
                      ncol=1, rel_heights=c(0.05,0.95))
    return(list(cell=A_all, bulk=B_all))
}

filterOutput <- function(grudat_spread_melt, out_filt){
    printRed <- function(comment, red_list){
        cat(paste(comment, paste(red_list, collapse=", "), "\n"))
    }

    grudat_filt <- grudat_spread_melt
    printRed("Total Cell types:", unique(grudat_filt$Cell))
    if (!is.null(out_filt$cells)){
        grudat_filt <- grudat_filt[grudat_filt$Cell %in% out_filt$cells,]
        printRed(" - selecting:", out_filt$cells)
    }
    printRed("Total Factors:", unique(grudat_spread_melt$Factors))
    if (!is.null(out_filt$facts)){
        grudat_filt <- grudat_filt[grudat_filt$Factors %in% out_filt$facts,]
        printRed(" - selecting:", out_filt$facts)
    }
    return(grudat_filt)
}


results <- musicOnAll(files)

if (heat_grouped_p) {
    plotGroupedHeatmaps(results)
} else {
    plotAllIndividualHeatmaps(results)
}

summat = summarizedMatrix(results)
grudat = groupByDataset(summat)
grudat_spread_melt = mergeFactorsSpread(grudat$spread, flattenFactorList(results))

##save.image(file="/tmp/sesh.RData")

## The output filters ONLY apply to boxplots, since these take
do.factors = (length(unique(grudat_spread_melt[["Factors"]])) > 1)

grudat_spread_melt_filt <- filterOutput(grudat_spread_melt, out_filt)

heat_maps = summarizeHeatmaps(grudat_spread_melt_filt, do.factors)
box_plots = summarizeBoxPlots(grudat_spread_melt_filt, do.factors)

pdf(out_heatsumm_pdf, width=14, height=14)
print(heat_maps)
print(box_plots)
dev.off()

## -- DEBUG --
## saveRDS(files, "/tmp/files.rds")
## saveRDS(results, "/tmp/results.rds")
## saveRDS(summat, "/tmp/summat.rds")
##saveRDS(grudat, "/tmp/grudat.rds")
##saveRDS(grudat_mod, "/tmp/grudatmod.rds")
##files = readRDS("/tmp/files.rds")
## results = readRDS("/tmp/results.rds")
## summat = readRDS("/tmp/summat.rds")
## grudat = readRDS("/tmp/grudat.rds")
