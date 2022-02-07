suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to estimate the effectiveness of proportion methods

## Load Conf
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

scale_yaxes <- function(gplot, value) {
    if (is.na(value)) {
        gplot
    } else {
        gplot + scale_y_continuous(lim = c(0, value))
    }
}

method_key <- list("MuSiC" = "est_music",
                   "NNLS" = "est_nnls")[[est_method]]


## Due to limiting sizes, we need to load and unload
## possibly very large datasets.
processPair <- function(sc_data, bulk_data,
                        ctypes_label, samples_label, ctypes){
    ## - Generate
    est_prop <- music_prop(
        bulk.eset = bulk_data, sc.eset = sc_data,
        clusters = ctypes_label,
        samples = samples_label, select.ct = ctypes, verbose = T)
    ## -
    estimated_music_props <- est_prop$Est.prop.weighted
    estimated_nnls_props <- est_prop$Est.prop.allgene
    ## -
    ##estimated_music_props_flat <- melt(estimated_music_props)
    ##estimated_nnls_props_flat <- melt(estimated_nnls_props)
    ## -
    ##saveRDS(bulk_data, "/
    return(list(est_music = estimated_music_props,
                est_nnls = estimated_nnls_props,
                bulk_sample_totals = colSums(exprs(bulk_data))))
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
            res <- processPair(sc_est, bulk_est,
                               celltypes_label, samples_label,
                               celltypes)
            results[[sc_name]][[bulk_name]] <- res
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
                    ## We do not assume unique bulk sample names
                    ## across different samples.
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
                    ddff[cell, sample] <- mat_prop[id_sample,cell]
                    ddff_scale[cell, sample] <- mat_prop[id_sample,cell] * vec_counts[[id_sample]]
                }
            }
        }
    }
    return(list(prop=ddff, scaled=ddff_scale))
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

makeHeatmapByGroup_single <- function(dataset, title, USE.LOG=TRUE){
    ## Convert from matrix to long format
    dataset["CT"] = rownames(dataset)
    melted = melt(dataset, value.name="VALS", variable.name="Bulk")

    if (USE.LOG){
        title = paste0("[Log10+1] ", title)
        melted$VALS = log10(melted$VALS + 1)
    }
    
    return(ggplot(melted) +
           geom_tile(aes(y=CT, x=Bulk,
                         fill=VALS), colour="white") +
           scale_fill_gradient2(low="steelblue", high="red", mid="white", name=element_blank()) + 
           theme(axis.text.x = element_text(angle=-90)) +
           ggtitle(title) +
           ylab("Cell types across all scRNA Datasets") +
           xlab("All Bulk Datasets"))
}

summarizeHeatmapsByGroup <- function(grudat){
    pheat_scale <- makeHeatmapByGroup_single(
        grudat$scaled,
        "Cell Types Proportions (Scaled to Number of Reads)")
    pheat_prop <- makeHeatmapByGroup_single(
        grudat$prop,
        "Cell Types Proportions (Normalised by Sample)")
    
    pdf(out_heatsumm_pdf, width=8, height=8)
    print(pheat_scale)
    print(pheat_prop)
    dev.off()
}


re = melt(lapply(grudat2$spread$scale, function(mat){mat["ct"]=rownames(mat); return(mat)}))
## Cell type by sample
## Sample by Cell type
ggplot(re.prop) + geom_boxplot(aes(y=ct, x=value, color=L1))
ggplot(re.prop) + geom_boxplot(aes(y=L1, x=value, color=ct))

results <- musicOnAll(files)
if (heat_grouped_p) {
    plotGroupedHeatmaps(results)
} else {
    plotAllIndividualHeatmaps(results)
}
summat = summarizedMatrix(results)
grudat = groupByDataset(summat)
summarizeHeatmapsByGroup(grudat)



## saveRDS(files, "/tmp/files.rds")
## saveRDS(results, "/tmp/results.rds")
## saveRDS(summat, "/tmp/summat.rds")
## saveRDS(grudat, "/tmp/grudat.rds")
## files = readRDS("/tmp/files.rds")
## results = readRDS("/tmp/results.rds")
## summat = readRDS("/tmp/summat.rds")
## grudat = readRDS("/tmp/grudat.rds")
