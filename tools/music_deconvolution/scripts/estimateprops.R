suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to estimate the effectiveness of proportion methods

## Load Conf
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

## Estimate cell type proportions
est_prop <- music_prop(
    bulk.eset = bulk_eset, sc.eset = scrna_eset,
    clusters = celltypes_label,
    samples = samples_label, select.ct = celltypes, verbose = T)


estimated_music_props <- est_prop$Est.prop.weighted
estimated_nnls_props <- est_prop$Est.prop.allgene
##
estimated_music_props_flat <- melt(estimated_music_props)
estimated_nnls_props_flat <- melt(estimated_nnls_props)

scale_yaxes <- function(gplot, value) {
    if (is.na(value)) {
        gplot
    } else {
        gplot + scale_y_continuous(lim = c(0, value))
    }
}

sieve_data <- function(func, music_data, nnls_data) {
    if (func == "list") {
        res <- list(if ("MuSiC" %in% methods) music_data else NULL,
                    if ("NNLS" %in% methods) nnls_data else NULL)
        res[lengths(res) > 0] ## filter out NULL elements
    } else if (func == "rbind") {
        rbind(if ("MuSiC" %in% methods) music_data else NULL,
              if ("NNLS" %in% methods) nnls_data else NULL)
    } else if (func == "c") {
        c(if ("MuSiC" %in% methods) music_data else NULL,
          if ("NNLS" %in% methods) nnls_data else NULL)
    }
}


## Show different in estimation methods
## Jitter plot of estimated cell type proportions
jitter_fig <- scale_yaxes(Jitter_Est(
    sieve_data("list",
               data.matrix(estimated_music_props),
               data.matrix(estimated_nnls_props)),
    method.name = methods, title = "Jitter plot of Est Proportions",
    size = 2, alpha = 0.7) + theme_minimal(), maxyscale)

## Make a Plot
## A more sophisticated jitter plot is provided as below. We separated
## the T2D subjects and normal subjects by their disease factor levels.
m_prop <- sieve_data("rbind",
                     estimated_music_props_flat,
                     estimated_nnls_props_flat)
colnames(m_prop) <- c("Sub", "CellType", "Prop")

if (is.null(celltypes)) {
    celltypes <- levels(m_prop$CellType)
    message("No celltypes declared, using:")
    message(celltypes)
}

if (is.null(phenotype_factors)) {
    phenotype_factors <- colnames(pData(bulk_eset))
}
## filter out unwanted factors like "sampleID" and "subjectName"
phenotype_factors <- phenotype_factors[
    !(phenotype_factors %in% phenotype_factors_always_exclude)]
message("Phenotype Factors to use:")
message(paste0(phenotype_factors, collapse = ", "))

m_prop$CellType <- factor(m_prop$CellType, levels = celltypes) # nolint
m_prop$Method <- factor(rep(methods, each = nrow(estimated_music_props_flat)), # nolint
                        levels = methods)

if (use_disease_factor) {

    if (phenotype_target_threshold == -99) {
        phenotype_target_threshold <- -Inf
        message("phenotype target threshold set to -Inf")
    }
    ## the "2" here is to do with the sample groups, not number of methods
    m_prop$Disease_factor <- rep(bulk_eset[[phenotype_target]], 2 * length(celltypes)) # nolint
    m_prop <- m_prop[!is.na(m_prop$Disease_factor), ]
    ## Generate a TRUE/FALSE table of Normal == 1 and Disease == 2
    sample_groups <- c("Normal", sample_disease_group)
    m_prop$Disease <- factor(sample_groups[(m_prop$Disease_factor > phenotype_target_threshold) + 1], # nolint
                             levels = sample_groups)

    ## Binary to scale: e.g. TRUE / 5 = 0.2
    m_prop$D <- (m_prop$Disease ==   # nolint
                 sample_disease_group) / sample_disease_group_scale
    ## NA's are not included in the comparison below
    m_prop <- rbind(subset(m_prop, Disease != sample_disease_group),
                    subset(m_prop, Disease == sample_disease_group))

    jitter_new <- scale_yaxes(
        ggplot(m_prop, aes(Method, Prop)) +
        geom_point(aes(fill = Method, color = Disease,
                       stroke = D, shape = Disease),
                   size = 2, alpha = 0.7,
                   position = position_jitter(width = 0.25, height = 0)) +
        facet_wrap(~ CellType, scales = "free") +
        scale_colour_manual(values = c("white", "gray20")) +
        scale_shape_manual(values = c(21, 24)) + theme_minimal(), maxyscale)

}

if (use_disease_factor) {

    ## Plot to compare method effectiveness
    ## Create dataframe for beta cell proportions and Disease_factor levels
    ## - Ugly code. Essentially, doubles the cell type proportions for each
    ##   set of MuSiC and NNLS methods
    m_prop_ana <- data.frame(
        pData(bulk_eset)[rep(1:nrow(estimated_music_props), length(methods)), #nolint
                         phenotype_factors],
        ## get proportions of target cell type
        ct.prop = sieve_data("c",
                             estimated_music_props[, phenotype_scrna_target],
                             estimated_nnls_props[, phenotype_scrna_target]),
        ##
        Method = factor(rep(methods,
                            each = nrow(estimated_music_props)),
                        levels = methods))
    ## - fix headers
    colnames(m_prop_ana)[1:length(phenotype_factors)] <- phenotype_factors #nolint
    ## - drop NA for target phenotype (e.g. hba1c)
    m_prop_ana <- subset(m_prop_ana, !is.na(m_prop_ana[phenotype_target]))
    m_prop_ana$Disease <- factor(   # nolint
        ## - Here we set Normal/Disease assignments across the methods
        sample_groups[(
            m_prop_ana[phenotype_target] > phenotype_target_threshold) + 1
            ],
        sample_groups)
    ## - Then we scale this binary assignment to a plotable factor
    m_prop_ana$D <- (m_prop_ana$Disease ==        # nolint
                     sample_disease_group) / sample_disease_group_scale

    jitt_compare <- scale_yaxes(
        ggplot(m_prop_ana, aes_string(phenotype_target, "ct.prop")) +
        geom_smooth(method = "lm",  se = FALSE, col = "black", lwd = 0.25) +
        geom_point(aes(fill = Method, color = Disease,
                       stroke = D, shape = Disease),
                   size = 2, alpha = 0.7) +  facet_wrap(~ Method) +
        ggtitle(paste0(toupper(phenotype_target), " vs. ",
                       toupper(phenotype_scrna_target),
                       " Cell Type Proportion")) +
        theme_minimal() +
        ylab(paste0("Proportion of ",
                    phenotype_scrna_target, " cells")) +
        xlab(paste0("Level of bulk factor (", phenotype_target, ")")) +
        scale_colour_manual(values = c("white", "gray20")) +
        scale_shape_manual(values = c(21, 24)), maxyscale)
}

## BoxPlot
plot_box <- scale_yaxes(Boxplot_Est(
    sieve_data("list",
               data.matrix(estimated_music_props),
               data.matrix(estimated_nnls_props)),
    method.name = methods) +
    theme(axis.text.x = element_text(angle = -90),
          axis.text.y = element_text(size = 8)) +
    ggtitle(element_blank()) + theme_minimal(), maxyscale)

## Heatmap
plot_hmap <- Prop_heat_Est(
    sieve_data(
        "list",
        data.matrix(estimated_music_props),
        data.matrix(estimated_nnls_props)),
    method.name = methods) +
    theme(axis.text.x = element_text(angle = -90),
          axis.text.y = element_text(size = 6))

pdf(file = outfile_pdf, width = 8, height = 8)
if (length(celltypes) <= 8) {
    plot_grid(jitter_fig, plot_box, labels = "auto", ncol = 1, nrow = 2)
} else {
    print(jitter_fig)
    plot_box
}
if (use_disease_factor) {
    plot_grid(jitter_new, jitt_compare, labels = "auto", ncol = 1, nrow = 2)
}
plot_hmap
message(dev.off())

writable <- function(obj, prefix, title) {
    write.table(obj,
                file = paste0("report_data/", prefix, "_",
                              title, ".tabular"),
                quote = F, sep = "\t", col.names = NA)
}

## Output Proportions
if ("NNLS" %in% methods) {
    writable(est_prop$Est.prop.allgene, "prop",
             "NNLS Estimated Proportions of Cell Types")
}

if ("MuSiC" %in% methods) {
    writable(est_prop$Est.prop.weighted, "prop",
             "Music Estimated Proportions of Cell Types")
    writable(est_prop$Weight.gene, "weightgene",
             "Music Estimated Proportions of Cell Types (by Gene)")
    writable(est_prop$r.squared.full, "rsquared",
             "Music R-sqr Estimated Proportions of Each Subject")
    writable(est_prop$Var.prop, "varprop",
             "Matrix of Variance of MuSiC Estimates")
}


if (use_disease_factor) {
    ## Summary table of linear regressions of disease factors
    for (meth in methods) {
        ##lm_beta_meth = lm(ct.prop ~ age + bmi + hba1c + gender, data =
        sub_data <- subset(m_prop_ana, Method == meth)

        ## We can only do regression where there are more than 1 factors
        ## so we must find and exclude the ones which are not
        gt1_facts <- sapply(phenotype_factors, function(facname) {
            return(length(unique(sort(sub_data[[facname]]))) == 1)
        })
        form_factors <- phenotype_factors
        exclude_facts <- names(gt1_facts)[gt1_facts]
        if (length(exclude_facts) > 0) {
            message("Factors with only one level will be excluded:")
            message(exclude_facts)
            form_factors <- phenotype_factors[
                !(phenotype_factors %in% exclude_facts)]
        }
        lm_beta_meth <- lm(as.formula(
            paste("ct.prop", paste(form_factors, collapse = " + "),
                  sep = " ~ ")), data = sub_data)
        message(paste0("Summary: ", meth))
        capture.output(summary(lm_beta_meth),
                       file = paste0("report_data/summ_Log of ",
                                     meth,
                                     " fitting.txt"))
    }
}
