options(
  show.error.messages = F,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
  }
)

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(rtracklayer)
  library(optparse)
  library(ggupset)
})

option_list <- list(
  make_option(c("-i", "--infile"), type = "character", help = "Peaks file to be annotated"),
  make_option(c("-H", "--header"), type = "logical", help = "Peaks file contains header row"),
  make_option(c("-G", "--gtf"), type = "character", help = "GTF to create TxDb."),
  make_option(c("-u", "--upstream"), type = "integer", help = "TSS upstream region"),
  make_option(c("-d", "--downstream"), type = "integer", help = "TSS downstream region"),
  make_option(c("-F", "--flankgeneinfo"), type = "logical", help = "Add flanking gene info"),
  make_option(c("-D", "--flankgenedist"), type = "integer", help = "Flanking gene distance"),
  make_option(c("-j", "--ignore_upstream"), type = "logical", help = "Ignore upstream"),
  make_option(
    c("-k", "--ignore_downstream"),
    type = "logical",
    help = "Ignore downstream"
  ),
  make_option(c("-f", "--format"), type = "character", help = "Output format (interval or tabular)."),
  make_option(c("-p", "--plots"), type = "logical", help = "PDF of plots."),
  make_option(c("-r", "--rdata"), type = "logical", help = "Output RData file.")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

peaks <- args$infile
gtf <- args$gtf
up <- args$upstream
down <- args$downstream
format <- args$format

if (!is.null(args$flankgeneinfo)) {
  flankgeneinfo <- TRUE
} else {
  flankgeneinfo <- FALSE
}

if (!is.null(args$ignore_upstream)) {
  ignore_upstream <- TRUE
} else {
  ignore_upstream <- FALSE
}

if (!is.null(args$ignore_downstream)) {
  ignore_downstream <- TRUE
} else {
  ignore_downstream <- FALSE
}

if (!is.null(args$header)) {
  header <- TRUE
} else {
  header <- FALSE
}

peaks <- readPeakFile(peaks, header = header)

# Make TxDb from GTF
txdb <- makeTxDbFromGFF(gtf, format = "gtf")

# Annotate peaks
peak_anno <-  annotatePeak(
  peaks,
  TxDb = txdb,
  tssRegion = c(-up, down),
  addFlankGeneInfo = flankgeneinfo,
  flankDistance = args$flankgenedist,
  ignoreUpstream = ignore_upstream,
  ignoreDownstream = ignore_downstream
)

# Add gene name
features <- import(gtf, format = "gtf")
ann <- unique(mcols(features)[, c("gene_id", "gene_name")])
res <- as.data.frame(peak_anno)
res <- merge(res, ann, by.x = "geneId", by.y = "gene_id")
names(res)[names(res) == "gene_name"] <- "geneName"

#Extract metadata cols, 1st is geneId, rest should be from col 7 to end
metacols <- res[, c(7:ncol(res), 1)]
# Convert from 1-based to 0-based format
if (format == "interval") {
  metacols <-
    apply(as.data.frame(metacols), 1, function(col)
      paste(col, collapse = "|"))
  resout <- data.frame(res$seqnames,
                       res$start - 1,
                       res$end,
                       metacols)
  colnames(resout)[1:4] <- c("Chrom", "Start", "End", "Comment")
} else {
  resout <- data.frame(res$seqnames,
                       res$start - 1,
                       res$end,
                       metacols)
  colnames(resout)[1:3] <- c("Chrom", "Start", "End")
}
write.table(
  resout,
  file = "out.tab",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (!is.null(args$plots)) {
  pdf("out.pdf", width = 14)
  plotAnnoPie(peak_anno)
  p1 <- plotAnnoBar(peak_anno)
  print(p1)
  vennpie(peak_anno)
  upsetplot(peak_anno)
  p2 <-
    plotDistToTSS(peak_anno, title = "Distribution of transcription factor-binding loci\nrelative to TSS")
  print(p2)
  dev.off()
  rm(p1, p2)
}

## Output RData file

if (!is.null(args$rdata)) {
  save.image(file = "ChIPseeker_analysis.RData")
}
