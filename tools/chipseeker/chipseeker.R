options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
    library(optparse)
})

option_list <- list(
    make_option(c("-i","--infile"), type="character", help="Peaks file to be annotated"),
    make_option(c("-G","--gtf"), type="character", help="GTF to create TxDb."),
    make_option(c("-f","--format"), type="character", help="Output format (interval or tabular)."),
    make_option(c("-p","--plots"), type="character", help="PDF of plots.")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

peaks = args$infile
gtf = args$gtf
format = args$format
plots = args$plots

peaks <- readPeakFile(peaks)

# Make TxDb from GTF
txdb <- makeTxDbFromGFF(gtf, format="gtf")
peakAnno <-  annotatePeak(peaks, TxDb=txdb)

# Convert from 1-based to 0-based format
res <- as.GRanges(peakAnno)
metacols <- mcols(res)
if (format == "interval") {
    metacols <- apply(as.data.frame(metacols), 1, function(col) paste(col, collapse="|"))
    resout  <- data.frame(Chrom=seqnames(res),
                    Start=start(res) - 1,
                    End=end(res),
                    Comment=metacols)
} else {
    resout <- data.frame(Chrom=seqnames(res),
                    Start=start(res) - 1,
                    End=end(res),
                    metacols)
}

write.table(resout, file="out.tab", sep="\t", row.names=FALSE, quote=FALSE)

if (!is.null(plots)) {
    pdf("out.pdf", width=14)
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
    dev.off()
}