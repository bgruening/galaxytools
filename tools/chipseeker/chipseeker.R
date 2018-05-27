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
    make_option(c("-t","--genome"), type="character", help="Source of reference genome to create TxDb."),
    make_option(c("-p","--plots"), type="character", help="PDF of plots.")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

peaks = args$infile
genome = args$genome
plots = args$plots

peaks <- readPeakFile(peaks)

if (genome %in% c("hg38","hg19","mm10")) {
    # Use TxDb from UCSC
    if (genome == "hg38") {
        suppressPackageStartupMessages({
            library(TxDb.Hsapiens.UCSC.hg38.knownGene)
            library(org.Hs.eg.db)
        })
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        annodb <- "org.Hs.eg.db"
    } else if (genome == "hg19") {
        suppressPackageStartupMessages({
            library(TxDb.Hsapiens.UCSC.hg19.knownGene)
            library(org.Hs.eg.db)
        })
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        annodb <- "org.Hs.eg.db"
    } else if (genome == "mm10") {
        suppressPackageStartupMessages({
            library(TxDb.Mmusculus.UCSC.mm10.knownGene)
            library(org.Mm.eg.db)
        })
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        annodb <- "org.Mm.eg.db"
    }
    
    peakAnno <-  annotatePeak(peaks, TxDb=txdb, annoDb=annodb)

} else {
    # Make TxDb from GTF
    txdb <- makeTxDbFromGFF(genome, format="gtf")
    peakAnno <-  annotatePeak(peaks, TxDb=txdb)

}

# Convert from 1-based to Interval 0-based format
res <- as.GRanges(peakAnno)
metacols <- mcols(res)
metacols <- apply(as.data.frame(metacols), 1, function(col) paste(col, collapse="|"))
resout  <- data.frame(Chrom=seqnames(res),
                Start=start(res) - 1,
                End=end(res),
                Comment=metacols)
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