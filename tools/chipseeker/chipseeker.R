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
    make_option(c("-u","--upstream"), type="integer", help="TSS upstream region"),
    make_option(c("-d","--downstream"), type="integer", help="TSS downstream region"),
    make_option(c("-F","--flankgeneinfo"), type="logical", help="Add flanking gene info"),
    make_option(c("-D","--flankgenedist"), type="integer", help="Flanking gene distance"),
    make_option(c("-f","--format"), type="character", help="Output format (interval or tabular)."),
    make_option(c("-p","--plots"), type="logical", help="PDF of plots."),
    make_option(c("-r","--rdata"), type="logical", help="Output RData file.")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

peaks = args$infile
gtf = args$gtf
up = args$upstream
down = args$downstream
format = args$format

peaks <- readPeakFile(peaks)

# Make TxDb from GTF
txdb <- makeTxDbFromGFF(gtf, format="gtf")
if (!is.null(args$flankgeneinfo)) {
    peakAnno <-  annotatePeak(peaks, TxDb=txdb, tssRegion=c(-up, down), addFlankGeneInfo=args$flankgeneinfo, flankDistance=args$flankgenedist)
} else {
    peakAnno <-  annotatePeak(peaks, TxDb=txdb, tssRegion=c(-up, down))
}

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

if (!is.null(args$plots)) {
    pdf("out.pdf", width=14)
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
    dev.off()
}

## Output RData file

if (!is.null(args$rdata)) {
  save.image(file = "ChIPseeker_analysis.RData")
}