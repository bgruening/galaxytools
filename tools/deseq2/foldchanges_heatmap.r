## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

args <- commandArgs(trailingOnly = TRUE)

library("gplots")
library("RColorBrewer")

tables <- {}
labels <- {}

labels <- unlist(strsplit(args[1], ","))
tables <- unlist(strsplit(args[2], ","))
ids_file=args[3]
output_pdf_file=args[4]

reds <- colorRampPalette(brewer.pal(9,"Reds"))(100)
reds <- rev(reds)
greens <- colorRampPalette(brewer.pal(9,"Greens"))(100)
cols <- c(reds, greens)

ids <- read.table(ids_file)

mat <- c()

for(i in 1:length(tables)) {
 # get data frame
    curr_table <- read.table(tables[[i]], sep="\t", header=FALSE)
    log2FCvect <- c()
    for(j in 1:length(ids$V1)) {
        log2FCvect <- c(log2FCvect, curr_table$V3[which(curr_table$V1 %in% ids$V1[j])])
    }
    # build foldChange data frame for heatmap
    mat <- cbind(mat, log2FCvect)
}
pdf(output_pdf_file) 
hm <- heatmap.2(mat, col = cols, trace="none", labCol=labels, labRow=NULL,scale="none",symm=F,symkey=T,symbreaks=T, cexCol=0.8, cexRow=0.8)
dev.off()  




