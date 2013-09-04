## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library('getopt');
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical",
  'outputfile' , 'o', 1, "character",
  'plots' , 'p', 2, "character",
  'input' , 'i', 1, "character",
  'samples', 's', 1, "character",
  'factors', 'f', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

htseqCountTable = read.table(opt$input, sep="\t", comment="", as.is=T)

colnames(htseqCountTable) <- htseqCountTable[2,]
names(htseqCountTable)<-colnames(htseqCountTable)
conditions <- htseqCountTable[1,]
conditions <- unlist(conditions[,-1])
tagnames <- htseqCountTable[-(1:2), 1]
htseqCountTable <- htseqCountTable[-(1:2),-1]
htseqCountTable[,1] <- gsub(",", ".", htseqCountTable[,1])

l <- unique(c(conditions))

library( "DESeq2" )

if ( !is.null(opt$plots) ) {
    pdf(opt$plots)
}
if (opt$samples=="all_vs_all"){
  # all versus all
  for (i in seq(1, length(l), by=1)) {
    k=i+1
    if(k<=length(l)){
      for (j in seq(k, length(l), by=1)) {
        currentColmuns <- which(conditions %in% c(l[[i]], l[[j]]))
            sampleNames <- names(htseqCountTable[currentColmuns])

        currentPairTable <- htseqCountTable[currentColmuns]
        ncondition1cols <- length(which(conditions %in% c(l[[i]])))
        ncondition2cols <- length(which(conditions %in% c(l[[i]])))

        ntabcols <- length(currentPairTable)
            currentPairTable <- as.integer(unlist(currentPairTable))
        currentPairTable <- matrix(unlist(currentPairTable), ncol = ntabcols, byrow = FALSE)
            colnames(currentPairTable) <- sampleNames

        pdata = data.frame(condition=factor(conditions[currentColmuns]),row.names=colnames(currentPairTable))

        dds = DESeqDataSetFromMatrix(countData = currentPairTable,  colData = pdata, design = ~ condition)

        dds <- DESeq(dds)
        sizeFactors(dds)
        if ( !is.null(opt$plots) ) {
            plotDispEsts(dds)
            plotMA(dds)
        }
        res <- results(dds)
        resCols <- colnames(res)

        cond <- c(rep(paste(unique(conditions[currentColmuns]),collapse="_"), nrow(res)))
        res[, "condition"] <- cond
        res[, "geneIds"] <- tagnames
        
        resSorted <- res[order(res$padj),]
        write.table(as.data.frame(resSorted[,c("condition", "geneIds", resCols)]), file = opt$outputfile, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}else{

    sampleSets <- unlist(strsplit(opt$samples, " "))
    samplesControl <- {}
    samplesExperiment <- {}
    samplesControl <- unlist(strsplit(sampleSets[1], ","))
    samplesExperiment <- unlist(strsplit(sampleSets[2], ","))

    # the minus one is needed because we get column indizes including the first column
    sampleColumns <- c()
    for (i in samplesControl) sampleColumns[[length(sampleColumns) + 1]] <- as.integer(i) - 1
    for (i in samplesExperiment) sampleColumns[[length(sampleColumns) + 1]] <- as.integer(i) - 1
    htseqSubCountTable <- htseqCountTable[,sampleColumns]

    ntabcols <- length(htseqSubCountTable)
    htseqSubCountTable <- as.integer(unlist(htseqSubCountTable))
    htseqSubCountTable <- matrix(unlist(htseqSubCountTable), ncol = ntabcols, byrow = FALSE)
    colnames(htseqSubCountTable) <- names(htseqCountTable)[sampleColumns]

    pdata = data.frame(condition=factor(conditions[sampleColumns]), row.names=colnames(htseqSubCountTable))

    if(!is.null(opt$factors)){
        factorSets <- unlist(strsplit(opt$factors, " "))
        for (factorSet in factorSets) {
            currentFactor <- unlist(strsplit(factorSet, ":"))
            for (factor in currentFactor[2]) {
                factoredCols <- unlist(strsplit(factor, ","))
                factoredCols <- as.numeric(factoredCols)
                factors <- c()
                for (i in sampleColumns){
                    j=i+1
                    if(j %in% factoredCols)
                        factors[[length(factors) + 1]] <- "yes"
                    else
                        factors[[length(factors) + 1]] <- "no"
                }
                # only add factor if factor is applied to the selected samples
                if((length(unique(factors))) >= 2){
                        pdata[[currentFactor[1]]]<-factor(factors)
                }else{
                    write("Input-Error 01: You can only apply factors to selected samples.", stderr())
                }
            }
        }
    }

    form <- as.formula(paste("", paste(names(pdata), collapse=" + "), sep=" ~ "))
    dds = DESeqDataSetFromMatrix(countData = htseqSubCountTable,  colData = pdata, design = form)

    dds <- DESeq(dds)
    sizeFactors(dds)
    if ( !is.null(opt$plots) ) {
        plotDispEsts(dds)
        plotMA(dds)
    }

    res <- results(dds)
    resCols <- colnames(res)

    cond <- c(rep(paste(unique(conditions[sampleColumns]),collapse="_"), nrow(res)))
    res[, "condition"] <- cond
    res[, "geneIds"] <- tagnames

    resSorted <- res[order(res$padj),]
    write.table(as.data.frame(resSorted[,c("condition", "geneIds", resCols)]), file = opt$outputfile, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
dev.off()

