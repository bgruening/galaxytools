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
	'outfile' , 'o', 1, "character",
	'outfilefiltered' , 'f', 1, "character",
	'plots' , 'p', 2, "character",
	'input' , 'i', 1, "character",
	'factors', 'm', 2, "character",
	'fittype', 't', 2, "character",
	'threshold', 'c', 2, "double"
#	'organism', 'g', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}

if( is.null(opt$fittype))
	opt$fittype = "parametric"

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

htseqCountTable = read.table(opt$input, sep="\t", comment="", as.is=T)

colnames(htseqCountTable) <- htseqCountTable[2,]
names(htseqCountTable)<-colnames(htseqCountTable)
conditions <- htseqCountTable[1,]
conditions <- unlist(conditions[,-1])
featurenames <- htseqCountTable[-(1:2), 1]
htseqCountTable <- htseqCountTable[-(1:2),-1]
htseqCountTable[,1] <- gsub(",", ".", htseqCountTable[,1])

l <- unique(c(conditions))

library('rjson')
library('DESeq2')
#library('biomaRt')

if ( !is.null(opt$plots) ) {
	pdf(opt$plots)
}

## The following function takes deseq data object, computes, writes and plots the results.
computeAndWriteResults <- function(dds, sampleCols, outputcsv, featurenames_filtered) {
	dds <- DESeq(dds, fitType= opt$fittype)
	sizeFactors(dds)
	res <- results(dds)
	resCols <- colnames(res)
	cond <- c(rep(paste(unique(conditions[sampleCols]),collapse="_"), nrow(res)))
	if(missing(featurenames_filtered)){
		res[, "geneIds"] <- featurenames
#		rownames(res) <- featurenames
		title_prefix = "Complete: "
	}else{
		res[, "geneIds"] <- featurenames_filtered
#		rownames(res) <- featurenames
		title_prefix = "Filtered: "
	}
	print(sum(res$padj < .1, na.rm=TRUE))
	print(opt$organism)
#	if(opt$organism != "other"){
#		dataset = ""
#		if(opt$organism == "mouse")
#			dataset = "mmusculus_gene_ensembl"
#		if(opt$organism == "human")
#			dataset = "hsapiens_gene_ensembl"
#		if(opt$organism == "fly")
#			dataset = "dmelanogaster_gene_ensembl"
#		ensembldb = useMart("ensembl",dataset=dataset)
#
#		annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_id","description"),
#					filters = "ensembl_gene_id",
#					values=res[, "geneIds"],
#					mart=ensembldb)

#		res <- merge(res, annot,
#					by.x = "geneIds",
#					by.y = "ensembl_gene_id",
#					all.x=TRUE)
		
#		resCols <- colnames(res)
#		resSorted <- res[order(res$padj),]
#		write.table(as.data.frame(resSorted[,c(resCols)]), file = outputcsv, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE, col.names = FALSE)					
#	}else{
		resSorted <- res[order(res$padj),]
		write.table(as.data.frame(resSorted[,c("geneIds", resCols)]), file = outputcsv, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE, col.names = FALSE)
#	}
	
	if ( !is.null(opt$plots) ) {
		plotDispEsts(dds, main= paste(title_prefix, "Dispersion estimate plot"))
		plotMA(dds, main= paste(title_prefix, "MA-plot"))

		library("RColorBrewer")
		library("gplots")
		select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
		hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

		rld <- rlogTransformation(dds, blind=TRUE)
		vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
		distsRL <- dist(t(assay(rld)))
		mat <- as.matrix(distsRL)
		heatmap.2(mat, trace="none", col = rev(hmcol), main = paste(title_prefix, "Sample-to-sample distances"))
	}
	return(res)
}

## This functions filters out the genes with low counts and plots p-value distribution
independentFilter <- function(deseqRes, countTable, sampleColumns, designFormula){
#	use = (deseqRes$baseMean > quantile(deseqRes$baseMean, probs = opt$threshold) & (!is.na(deseqRes$pvalue)))
	use = (deseqRes$baseMean > opt$threshold & !is.na(deseqRes$pvalue))
	ddsFilt <- DESeqDataSetFromMatrix(countData = countTable[use, ], colData = sampleTable, design = designFormula)
	computeAndWriteResults(ddsFilt, sampleColumns, opt$outfilefiltered, featurenames[use])
	
	## p-value distribution
	h1 <- hist(deseqRes$pvalue[!use], breaks=50, plot=FALSE)
	h2 <- hist(deseqRes$pvalue[use], breaks=50, plot=FALSE)
	colori <- c(`filtered out`="khaki", `complete`="powderblue")
	barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
			col = colori, space = 0, main = "Distribution of p-values", ylab="frequency")
	text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
	legend("topright", fill=rev(colori), legend=rev(names(colori)))	
}

parser <- newJSONParser()
parser$addData( opt$factors )
factorsList <- parser$getObject()

sampleColumns <- c()

for (j in 2:length(factorsList[[1]])){
    for (k in 1:length(names(factorsList[[1]][[j]]))){
        for (l in 1:length(factorsList[[1]][[j]][[k]])){
            sampleColumns[[length(sampleColumns) + 1]] <- as.integer(factorsList[[1]][[j]][[k]][l]) - 1
        }
    }
}

sampleColumns<-sort(unique(sampleColumns))
htseqSubCountTable <- htseqCountTable[,sampleColumns]
ntabcols <- length(htseqSubCountTable)
htseqSubCountTable <- as.integer(unlist(htseqSubCountTable))
htseqSubCountTable <- matrix(unlist(htseqSubCountTable), ncol = ntabcols, byrow = FALSE)
colnames(htseqSubCountTable) <- names(htseqCountTable)[sampleColumns]
sampleTable = data.frame(row.names=colnames(htseqSubCountTable))
 
for(i in 1:length(factorsList)){
    factorName<-factorsList[[i]][[1]]
    for (j in 2:length(factorsList[[i]])){
        effected_cols <- c()
        for (k in 1:length(names(factorsList[[i]][[j]]))){
            for (l in 1:length(factorsList[[i]][[j]][[k]])){
                if(colnames(htseqCountTable)[factorsList[[i]][[j]][[k]][l]-1] %in% rownames(sampleTable)){
                    sampleTable[colnames(htseqCountTable)[factorsList[[i]][[j]][[k]][l]-1],factorName] <- names(factorsList[[i]][[j]])[k]
                }
            }
        }
    }
}
sampleTable[is.na(sampleTable)] <- "undefined"
sampleTable

factorNames <- c(names(sampleTable)[-1], names(sampleTable)[1])
designFormula <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))
dds = DESeqDataSetFromMatrix(countData = htseqSubCountTable,  colData = sampleTable, design = designFormula)
deseqres <- computeAndWriteResults(dds, sampleColumns, opt$outfile) 

## independent filtering
independentFilter(deseqres, htseqSubCountTable, sampleColumns, designFormula)

dev.off()
sessionInfo()
