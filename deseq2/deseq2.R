library('getopt');
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical",
  'outputfile' , 'o', 1, "character",
  'plots' , 'p', 1, "character",
  'input' , 'i', 1, "character",
  'samples', 's', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples<-trim(opt$samples)


htseqCountTable = read.table(opt$input, sep="\t", comment="", as.is=T)

l<-unique(c(htseqCountTable[1,-1]))
#print(l)

colnames(htseqCountTable) <- htseqCountTable[1,]
names(htseqCountTable)<-colnames(htseqCountTable)
tagnames <- htseqCountTable[-(1:2), 1]
htseqCountTable <- htseqCountTable[-(1:2),-1]
htseqCountTable[,1] <- gsub(",", ".", htseqCountTable[,1])


library( "DESeq2" )

opt$samples

print('----------------')

pdf(opt$plots)

if (opt$samples=="all_vs_all"){
  # all versus all
  for (i in seq(1, length(l), by=1)) {
    k=i+1
    if(k<=length(l)){
      for (j in seq(k, length(l), by=1)) {
        
        colData<-names(htseqCountTable[, which(gsub("[.][0-9]+","",names(htseqCountTable)) %in% c(l[[i]], l[[j]]))])
        currentPairTable<-htseqCountTable[ , which(gsub("[.][0-9]+","",names(htseqCountTable)) %in% c(l[[i]], l[[j]]))]
        #      rownames(currentPairTable)<-tagnames
        
        write.table(currentPairTable, file=paste(l[[i]],"_", l[[j]],".txt"))
        currentPairTable<-read.table(paste(l[[i]],"_", l[[j]],".txt"),row.names=1,header=T)
        
        currentPairTable <- as.data.frame(currentPairTable,stringsAsFactors=F)
        #print(currentPairTable)
        pdata = data.frame(condition=factor(c(rep(l[[i]], 2), rep(l[[j]], 2))),row.names=colnames(currentPairTable))
        dds = DESeqDataSetFromMatrix(countData = currentPairTable,  colData = pdata, design = ~ condition)
        
        #      colData(dds)$condition
        dds <- DESeq(dds)
        sizeFactors(dds)
        plotDispEsts(dds)
        plotMA(dds)
        
        res <- results(dds)
        rownames(res)<-tagnames
        resSorted <- res[order(res$padj),];
        #column.names=FALSE
        write.csv( as.data.frame(resSorted), file = opt$outputfile, sep="\t")
      }
    }
  }
}else{

  sampleSets<-unlist(strsplit(opt$samples, " "))
  sampleSets
  samplesControl <- {}
  samplesExperiment <- {}
  sampleNames <- {}

  samplesControl<-unlist(strsplit(sampleSets[1], ","))
  samplesExperiment<-unlist(strsplit(sampleSets[2], ","))

  print(samplesControl)
  print(samplesExperiment)

  # the minus one is needed because we get column indizes including the first column
  samplecolumns <- c()
  for (i in samplesControl) samplecolumns[[length(samplecolumns) + 1]] <- as.integer(i) - 1
  for (i in samplesExperiment) samplecolumns[[length(samplecolumns) + 1]] <- as.integer(i) - 1
  
print(samplecolumns)
  
  htseqSubCountTable <- htseqCountTable[,samplecolumns]

  #TODO: pavan hack :)
  #write.table(htseqSubCountTable, file=paste(gsub(" ","_",opt$samples),".txt"))
  #htseqSubCountTable<-read.table(paste(gsub(" ","_",opt$samples),".txt"),row.names=1,header=T)

  pdata = data.frame(condition=factor(c(names( htseqSubCountTable )[ samplecolumns ])),row.names=colnames(htseqSubCountTable))
  dds = DESeqDataSetFromMatrix(countData = bla,  colData = pdata, design = ~ condition)
  
  dds <- DESeq(dds)
  sizeFactors(dds)
  plotDispEsts(dds)
  plotMA(dds)
  
  res <- results(dds)
  rownames(res)<-tagnames
  resSorted <- res[order(res$padj),];
  #column.names=FALSE
  write.csv( as.data.frame(resSorted), file = opt$outputfile, sep="\t")
}
dev.off()
