##given the list of MSIGDB gene sets containing ENTREZ identifiers and orthologous groups from EGGNOG; map human genes from MSIGDB into mouse (mm9), fly(dm3), and worm(ce10) via orthologouse relationships. 
#Supported annotations 
#human (hg19)=  #host='feb2014.archive.ensembl.org', dataset = "hsapiens_gene_ensembl") #ENSEMBL75; GRCh37 release
#worm (ce10) #host='jan2013.archive.ensembl.org', dataset = "celegans_gene_ensembl" # ENSEMBL70; WormBase v. WS220 WBcel215; CE10
#fly (dm3) =  #host='dec2014.archive.ensembl.org', dataset = "dmelanogaster_gene_ensembl") #ENSEMBL78; #BDGP-5 release dm3
#mouse (mm9) =  #host='may2012.archive.ensembl.org', dataset = "mmusculus_gene_ensembl") #ENSEMBL67 #NCBI 37 mm9

#install necessary packages 
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
suppressWarnings(suppressMessages(library(biomaRt)))
suppressWarnings(suppressMessages(library(data.table)))
#suppressWarnings(suppressMessages(library('org.Hs.eg.db'))) #mouse: org.Mm.eg.db

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

help_command = "
map_msigdb_orthologs.R: A script to map MSIGDB human gene set annotations into mouse (mm9), fly (dm3), and worm (ce10) genomes

Arguments:
--msigdb=<path to human MSIGDB gene set annotations in ENTREZ format>  - Download from http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.entrez.gmt
--help              - print this text

Example:
Rscript map_msigdb_orthologs.R --msigdb=/path/to/c2.cp.v5.0.entrez.gmt"

## Help section
if("--help" %in% args) {
  cat(help_command, "\n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## Arg1 default
if(!("msigdb" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to MSIGDB gene set annotation file")
}

msigdb = argsL$msigdb

###Define functions###
#1. parse gene sets from MSIGDB annotation file
#4. Look for enriched MSigDB gene sets for the genes found in the anot.bed compared to the background genes (found in gtf)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c(org.Hs.eg.db"))

parse_msigdb = function(msigdb){
  data = readLines(msigdb)
  
  gene_lists = list()
  mynames = c()
  for (i in 1:length(data)){
    info = unlist(strsplit(data[i], "\t"))
    gene_lists = c(gene_lists, list(paste(info[3:length(info)], sep= "\t")))
    mynames = c(mynames, info[1])
  }
  names(gene_lists) = mynames	
  
  return(gene_lists)
}

#given two biomart connections and a set of entrez gene identifiers; retrieve orthologs between mart1 and mart2 for the given list of genes 
retrieve_orthologs = function(mart1, mart2, gene_set){
  getLDS(attributes = c("entrezgene"),
         filters = "entrezgene", values = gene_set, mart = mart1,
         attributesL = c("entrezgene"), martL = mart2)
}

create_msigdb_dataset = function(hg19_genesets, orthologs){
  
  dataset = list()
  counts1 = c()
  counts2 = c()
  for (i in 1:length(hg19_genesets)){
    my_set = unique(as.numeric(paste(unlist(hg19_genesets[i]))))
    my_orth = unique(orthologs[orthologs$EntrezGene.ID %in% my_set,]$EntrezGene.ID.1)
    my_name = names(hg19_genesets)[i]
    #cat(my_name,'\n',orthologs[orthologs$EntrezGene.ID %in% my_set,]$EntrezGene.ID,'\n',my_orth,'\n')
    dataset[[my_name]] = my_orth
    counts1 = c(counts1, length(my_set))
    counts2 = c(counts2, length(my_orth))
    #if(length(my_set) > 800){cat(my_name, "\n")}
    #if((length(my_set) / length(my_orth)) < 0.5){ cat(my_name, length(my_set), length(my_orth), "\n")}
  }
  plot(counts1, counts2) #diagnostic plot to see if number of found orthologs are correlated to number of input genes
  return(dataset)
}

print_msigdb_dataset = function(dataset, output_filename){
  if (file.exists(output_filename)){
    cat('A filename with the same name exists. Delete the file first or choose a different file')
    stop()
  }
  sink(file = output_filename, append = TRUE)
  for (i in 1:length(dataset)){
    my_name = names(dataset)[i]
    my_set = paste(unique(paste(unlist(dataset[i]))), collapse = '\t')
    cat(my_name, 'no_http', my_set, sep='\t', "\n")
  }
  sink()
}

#parse lists of genes from MSigDB 
hg19_genesets = parse_msigdb(msigdb)
cat("got the gene lists from",msigdb,"\n")

#get homologs of human genes from the corresponding species using ENSEMBL biomart
#define the database connections.
hg19 = useMart(biomart='ENSEMBL_MART_ENSEMBL', host='feb2014.archive.ensembl.org', dataset = "hsapiens_gene_ensembl")
mm9 = useMart(biomart='ENSEMBL_MART_ENSEMBL', host='may2012.archive.ensembl.org', dataset = "mmusculus_gene_ensembl")
dm3 = useMart(biomart='ENSEMBL_MART_ENSEMBL', host='dec2014.archive.ensembl.org', dataset = "dmelanogaster_gene_ensembl")
ce10 = useMart(biomart='ENSEMBL_MART_ENSEMBL', host='jan2013.archive.ensembl.org', dataset = "celegans_gene_ensembl")
cat('Created the database connections to biomart\n')

#define all available human genes in all the given gene sets 
hg19_genes = unique(paste(unlist(hg19_genesets)))

#retrieve orthologs lists for all human genes found in the MSIGDB gene sets 
mm9_orth = retrieve_orthologs(mart1 = hg19, mart2 = mm9, hg19_genes)
cat('retrieved mm9 orthologs\n')
dm3_orth = retrieve_orthologs(mart1 = hg19, mart2 = dm3, hg19_genes)
cat('retrieved dm3 orthologs\n')
ce10_orth = retrieve_orthologs(mart1 = hg19, mart2 = ce10, hg19_genes)
cat('retrieved ce10 orthologs\n')


#create MSIGDB datasets for each species (worm, fly, and mouse) by using orthologous relationships to human
library('tools')
ext = file_ext(msigdb)
base = gsub(paste0('.', ext, '$'), '', msigdb)

cat('creating msigdb datasets...\n')

mm9_msigdb = create_msigdb_dataset(hg19_genesets = hg19_genesets, orthologs = mm9_orth)
print_msigdb_dataset(dataset = mm9_msigdb, output_filename = paste(base, 'mm9', ext, sep='.'))

dm3_msigdb = create_msigdb_dataset(hg19_genesets = hg19_genesets, orthologs = dm3_orth)
print_msigdb_dataset(dataset = mm9_msigdb, output_filename = paste(base, 'dm3', ext, sep='.'))

ce10_msigdb = create_msigdb_dataset(hg19_genesets = hg19_genesets, orthologs = ce10_orth)
print_msigdb_dataset(dataset = mm9_msigdb, output_filename = paste(base, 'ce10', ext, sep='.'))
