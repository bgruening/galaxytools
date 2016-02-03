##Download curated pathways from MSigDB: 
##http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/c2.cp.v5.0.symbols.gmt

suppressWarnings(suppressMessages(library('org.Hs.eg.db'))) #mouse: org.Mm.eg.db
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library('org.Ce.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Mm.eg.db'))) 
suppressWarnings(suppressMessages(library('org.Dm.eg.db'))) 
suppressWarnings(suppressMessages(library('rtracklayer')))
  
#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}

help_command = "
rcas.msigdb.R: Gene Set Enrichment Analysis module of RCAS based on Gene Sets from MSigDB

Arguments:
--gmt=<path to MSigDB gmt file>  - e.g /home/buyar/projecst/RCAS/c2.cp.v5.0.symbols.gmt
--gff3=<path to GENCODE Annotation File>  - e.g /data/akalin/Base/Annotation/GenomeAnnotation/hg19/gencode/gencode.v19.annotation.gff3
--anot=<path to parse_anot.py output>  - e.g /home/buyar/projects/RCAS/test/PARCLIP_AGO1234_Hafner2010a_hg19_xaa.anot.tsv
--out=<output prefix>   -e.g Hafner2010.hg19
--species=<species name>    -Choose between (human, fly, worm, mouse)
--help              - print this text

Example:
Rscript rcas.msigdb.R --gmt=c2.cp.v5.0.entrez.gmt --gff3=gencode.v19.annotation.gff3 --anot=PARCLIP_AGO1234_Hafner2010a_hg19_xaa.anot.tsv --out=myproject.hg19 --species=human"

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

if(!("gmt" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to gmt file")
}

if(!("gff3" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to gff3 file")
}

if(!("anot" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the path to anot.bed file")
}

if(!("out" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("provide the output prefix to anot.bed file")
}

## Arg4
if(!("species" %in% argsDF$V1)) {
  cat(help_command, "\n")
  stop("Choose which species' annotation to retrieve: human, fly, worm, or mouse")
}

msigdb = argsL$gmt
gff3_file = argsL$gff3
anot_file = argsL$anot
out_file = argsL$out
species = argsL$species

if (!(species %in% c('human', 'fly', 'worm', 'mouse'))){
  cat(help_command, "\n")
  cat('selected species \"',species, '\" is not supported','\n')
  stop("Choose which species' annotation to retrieve. Supported species are: human, fly, worm, or mouse")
}

########################################################################################################################################
#2. Get the list of all protein-coding genes from the GTF file 

if(file.exists(paste0(gff3_file, ".rds")))
{
  gff = readRDS(paste0(gff3_file, ".rds"))
}else
{
  gff = import.gff3(gff3_file)
  saveRDS(gff, file=paste0(gff3_file, ".rds"))
}
idx <- mcols(gff)$gene_type == "protein_coding" & mcols(gff)$type=="gene"
all_gene_ids = unique(gff[idx]$gene_id)

cat("Read ",length(all_gene_ids),"genes from ",gff3_file,"\n")

#3. Get the list of protein-coding genes from the anot.bed file => those that are found to be targeted in the PAR-CLIP experiments
ann = fread(anot_file)
ann.gene_ids = unique(ann$gene_id)
cat("Read ",length(ann.gene_ids),"genes from ",anot_file, "\n")

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

count_associations = function(treatment, background, gene_lists){
  i = 0
  mynames = names(gene_lists)
  t_size = length(treatment)
  b_size = length(background)
  t_counts = c()
  b_counts = c()
  exp_vals = c()
  enrichment_status = c()
  pval_calc = c()
  for (l in gene_lists){
    i = i + 1
    name = mynames[i]
    t = sum(treatment %in% l)
    b = sum(background %in% l)
    exp = t_size * (b/b_size)
    comparison =  matrix(c(t, b, t_size - t, b_size - b), nrow = 2, dimnames = list(c("treatment", "background"), c("found", "not_found")))
    pval  = fisher.test(comparison, alternative = "two.sided")$p.value
    t_counts = c(t_counts, t)
    b_counts = c(b_counts, b)
    exp_vals = c(exp_vals, exp)
    if (t >= exp){enrichment_status = c(enrichment_status, "enriched")}
    if (t < exp){enrichment_status = c(enrichment_status, "depleted")}
    pval_calc = c(pval_calc, pval)
  }
  results = cbind.data.frame(mynames, t_counts, rep(t_size, length(mynames)), b_counts, rep(b_size, length(mynames)), exp_vals, enrichment_status, pval_calc)
  colnames(results) = c("list_name", "treatment_count", "treatment_size", "background_count", "background_size", "expected_in_treatment", "enrichment_status", "pval")
  results$bonferroni = p.adjust(results$pval, method = "bonferroni")
  results$BH = p.adjust(results$pval, method = "BH")
  return(data.table(results))
}

get_unique_items_from_list = function(mylist){
  
  items = c()
  for (l in mylist)
  {
    items = c(items, unique(l))
  }
  return(items)
}

#parse lists of genes from MSigDB 
gene_lists = parse_msigdb(msigdb)
cat("got the gene lists from",msigdb,"\n")

#map ENSEMBL gene ids to ENtrez Gene Ids
if (species == 'human'){
  mapping_dict = org.Hs.egENSEMBL2EG
}else if (species == 'mouse'){
  mapping_dict = org.Mm.egENSEMBL2EG
}else if (species == 'fly'){
  mapping_dict = org.Dm.egENSEMBL2EG
}else if (species == 'worm'){
  mapping_dict = org.Ce.egENSEMBL2EG
}
ens2eg <- as.list(mapping_dict) #mapping dictionary from ENSEMBL to ENTREZ 
treatment = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", ann.gene_ids)])
background = get_unique_items_from_list(ens2eg[gsub("\\..*$", "", all_gene_ids)])

cat("mapped ",length(ann.gene_ids),"from ensembl to",length(treatment),"gene ids from entrez for treatment", "\n")
cat("mapped ",length(all_gene_ids)," from ensembl to",length(background),"gene ids from entrez for background", "\n")

#count the number of times each pathway is associated to the gene list of interest
#compare the frequency of associations of different pathways to the treatment set versus the background set. 
#find pathways that are enriched in the treatment set. 
results = count_associations(treatment = treatment, background = background, gene_lists)
results = results[order(pval)]

#Write pathways with FDR < 0.001 to file
#write.table(results[bonferroni < 0.001], file = out_file, quote = FALSE, row.names = FALSE, sep="\t")
write.table(results, file = out_file, quote = FALSE, row.names = FALSE, sep="\t")




