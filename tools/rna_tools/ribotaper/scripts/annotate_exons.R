#!/usr/bin/Rscript

###################################################################
#    This file is part of RiboTaper.
#    RiboTaper is a method for defining traslated ORFs using
#    Ribosome Profiling data.
#   
#    Copyright (C) 2015  Lorenzo Calviello
#
#    RiboTaper is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    RiboTaper is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RiboTaper.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Lorenzo.Calviello@mdc-berlin.de
#######################################################################


###script for annotating exons, it takes as arguments the annotation directory, the RiboTaper scripts directory, the n of cores

args <- commandArgs(trailingOnly = TRUE)


print(paste("--- annotating exons","---",date(),sep=" "))


suppressMessages(source(paste(args[2],"functions.R",sep = "/")))


genes<-paste(args[1],"gene_annot_names",sep = "/")

genes_annot<-read.table(genes,stringsAsFactors=F,header=F)

colnames(genes_annot)<-c("gene_id","annotation","gene_symbol")

nonccds_res<-read.table("results_nonccds",header=T,stringsAsFactors=F)

ids_nonccds<-sapply(nonccds_res$exon_id,FUN=function(x){strsplit(x,split="_")})
nonccds_res$gene_id<-as.character(lapply(ids_nonccds,"[[",5))

nonccds_res<-merge(nonccds_res,genes_annot,by="gene_id")
nonccds_res$type<-"non_ccds_exon"



ccds_res<-read.table("results_ccds",header=T,stringsAsFactors=F)
ids_ccds<-sapply(ccds_res$exon_id,FUN=function(x){strsplit(x,split="_")})
ccds_res$gene_id<-as.character(lapply(ids_ccds,"[[",5))

ccds_res<-merge(ccds_res,genes_annot,by="gene_id")
ccds_res$type<-"ccds"


exons_ccds_res<-read.table("results_exonsccds",header=T,stringsAsFactors=F)
ids_exons_ccds<-sapply(exons_ccds_res$exon_id,FUN=function(x){strsplit(x,split="_")})
exons_ccds_res$gene_id<-as.character(lapply(ids_exons_ccds,"[[",5))

exons_ccds_res<-merge(exons_ccds_res,genes_annot,by="gene_id")
exons_ccds_res$type<-"exon"


all<-rbind(ccds_res,exons_ccds_res)

coords<-matrix(nrow=dim(all)[1],ncol=1)
for(i in seq(1,dim(all)[1])){
        coords[i,1]<-paste(strsplit(all$exon_id[i],split="_")[[1]][1:3],collapse="_")
}
all$coords<-coords

all <-all[order(all$coords,all$type,decreasing=F),]


coords2<-matrix(nrow=dim(all)[1],ncol=3)
for(i in seq(1,dim(all)[1])){
        coords2[i,1]<-strsplit(all$exon_id[i],split="_")[[1]][1]
        coords2[i,2]<-strsplit(all$exon_id[i],split="_")[[1]][2]
        coords2[i,3]<-strsplit(all$exon_id[i],split="_")[[1]][3]
}


all$chr<-coords2[,1]
all$start<-as.integer(coords2[,2])
all$end<-as.integer(coords2[,3])

all$nt_more<-NA
all$nt_more_ribocovered<-NA
all$nt_more_P_sites<-NA
all$nt_more_rnacovered<-NA
all$nt_more_cent_sites<-NA
all$overlapping_ccds_start<-NA
all$overlapping_ccds_end<-NA


list_genes_exon_ccds<-split.data.frame(all,f=all$gene_id,drop=T)


list_genes_exon_ccds_annot<-list()

list_genes_exon_ccds_annot<-mclapply(X=list_genes_exon_ccds,FUN=annotate_exons,mc.cores=args[3],mc.preschedule = TRUE)


all_annot<-do.call(rbind.data.frame,list_genes_exon_ccds_annot)



write.table(file="results_nonccds_annot",sep="\t",nonccds_res,quote=F,row.names=F)


write.table(file="all_calculations_ccdsgenes_annot_new",sep="\t",all_annot,quote=F,row.names=F)
print(paste("--- annotating exons, Done!","---",date(),sep=" "))
