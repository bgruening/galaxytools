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


###script for grouping and filtering identified ORFs, writing BED files and creating a protein FASTA file.


suppressMessages(library("seqinr"))

print(paste("--- create protein db and output final ORFs ---",date(),sep=" "))


ORFs_new_more<-read.table("ORFs_CCDS/more_tapers/ORFs_sign_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

ORFs_new_max<-read.table("ORFs_CCDS/max_P_sites/ORFs_sign_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

ORFs_new_best<-read.table("ORFs_CCDS/best_periodicity/ORFs_sign_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

nonccdsORFS_new_more<-read.table("ORFs_NONCCDS/more_tapers/ORFs_sign_nocds_nofilter",stringsAsFactors=F,header=T,quote = "")

nonccdsORFS_new_max<-read.table("ORFs_NONCCDS/max_P_sites/ORFs_sign_nocds_nofilter",stringsAsFactors=F,header=T,quote = "")

nonccdsORFS_new_best<-read.table("ORFs_NONCCDS/best_periodicity/ORFs_sign_nocds_nofilter",stringsAsFactors=F,header=T,quote = "")


sORFS_new_more<-read.table("ORFs_CCDS/more_tapers/sORFs_sign_filtered_cds",stringsAsFactors=F,header=T,quote = "")

sORFS_new_max<-read.table("ORFs_CCDS/max_P_sites/sORFs_sign_filtered_cds",stringsAsFactors=F,header=T,quote = "")

sORFS_new_best<-read.table("ORFs_CCDS/best_periodicity/sORFs_sign_filtered_cds",stringsAsFactors=F,header=T,quote = "")


ncORFS_new_more<-read.table("ORFs_NONCCDS/more_tapers/ORFs_sign_cds_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

ncORFS_new_max<-read.table("ORFs_NONCCDS/max_P_sites/ORFs_sign_cds_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

ncORFS_new_best<-read.table("ORFs_NONCCDS/best_periodicity/ORFs_sign_cds_notfiltered_multi",stringsAsFactors=F,header=T,quote = "")

if(dim(ORFs_new_more)[1]>0 & dim(ORFs_new_max)[1]>0 & dim(ORFs_new_best)[1]>0){
        
        ORFs_new<-rbind(ORFs_new_more,ORFs_new_max[!ORFs_new_max[,"gene_id"]%in%ORFs_new_more[,"gene_id"],])
        ORFs_new<-rbind(ORFs_new,ORFs_new_best[!ORFs_new_best[,"gene_id"]%in%ORFs_new[,"gene_id"],])
        ORFs_new$category<-"ORFs_ccds"
        ORFs_new$annotation<-"protein_coding"
        ORFs_new$header_tofasta<-paste(ORFs_new$ORF_id_tr,ORFs_new$gene_id,ORFs_new$Method,ORFs_new$annotation,ORFs_new$category,ORFs_new$ORF_P_sites,ORFs_new$ORF_spec3_spec_ribo,ORFs_new$ORF_spec_multi_ribo,sep=":")
}
if(dim(ncORFS_new_more)[1]>0 & dim(ncORFS_new_max)[1]>0 & dim(ncORFS_new_best)[1]>0){
        
        ncORFS_new<-rbind(ncORFS_new_more,ncORFS_new_max[!ncORFS_new_max[,"gene_id"]%in%ncORFS_new_more[,"gene_id"],])
        ncORFS_new<-rbind(ncORFS_new,ncORFS_new_best[!ncORFS_new_best[,"gene_id"]%in%ncORFS_new[,"gene_id"],])
        ncORFS_new$category<-"ncORFS"
        ncORFS_new$header_tofasta<-paste(ncORFS_new$ORF_id_tr,ncORFS_new$gene_id,ncORFS_new$Method,ncORFS_new$annotation,ncORFS_new$category,ncORFS_new$ORF_P_sites,ncORFS_new$ORF_spec3_spec_ribo,ncORFS_new$ORF_spec_multi_ribo,sep=":")
        ncORFS_new[,c("annotated_start","annotated_stop","ORF_id_tr_annotated")]<-NA
}
if(dim(nonccdsORFS_new_more)[1]>0 & dim(nonccdsORFS_new_max)[1]>0 & dim(nonccdsORFS_new_best)[1]>0){
        nonccdsORFS_new<-rbind(nonccdsORFS_new_more,nonccdsORFS_new_max[!nonccdsORFS_new_max[,"gene_id"]%in%nonccdsORFS_new_more[,"gene_id"],])
        nonccdsORFS_new<-rbind(nonccdsORFS_new,nonccdsORFS_new_best[!nonccdsORFS_new_best[,"gene_id"]%in%nonccdsORFS_new[,"gene_id"],])
        nonccdsORFS_new$category<-"nonccds_coding_ORFs"
        nonccdsORFS_new$header_tofasta<-paste(nonccdsORFS_new$ORF_id_tr,nonccdsORFS_new$gene_id,nonccdsORFS_new$Method,nonccdsORFS_new$annotation,nonccdsORFS_new$category,nonccdsORFS_new$ORF_P_sites,nonccdsORFS_new$ORF_spec3_spec_ribo,nonccdsORFS_new$ORF_spec_multi_ribo,sep=":")
        nonccdsORFS_new[,c("annotated_start","annotated_stop","ORF_id_tr_annotated")]<-NA
}

if(dim(sORFS_new_more)[1]>0 & dim(sORFS_new_max)[1]>0 & dim(sORFS_new_best)[1]>0){
        
        sORFS_new<-rbind(sORFS_new_more,sORFS_new_max[!sORFS_new_max[,"gene_id"]%in%sORFS_new_more[,"gene_id"],])
        sORFS_new<-rbind(sORFS_new,sORFS_new_best[!sORFS_new_best[,"gene_id"]%in%sORFS_new[,"gene_id"],])
        sORFS_new$category<-sORFS_new$type
        sORFS_new$type<-NULL
        sORFS_new$annotation<-"protein_coding"
        sORFS_new$header_tofasta<-paste(sORFS_new$ORF_id_tr,sORFS_new$gene_id,sORFS_new$Method,sORFS_new$annotation,sORFS_new$category,sORFS_new$ORF_P_sites,sORFS_new$ORF_spec3_spec_ribo,sORFS_new$ORF_spec_multi_ribo,sep=":")
}
if(is.null(sORFS_new$annotated_start)){sORFS_new$annotated_start<-NA}
if(is.null(sORFS_new$annotated_stop)){sORFS_new$annotated_stop<-NA}
if(is.null(sORFS_new$ORF_id_tr_annotated)){sORFS_new$ORF_id_tr_annotated<-NA}

if(is.null(ORFs_new$annotated_start)){ORFs_new$annotated_start<-NA}
if(is.null(ORFs_new$annotated_stop)){ORFs_new$annotated_stop<-NA}
if(is.null(ORFs_new$ORF_id_tr_annotated)){ORFs_new$ORF_id_tr_annotated<-NA}


cat_obj<-c("ORFs_new","ncORFS_new","nonccdsORFS_new","sORFS_new")
present<-c()
for(q in 1:length(cat_obj)){
        present[q]<-exists(cat_obj[q])
}

ORFs_ALL<-do.call(rbind.data.frame,mget(cat_obj[present]))

ORFs_ALL<-ORFs_ALL[(ORFs_ALL$ORF_pval_multi_ribo<0.05),]

ORFs_ALL<-ORFs_ALL[!is.na(ORFs_ALL$ORF_pept),]
names(ORFs_ALL)<-gsub(x=names(ORFs_ALL),pattern="st2vect",replacement="stop_pos")

ORFs_ALL_filt<-ORFs_ALL[which((ORFs_ALL$pct_covered_onlymulti_ribo/ORFs_ALL$pct_region_covered_ribo)<0.3),]
rem<-which(ORFs_ALL_filt[,"category"]=="nonccds_coding_ORFs" & ORFs_ALL_filt[,"annotation"]!="protein_coding")
if(length(rem)>0){
        ORFs_ALL_filt<-ORFs_ALL_filt[-rem,]
}
ORFs_ALL_filt<-ORFs_ALL_filt[!is.na(ORFs_ALL_filt$ORF_pept),]



ORFs_ALL<-ORFs_ALL[,c("gene_id","gene_symbol","transcript_id","annotation",
                      "length","strand", "n_exons", "P_sites_sum", "RNA_sites", "Ribo_cov_aver", 
                      "RNA_cov_aver","category","ORF_id_tr", "start_pos","stop_pos", "annotated_start", "annotated_stop", "ORF_id_gen", 
                      "ORF_length", "reads_ribo", "reads_rna", "ORF_P_sites","ORF_Psit_pct_in_frame", 
                      "ORF_RNA_sites", "ORF_RNAsit_pct_in_frame", "ORF_pval_multi_ribo", 
                      "ORF_pval_multi_rna","ORF_spec_multi_ribo","ORF_spec_multi_rna", "ORF_id_tr_annotated", "n_exons_ORF","pct_region_covered_ribo", "pct_covered_onlymulti_ribo", "pct_region_covered_rna",
                      "pct_covered_onlymulti_rna", "Method", "header_tofasta", "ORF_pept")
                   ]

ORFs_ALL_filt<-ORFs_ALL_filt[,c("gene_id","gene_symbol","transcript_id","annotation",
                                "length","strand", "n_exons", "P_sites_sum", "RNA_sites", "Ribo_cov_aver", 
                                "RNA_cov_aver","category","ORF_id_tr", "start_pos","stop_pos", "annotated_start", "annotated_stop", "ORF_id_gen", 
                                "ORF_length","reads_ribo","reads_rna", "ORF_P_sites", "ORF_Psit_pct_in_frame", 
                                "ORF_RNA_sites", "ORF_RNAsit_pct_in_frame", "ORF_pval_multi_ribo", 
                                "ORF_pval_multi_rna","ORF_spec_multi_ribo","ORF_spec_multi_rna", "ORF_id_tr_annotated", "n_exons_ORF","pct_region_covered_ribo", "pct_covered_onlymulti_ribo", "pct_region_covered_rna",
                                "pct_covered_onlymulti_rna", "Method", "header_tofasta", "ORF_pept")
                             ]

names(ORFs_ALL)[which(names(ORFs_ALL)=="reads_ribo")]<-"ORF_reads_ribo"
names(ORFs_ALL)[which(names(ORFs_ALL)=="reads_rna")]<-"ORF_reads_rna"
names(ORFs_ALL_filt)[which(names(ORFs_ALL_filt)=="reads_ribo")]<-"ORF_reads_ribo"
names(ORFs_ALL_filt)[which(names(ORFs_ALL_filt)=="reads_rna")]<-"ORF_reads_rna"

#write.table(ORFs_ALL_filt,file="ORFs_more_filt",quote=F,col.names=T,row.names=F,sep="\t")
#write.table(ORFs_ALL,file="ORFs_more",quote=F,col.names=T,row.names=F,sep="\t")
#write.fasta(sequences=as.list(ORFs_ALL$ORF_pept),names=ORFs_ALL$header_tofasta,file.out="protein_db_more.fasta")

if(dim(ORFs_new_more)[1]>0 & dim(ORFs_new_max)[1]>0 & dim(ORFs_new_best)[1]>0){
        ORFs_new<-rbind(ORFs_new_max,ORFs_new_more[!ORFs_new_more[,"gene_id"]%in%ORFs_new_max[,"gene_id"],])
        ORFs_new<-rbind(ORFs_new,ORFs_new_best[!ORFs_new_best[,"gene_id"]%in%ORFs_new[,"gene_id"],])
        ORFs_new$category<-"ORFs_ccds"
        ORFs_new$annotation<-"protein_coding"
        ORFs_new$header_tofasta<-paste(ORFs_new$ORF_id_tr,ORFs_new$gene_id,ORFs_new$Method,ORFs_new$annotation,ORFs_new$category,ORFs_new$ORF_P_sites,ORFs_new$ORF_spec3_spec_ribo,ORFs_new$ORF_spec_multi_ribo,sep=":")
}
if(dim(ncORFS_new_more)[1]>0 & dim(ncORFS_new_max)[1]>0 & dim(ncORFS_new_best)[1]>0){
        
        ncORFS_new<-rbind(ncORFS_new_max,ncORFS_new_more[!ncORFS_new_max[,"gene_id"]%in%ncORFS_new_max[,"gene_id"],])
        ncORFS_new<-rbind(ncORFS_new,ncORFS_new_best[!ncORFS_new_best[,"gene_id"]%in%ncORFS_new[,"gene_id"],])
        ncORFS_new$category<-"ncORFS"
        ncORFS_new$header_tofasta<-paste(ncORFS_new$ORF_id_tr,ncORFS_new$gene_id,ncORFS_new$Method,ncORFS_new$annotation,ncORFS_new$category,ncORFS_new$ORF_P_sites,ncORFS_new$ORF_spec3_spec_ribo,ncORFS_new$ORF_spec_multi_ribo,sep=":")
        ncORFS_new[,c("annotated_start","annotated_stop","ORF_id_tr_annotated")]<-NA
}
if(dim(nonccdsORFS_new_more)[1]>0 & dim(nonccdsORFS_new_max)[1]>0 & dim(nonccdsORFS_new_best)[1]>0){
        
        nonccdsORFS_new<-rbind(nonccdsORFS_new_max,nonccdsORFS_new_more[!nonccdsORFS_new_max[,"gene_id"]%in%nonccdsORFS_new_max[,"gene_id"],])
        nonccdsORFS_new<-rbind(nonccdsORFS_new,nonccdsORFS_new_best[!nonccdsORFS_new_best[,"gene_id"]%in%nonccdsORFS_new[,"gene_id"],])
        nonccdsORFS_new$category<-"nonccds_coding_ORFs"
        nonccdsORFS_new$header_tofasta<-paste(nonccdsORFS_new$ORF_id_tr,nonccdsORFS_new$gene_id,nonccdsORFS_new$Method,nonccdsORFS_new$annotation,nonccdsORFS_new$category,nonccdsORFS_new$ORF_P_sites,nonccdsORFS_new$ORF_spec3_spec_ribo,nonccdsORFS_new$ORF_spec_multi_ribo,sep=":")
        nonccdsORFS_new[,c("annotated_start","annotated_stop","ORF_id_tr_annotated")]<-NA
}
if(dim(sORFS_new_more)[1]>0 & dim(sORFS_new_max)[1]>0 & dim(sORFS_new_best)[1]>0){
        
        sORFS_new<-rbind(sORFS_new_max,sORFS_new_more[!sORFS_new_max[,"gene_id"]%in%sORFS_new_max[,"gene_id"],])
        sORFS_new<-rbind(sORFS_new,sORFS_new_best[!sORFS_new_best[,"gene_id"]%in%sORFS_new[,"gene_id"],])
        sORFS_new$category<-sORFS_new$type
        sORFS_new$type<-NULL
        sORFS_new$annotation<-"protein_coding"
        sORFS_new$header_tofasta<-paste(sORFS_new$ORF_id_tr,sORFS_new$gene_id,sORFS_new$Method,sORFS_new$annotation,sORFS_new$category,sORFS_new$ORF_P_sites,sORFS_new$ORF_spec3_spec_ribo,sORFS_new$ORF_spec_multi_ribo,sep=":")
}
if(is.null(sORFS_new$annotated_start)){sORFS_new$annotated_start<-NA}
if(is.null(sORFS_new$annotated_stop)){sORFS_new$annotated_stop<-NA}
if(is.null(sORFS_new$ORF_id_tr_annotated)){sORFS_new$ORF_id_tr_annotated<-NA}

if(is.null(ORFs_new$annotated_start)){ORFs_new$annotated_start<-NA}
if(is.null(ORFs_new$annotated_stop)){ORFs_new$annotated_stop<-NA}
if(is.null(ORFs_new$ORF_id_tr_annotated)){ORFs_new$ORF_id_tr_annotated<-NA}


cat_obj<-c("ORFs_new","ncORFS_new","nonccdsORFS_new","sORFS_new")
present<-c()
for(q in 1:length(cat_obj)){
        present[q]<-exists(cat_obj[q])
}

ORFs_ALL<-do.call(rbind.data.frame,mget(cat_obj[present]))
ORFs_ALL<-ORFs_ALL[(ORFs_ALL$ORF_pval_multi_ribo<0.05),]

ORFs_ALL<-ORFs_ALL[!is.na(ORFs_ALL$ORF_pept),]
names(ORFs_ALL)<-gsub(x=names(ORFs_ALL),pattern="st2vect",replacement="stop_pos")

ORFs_ALL_filt<-ORFs_ALL[which((ORFs_ALL$pct_covered_onlymulti_ribo/ORFs_ALL$pct_region_covered_ribo)<0.3),]
rem<-which(ORFs_ALL_filt[,"category"]=="nonccds_coding_ORFs" & ORFs_ALL_filt[,"annotation"]!="protein_coding")
if(length(rem)>0){
        ORFs_ALL_filt<-ORFs_ALL_filt[-rem,]
}

ORFs_ALL_filt<-ORFs_ALL_filt[!is.na(ORFs_ALL_filt$ORF_pept),]


list_coords_bed<-list()

for(i in 1:dim(ORFs_ALL)[1]){
        orf<-ORFs_ALL[i,]
        strand<-orf$strand
        P_sites_sum<-orf$ORF_P_sites
        orf_id<-orf$ORF_id_tr
        orf_category<-orf$category
        orf_annotation<-orf$annotation
        
        all_ex<-strsplit(orf$to_check_ALL,split=";")[[1]]
        all_ex<-all_ex[all_ex!="NA"]
        list_exs<-list()
        for(j in 1:length(all_ex)){
                ex<-strsplit(all_ex[j],split="_")[[1]]
                bed<-data.frame(chr=ex[1],start=ex[2],end=ex[3],orf_name=paste(orf_id,orf_category,orf_annotation,sep=";"),P_sites=P_sites_sum,strand_bed=strand,stringsAsFactors=F)
                list_exs[[j]]<-bed
        }
        exs<-do.call(args=list_exs,what=rbind.data.frame)
        list_coords_bed[[i]]<-exs
        
}

coords_bed<-do.call(args=list_coords_bed,what=rbind.data.frame)


write.table(file="translated_ORFs.bed",x=coords_bed,col.names=F,row.names=F,quote=F,sep="\t")

system("sort -k1,1 -k2,2n translated_ORFs.bed > translated_ORFs_sorted.bed")
system("rm translated_ORFs.bed")


list_coords_bed<-list()

for(i in 1:dim(ORFs_ALL_filt)[1]){
        orf<-ORFs_ALL_filt[i,]
        strand<-orf$strand
        P_sites_sum<-orf$ORF_P_sites
        orf_id<-orf$ORF_id_tr
        orf_category<-orf$category
        orf_annotation<-orf$annotation
        
        all_ex<-strsplit(orf$to_check_ALL,split=";")[[1]]
        all_ex<-all_ex[all_ex!="NA"]
        list_exs<-list()
        for(j in 1:length(all_ex)){
                ex<-strsplit(all_ex[j],split="_")[[1]]
                bed<-data.frame(chr=ex[1],start=ex[2],end=ex[3],orf_name=paste(orf_id,orf_category,orf_annotation,sep=";"),P_sites=P_sites_sum,strand_bed=strand,stringsAsFactors=F)
                list_exs[[j]]<-bed
        }
        exs<-do.call(args=list_exs,what=rbind.data.frame)
        list_coords_bed[[i]]<-exs
        
}

coords_bed<-do.call(args=list_coords_bed,what=rbind.data.frame)


write.table(file="translated_ORFs_filtered.bed",x=coords_bed,col.names=F,row.names=F,quote=F,sep="\t")

system("sort -k1,1 -k2,2n translated_ORFs_filtered.bed > translated_ORFs_filtered_sorted.bed")
system("rm translated_ORFs_filtered.bed")




ORFs_ALL<-ORFs_ALL[,c("gene_id","gene_symbol","transcript_id","annotation",
                      "length","strand", "n_exons", "P_sites_sum", "RNA_sites", "Ribo_cov_aver", 
                      "RNA_cov_aver","category","ORF_id_tr", "start_pos","stop_pos", "annotated_start", "annotated_stop", "ORF_id_gen", 
                      "ORF_length", "reads_ribo", "reads_rna", "ORF_P_sites","ORF_Psit_pct_in_frame", 
                      "ORF_RNA_sites", "ORF_RNAsit_pct_in_frame", "ORF_pval_multi_ribo", 
                      "ORF_pval_multi_rna","ORF_spec_multi_ribo","ORF_spec_multi_rna", "ORF_id_tr_annotated", "n_exons_ORF","pct_region_covered_ribo", "pct_covered_onlymulti_ribo", "pct_region_covered_rna",
                      "pct_covered_onlymulti_rna", "Method", "header_tofasta", "ORF_pept")
                   ]

ORFs_ALL_filt<-ORFs_ALL_filt[,c("gene_id","gene_symbol","transcript_id","annotation",
                                "length","strand", "n_exons", "P_sites_sum", "RNA_sites", "Ribo_cov_aver", 
                                "RNA_cov_aver","category","ORF_id_tr", "start_pos","stop_pos", "annotated_start", "annotated_stop", "ORF_id_gen", 
                                "ORF_length","reads_ribo", "reads_rna", "ORF_P_sites", "ORF_Psit_pct_in_frame", 
                                "ORF_RNA_sites", "ORF_RNAsit_pct_in_frame", "ORF_pval_multi_ribo", 
                                "ORF_pval_multi_rna","ORF_spec_multi_ribo","ORF_spec_multi_rna", "ORF_id_tr_annotated", "n_exons_ORF","pct_region_covered_ribo", "pct_covered_onlymulti_ribo", "pct_region_covered_rna",
                                "pct_covered_onlymulti_rna", "Method", "header_tofasta", "ORF_pept")
                             ]


names(ORFs_ALL)[which(names(ORFs_ALL)=="reads_ribo")]<-"ORF_reads_ribo"
names(ORFs_ALL)[which(names(ORFs_ALL)=="reads_rna")]<-"ORF_reads_rna"
names(ORFs_ALL_filt)[which(names(ORFs_ALL_filt)=="reads_ribo")]<-"ORF_reads_ribo"
names(ORFs_ALL_filt)[which(names(ORFs_ALL_filt)=="reads_rna")]<-"ORF_reads_rna"


write.table(ORFs_ALL_filt,file="ORFs_max_filt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(ORFs_ALL,file="ORFs_max",quote=F,col.names=T,row.names=F,sep="\t")
write.fasta(sequences=as.list(ORFs_ALL$ORF_pept),names=ORFs_ALL$header_tofasta,file.out="protein_db_max.fasta")

print(paste("--- protein db and output final ORFs, Done! ---",date(),sep=" "))

