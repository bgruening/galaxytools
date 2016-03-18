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


###script for analyzing data tracks, it takes as arguments the name tag of the track (i.e. "ccds"), the RiboTaper scripts directory, the n of cores

args <- commandArgs(trailingOnly = TRUE)

print(paste("--- analyzing",args[1],"exonic tracks","---",date(),sep=" "))


suppressMessages(source(paste(args[2],"functions.R",sep = "/")))


registerDoMC(args[3])

sink(file=NULL,type="message")


all_tracks<-readBigText(paste("data_tracks/Psit_Ribo_Rna_Cent_tracks",as.character(args[1]),sep="_"))

all_index<-read.table(paste("data_tracks/index_tracks",as.character(args[1]),sep="_"),stringsAsFactors=F,header=F)

index_chrs<-sapply(strsplit(all_index$V1,split="_"),"[[",1)


regions<-unique(all_index)

results_regions<-list()



results_regions<-foreach(s=1:dim(regions)[1],.combine=rbind,.multicombine=T) %dopar%{
        tryCatch({
                chr_reg<-sapply(strsplit(regions[s,1],split="_"),"[[",1)
                all_tr<-all_tracks[which(index_chrs==chr_reg)]
                ind_tr<-subset(all_index,index_chrs==chr_reg)
                x<-all_tr[ind_tr==regions[s,1]]
                x<-t(data.frame(strsplit(x,split=" ")))
                return(make_analysis_exons(x))
                
        }, error=function(x){
                
                return("error")
        }
        )
}




names_covbeds<-c("exon_id","strand","reads","bases_covered","total_bases","pct_region_covered")


RIBO_best<-read.table(paste("RIBO_best_counts",as.character(args[1]),sep="_"),stringsAsFactors=F,header=F)

colnames(RIBO_best)<-names_covbeds


RIBO_unique<-read.table(paste("RIBO_unique_counts",as.character(args[1]),sep="_"),stringsAsFactors=F,header=F)

colnames(RIBO_unique)<-names_covbeds


RNA_best<-read.table(paste("RNA_best_counts",as.character(args[1]),sep="_"),stringsAsFactors=F,header=F)

colnames(RNA_best)<-names_covbeds


RNA_unique<-read.table(paste("RNA_unique_counts",as.character(args[1]),sep="_"),stringsAsFactors=F,header=F)

colnames(RNA_unique)<-names_covbeds




multi_table_RIBO<-merge(RIBO_best,RIBO_unique,by="exon_id")
multi_table_RNA<-merge(RNA_best,RNA_unique,by="exon_id")
multi_table_RIBO$pct_covered_onlymulti<-multi_table_RIBO$pct_region_covered.x-multi_table_RIBO$pct_region_covered.y
multi_table_RIBO$reads_multi<-multi_table_RIBO$reads.x-multi_table_RIBO$reads.y

multi_table_RNA$pct_covered_onlymulti<-multi_table_RNA$pct_region_covered.x-multi_table_RNA$pct_region_covered.y
multi_table_RNA$reads_multi<-multi_table_RNA$reads.x-multi_table_RNA$reads.y


multi_table_RIBO<-multi_table_RIBO[,c(1,2,5,3,13,6,12)]
names(multi_table_RIBO)<-c("exon_id", "strand", "length", "reads_ribo", "reads_multi_ribo","pct_region_covered_ribo", 
                           "pct_covered_onlymulti_ribo")



multi_table_RNA<-multi_table_RNA[,c(1,3,13,6,12)]
names(multi_table_RNA)<-c("exon_id", "reads_rna", "reads_multi_rna", "pct_region_covered_rna", 
                          "pct_covered_onlymulti_rna")

multi_table<-merge(multi_table_RIBO,multi_table_RNA,by="exon_id")

RESULTS<-merge(results_regions,multi_table,by=1)

write.table(RESULTS,file=paste("results",args[1],sep="_"),quote=F,sep="\t",row.names=F)


print(paste("--- track_analysis Done!","---",date(),sep=" "))

