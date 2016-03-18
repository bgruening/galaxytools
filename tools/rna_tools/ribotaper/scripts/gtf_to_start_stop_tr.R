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


###script for creating transcript-level coordinates of CDS positions (start and stop) from a .gtf file

print(paste("--- extracting transcript-level CDS cordinates","---",date(),sep=" "))

exons_cds_all<-read.table("exons_cds_all",stringsAsFactors=F,header=F)
colnames(exons_cds_all)<-c("chr","type","start","end","strand","transcript_id")
tr_cds<-unique(exons_cds_all[exons_cds_all[,"type"]=="CDS","transcript_id"])

exons_cds_all2<-exons_cds_all[exons_cds_all[,"transcript_id"]%in%tr_cds,]
exons_cds_all2$length<-1+(exons_cds_all2$end-exons_cds_all2$start)
list_exons_cds_tr<-split.data.frame(x=exons_cds_all2,f=exons_cds_all2$transcript_id,drop=T)

list_coords<-list()
for(i in 1:length(list_exons_cds_tr)){
        transcr<-tr_cds[i]
        trascr_data<-list_exons_cds_tr[[transcr]]
        
        strand<-trascr_data$strand[1]
        
        exons_in_transcr<-trascr_data[trascr_data[,"type"]=="exon",]
        if(strand=="-"){exons_in_transcr<-exons_in_transcr[dim(exons_in_transcr)[1]:1,]}
        
        
        
        cds_in_transcr<-trascr_data[trascr_data[,"type"]=="CDS",]
        if(strand=="-"){cds_in_transcr<-cds_in_transcr[dim(cds_in_transcr)[1]:1,]}
        
        
        cumsumexons<-cumsum(exons_in_transcr$length)
        revcumsumexons<-cumsum(rev(exons_in_transcr$length))
        
        cumsumcds<-cumsum(cds_in_transcr$length)
        
        st_cod<-cds_in_transcr[1,"start"]
        if(strand=="-"){st_cod<-cds_in_transcr[1,"end"]}
        
        end_cod<-(cds_in_transcr[dim(cds_in_transcr)[1],"end"])
        if(strand=="-"){end_cod<-(cds_in_transcr[dim(cds_in_transcr)[1],"start"])}
        
        st_ex<-which((st_cod>=exons_in_transcr$start & st_cod<=exons_in_transcr$end))
        
        end_ex<-which((end_cod>=exons_in_transcr$start & end_cod<=exons_in_transcr$end))
        
        nt_dist_start<-st_cod-exons_in_transcr[st_ex,"start"]
        if(strand=="-"){nt_dist_start<-exons_in_transcr[st_ex,"end"]-st_cod}
        
        if(st_ex>1){nt_dist_start<-nt_dist_start+cumsumexons[st_ex-1]}
        
        nt_dist_stop<-exons_in_transcr[end_ex,"end"]-end_cod
        if(strand=="-"){nt_dist_stop<-end_cod-exons_in_transcr[end_ex,"start"]}
        
        if(end_ex<dim(exons_in_transcr)[1]){nt_dist_stop<-nt_dist_stop+revcumsumexons[dim(exons_in_transcr)[1]-end_ex]}
        
        tr_len<-sum(exons_in_transcr$length)
        start_coord<-nt_dist_start+1
        stop_coord<-(tr_len-nt_dist_stop)+3
        if(nt_dist_stop==0){stop_coord<-tr_len}
        x<-data.frame(transcript_id<-transcr,start_tr<-start_coord,stop_tr<-stop_coord)
        list_coords[[i]]<-x
}

coords<-do.call(what=rbind.data.frame,args=list_coords)
colnames(coords)<-c("transcript_id","start_tx","stop_tx")
write.table(coords,file="cds_coords_transcripts",row.names=F,col.names=F,quote=F,sep="\t")

print(paste("--- extracting transcript-level CDS cordinates, Done!","---",date(),sep=" "))

