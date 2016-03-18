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


###script for making aggregate plots around start-stop codons, takes as argument the bed file from create_metaplots.bash

print(paste("--- plotting aggreate start-stop profiles","---",date(),sep=" "))

args <- commandArgs(trailingOnly = TRUE)

reads<-read.table(args[1],stringsAsFactors=F,header=F,sep="\t",comment.char="")
colnames(reads)<-c("chr","start","end","read_id","map_quality","strand",".1",".2",".3","spanning_exons","length_per_exon","length_introns","chr_stst","start_stst","end_stst","type_stst","gene_id_stst","strand_stst")


reads_simpl<-reads[reads[,"length_introns"]=="0",]

reads_simpl$count<-1

list_str<-split.data.frame(reads_simpl,f=reads_simpl[,"strand"])
list_str[["+"]]$distance<-list_str[["+"]][,"start"]-list_str[["+"]][,"start_stst"]
list_str[["-"]]$distance<-list_str[["-"]][,"end_stst"]-list_str[["-"]][,"end"]

reads_simpl<-do.call(rbind.data.frame,list_str)

dists_all<-with(reads_simpl,aggregate(count,by=list(type_stst,length_per_exon,distance),FUN=sum))
colnames(dists_all)<-c("type","length","distance","counts")
dists_all<-dists_all[c(which(dists_all$type=="start_codon" & dists_all$distance%in%c(-23:17,16:56)),which(dists_all$type=="stop_codon" & dists_all$distance%in%c(-68:-28,-32:10))),]

lw<-3
lengths<-as.numeric(sort(unique(dists_all$length),decreasing=F))
for(i in lengths){
        names<-paste(args,"_",as.character(i),".pdf",sep="")        
        pdf(file=names,width=31.5,height=20.5,fonts="Helvetica")
        par(mfrow=c(2,2),cex=2.2,mar=c(4.5,4,3,1))        
        starts_ok<-dists_all[dists_all[,"length"]==i & dists_all[,"type"]=="start_codon",]
        stops_ok<-dists_all[dists_all[,"length"]==i & dists_all[,"type"]=="stop_codon",]
        
        starts<-starts_ok[starts_ok[,"distance"]%in%c(-23:17),]
        if(dim(starts)[1]==0){plot(1,1,type="n")}
        if(dim(starts)[1]>0){
                plotto<-as.data.frame(t(t(-23:17)),stringsAsFactors=F)
                colnames(plotto)<-"distance"
                plotto$counts<-0
                for(g in 1:dim(plotto)[1]){
                        dis<-plotto$distance[g]
                        if(sum(starts$distance==dis)>0){
                                plotto[g,"counts"]<-starts$counts[starts$distance==dis]
                        }
                }
                plot(plotto$counts,col=c("red","blue","green"),type="h",xlab="distance (nt)",ylab="Alignments",xaxt="n",main=paste("Distance 5' - start codons\n",as.character(starts$length[1]),"nt reads",sep=" "),lwd=lw)
                axis(1, at=seq(1,length(plotto$counts),by=1), labels=as.character(seq(min(plotto$distance),max(plotto$distance),by=1)),xaxp = c(-40,40,80),las=2)
        }
        starts<-starts_ok[starts_ok[,"distance"]%in%c(16:56),]
        
        if(dim(starts)[1]==0){plot(1,1,type="n")}
        if(dim(starts)[1]>0){
                plotto<-as.data.frame(t(t(16:56)),stringsAsFactors=F)
                colnames(plotto)<-"distance"
                plotto$counts<-0
                for(g in 1:dim(plotto)[1]){
                        dis<-plotto$distance[g]
                        if(sum(starts$distance==dis)>0){
                                plotto[g,"counts"]<-starts$counts[starts$distance==dis]
                        }
                }
                plot(plotto$counts,col=c("red","blue","green"),type="h",xlab="distance (nt)",ylab="Alignments",xaxt="n",main=paste("Distance 5' - start codons\n",as.character(starts$length[1]),"nt reads",sep=" "),lwd=lw)
                axis(1, at=seq(1,length(plotto$counts),by=1), labels=as.character(seq(min(plotto$distance),max(plotto$distance),by=1)),xaxp = c(-40,40,80),las=2)
        }
        
        stops<-stops_ok[stops_ok[,"distance"]%in%c(-68:-28),]
        
        if(dim(stops)[1]==0){plot(1,1,type="n")}
        if(dim(stops)[1]>0){
                
                plotto<-as.data.frame(t(t(-68:-28)),stringsAsFactors=F)
                colnames(plotto)<-"distance"
                plotto$counts<-0
                for(g in 1:dim(plotto)[1]){
                        dis<-plotto$distance[g]
                        if(sum(stops$distance==dis)>0){
                                plotto[g,"counts"]<-stops$counts[stops$distance==dis]
                        }
                }
                plot(plotto$counts,col=c("red","blue","green"),type="h",xlab="distance (nt)",ylab="Alignments",xaxt="n",main=paste("Distance 5' - stops codons\n",as.character(stops$length[1]),"nt reads",sep=" "),lwd=lw)
                axis(1, at=seq(1,length(plotto$counts),by=1), labels=as.character(seq(min(plotto$distance),max(plotto$distance),by=1)),xaxp = c(-40,40,80),las=2)
        }
        
        
        stops<-stops_ok[stops_ok[,"distance"]%in%c(-32:10),]
        
        if(dim(stops)[1]==0){plot(1,1,type="n")}
        if(dim(stops)[1]>0){
                plotto<-as.data.frame(t(t(-32:10)),stringsAsFactors=F)
                colnames(plotto)<-"distance"
                plotto$counts<-0
                for(g in 1:dim(plotto)[1]){
                        dis<-plotto$distance[g]
                        if(sum(stops$distance==dis)>0){
                                plotto[g,"counts"]<-stops$counts[stops$distance==dis]
                        }
                }
                plot(plotto$counts,col=c("red","blue","green"),type="h",xlab="distance (nt)",ylab="Alignments",xaxt="n",main=paste("Distance 5' - stop codons\n",as.character(stops$length[1]),"nt reads",sep=" "),lwd=lw)
                axis(1, at=seq(1,length(plotto$counts),by=1), labels=as.character(seq(min(plotto$distance),max(plotto$distance),by=1)),xaxp = c(-40,40,80),las=2)
        }
        
        dev.off()
}

print(paste("--- aggregate start-stop plots, Done!","---",date(),sep=" "))

