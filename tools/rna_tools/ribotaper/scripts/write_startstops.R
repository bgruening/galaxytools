#!/usr/bin/Rscript

###################################################################
#    This file is part of RiboTaper.
#    RiboTaper is a method for defining traslated ORFs using
#    Ribosome Profiling data.
#   
#    Copyright (C) 2016  Lorenzo Calviello
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

command<-paste("wc -l","start_stops_FAR.bed")
lines_in_file<-system(command,intern=T)
lines_in_file<-as.numeric(strsplit(lines_in_file,split=" ")[[1]][1])

if(lines_in_file<2){
        
        print("--------creating start_stops_FAR.bed from CDS regions")
        cdss<-read.table("all_cds.bed",stringsAsFactors=F,header=F)
        spls<-split.data.frame(cdss,f=cdss$V5)
        
        stst<-lapply(spls,FUN=function(x){
                str<-x[1,6]
                dims<-dim(x)[1]
                y<-x[1,]
                z<-x[dims,]
                if(str=="+"){
                        y[3]<-as.numeric(y[2])+1
                        y[4]<-"start_codon"
                        z[2]<-as.numeric(z[3])
                        z[3]<-as.numeric(z[2])+1
                        z[4]<-"stop_codon"
                }
                if(str=="-"){
                        y[2]<-as.numeric(y[2])-1
                        y[3]<-as.numeric(y[2])+1
                        y[4]<-"stop_codon"
                        z[2]<-as.numeric(z[3])-1
                        z[3]<-as.numeric(z[2])+1
                        z[4]<-"start_codon"
                }
                all<-rbind.data.frame(y,z)
                return(all)
        })
        stst<-do.call(rbind.data.frame,args=stst)
        
        stst<-stst[order(stst$V1,stst$V2),]
        diffs<-diff(stst$V2)
        stst<-stst[-which(abs(diffs)<100),]
        
        write.table(stst,file="start_stops_FAR.bed",quote=F,col.names=F,row.names=F,sep="\t")
}



