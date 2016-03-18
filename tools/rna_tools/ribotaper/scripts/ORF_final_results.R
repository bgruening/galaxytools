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


###script for plotting general results about the identified ORFs, takes no arguments


print(paste("--- plotting ORF finding results ---",date(),sep= " "))
ORFs<-read.table("ORFs_max_filt",stringsAsFactors=F,header=T,quote = "")
ORFs_all<-read.table("ORFs_max",stringsAsFactors=F,header=T,quote = "")

df<-(data.frame(table(ORFs_all$category,ORFs_all$annotation)))
names(df)<-c("category","annotation","n_ORFs")
df<-df[df[,"n_ORFs"]>0,]

ORFs_genes<-unique(ORFs_all[,c("category","annotation","gene_id")])
df_genes<-(data.frame(table(ORFs_genes$category,ORFs_genes$annotation)))
names(df_genes)<-c("category","annotation","n_genes")
df_genes<-df_genes[df_genes[,"n_genes"]>0,]

df_filt<-(data.frame(table(ORFs$category,ORFs$annotation)))
names(df_filt)<-c("category","annotation","n_ORFs_filtered")
df_filt<-df_filt[df_filt[,"n_ORFs_filtered"]>0,]

ORFs_genes_filt<-unique(ORFs[,c("category","annotation","gene_id")])
df_genes_filt<-(data.frame(table(ORFs_genes_filt$category,ORFs_genes_filt$annotation)))
names(df_genes_filt)<-c("category","annotation","n_genes_filtered")
df_genes_filt<-df_genes_filt[df_genes_filt[,"n_genes_filtered"]>0,]



df_new<-merge(df,df_filt,by=c("category","annotation"),all.x=T)

df_new<-merge(df_new,df_genes,by=c("category","annotation"),all.x=T)
df_new<-merge(df_new,df_genes_filt,by=c("category","annotation"),all.x=T)
df_new<-df_new[order(df_new$n_ORFs,decreasing=T),]

write.table("ORFs_genes_found",x=df_new,quote=F,row.names=F,col.names=T,sep="\t")

ORFs_coding<-ORFs[ORFs[,"annotation"]=="protein_coding",]
if(dim(ORFs_coding)[1]>0){
        tb<-as.data.frame(table(ORFs_coding$category),stringsAsFactors=F)
        names(tb)<-c("category","counts")
        tb<-tb[order(tb$counts,decreasing=T),]
        if(dim(tb)[1]>4){
                tb_ok<-tb[1:4,]
                tb_more<-tb[5:dim(tb)[1],]
                tb_other<-data.frame(counts=sum(tb_more$counts),category="other_coding",stringsAsFactors=F)
                tb<-rbind.data.frame(tb_ok,tb_other)
                ORFs_coding[ORFs_coding[,"category"]%in%tb_more$category,"category"]<-"other_coding"
                
                ORFs_coding$category<-factor(ORFs_coding$category,levels=tb$category)
                
                
        }
        if(dim(tb)[1]>=4){
                ORFs_coding$category<-factor(ORFs_coding$category,levels=tb$category)
        }
        
        
}
ncORFs<-ORFs[ORFs[,"category"]=="ncORFS",]

if(dim(ncORFs)[1]>0){
        tb<-as.data.frame(table(ncORFs$annotation),stringsAsFactors=F)
        names(tb)<-c("annotation","counts")
        tb<-tb[order(tb$counts,decreasing=T),]
        if(dim(tb)[1]>4){
                tb_ok<-tb[1:4,]
                tb_more<-tb[5:dim(tb)[1],]
                tb_other<-data.frame(counts=sum(tb_more$counts),annotation="other_ncORFs",stringsAsFactors=F)
                tb<-rbind.data.frame(tb_ok,tb_other)
                ncORFs[ncORFs[,"annotation"]%in%tb_more$annotation,"annotation"]<-"other_ncORFs"
                
                ncORFs$category<-factor(ncORFs$annotation,levels=tb$annotation)
                
                
        }
        if(dim(tb)[1]>=4){
                ncORFs$category<-factor(ncORFs$annotation,levels=tb$annotation)
        }
        
        
}


all<-rbind.data.frame(ORFs_coding[,c("category","ORF_length","ORF_P_sites")],ncORFs[,c("category","ORF_length","ORF_P_sites")])
all_filt<-all


ORFs_coding<-ORFs_all[ORFs_all[,"annotation"]=="protein_coding",]
if(dim(ORFs_coding)[1]>0){
        tb<-as.data.frame(table(ORFs_coding$category),stringsAsFactors=F)
        names(tb)<-c("category","counts")
        tb<-tb[order(tb$counts,decreasing=T),]
        if(dim(tb)[1]>4){
                tb_ok<-tb[1:4,]
                tb_more<-tb[5:dim(tb)[1],]
                tb_other<-data.frame(counts=sum(tb_more$counts),category="other_coding",stringsAsFactors=F)
                tb<-rbind.data.frame(tb_ok,tb_other)
                ORFs_coding[ORFs_coding[,"category"]%in%tb_more$category,"category"]<-"other_coding"
                
                ORFs_coding$category<-factor(ORFs_coding$category,levels=tb$category)
                
                
        }
        if(dim(tb)[1]>=4){
                ORFs_coding$category<-factor(ORFs_coding$category,levels=tb$category)
        }
        
        
}
ncORFs<-ORFs_all[ORFs_all[,"category"]=="ncORFS",]

if(dim(ncORFs)[1]>0){
        tb<-as.data.frame(table(ncORFs$annotation),stringsAsFactors=F)
        names(tb)<-c("annotation","counts")
        tb<-tb[order(tb$counts,decreasing=T),]
        if(dim(tb)[1]>4){
                tb_ok<-tb[1:4,]
                tb_more<-tb[5:dim(tb)[1],]
                tb_other<-data.frame(counts=sum(tb_more$counts),annotation="other_ncORFs",stringsAsFactors=F)
                tb<-rbind.data.frame(tb_ok,tb_other)
                ncORFs[ncORFs[,"annotation"]%in%tb_more$annotation,"annotation"]<-"other_ncORFs"
                
                ncORFs$category<-factor(ncORFs$annotation,levels=tb$annotation)
                
                
        }
        if(dim(tb)[1]>=4){
                ncORFs$category<-factor(ncORFs$annotation,levels=tb$annotation)
        }
        
        
}


all<-rbind.data.frame(ORFs_coding[,c("category","ORF_length","ORF_P_sites")],ncORFs[,c("category","ORF_length","ORF_P_sites")])



pdf(file="Final_ORF_results.pdf",width=7,height=10,onefile=T,title="ORFs_results")
par(mar=c(10,4,4,4))
par(mfrow=c(2,2))
barplot(table(all_filt$category),col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORFs_filtered",las=2)
grid(lwd=1.2,col="black")
barplot(log10(table(all_filt$category)),col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORFs_filtered(logscale)",yaxt="n",las=2)
axis(side=2,at=0:4,labels=10^(0:4))
grid(lwd=1.2,col="black")

barplot(table(all$category),col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORFs_all",las=2)
grid(lwd=1.2,col="black")
barplot(log10(table(all$category)),col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORFs_all(logscale)",yaxt="n",las=2)
axis(side=2,at=0:4,labels=10^(0:4))
grid(lwd=1.2,col="black")



par(mfrow=c(3,1))

boxplot(log10(all$ORF_P_sites)~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_P_sites",yaxt="n",main="ORFs_filtered",las=2)
axis(side=2,at=1:4,labels=10^(1:4),las=2)
grid(lwd=1.2,col="black")
boxplot(log10(all$ORF_length)~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_length",yaxt="n",las=2)
axis(side=2,at=1:4,labels=10^(1:4),las=2)
grid(lwd=1.2,col="black")
boxplot((all$ORF_P_sites/(all$ORF_length/3))~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_P_sites_per_codon",ylim=c(0,6),las=2)
grid(lwd=1.2,col="black")

par(mfrow=c(3,1))

boxplot(log10(all$ORF_P_sites)~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_P_sites",yaxt="n",main="ORFs_all",las=2)
axis(side=2,at=1:4,labels=10^(1:4),las=2)
grid(lwd=1.2,col="black")
boxplot(log10(all$ORF_length)~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_length",yaxt="n",las=2)
axis(side=2,at=1:4,labels=10^(1:4),las=2)
grid(lwd=1.2,col="black")
boxplot((all$ORF_P_sites/(all$ORF_length/3))~all$category,col=c("red","dark red","yellow","orange","grey","dark blue","blue","cornflowerblue","cyan4","grey"),ylab="ORF_P_sites_per_codon",ylim=c(0,6),las=2)
grid(lwd=1.2,col="black")


dev.off()


print(paste("--- ORF finding results Done! ---",date(),sep= " "))

