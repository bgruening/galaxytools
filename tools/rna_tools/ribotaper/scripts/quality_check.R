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


###script for plotting general results about the exon-level analysis as a QC step, takes as arguments the annotation directory


print(paste("--- create QC plots ---",date(),sep=" "))




args <- commandArgs(trailingOnly = TRUE)

ribo_best<-system("samtools view -c RIBO_best.bam",intern=T)
ribo_unique<-system("samtools view -c RIBO_unique.bam",intern=T)
ribo_psit<-system("wc -l P_sites_all ",intern=T)
ribo_psit<-strsplit(ribo_psit,split=" ")[[1]][1]

rna_best<-system("samtools view -c RNA_best.bam",intern=T)
rna_unique<-system("samtools view -c RNA_unique.bam",intern=T)
rna_psit<-system("wc -l Centered_RNA ",intern=T)
rna_psit<-strsplit(rna_psit,split=" ")[[1]][1]

all_annot<-read.table("all_calculations_ccdsgenes_annot_new",header=T,stringsAsFactors=F,quote = "")

ccds<-all_annot[all_annot[,"type"]=="ccds" & all_annot[,"P_sites_sum"]>5 & !is.na(all_annot$pval_multit_3nt_ribo),]

nonccds<-read.table("results_nonccds_annot",header=T,stringsAsFactors=F,quote = "")
noncoding<-nonccds[nonccds[,"annotation"]!="protein_coding",]
noncoding<-noncoding[noncoding[,"P_sites_sum"]>5,]
noncoding<-noncoding[noncoding[,"length.x"]>5,]
utrs<-all_annot[all_annot[,"type"]%in%c("3_utrs_ex","5_utrs_ex"),]
utrs<-utrs[utrs[,"P_sites_sum"]>5,]

fra<-paste(args[1],"frames_ccds",sep = "/")

frames<-read.table(fra,stringsAsFactors=F,header=F)

colnames(frames)<-c("exon_id","frame_start_annot","strand","length")
ccds_frames<-merge(ccds,frames,by="exon_id")
ccds_frames$ok_annot<-FALSE

ccds_frames[ccds_frames[,"strand"]=="+" & ccds_frames[,"frame_start_annot"]==ccds_frames[,"frame_start_pred"],"ok_annot"]<-TRUE
ccds_frames[ccds_frames[,"strand"]=="-" & ccds_frames[,"frame_start_annot"]==ccds_frames[,"frame_end_pred"],"ok_annot"]<-TRUE


lib_size<-as.numeric(ribo_psit)

all_ccds<-ccds

all_ccds$RPKM_ribo<-apply(X=all_ccds,MARGIN=1,FUN=function(x){(10^9 * as.numeric(x["P_sites_sum"]))/(lib_size * as.numeric(x["length.x"]))})

quantiles_RPKM_ribo<-quantile(all_ccds$RPKM_ribo,probs=seq(0,1,length.out=8))
all_ccds$quant_RPKM_ribo<-cut(x=all_ccds$RPKM_ribo,breaks=quantiles_RPKM_ribo,labels=as.character(1:7))
quantiles_length<-quantile(all_ccds$length.x,probs=seq(0,1,length.out=8))
all_ccds$quant_length<-cut(x=all_ccds$length.x,breaks=quantiles_length,labels=as.character(1:7))

length_rpkm<-rbind(c(2,2),c(2,4),c(2,6),c(4,2),c(4,4),c(4,6),c(6,2),c(6,4),c(6,6))
rownames(length_rpkm)<-c("short_low","short_med","short_high","medium_low","medium_med","medium_high","long_low","long_med","long_high")
colnames(length_rpkm)<-c("length","rpkm")
results<-list()
for(i in 1:dim(length_rpkm)[1]){
        combin<-length_rpkm[i,]
        name<-rownames(length_rpkm)[i]
        exons_all<-all_ccds[all_ccds[,"quant_length"]==combin["length"] & all_ccds[,"quant_RPKM_ribo"]==combin["rpkm"],]
        if(dim(exons_all)[1]>10){
                res<-as.data.frame(t(as.matrix(table(exons_all[,"pval_multit_3nt_ribo"]<0.05)/dim(exons_all)[1])),stringsAsFactors=F)
                if(dim(res)[2]==2){
                        colnames(res)<-c("non-periodic","periodic")
                }
                if(dim(res)[2]==1){
                        if(names(res)[1]==FALSE){
                                res[,2]<-0
                                colnames(res)<-c("non-periodic","periodic")
                        }
                        if(names(res)[1]==TRUE){
                                res[,2]<-res[,1]
                                res[,1]<-0
                                colnames(res)<-c("non-periodic","periodic")
                                
                        }
                }
                res_rna<-as.data.frame(t(as.matrix(table(exons_all[,"pval_multit_3nt_rna"]<0.05)/sum(!is.na(exons_all$pval_multit_3nt_rna)))),stringsAsFactors=F)
                
                if(dim(res_rna)[2]==2){
                        colnames(res_rna)<-c("non-periodic","periodic")
                }
                if(dim(res_rna)[2]==1){
                        if(names(res)[1]==TRUE){
                                res_rna[,2]<-0
                                colnames(res_rna)<-c("non-periodic","periodic")
                        }
                        if(names(res)[1]==FALSE){
                                res_rna[,2]<-res_rna[,1]
                                res_rna[,1]<-0
                                colnames(res_rna)<-c("non-periodic","periodic")
                                
                        }
                }
                res<-cbind(res,res_rna)
                res[,"n_exons"]<-dim(exons_all)[1]
                res[,"RPKM"]<-paste(paste(round(quantiles_RPKM_ribo[combin[2]],digits=1),round(quantiles_RPKM_ribo[combin[2]+1],digits=1),sep="-"),"RPKM")
                res[,"length"]<-paste(paste(round(quantiles_length[combin[1]],digits=1),round(quantiles_length[combin[1]+1],digits=1),sep="-"),"nt")
                res[,"category"]<-name
        }
        if(dim(exons_all)[1]<=10){
                res<-NULL
        }
        results[[i]]<-res
        
}
results<-do.call(rbind.data.frame,args=results)

results$length<-factor(results$length, levels=unique(results$length))
results$RPKM<-factor(results$RPKM, levels=unique(results$RPKM))

###

pdf(file="quality_check_plots.pdf",width=35,height=25,onefile=T,title="")


lefts<-c(0,.25,.5,.75, 0,.25,.5,.75, 0,.5, 0,.5, 0,.5)

rights<-c(.25,.5,.75,1, .25,.5,.75,1, .5,1, .5,1, .5,1)

bottoms<-c(.75,.75,.75,.75, .5,.5,.5,.5, .33,.33,.166,.166,0,0)

tops<-c(1,1,1,1, .75,.75,.75,.75, .5,.5,.33,.33,.166,.166)


matfig<-(cbind(lefts,rights,bottoms,tops))

close.screen(a=T)
#par(mgp=c(13, 1, 0))
n_ccds<-length(which(all_annot$type=="ccds"))
n_ccds_5nt<-dim(all_ccds)[1]
n_ccds_5nt_rna<-length(which(!is.na(ccds$chisq_rna)))


n_ccds_period<-length(which(all_ccds$pval_multit_3nt_ribo<0.05))
n_ccds_chisq<-length(which(all_ccds$chisq_ribo<0.05))
n_ccds_period_rna<-length(which(all_ccds$pval_multit_3nt_rna<0.05))
n_ccds_chisq_rna<-length(which(all_ccds$chisq_rna<0.05))

split.screen(matfig)

screen(1)
par(mar=c(6.1,8,2,2))

barp<-barplot(as.numeric(c(ribo_best,ribo_unique,ribo_psit)),beside=T,col=c("orange","red","dark red"),names.arg=c(""),cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.8,mgp=c(13, 1, 0),main="")
axis(side=1,labels=c("Ribo\naligned reads","Ribo\nunique reads","P-sites\npositions"),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.1,0))
screen(2)
par(mar=c(6.1,8,2,2))

barp<-barplot(as.numeric(c(rna_best,rna_unique,rna_psit)),beside=T,col=c("white","grey","dark grey"),names.arg=c(""),cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.8,main="")
axis(side=1,labels=c("RNA\naligned reads","RNA\nunique reads","RNA-sites\npositions"),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.1,0))
screen(3)
par(mar=c(6.1,8,2,2))

barp<-barplot(c(n_ccds,n_ccds_5nt,n_ccds_5nt_rna),names.arg=c(""),col=c("white","indianred2","red"),cex.main=1.8,cex.axis=1.8,cex.lab=1.8,cex.names=1.8,main="")
axis(side=1,labels=c("all ccds\nexons","ccds exons\n>5 P-sites","ccds exons\n>5 RNA & P-sit"),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.1,0))
screen(4)
par(mar=c(6.1,8,2,2))

rna_ok<-length(which(!is.na(ccds$pval_multit_3nt_rna)))
rna_ok2<-length(which(!is.na(ccds$chisq_rna)))
m_r<-table(ccds$pval_multit_3nt_ribo<0.05)/dim(ccds)[1]
m_rn<-table(ccds$pval_multit_3nt_rna<0.05)/rna_ok
# 
c_r<-table(ccds$chisq_ribo<0.05)/dim(ccds)[1]
c_rn<-table(ccds$chisq_rna<0.05)/rna_ok2


barp<-barplot(c(m_r[2],c_r[2],m_rn[1],c_rn[1]),xpd=F,col=c("red","red","grey","grey"),space=0.1,names.arg="",ylab="% CCDS exons",main="",cex.main=1.8,cex.axis=1.8,cex.lab=1.8,cex.names=1.8)
axis(side=1,labels=c("Multitap\nribo","Chi-sq\nribo","Multitap\nrna","Chi-sq\nrna"),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.1,0))

#barplot(c(m_rn[1],c_rn[1]),xpd=F,col="grey",space=0.1,names.arg=c("Multi-taper test","Chi-squared test"),ylab="% CCDS exons",main="Negative exons, (P-value > 0.05) \n This_study RNA-seq",cex.main=1.8,cex.axis=1.8,cex.lab=1.8,cex.names=1.8)
screen(5)
par(mar=c(6.1,8,2,2))

hist(ccds$pval_multit_3nt_rna,breaks=50,col="grey",main="",cex.main=1.8,cex.axis=1.8,cex.lab=1.8,xlab="")
axis(side=1,labels="P-values multitaper test \n RNA-seq",at=.5,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,5,0))
screen(6)
par(mar=c(6.1,8,2,2))

hist(ccds$chisq_rna,breaks=50,col="grey",main="",xlab="",cex.main=1.8,cex.axis=1.8,cex.lab=1.8)
axis(side=1,labels="P-values Chi-squared test \n RNA-seq",at=.5,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,5,0))
screen(7)
par(mar=c(6.1,8,2,2))

gino<-density(apply(ccds[,c("pctPhase_frame","pctPhase_frame_1","pctPhase_frame_2")],FUN=max,1),from=0,to=1)
plot(gino,col="violet",main="",xlab="% of P-sites on the max frame",cex.main=2,cex.axis=2,lwd=5,cex.lab=2)
if(dim(utrs)[1]>0){
        gino<-density(apply(utrs[,c("pctPhase_frame","pctPhase_frame_1","pctPhase_frame_2")],FUN=max,1),from=0,to=1)
        lines(gino,col="orange",lwd=5)
}
gino<-density(apply(noncoding[,c("pctPhase_frame","pctPhase_frame_1","pctPhase_frame_2")],FUN=max,1),from=0,to=1)
lines(gino,col="dark grey",lwd=5)
legend("topleft",c("CCDS exons","UTRs","non-coding"),lty=c(1,1,1),col=c("violet","orange","dark grey"),lwd=c(2.2,2.2,2.2),cex=1.8)


screen(8)
par(mar=c(6.1,1,2,2))

same_as_annot<-table(ccds_frames$ok_annot)/dim(ccds_frames)[1]
names(same_as_annot)<-c("diff_annot","same_annot")
same_as_annot<-round(same_as_annot,digits=3)
pie(same_as_annot,labels=paste(names(same_as_annot),":\n",as.character(100 * same_as_annot),"%"),col=c("dark grey","violet"),cex=1.8,cex.main=2,main="",init.angle=270)



for(j in 1:3){
        to_barpl<-split.data.frame(results,f=results$length)
        barpl<-(to_barpl[[j]])
        if(j==1){screen(9)}
        if(j==2){screen(11)}
        if(j==3){screen(13)}
        par(mar=c(6.1,8,2,2))
        
        barp<-barplot(t(100*as.matrix(barpl[,2:1])),ylim=c(0,100),col=c("dark red","grey"),names.arg=rep("",dim(barpl)[1]),ylab="% periodic CCDS exons\nRibo-seq",main=paste(barpl$length[1],"exons"),cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.8)
        axis(side=1,labels=paste(barpl$RPKM,"\n n_of exons=",barpl$n_exons),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.7,0))
        
        #         if(j==1){
        #                 mtext(side=3,"3nt periodicity in ccds exons, Ribo-seq",line=+2.6)
        #         }
        if(j==1){screen(10)}
        if(j==2){screen(12)}
        if(j==3){screen(14)} 
        par(mar=c(6.1,8,2,2))
        
        barp<-barplot(t(100*as.matrix(barpl[,4:3])),ylim=c(0,100),col=c("dark red","dark grey"),names.arg=rep("",dim(barpl)[1]),ylab="% periodic CCDS exons\nRNA-seq",main=paste(barpl$length[1],"exons"),cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.8)
        axis(side=1,labels=paste(barpl$RPKM,"\n n_of exons=",barpl$n_exons),at=barp,cex.axis=1.8,cex.main=1.8,cex=1.8,mgp=c(3,2.7,0))
        
        #         if(j==1){
        #                 mtext(side=3,"3nt periodicity in ccds exons, RNA-seq",line=+2.6)
        #         }
}
###


dev.off()

print(paste("--- QC plots Done! ---",date(),sep=" "))
