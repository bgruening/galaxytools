#!/usr/bin/Rscript

all_ex<-read.table("all_exons.bed",stringsAsFactors=F,header=F)
spli<-split.data.frame(all_ex,f=all_ex$V5)

spli2<-lapply(spli,FUN=function(x){
        minc<-min(x[,2])
        maxc<-max(x[,3])
        data.frame(chr=x[1,1],start=minc,end=maxc,le=maxc-minc,gene_id=x[1,5],strand=x[1,6],stringsAsFactors=F)
})
spli3<-do.call(what=rbind.data.frame,args=spli2)
write.table(file="genes_start_end",x=spli3,col.names=F,row.names=F,quote=F,sep="\t")
system("sort -k1,1 -k2,2n genes_start_end > genes_start_end.bed ")
system("rm genes_start_end")