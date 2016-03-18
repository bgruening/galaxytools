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


###script for CCDS genes ORF Finding, takes as arguments annotation dir, RiboTaper scripts dir, bedtools dir, n of cores



args <- commandArgs(trailingOnly = TRUE)

print(paste("--- CCDS ORF finding","---",date(),sep=" "))
###loads functions

suppressMessages(source(paste(args[2],"functions.R",sep = "/")))

###takes n of cores

registerDoMC(args[4])

###loads annotation files

annot<-paste(args[1],"cds_coords_transcripts",sep = "/")
cdss_transcripts<-read.table(annot,stringsAsFactors=F,header=F)
colnames(cdss_transcripts)<-c("transcript_id","annotated_start","annotated_stop")

###loads exonic results


results_ccds_ORFs<-read.table("all_calculations_ccdsgenes_annot_new",sep="\t",quote = "",stringsAsFactors=F,header=T)
all_annot_notCCDS<-results_ccds_ORFs[results_ccds_ORFs[,"type"]!="ccds",]

###calculates coordinates for sequence search



fives_threes<-all_annot_notCCDS[all_annot_notCCDS[,"type"]%in%c("3_utrs_ex","3_utrs_st","5_utrs_ex","5_utrs_st"),]
fives_threes_nonov<-fives_threes[is.na(fives_threes["overlapping_ccds_start"]),]
fives_threes_ov<-fives_threes[!is.na(fives_threes["overlapping_ccds_start"]),]

fives_threes_ok<-rbind(fives_threes_nonov,fives_threes_ov[fives_threes_ov[,"type"]%in%c("3_utrs_st","5_utrs_st"),])


###reads data tracks


all_tracks_ccds<-readBigText("data_tracks/Psit_Ribo_Rna_Cent_tracks_ccds")

index_ccds<-read.table("data_tracks/index_tracks_ccds",stringsAsFactors=F,header=F)
colnames(index_ccds)<-"exon_id"


all_tracks_exonsccds<-readBigText("data_tracks/Psit_Ribo_Rna_Cent_tracks_exonsccds")

index_exonsccds<-read.table("data_tracks/index_tracks_exonsccds",stringsAsFactors=F,header=F)
colnames(index_exonsccds)<-"exon_id"

tr_ex<-paste(args[1],"transcr_exons_ccds.bed",sep = "/")

transcr_ccds<-read.table(tr_ex,stringsAsFactors=F,header=F)
colnames(transcr_ccds)<-c("chr","start","end","transcript_id","gene_id","strand")
transcr_ccds$coords_id<-paste(transcr_ccds[,1],transcr_ccds[,2],transcr_ccds[,3],sep="_")
transcr_ccds$coords_ok<-paste(transcr_ccds[,1],transcr_ccds[,2],transcr_ccds[,3],transcr_ccds[,6],sep="_")
results_ccds_ORFs$coords_ok<-paste(results_ccds_ORFs$coords,results_ccds_ORFs$strand.x,sep="_")

###selects transcript with >2 reads

transcr_ccds<-unique(merge(transcr_ccds,results_ccds_ORFs[,c("coords_ok","P_sites_sum")],by="coords_ok",all.x=T))
transcr_sites<-aggregate(transcr_ccds$P_sites_sum,by=list(transcr_ccds$transcript_id),FUN=sum)
colnames(transcr_sites)<-c("transcript_id","n_P_sites")
transcr_sites<-transcr_sites[transcr_sites[,"n_P_sites"]>2,]
transcr_sites<-unique(transcr_sites$transcript_id)
###checks  for CCDS transcripts (if available)


if(sum(list.files(path =args[1])=="transcr_exons_ccds_ccdsid.bed")>0){
        tr_cc_app_ccdsid<-paste(args[1],"transcr_exons_ccds_ccdsid.bed",sep = "/")
        transcr_ccds_ccdsid<-read.table(tr_cc_app_ccdsid,stringsAsFactors=F,header=F)
        colnames(transcr_ccds_ccdsid)<-c("chr","start","end","transcript_id","gene_id","strand")
        transcr_ccds_ccdsid$coords_id<-paste(transcr_ccds_ccdsid[,1],transcr_ccds_ccdsid[,2],transcr_ccds_ccdsid[,3],sep="_")
        transcr_ccds_ccdsid$exon_id<-paste(transcr_ccds_ccdsid$coords_id,"EXONCCDS",transcr_ccds_ccdsid[,5],sep="_")
        transcript_ccds_transl_uORF<-unique(transcr_ccds_ccdsid[transcr_ccds_ccdsid[,"exon_id"]%in%fives_threes_ok[,"exon_id"],"transcript_id"])        
        transcript_ccds_transl_uORF<-transcript_ccds_transl_uORF[!is.na(transcript_ccds_transl_uORF)]
        if(sum(list.files(path =args[1])=="transcr_exons_ccds_appris.bed")==0){
                transcr_sites<-unique(c(transcript_ccds_transl_uORF))
                
        }
        
}

###checks  for APPRIS transcripts (if available)


if(sum(list.files(path =args[1])=="transcr_exons_ccds_appris.bed")>0){
        tr_cc_app<-paste(args[1],"transcr_exons_ccds_appris.bed",sep = "/")
        transcr_ccds_appr<-read.table(tr_cc_app,stringsAsFactors=F,header=F)
        colnames(transcr_ccds_appr)<-c("chr","start","end","transcript_id","gene_id","strand")
        transcr_ccds_appr$coords_id<-paste(transcr_ccds_appr[,1],transcr_ccds_appr[,2],transcr_ccds_appr[,3],sep="_")
        transcr_ccds_appr$exon_id<-paste(transcr_ccds_appr$coords_id,"EXONnonCCDS",transcr_ccds_appr[,5],sep="_")
        transcr_ccds_appr$coords2<-paste(transcr_ccds_appr$chr,":",transcr_ccds_appr$start,"-",transcr_ccds_appr$end,"(",transcr_ccds_appr$strand,")",sep="")
        #see prev versions to change this
        #transcript_ccds_transl<-results_ccds_ORFs[results_ccds_ORFs[,"pval_multit_3nt_ribo"]<0.05 & results_ccds_ORFs[,"P_sites_sum"]>5 ,]
        transcript_ccds_transl<-transcr_ccds_appr[!is.na(transcr_ccds_appr[,"gene_id"]),]
        transcr_sites<-unique(transcript_ccds_transl$transcript_id)
        if(sum(list.files(path =args[1])=="transcr_exons_ccds_ccdsid.bed")>0){
                transcript_ccds_transl<-results_ccds_ORFs[results_ccds_ORFs[,"pval_multit_3nt_ribo"]<0.05 & results_ccds_ORFs[,"P_sites_sum"]>5 ,]
                transcript_ccds_transl<-transcript_ccds_transl[!is.na(transcript_ccds_transl[,"gene_id"]),]
                transcript_ccds_transl<-unique(transcr_ccds_appr[transcr_ccds_appr[,"coords_id"]%in%transcript_ccds_transl[,"coords"],"transcript_id"])        
                transcript_ccds_transl<-transcript_ccds_transl[!is.na(transcript_ccds_transl)]
                transcr_sites<-unique(c(transcript_ccds_transl,transcript_ccds_transl_uORF))
                
        }
        
}


#reduce the search space to enhance speed



index_coords_ccds<-sapply(strsplit(index_ccds$exon_id,split="_"),function(x){paste(x[1],x[2],x[3],sep="_")})
index_coords_exonsccds<-sapply(strsplit(index_exonsccds$exon_id,split="_"),function(x){paste(x[1],x[2],x[3],sep="_")})


#index_coords_ccds<-index_coords_ccds[index_coords_ccds%in%transcr_ccds_fin_ids]

#index_coords_exonsccds<-index_coords_exonsccds[index_coords_exonsccds%in%transcr_ccds_fin_ids]

if(sum(list.files(path =args[1])=="transcr_exons_ccds_ccdsid.bed")==1){
        transcr_sites<-unique(c(transcript_ccds_transl_uORF,transcr_sites))
}

transcr_sites<-transcr_sites[!is.na(transcr_sites)]

transcr_ccds_fin<-transcr_ccds[transcr_ccds[,"transcript_id"]%in%transcr_sites,]
transcr_ccds_fin_ids<-unique(unlist(transcr_ccds_fin[,"coords_id"]))
# 
# all_tracks_ccds<-all_tracks_ccds[(index_coords_ccds%in%transcr_ccds_fin_ids)]
# index_ccds<-subset(index_ccds,index_coords_ccds%in%transcr_ccds_fin_ids)
# 
# all_tracks_exonsccds<-all_tracks_exonsccds[(index_coords_exonsccds%in%transcr_ccds_fin_ids)]
# index_exonsccds<-subset(index_exonsccds,index_coords_exonsccds%in%transcr_ccds_fin_ids)
# 

st_st_NA<-data.frame(start_pos=NA,st2vect=NA)
st_st_NA$ORF_frame<-NA
st_st_NA$ORF_length<-NA
st_st_NA$ORF_P_sites<-NA
st_st_NA$ORF_Psit_pct_in_frame<-NA
st_st_NA$ORF_RNA_sites<-NA
st_st_NA$ORF_RNAsit_pct_in_frame<-NA
st_st_NA$ORF_freq_multi_ribo<-NA
st_st_NA$ORF_pval_multi_ribo<-NA
st_st_NA$ORF_spec_multi_ribo<-NA
st_st_NA$ORF_freq_multi_rna<-NA
st_st_NA$ORF_pval_multi_rna<-NA
st_st_NA$ORF_spec_multi_rna<-NA
st_st_NA$ORF_freq3_fft_ribo<-NA
st_st_NA$ORF_spec3_fft_ribo<-NA
st_st_NA$ORF_freq3_spec_ribo<-NA
st_st_NA$ORF_spec3_spec_ribo<-NA
st_st_NA$ORF_freq3_fft_rna<-NA
st_st_NA$ORF_spec3_fft_rna<-NA
st_st_NA$ORF_freq3_spec_rna<-NA
st_st_NA$ORF_spec3_spec_rna<-NA
st_st_NA$ORF_ORF_score_ribo<-NA
st_st_NA$ORF_ORF_score_rna<-NA
st_st_NA$ORF_chisq_ribo<-NA
st_st_NA$ORF_chisq_rna<-NA
st_st_NA$ORF_Ribo_cov_aver<-NA
st_st_NA$ORF_RNA_cov_aver<-NA
st_st_NA$ORF_pept<-NA
st_st_NA$nt_tocheck_next_start<-0
st_st_NA$pval_next_start<-1
st_st_NA$P_sites_next_start<-0
st_st_NA$pct_P_sites_inframe_next_start<-0
st_st_NA$Method<-NA
st_st_NA$to_check<-NA
st_st_NA$to_check_rem<-NA
st_st_NA$ORF_id_tr<-NA
st_st_NA$ORF_id_gen<-NA
st_st_NA$to_check_ALL<-NA

CCDS_orfs<-foreach(j=1:length(transcr_sites),.combine=rbind,.multicombine=T) %dopar%{
        transcript<-transcr_sites[j]
        
        ###assembles transcript
        
        exons_in_transcr<-transcr_ccds[transcr_ccds[,"transcript_id"]==transcript,]
        #order exons
        exons_in_transcr<-exons_in_transcr[order(exons_in_transcr$start,decreasing=F),]
        list_exons_transcr<-list()
        list_exons_seqs<-list()
        
        for(k in seq(1,dim(exons_in_transcr)[1])){
                exon_track<-c()
                subs_ccds<-index_coords_ccds==exons_in_transcr[k,"coords_id"]
                if(sum(subs_ccds)>0){
                        if(sum(subs_ccds)==5){
                                exon_track<-all_tracks_ccds[subs_ccds]
                        }
                        if(sum(subs_ccds)>5){
                                exon_track<-all_tracks_ccds[which(subs_ccds)[1:5]]
                        }
                }
                if(length(exon_track)==0){
                        subs_exonsccds<-index_coords_exonsccds==exons_in_transcr[k,"coords_id"]
                        if(sum(subs_exonsccds)==5){
                                
                                exon_track<-all_tracks_exonsccds[subs_exonsccds]
                        }
                        if(sum(subs_exonsccds)>5){
                                exon_track<-all_tracks_exonsccds[which(subs_exonsccds)[1:5]]
                                
                        }
                        
                        
                }
                
                withsep<-strsplit(exon_track,split=" ")
                x<-t(data.frame(withsep))
                
                strand<-x[1,2]
                tracks<-t(x[,-c(1:2)])
                
                colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent","Seq")
                seq<-tracks[,5]
                tracks<-tracks[,1:4]
                mode(tracks)<-"numeric"
                length<-dim(tracks)[1]
                list_exons_transcr[[k]]<-tracks
                list_exons_seqs[[k]]<-seq
                
        }
        
        merged_tracks<-do.call(what=rbind,list_exons_transcr)
        
        if(strand=="-"){
                merged_tracks<-cbind(rev(merged_tracks[,1]),rev(merged_tracks[,2]),rev(merged_tracks[,3]),rev(merged_tracks[,4]))
        }
        
        tracks<-merged_tracks
        length<-dim(tracks)[1]
        
        if(strand=="+"){
                seq_transcr<-unlist(list_exons_seqs)
        }
        if(strand=="-"){
                
                seq_transcr<-unlist(list_exons_seqs)
                seq_transcr<-comp(rev((seq_transcr)),forceToLower=F)
        }
        transcr_data<-data.frame(transcript_id=transcript,stringsAsFactors=F)
        transcr_data$gene_id<-unique(transcr_ccds[transcr_ccds[,"transcript_id"]==transcript,"gene_id"])[1]
        transcr_data$annotation<-unique(results_ccds_ORFs[results_ccds_ORFs[,"gene_id"]==transcr_data$gene_id,"annotation",])[1]
        transcr_data$gene_symbol<-unique(results_ccds_ORFs[results_ccds_ORFs[,"gene_id"]==transcr_data$gene_id,"gene_symbol",])[1]
        
        P_sites_sum<-sum(tracks[,1])
        RNA_sites_sum<-sum(tracks[,4])
        transcr_data$strand<-strand
        transcr_data$length<-length
        transcr_data$n_exons<-dim(exons_in_transcr)[1]
        transcr_data$P_sites_sum<-P_sites_sum
        transcr_data$RNA_sites<-RNA_sites_sum
        transcr_data$Ribo_cov_aver<-mean(tracks[,2])
        transcr_data$RNA_cov_aver<-mean(tracks[,3])
        
        transcr_data$freq_multit_3nt<-NA
        transcr_data$pval_multit_3nt<-NA
        transcr_data$spec_multit_3nt<-NA
        if(P_sites_sum>2 & length>5){
                if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                vals_mtm<-take_freqs_Fvalues_all_around_3nt_spec(n_tapers=24,time_bw=12,tracks[,1],slepians_values=slepians)[c(1,6,7)]
                transcr_data$freq_multit_3nt<-vals_mtm[1]
                transcr_data$pval_multit_3nt<-vals_mtm[2]
                transcr_data$spec_multit_3nt<-vals_mtm[3]
                
        }
        
        Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
        Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
        Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
        
        transcr_data$chisq_noccds_psit<-NA
        if(P_sites_sum>15){
                transcr_data$chisq_noccds_psit<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
        if(P_sites_sum<16 & P_sites_sum>0){
                transcr_data$chisq_noccds_psit<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }    
        pctPhase_frame<-Phase_P_sites_frame/P_sites_sum
        pctPhase_frame_1<-Phase_P_sites_frame_1/P_sites_sum
        pctPhase_frame_2<-Phase_P_sites_frame_2/P_sites_sum
        
        Centered_sites_sum<-round(sum(tracks[,4]),digits=6)
        
        Phase_Centered_sites_frame<-sum(tracks[seq(1,length,by=3),4])
        Phase_Centered_sites_frame_1<-sum(tracks[seq(2,length,by=3),4])
        Phase_Centered_sites_frame_2<-sum(tracks[seq(3,length,by=3),4])
        
        pctPhaseCentered_frame<-Phase_Centered_sites_frame/Centered_sites_sum
        pctPhaseCentered_frame_1<-Phase_Centered_sites_frame_1/Centered_sites_sum
        pctPhaseCentered_frame_2<-Phase_Centered_sites_frame_2/Centered_sites_sum
        
        transcr_data$chisq_noccds_rna<-NA
        if(Centered_sites_sum>15){
                chisq_rna<-chisq.test(as.table(c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2)))$p.value}
        if(Centered_sites_sum<16 & Centered_sites_sum>0){
                chisq_rna<-xmulti(obs=c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }
        
        
        MAXPhase_frame<-max(c(pctPhase_frame,pctPhase_frame_1,pctPhase_frame_2))
        FRAME_MAX_phase<-max.col(t(c(pctPhase_frame,pctPhase_frame_1,pctPhase_frame_2)))-1
        
        MAXPhaseCentered_frame<-max(c(pctPhaseCentered_frame,pctPhaseCentered_frame_1,pctPhaseCentered_frame_2))
        FRAME_MAX_phaseCentered<-max.col(t(c(pctPhaseCentered_frame,pctPhaseCentered_frame_1,pctPhaseCentered_frame_2)))-1
        
        frame_start_pred<-FRAME_MAX_phase
        frame_end_pred<-(length-(FRAME_MAX_phase+1))%%3
        
        ###Finds ORFs on the 3 different frames
        
        all_sign_frames<-list()
        for(u in 0:2){
                
                pept<-NA
                pept<-unlist(getTrans(seq_transcr,sens="F",frame=u))
                
                starts<-pept=="M"
                
                stops<-pept=="*"
                transcr_data$orf_position<-"undetected"
                
                start_pos<-((1:length(pept))[starts])*3
                if(length(start_pos)>0){
                        start_pos<-start_pos+u-2
                } else {start_pos<-NA}
                
                stop_pos<-((1:length(pept))[stops])*3
                if(length(stop_pos)>0){
                        stop_pos<-stop_pos+u-2
                } else {stop_pos<-NA}
                
                #NAs
                if(sum(!is.na(start_pos))==0 | sum(!is.na(stop_pos))==0){
                        st_st<-st_st_NA
                        transcr_data_fr_sORFs<-cbind(transcr_data,st_st_NA)
                }
                
                if(sum(!is.na(start_pos))>0 & sum(!is.na(stop_pos))>0){
                        st2vect<-c()
                        for(h in 1:length(start_pos)){
                                st1<-start_pos[h]
                                diff<-stop_pos-st1
                                diff<-diff[diff>0]
                                if(length(diff)>0){st2<-st1+min(diff)}
                                if(length(diff)==0){st2<-NA}
                                st2vect[h]<-st2
                                
                        }
                        st_st<-data.frame(cbind(start_pos,st2vect))
                        
                        st_st<-st_st[!is.na(st_st[,"st2vect"]),]
                        if(dim(st_st)[1]>0){
                                if(dim(st_st)[1]==1){
                                        list_coords=list()
                                        list_coords[[1]]<-st_st[,1]:st_st[,2]
                                }
                                if(dim(st_st)[1]>1){
                                        list_coords<-apply(st_st,FUN=function(x){x[1]:x[2]},1)
                                }
                                
                                max_period<-NA
                                start_pos<-NA
                                stop_pos<-NA
                                pval_max_period<-NA
                        }
                        if(dim(st_st)[1]>0){
                                st_st$ORF_frame<-u
                                st_st$ORF_length<-NA
                                st_st$ORF_P_sites<-NA
                                st_st$ORF_Psit_pct_in_frame<-NA
                                st_st$ORF_RNA_sites<-NA
                                st_st$ORF_RNAsit_pct_in_frame<-NA
                                st_st$ORF_freq_multi_ribo<-NA
                                st_st$ORF_pval_multi_ribo<-NA
                                st_st$ORF_spec_multi_ribo<-NA
                                st_st$ORF_freq_multi_rna<-NA
                                st_st$ORF_pval_multi_rna<-NA
                                st_st$ORF_spec_multi_rna<-NA
                                
                                st_st$ORF_freq3_fft_ribo<-NA
                                st_st$ORF_spec3_fft_ribo<-NA
                                st_st$ORF_freq3_spec_ribo<-NA
                                st_st$ORF_spec3_spec_ribo<-NA
                                st_st$ORF_freq3_fft_rna<-NA
                                st_st$ORF_spec3_fft_rna<-NA
                                st_st$ORF_freq3_spec_rna<-NA
                                st_st$ORF_spec3_spec_rna<-NA
                                st_st$ORF_ORF_score_ribo<-NA
                                st_st$ORF_ORF_score_rna<-NA
                                st_st$ORF_chisq_ribo<-NA
                                st_st$ORF_chisq_rna<-NA
                                st_st$ORF_Ribo_cov_aver<-NA
                                st_st$ORF_RNA_cov_aver<-NA
                                st_st$ORF_pept<-NA
                                st_st$Method<-NA
                                st_st$to_check<-NA
                                st_st$to_check_rem<-NA
                                st_st$ORF_id_tr<-NA
                                st_st$ORF_id_gen<-NA
                                st_st$to_check_ALL<-NA
                                for(r in 1:dim(st_st)[1]){
                                        tracks_stst<-tracks[st_st[r,1]:st_st[r,2],]
                                        length<-dim(tracks_stst)[1]
                                        P_sites_sum<-sum(tracks_stst[,1])
                                        RNA_sites_sum<-sum(tracks_stst[,4])
                                        st_st[r,"ORF_length"]<-length-1
                                        st_st[r,"ORF_P_sites"]<-P_sites_sum
                                        st_st[r,"ORF_RNA_sites"]<-RNA_sites_sum
                                        st_st[r,"ORF_Ribo_cov_aver"]<-mean(tracks_stst[,2])
                                        st_st[r,"ORF_RNA_cov_aver"]<-mean(tracks_stst[,3])
                                        if(P_sites_sum>5 & length>5){
                                                Phase_P_sites_frame<-sum(tracks_stst[seq(1,length,by=3),1])
                                                Phase_P_sites_frame_1<-sum(tracks_stst[seq(2,length,by=3),1])
                                                Phase_P_sites_frame_2<-sum(tracks_stst[seq(3,length,by=3),1])
                                                st_st[r,"ORF_Psit_pct_in_frame"]<-Phase_P_sites_frame/P_sites_sum
                                                if((Phase_P_sites_frame/P_sites_sum)>0.5){
                                                        score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                                                        score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                                                        score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                                                        
                                                        orfsc<-log2(score1+score2+score3+1)
                                                        st_st[r,"ORF_ORF_score_ribo"]<-orfsc
                                                        if(Phase_P_sites_frame<=Phase_P_sites_frame_1 | Phase_P_sites_frame<=Phase_P_sites_frame_2){
                                                                st_st[r,"ORF_ORF_score_ribo"]<--orfsc
                                                        }
                                                        
                                                        if(max(tracks_stst[,1])>(P_sites_sum*.7)){
                                                                new_track<-tracks_stst
                                                                new_track[which(new_track[,1]==max(new_track[,1]))]<-0
                                                                st_st[r,"ORF_ORF_score_ribo"]<-NA
                                                                if(sum(new_track[,1])>2){
                                                                        Phase_P_sites_frame_corr<-sum(new_track[seq(1,length,by=3),1])
                                                                        Phase_P_sites_frame_1_corr<-sum(new_track[seq(2,length,by=3),1])
                                                                        Phase_P_sites_frame_2_corr<-sum(new_track[seq(3,length,by=3),1])
                                                                        score1<-((Phase_P_sites_frame_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                                                        score2<-((Phase_P_sites_frame_1_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                                                        score3<-((Phase_P_sites_frame_2_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                                                        st_st[r,"ORF_ORF_score_ribo"]<-log2(score1+score2+score3+1)
                                                                        if(Phase_P_sites_frame_corr<=Phase_P_sites_frame_1_corr | Phase_P_sites_frame<=Phase_P_sites_frame_2_corr){
                                                                                st_st[r,"ORF_ORF_score_ribo"]<--log2(score1+score2+score3+1)
                                                                        }                                                                }
                                                        }
                                                        
                                                        if(P_sites_sum>15){
                                                                st_st[r,"ORF_chisq_ribo"]<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
                                                        if(P_sites_sum<16 & P_sites_sum>0){
                                                                st_st[r,"ORF_chisq_ribo"]<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
                                                        }  
                                                        if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                                                        if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                                                        values_mtm_orf<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks_stst[,1],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                                                        st_st[r,"ORF_freq_multi_ribo"]<-values_mtm_orf[1]
                                                        
                                                        st_st[r,"ORF_pval_multi_ribo"]<-values_mtm_orf[2]
                                                        st_st[r,"ORF_spec_multi_ribo"]<-values_mtm_orf[3]
                                                        fft_sp<-take_maxfreq_and_power_FFT_Spec(tracks_stst[,1])
                                                        st_st[,"ORF_freq3_fft_ribo"]<-fft_sp[1]
                                                        st_st[,"ORF_spec3_fft_ribo"]<-fft_sp[2]
                                                        st_st[,"ORF_freq3_spec_ribo"]<-fft_sp[3]
                                                        st_st[,"ORF_spec3_spec_ribo"]<-fft_sp[4]
                                                        
                                                        pept<-unlist(getTrans(seq_transcr[st_st[r,1]:st_st[r,2]],sens="F"))
                                                        st_st[r,"ORF_pept"]<-paste(pept,sep="",collapse="")
                                                }
                                                if(RNA_sites_sum>5 & (Phase_P_sites_frame/P_sites_sum)>0.5){
                                                        if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                                                        if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                                                        values_mtm_orf_rna<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks_stst[,4],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                                                        
                                                        st_st[r,"ORF_freq_multi_rna"]<-values_mtm_orf_rna[1]
                                                        st_st[r,"ORF_pval_multi_rna"]<-values_mtm_orf_rna[2]
                                                        st_st[r,"ORF_spec_multi_rna"]<-values_mtm_orf_rna[3]
                                                        
                                                        fft_sp<-take_maxfreq_and_power_FFT_Spec(tracks_stst[,4])
                                                        st_st[,"ORF_freq3_fft_rna"]<-fft_sp[1]
                                                        st_st[,"ORF_spec3_fft_rna"]<-fft_sp[2]
                                                        st_st[,"ORF_freq3_spec_rna"]<-fft_sp[3]
                                                        st_st[,"ORF_spec3_spec_rna"]<-fft_sp[4]
                                                        
                                                        Phase_Centered_sites_frame<-sum(tracks_stst[seq(1,length,by=3),4])
                                                        Phase_Centered_sites_frame_1<-sum(tracks_stst[seq(2,length,by=3),4])
                                                        Phase_Centered_sites_frame_2<-sum(tracks_stst[seq(3,length,by=3),4])
                                                        st_st[r,"ORF_RNAsit_pct_in_frame"]<-Phase_Centered_sites_frame/RNA_sites_sum
                                                        score1<-((Phase_Centered_sites_frame-P_sites_sum/3)^2)/(RNA_sites_sum/3)
                                                        score2<-((Phase_Centered_sites_frame_1-P_sites_sum/3)^2)/(RNA_sites_sum/3)
                                                        score3<-((Phase_Centered_sites_frame_2-P_sites_sum/3)^2)/(RNA_sites_sum/3)
                                                        
                                                        orfsc<-log2(score1+score2+score3+1)
                                                        st_st[r,"ORF_ORF_score_rna"]<-orfsc
                                                        if(Phase_Centered_sites_frame<=Phase_Centered_sites_frame_1 | Phase_Centered_sites_frame<=Phase_Centered_sites_frame_2){
                                                                st_st[r,"ORF_ORF_score_rna"]<--orfsc
                                                        }
                                                        
                                                        if(max(tracks_stst[,4])>(RNA_sites_sum*.7)){
                                                                new_track<-tracks_stst
                                                                new_track[which(new_track[,4]==max(new_track[,4]))]<-0
                                                                st_st[r,"ORF_ORF_score_rna"]<-NA
                                                                if(sum(new_track[,4])>2){
                                                                        Phase_Centered_sites_frame_corr<-sum(new_track[seq(1,length,by=3),4])
                                                                        Phase_Centered_sites_frame_1_corr<-sum(new_track[seq(2,length,by=3),4])
                                                                        Phase_Centered_sites_frame_2_corr<-sum(new_track[seq(3,length,by=3),4])
                                                                        score1<-((Phase_Centered_sites_frame_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                                                        score2<-((Phase_Centered_sites_frame_1_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                                                        score3<-((Phase_Centered_sites_frame_2_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                                                        st_st[r,"ORF_ORF_score_rna"]<-log2(score1+score2+score3+1)
                                                                        if(Phase_Centered_sites_frame_corr<=Phase_Centered_sites_frame_1_corr | Phase_Centered_sites_frame_corr<=Phase_Centered_sites_frame_2_corr){
                                                                                st_st[r,"ORF_ORF_score_rna"]<--log2(score1+score2+score3+1)
                                                                        }
                                                                }
                                                        }
                                                        
                                                        if(RNA_sites_sum>15){
                                                                st_st[r,"ORF_chisq_rna"]<-chisq.test(as.table(c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2)))$p.value}
                                                        if(RNA_sites_sum<16 & RNA_sites_sum>0){
                                                                st_st[r,"ORF_chisq_rna"]<-xmulti(obs=c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
                                                        }
                                                        
                                                }
                                        }
                                }
                                if(dim(st_st)[1]>0){st_st<-st_st[!is.na(st_st[,"ORF_pval_multi_ribo"]),]}
                                if(dim(st_st)[1]>0){st_st<-st_st[st_st[,"ORF_Psit_pct_in_frame"]>0.5,]}
                                if(dim(st_st)[1]>0){
                                        st_st$nt_tocheck_next_start<-0
                                        st_st$pval_next_start<-1
                                        st_st$P_sites_next_start<-0
                                        st_st$pct_P_sites_inframe_next_start<-0
                                        #find starts per each stop codon
                                        list_stopsorfs<-split.data.frame(x=st_st,f=st_st[,2],drop=T)
                                        
                                        transcr_data_fr<-transcr_data
                                        
                                        list_sORFs_frame_moretap<-list()
                                        list_sORFs_frame_bestperiod<-list()
                                        list_sORFs_frame_maxsit<-list()
                                        
                                        for(g in 1:length(list_stopsorfs)){
                                                
                                                stoplist<-list_stopsorfs[[g]]
                                                max_period<-stoplist[stoplist[,"ORF_pval_multi_ribo"]==min(stoplist[,"ORF_pval_multi_ribo"]),]
                                                list_sORFs_frame_bestperiod[[g]]<-max_period 
                                                stoplists_period<-stoplist[stoplist[,"ORF_pval_multi_ribo"]<0.05,]
                                                if(dim(stoplists_period)[1]>0){
                                                        stoplists_period<-stoplists_period[!is.na(stoplists_period[,"ORF_pval_multi_ribo"]),]
                                                }
                                                if(dim(stoplists_period)[1]>1){
                                                        
                                                        for(b in 1:(dim(stoplists_period)[1]-1)){
                                                                stoplists_period[b,"nt_tocheck_next_start"]<-stoplists_period[b+1,"start_pos"]-stoplist[b,"start_pos"]
                                                                tracks_stst<-tracks[stoplists_period[b,"start_pos"]:stoplists_period[b+1,"start_pos"],]
                                                                length<-dim(tracks_stst)[1]
                                                                P_sites_sum<-sum(tracks_stst[,1])
                                                                pval_to_next<-1                                                              
                                                                
                                                                Phase_P_sites_frame<-sum(tracks_stst[seq(1,length,by=3),1])
                                                                Phase_P_sites_frame_1<-sum(tracks_stst[seq(2,length,by=3),1])
                                                                Phase_P_sites_frame_2<-sum(tracks_stst[seq(3,length,by=3),1])
                                                                
                                                                pctPhase_frame<-Phase_P_sites_frame/P_sites_sum
                                                                pctPhase_frame_1<-Phase_P_sites_frame_1/P_sites_sum
                                                                pctPhase_frame_2<-Phase_P_sites_frame_2/P_sites_sum
                                                                
                                                                if(P_sites_sum>5){
                                                                        if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                                                                        if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                                                                        
                                                                        pval_to_next<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks_stst[,1],n_tapers=24,time_bw=12,slepians_values=slepians)[6]
                                                                }
                                                                stoplists_period[b,"P_sites_next_start"]<-P_sites_sum
                                                                
                                                                stoplists_period[b,"pct_P_sites_inframe_next_start"]<-pctPhase_frame
                                                                
                                                                stoplists_period[b,"pval_next_start"]<-pval_to_next
                                                        }
                                                        
                                                        max_sit<-stoplists_period[which(stoplists_period[,"P_sites_next_start"]>5 & stoplists_period[,"pct_P_sites_inframe_next_start"]>0.5)[1],]
                                                        max_sit<-max_sit[!is.na(max_sit[,"ORF_length"]),]
                                                        
                                                        if(dim(max_sit)[1]==0){
                                                                max_sit<-max_period
                                                        }
                                                        list_sORFs_frame_maxsit[[g]]<-max_sit   
                                                        
                                                        more_tap<-stoplists_period[which(stoplists_period[,"pval_next_start"]<0.05)[1],]
                                                        more_tap<-more_tap[!is.na(more_tap[,"ORF_length"]),]
                                                        if(dim(more_tap)[1]==0){
                                                                more_tap<-max_period
                                                        }
                                                        list_sORFs_frame_moretap[[g]]<-more_tap
                                                        
                                                }
                                                if(dim(stoplists_period)[1]<2){
                                                        list_sORFs_frame_maxsit[[g]]<-max_period
                                                        list_sORFs_frame_moretap[[g]]<-max_period
                                                        
                                                }
                                                
                                                
                                        }
                                        sORFs_frame_moretap<-do.call(what=rbind.data.frame,args=list_sORFs_frame_moretap)
                                        sORFs_frame_moretap$Method<-"more_tapers"
                                        sORFs_frame_maxsit<-do.call(what=rbind.data.frame,args=list_sORFs_frame_maxsit)
                                        sORFs_frame_maxsit$Method<-"max_P_sites"
                                        sORFs_frame_bestperiod<-do.call(what=rbind.data.frame,args=list_sORFs_frame_bestperiod)
                                        sORFs_frame_bestperiod$Method<-"best_periodicity"
                                        sORFs_frames<-rbind(sORFs_frame_moretap,sORFs_frame_maxsit,sORFs_frame_bestperiod)
                                        
                                        for(w in 1:dim(sORFs_frames)[1]){
                                                transcr_data_fr[w,]<-transcr_data_fr[1,]
                                        }
                                        
                                        transcr_data_fr_sORFs<-cbind(transcr_data_fr,sORFs_frames)
                                        transcr_data_fr_sORFs$orf_position<-"detected"
                                }
                        }
                        
                        if(dim(st_st)[1]==0){
                                st_st<-st_st_NA
                                transcr_data_fr_sORFs<-cbind(transcr_data,st_st_NA)
                        }
                }
                
                
                all_sign_frames[[u+1]]<-transcr_data_fr_sORFs
        }
        all_sign_frames<-do.call(what=rbind.data.frame,args=all_sign_frames)
        transcr_all_frames_res<-unique(all_sign_frames)
        transcr_all_frames_res$ORF_id_tr<-paste(transcr_all_frames_res$transcript_id,transcr_all_frames_res$start_pos,transcr_all_frames_res$st2vect,sep="_")
        transcr_all_frames_ok<-transcr_all_frames_res[!is.na(transcr_all_frames_res$ORF_pept),]
        if(dim(transcr_all_frames_ok)[1]>0){
                all_orfs<-unique(transcr_all_frames_ok[,c("transcript_id","length","strand","start_pos","st2vect","ORF_length","gene_id")])
                transcr<-all_orfs$transcript_id[1]
                trascr_length<-all_orfs$length[1]
                orf_strand<-all_orfs$strand[1]                
                ex_intr_coords<-exons_in_transcr$coords_id
                if(orf_strand=="-"){ex_intr_coords<-rev(ex_intr_coords)}
                
                exons_in_transcr_data<-results_ccds_ORFs[results_ccds_ORFs[,"coords"]%in%ex_intr_coords,]
                exons_in_transcr_data<-exons_in_transcr_data[match(ex_intr_coords,exons_in_transcr_data$coords),]
                cumsumexons<-cumsum(exons_in_transcr_data$length.x)
                
                list_orfas<-list()
                for(z in 1:dim(all_orfs)[1]){
                        orfa<-all_orfs[z,]
                        
                        transcr_data<-data.frame(transcript_id=transcr)
                        
                        orf_start<-orfa$start_pos
                        orf_end<-orfa$st2vect
                        
                        st_ex<-which((cumsumexons-orf_start)==min(cumsumexons[cumsumexons>orf_start]-orf_start))
                        end_ex<-which((cumsumexons-orf_end)==min(cumsumexons[cumsumexons>=orf_end]-orf_end))
                        in_betw_ex<-st_ex:end_ex
                        in_betw_ex<-in_betw_ex[!in_betw_ex%in%c(st_ex,end_ex)>0]
                        exon_inbetween_data<-exons_in_transcr_data[in_betw_ex,]
                        
                        
                        coord_start<-NA
                        coord_end<-NA
                        nt_to_rem<-NA
                        rem_len<-0
                        if(st_ex>1){rem_len<-cumsumexons[st_ex-1]}
                        if(orfa$strand=="+"){coord_start<-exons_in_transcr_data[st_ex,"start"] + (orf_start-rem_len)}
                        if(orfa$strand=="-"){coord_start<-exons_in_transcr_data[st_ex,"end"] - (orf_start-rem_len)}
                        
                        if(length(in_betw_ex)==0){
                                if(st_ex==end_ex){nt_to_rem<-0}
                                if(st_ex!=end_ex){if(orfa$strand=="+"){
                                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                                }
                                                  if(orfa$strand=="-"){
                                                          nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                                                  }
                                }
                        }
                        
                        if(length(in_betw_ex)>0){
                                nt_in_betw<-sum(exons_in_transcr_data[in_betw_ex,"length.x"])
                                if(orfa$strand=="+"){
                                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                                }
                                if(orfa$strand=="-"){
                                        nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                                }
                                nt_to_rem<-nt_to_rem+nt_in_betw
                        }
                        
                        if(st_ex==end_ex & orfa$strand=="+"){coord_end<-coord_start+orfa$ORF_length+1}
                        if(st_ex==end_ex & orfa$strand=="-"){coord_end<-coord_start-orfa$ORF_length+1}
                        
                        if(st_ex!=end_ex & orfa$strand=="+"){coord_end<-exons_in_transcr_data[end_ex,"start"] + (orfa$ORF_length-nt_to_rem)+1}
                        if(st_ex!=end_ex & orfa$strand=="-"){coord_end<-exons_in_transcr_data[end_ex,"end"] - (orfa$ORF_length-nt_to_rem)+1}
                        
                        if(orfa$strand=="-"){
                                coord_start2<-coord_start
                                coord_start<-coord_end
                                coord_end<-coord_start2
                        }
                        
                        
                        if(st_ex!=end_ex & orfa$strand=="+"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,exons_in_transcr_data[st_ex,"end"],"CCDS",orfa$gene_id,orfa$strand,sep="_")
                                                             to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],exons_in_transcr_data[end_ex,"start"],coord_end,"CCDS",orfa$gene_id,orfa$strand,sep="_")
                                                             to_check<-paste(to_check_st,to_check_end,sep=";")
                                                             
                        }
                        if(st_ex!=end_ex & orfa$strand=="-"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],exons_in_transcr_data[st_ex,"start"],coord_end,"CCDS",orfa$gene_id,orfa$strand,sep="_")
                                                             to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],coord_start,exons_in_transcr_data[end_ex,"end"],"CCDS",orfa$gene_id,orfa$strand,sep="_")
                                                             to_check<-paste(to_check_st,to_check_end,sep=";")
                        }
                        
                        if(st_ex==end_ex){to_check<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,"CCDS",orfa$gene_id,orfa$strand,sep="_")}
                        orfa$to_check<-to_check
                        orfa$to_check_rem<-NA
                        if(length(in_betw_ex)>0){
                                orfa$to_check_rem<-paste(exon_inbetween_data$exon_id,collapse=";")
                                
                        }
                        orfa$ORF_id_tr<-paste(transcr_data$transcript_id,orf_start,orf_end,sep="_")
                        orfa$ORF_id_gen<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,sep="_")
                        orfa$to_check_ALL<-paste(orfa$to_check,orfa$to_check_rem,sep=";")
                        list_orfas[[z]]<-orfa
                        
                        
                }
                list_orfas<-do.call(rbind.data.frame,args=list_orfas)
                transcr_all_frames_ok$ORF_id_gen<-NULL
                transcr_all_frames_ok$to_check<-NULL
                transcr_all_frames_ok$to_check_rem<-NULL
                transcr_all_frames_ok$to_check_ALL<-NULL
                
                transcr_all_frames_ok<-merge(transcr_all_frames_ok,list_orfas[,c("ORF_id_tr","ORF_id_gen","to_check","to_check_rem","to_check_ALL")],by="ORF_id_tr")
                #reconcile and maybe add the rest
                return(transcr_all_frames_ok)
        }
        if(dim(transcr_all_frames_ok)[1]==0){return(transcr_all_frames_res)}
        
        
}
CCDS_orfs_found<-CCDS_orfs[!is.na(CCDS_orfs[,"ORF_pept"]),]


CCDS_orfs<-merge(CCDS_orfs_found,cdss_transcripts,by="transcript_id",all.x=T)

write.table(CCDS_orfs,file="orfs_found",quote=F,row.names=F,sep="\t",col.names=T)

options(scipen=999)

CCDS_orfs$ORF_id_tr_minus2<-paste(CCDS_orfs$transcript_id,CCDS_orfs$start_pos,CCDS_orfs$st2vect+2,sep="_")
CCDS_orfs$ORF_id_tr_annotated<-paste(CCDS_orfs$transcript_id,CCDS_orfs$annotated_start,CCDS_orfs$annotated_stop,sep="_")


#nonccds_res<-results_ccds_ORFs
CCDS_orfs_uniq<-CCDS_orfs

print(paste("--- checking CCDS ORF coverage and multi-mapping ratio,",date(),sep=" "))

all_sORFs_CCDS_multi<-CCDS_orfs_uniq


ex_to_check<-strsplit(all_sORFs_CCDS_multi$to_check,split=";")

ex_to_check<-unique(unlist(ex_to_check))

ex_to_check_spl<-strsplit(ex_to_check,split="_")

bedfiles_to_check<-data.frame(chr=NA,start=NA,end=NA,type=NA,gene_id=NA,strand=NA)
for(h in 1:length(ex_to_check_spl)){
        to_bed<-ex_to_check_spl[[h]]
        bedfiles_to_check[h,"chr"]<-to_bed[1]
        bedfiles_to_check[h,"start"]<-as.numeric(to_bed[2])
        bedfiles_to_check[h,"end"]<-as.numeric(to_bed[3])
        bedfiles_to_check[h,"type"]<-to_bed[4]
        bedfiles_to_check[h,"gene_id"]<-to_bed[5]
        bedfiles_to_check[h,"strand"]<-to_bed[6]
        
}



write.table(bedfiles_to_check,file="bed_tocheck_ccds.bed",quote=F,row.names=F,sep="\t",col.names=F)

scr<-paste(args[2],"analyze_multi_clust.bash",sep="/")
syst_scr<-paste(scr,"bed_tocheck_ccds.bed bed_tocheck_ccds",args[3],sep = " ")
system(syst_scr)

scr<-paste(args[2],"include_multi_nomerge.R",sep="/")
syst_scr<-paste(scr,"bed_tocheck_ccds",sep = " ")

system(syst_scr)

res_to_check<-read.table(file="multi_table_bed_tocheck_ccds",header=T,stringsAsFactors=F)

dir.create("tmp_ccds", showWarnings = FALSE)

system("mv *tocheck_ccds* tmp_ccds/")

setwd("tmp_ccds")


ex_rem<-strsplit(as.character(all_sORFs_CCDS_multi$to_check_rem),split=";")

ex_rem<-unique(unlist(ex_rem))
ex_rem<-ex_rem[!is.na(ex_rem)]


res_ex_rem<-results_ccds_ORFs[results_ccds_ORFs[,"exon_id"]%in%ex_rem,c(c("exon_id","strand.x","length.y","reads_ribo","reads_multi_ribo","pct_region_covered_ribo","pct_covered_onlymulti_ribo","reads_rna","reads_multi_rna","pct_region_covered_rna","pct_covered_onlymulti_rna"))]
names(res_ex_rem)<-names(res_to_check)

res_all_multi<-rbind.data.frame(res_ex_rem,res_to_check)

res_all_multi$exon_id_2<-paste(res_all_multi$exon_id,res_all_multi$strand,sep="_")


all_sORFs_CCDS_multi_final<-foreach(g=1:(dim(all_sORFs_CCDS_multi)[1]),.combine=rbind,.multicombine=T) %dopar%{
        s<-all_sORFs_CCDS_multi[g,]
        list_ex<-strsplit(s$to_check_ALL,split=";")[[1]]
        with_exon2<-which(res_all_multi[,"exon_id_2"]%in%list_ex)
        with_exon1<-which(res_all_multi[,"exon_id"]%in%list_ex)
        to_take<-unique(c(with_exon2,with_exon1))
        res_multi<-res_all_multi[to_take,]
        res_multi$reads_ribo<-sum(res_multi$reads_ribo)
        res_multi$reads_multi_ribo<-sum(res_multi$reads_multi_ribo)
        res_multi$pct_region_covered_ribo_ALL<-res_multi$pct_region_covered_ribo*res_multi$length.y
        res_multi$pct_covered_onlymulti_ribo_ALL<-res_multi$pct_covered_onlymulti_ribo*res_multi$length.y
        res_multi$pct_region_covered_ribo<-sum(res_multi$pct_region_covered_ribo_ALL)/(sum(res_multi$length.y))
        res_multi$pct_covered_onlymulti_ribo<-sum(res_multi$pct_covered_onlymulti_ribo_ALL)/(sum(res_multi$length.y))
        res_multi$reads_rna<-sum(res_multi$reads_rna)
        res_multi$reads_multi_rna<-sum(res_multi$reads_multi_rna)
        res_multi$pct_region_covered_rna_ALL<-res_multi$pct_region_covered_rna*res_multi$length.y
        res_multi$pct_covered_onlymulti_rna_ALL<-res_multi$pct_covered_onlymulti_rna*res_multi$length.y
        res_multi$pct_region_covered_rna<-sum(res_multi$pct_region_covered_rna_ALL)/sum(res_multi$length.y)
        res_multi$pct_covered_onlymulti_rna<-sum(res_multi$pct_covered_onlymulti_rna_ALL)/sum(res_multi$length.y)
        
        s<-cbind(s,res_multi[1,])
        s
}

print(paste("--- Selecting best transcript per CCDS ORF,",date(),sep=" "))


write.table(all_sORFs_CCDS_multi_final,file="orfs_bef_ag",quote=F,row.names=F,sep="\t",col.names=T)


agg<-aggregate(x=all_sORFs_CCDS_multi_final[,"RNA_sites"],by=list(all_sORFs_CCDS_multi_final[,"gene_id"],all_sORFs_CCDS_multi_final[,"ORF_pept"],all_sORFs_CCDS_multi_final[,"Method"]),FUN=max)
names(agg)<-c("gene_id","ORF_pept","Method","RNA_sites")
agg2<-merge(x=all_sORFs_CCDS_multi_final[,c("ORF_id_tr_minus2","length","gene_id","ORF_pept","Method","RNA_sites")],agg,by=c("gene_id","ORF_pept","Method","RNA_sites"))

agg3<-aggregate(x=agg2[,"length"],by=list(agg2[,"gene_id"],agg2[,"ORF_pept"],agg2[,"Method"],agg2[,"RNA_sites"]),FUN=max)

names(agg3)<-c("gene_id","ORF_pept","Method","RNA_sites","length")
agg4<-merge(x=all_sORFs_CCDS_multi_final[,c("ORF_id_tr_minus2","length","gene_id","ORF_pept","Method","RNA_sites")],agg3,by=c("gene_id","ORF_pept","Method","length","RNA_sites"))
all_sORFs_CCDS_multi_final<-all_sORFs_CCDS_multi_final[all_sORFs_CCDS_multi_final[,"ORF_id_tr_minus2"]%in%agg4[,"ORF_id_tr_minus2"],]


all_sORFs_CCDS_periodic<-all_sORFs_CCDS_multi_final[all_sORFs_CCDS_multi_final[,"ORF_pval_multi_ribo"]<0.05,]
all_sORFs_CCDS_periodic<-all_sORFs_CCDS_multi_final[!is.na(all_sORFs_CCDS_multi_final[,"transcript_id"]),]



all_sORFs_CCDS_periodic$n_exons_ORF<-sapply(strsplit(all_sORFs_CCDS_periodic$to_check_ALL,split=";"),FUN=function(x){sum(x!="NA")})


print(paste("--- Checking CCDS ORFs intersections with annotated CDS regions,",date(),sep=" "))



ex_to_check<-strsplit(all_sORFs_CCDS_periodic$to_check_ALL,split=";")

ex_to_check_spl<-unique(unlist(ex_to_check))

ex_to_check_spl<-strsplit(ex_to_check_spl,split="_")

bedfiles_to_check<-data.frame(chr=NA,start=NA,end=NA,type=NA,gene_id=NA,strand=NA)
for(h in 1:length(ex_to_check_spl)){
        to_bed<-ex_to_check_spl[[h]]
        bedfiles_to_check[h,"chr"]<-to_bed[1]
        bedfiles_to_check[h,"start"]<-as.numeric(to_bed[2])
        bedfiles_to_check[h,"end"]<-as.numeric(to_bed[3])
        bedfiles_to_check[h,"type"]<-to_bed[4]
        bedfiles_to_check[h,"gene_id"]<-to_bed[5]
        bedfiles_to_check[h,"strand"]<-to_bed[6]
        
}
bedfiles_to_check<-bedfiles_to_check[!is.na(bedfiles_to_check[,"chr"]),]
bedfiles_to_check<-bedfiles_to_check[bedfiles_to_check[,"chr"]!="NA",]

write.table(bedfiles_to_check,file="sORFs_totest",quote=F,row.names=F,sep="\t",col.names=F)

system("sort -k1,1 -k2,2n sORFs_totest > sORFs_totest.bed")

bedfiles_to_check<-read.table("sORFs_totest.bed",stringsAsFactors=F,header=F)
colnames(bedfiles_to_check)<-c("chr","start","end","type","gene_id","strand")
bedfiles_to_check<-bedfiles_to_check[!is.na(bedfiles_to_check[,"chr"]),]

bedfiles_to_check[is.na(bedfiles_to_check["strand"]),"strand"]<-"+"

write.table(bedfiles_to_check,file="sORFs_totest.bed",quote=F,row.names=F,sep="\t",col.names=F)

fhalf_scr<-paste(args[3],"intersectBed -v -a sORFs_totest.bed -b",sep = "/")

shalf_scr<-paste(args[1],"all_cds.bed > sORFs_totest_nocds.bed",sep = "/")

system(paste(fhalf_scr,shalf_scr,sep = " "))


command<-paste("wc -l","sORFs_totest_nocds.bed")
lines_in_file<-system(command,intern=T)
lines_in_file<-as.numeric(strsplit(lines_in_file,split=" ")[[1]][1])

if(lines_in_file>0){
        results_nonoverlapcdss<-read.table("sORFs_totest_nocds.bed",stringsAsFactors=F,header=F)
        names(results_nonoverlapcdss)<-names(bedfiles_to_check)
        results_nonoverlapcdss[,"exon_id"]<-paste(results_nonoverlapcdss[,"chr"],results_nonoverlapcdss[,"start"],results_nonoverlapcdss[,"end"],results_nonoverlapcdss[,"type"],results_nonoverlapcdss[,"gene_id"],results_nonoverlapcdss[,"strand"],sep="_")
        NA_str<-which(is.na(results_nonoverlapcdss[,"strand"]))
        if(length(NA_str)>0){
                for(o in NA_str){
                        results_nonoverlapcdss[o,"exon_id"]<-paste(results_nonoverlapcdss[o,"chr"],results_nonoverlapcdss[o,"start"],results_nonoverlapcdss[o,"end"],results_nonoverlapcdss[o,"type"],results_nonoverlapcdss[o,"gene_id"],sep="_")
                        
                }
        }
}

if(lines_in_file==0){
        results_nonoverlapcdss<-data.frame(exon_id=NA,stringsAsFactors=F)
}

overl_cds<-c()
for(i in 1:length(ex_to_check)){
        a<-ex_to_check[[i]]
        a<-a[a!="NA"]
        overl_cds[i]<-sum(!a%in%results_nonoverlapcdss$exon_id)>0
        
}

all_sORFs_CCDS_periodic_nocds<-all_sORFs_CCDS_periodic[!overl_cds,]

all_sORFs_CCDS_periodic_nocds<-all_sORFs_CCDS_periodic_nocds[all_sORFs_CCDS_periodic_nocds[,"ORF_pval_multi_ribo"]<0.05,]
all_sORFs_CCDS_periodic_nocds<-all_sORFs_CCDS_periodic_nocds[!is.na(all_sORFs_CCDS_periodic_nocds[,"transcript_id"]),]
write.table(all_sORFs_CCDS_periodic_nocds,file="orfs_before_u_dorfs",quote=F,row.names=F,sep="\t",col.names=T)
if(dim(all_sORFs_CCDS_periodic_nocds)[1]>0){
        all_sORFs_CCDS_periodic_nocds$type<-NA
        
        for(r in 1:dim(all_sORFs_CCDS_periodic_nocds)[1]){
                
                x<-all_sORFs_CCDS_periodic_nocds[r,]
                type<-NA
                if(!is.na(as.numeric(x[,"annotated_start"])) &  !is.na(as.numeric(x[,"annotated_stop"]))){
                        if(as.numeric(x[,"start_pos"])<as.numeric(x[,"annotated_start"])){type<-"uORF"}
                        if(as.numeric(x[,"start_pos"])>as.numeric(x[,"annotated_stop"])){type<-"dORF"}
                        if(as.numeric(x[,"start_pos"])>as.numeric(x[,"annotated_start"]) & x[,"start_pos"]<as.numeric(x[,"annotated_stop"]) & as.numeric(x[,"st2vect"])>as.numeric(x[,"annotated_stop"])){type<-"Overl_dORF"}
                        if(as.numeric(x[,"start_pos"])<as.numeric(x[,"annotated_start"]) & x[,"st2vect"]>as.numeric(x[,"annotated_start"])){type<-"Overl_uORF"}
                        
                }
                all_sORFs_CCDS_periodic_nocds[r,"type"]<-type
        }
}

if(dim(all_sORFs_CCDS_periodic_nocds)[1]==0){
        print("Warning! No u/dORFs found ! all ORFs overlap annotated CDS exons")
        all_sORFs_CCDS_periodic_nocds[1,]<-NA
        all_sORFs_CCDS_periodic_nocds$type<-NA
        
}

all_sORFs_CCDS_periodic_nocds_filtered_multi<-all_sORFs_CCDS_periodic_nocds[(all_sORFs_CCDS_periodic_nocds$pct_covered_onlymulti_ribo/all_sORFs_CCDS_periodic_nocds$pct_region_covered_ribo)<0.3,]
all_sORFs_CCDS_periodic_nocds_filtered_multi<-all_sORFs_CCDS_periodic_nocds_filtered_multi[all_sORFs_CCDS_periodic_nocds_filtered_multi$pct_region_covered_ribo>0.3,]
all_sORFs_CCDS_periodic_nocds_filtered_multi<-all_sORFs_CCDS_periodic_nocds_filtered_multi[!is.na(all_sORFs_CCDS_periodic_nocds_filtered_multi[,"transcript_id"]),]

all_sORFs_CCDS_periodic<-all_sORFs_CCDS_periodic[overl_cds,]
all_sORFs_CCDS_periodic_nofilt<-all_sORFs_CCDS_periodic
all_sORFs_CCDS_periodic<-all_sORFs_CCDS_periodic[(all_sORFs_CCDS_periodic$pct_covered_onlymulti_ribo/all_sORFs_CCDS_periodic$pct_region_covered_ribo)<0.3,]
all_sORFs_CCDS_periodic<-all_sORFs_CCDS_periodic[!is.na(all_sORFs_CCDS_periodic[,"transcript_id"]),]


setwd("../")


dir.create("ORFs_CCDS", showWarnings = FALSE)
dir.create("ORFs_CCDS/best_periodicity", showWarnings = FALSE)
dir.create("ORFs_CCDS/max_P_sites", showWarnings = FALSE)
dir.create("ORFs_CCDS/more_tapers", showWarnings = FALSE)


sORFs_sign_filtered_cds<-all_sORFs_CCDS_periodic_nocds[all_sORFs_CCDS_periodic_nocds[,"Method"]=="best_periodicity",]
write.table(sORFs_sign_filtered_cds,file="ORFs_CCDS/best_periodicity/sORFs_sign_filtered_cds",quote=F,row.names=F,sep="\t",col.names=T)
sORFs_sign_filtered_cds_multi<-all_sORFs_CCDS_periodic_nocds_filtered_multi[all_sORFs_CCDS_periodic_nocds_filtered_multi[,"Method"]=="best_periodicity",]
write.table(sORFs_sign_filtered_cds_multi,file="ORFs_CCDS/best_periodicity/sORFs_sign_filtered_cds_multi",quote=F,row.names=F,sep="\t",col.names=T)
ORFs_sign_filtered_multi<-all_sORFs_CCDS_periodic[all_sORFs_CCDS_periodic[,"Method"]=="best_periodicity",]
ORFs_sign_notfiltered_multi<-all_sORFs_CCDS_periodic_nofilt[all_sORFs_CCDS_periodic_nofilt[,"Method"]=="best_periodicity",]

write.table(ORFs_sign_filtered_multi,file="ORFs_CCDS/best_periodicity/ORFs_sign_filtered_multi",quote=F,row.names=F,sep="\t",col.names=T)
write.table(ORFs_sign_notfiltered_multi,file="ORFs_CCDS/best_periodicity/ORFs_sign_notfiltered_multi",quote=F,row.names=F,sep="\t",col.names=T)

ORFs_all<-CCDS_orfs[CCDS_orfs[,"Method"]=="best_periodicity",]
ORFs_all<-ORFs_all[!is.na(ORFs_all[,"transcript_id"]),]
write.table(ORFs_all,file="ORFs_CCDS/best_periodicity/ORFs_all",quote=F,row.names=F,sep="\t",col.names=T)


sORFs_sign_filtered_cds<-all_sORFs_CCDS_periodic_nocds[all_sORFs_CCDS_periodic_nocds[,"Method"]=="max_P_sites",]
write.table(sORFs_sign_filtered_cds,file="ORFs_CCDS/max_P_sites/sORFs_sign_filtered_cds",quote=F,row.names=F,sep="\t",col.names=T)
sORFs_sign_filtered_cds_multi<-all_sORFs_CCDS_periodic_nocds_filtered_multi[all_sORFs_CCDS_periodic_nocds_filtered_multi[,"Method"]=="max_P_sites",]
write.table(sORFs_sign_filtered_cds_multi,file="ORFs_CCDS/max_P_sites/sORFs_sign_filtered_cds_multi",quote=F,row.names=F,sep="\t",col.names=T)
ORFs_sign_filtered_multi<-all_sORFs_CCDS_periodic[all_sORFs_CCDS_periodic[,"Method"]=="max_P_sites",]

ORFs_sign_notfiltered_multi<-all_sORFs_CCDS_periodic_nofilt[all_sORFs_CCDS_periodic_nofilt[,"Method"]=="max_P_sites",]

write.table(ORFs_sign_notfiltered_multi,file="ORFs_CCDS/max_P_sites/ORFs_sign_notfiltered_multi",quote=F,row.names=F,sep="\t",col.names=T)


write.table(ORFs_sign_filtered_multi,file="ORFs_CCDS/max_P_sites/ORFs_sign_filtered_multi",quote=F,row.names=F,sep="\t",col.names=T)
ORFs_all<-CCDS_orfs[CCDS_orfs[,"Method"]=="max_P_sites",]
ORFs_all<-ORFs_all[!is.na(ORFs_all[,"transcript_id"]),]

write.table(ORFs_all,file="ORFs_CCDS/max_P_sites/ORFs_all",quote=F,row.names=F,sep="\t",col.names=T)


sORFs_sign_filtered_cds<-all_sORFs_CCDS_periodic_nocds[all_sORFs_CCDS_periodic_nocds[,"Method"]=="more_tapers",]
write.table(sORFs_sign_filtered_cds,file="ORFs_CCDS/more_tapers/sORFs_sign_filtered_cds",quote=F,row.names=F,sep="\t",col.names=T)
sORFs_sign_filtered_cds_multi<-all_sORFs_CCDS_periodic_nocds_filtered_multi[all_sORFs_CCDS_periodic_nocds_filtered_multi[,"Method"]=="more_tapers",]
write.table(sORFs_sign_filtered_cds_multi,file="ORFs_CCDS/more_tapers/sORFs_sign_filtered_cds_multi",quote=F,row.names=F,sep="\t",col.names=T)
ORFs_sign_filtered_multi<-all_sORFs_CCDS_periodic[all_sORFs_CCDS_periodic[,"Method"]=="more_tapers",]

ORFs_sign_notfiltered_multi<-all_sORFs_CCDS_periodic_nofilt[all_sORFs_CCDS_periodic_nofilt[,"Method"]=="more_tapers",]

write.table(ORFs_sign_notfiltered_multi,file="ORFs_CCDS/more_tapers/ORFs_sign_notfiltered_multi",quote=F,row.names=F,sep="\t",col.names=T)


write.table(ORFs_sign_filtered_multi,file="ORFs_CCDS/more_tapers/ORFs_sign_filtered_multi",quote=F,row.names=F,sep="\t",col.names=T)
ORFs_all<-CCDS_orfs[CCDS_orfs[,"Method"]=="more_tapers",]
ORFs_all<-ORFs_all[!is.na(ORFs_all[,"transcript_id"]),]

write.table(ORFs_all,file="ORFs_CCDS/more_tapers/ORFs_all",quote=F,row.names=F,sep="\t",col.names=T)

print(paste("--- CCDS ORF finding Done!","---",date(),sep=" "))


