library("XNomial")
library("foreach")
library("doMC")

library("multitaper")
library("seqinr")

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


###This functions reads big text files efficiently: from http://www.r-bloggers.com/faster-files-in-r/
readBigText<-function(x){
        f=file(x,"rb")
        a=readChar(f,file.info(x)$size,useBytes=T);a<-strsplit(a,"\n",fixed=T,useBytes=T)[[1]]
        close(f)
        return(a)
}



### This function gets the FFT frequence and power values, from http://stackoverflow.com/questions/3485456/useful-little-functions-in-r


getFFTFreqs<-function(Nyq.Freq, data)
{
        if ((length(data) %% 2) == 1) # Odd number of samples
        {
                FFTFreqs <- c(seq(0, Nyq.Freq, length.out=(length(data)+1)/2), 
                              seq(-Nyq.Freq, 0, length.out=(length(data)-1)/2))
        }
        else # Even number
        {
                FFTFreqs <- c(seq(0, Nyq.Freq, length.out=length(data)/2), 
                              seq(-Nyq.Freq, 0, length.out=length(data)/2))
        }
        
        return (FFTFreqs)
}

### This function outputs the max FFT frequence and power values

take_maxfreq_and_power_FFT_Spec<-function(x){
        
        if(length(x)<10){x<-c(rep(0,3),x,rep(0,3))}
        gino<-getFFTFreqs(Nyq.Freq=0.5,data=x)
        modFFT <- Mod(fft(x))
        FFTdata <- cbind(gino, modFFT)
        
        freq3_fft<-abs(FFTdata[which(abs((abs(FFTdata[,1])-(1/3)))==min(abs((abs(FFTdata[,1])-(1/3))))),1])
        
        power3_fft<-FFTdata[which(abs((abs(FFTdata[,1])-(1/3)))==min(abs((abs(FFTdata[,1])-(1/3))))),2]
        
       
        
        spect_x<-spectrum(x,plot=FALSE)
        
        freq3_sp<-abs(spect_x$freq[which(abs(spect_x$freq-(1/3))==min(abs((spect_x$freq)-(1/3))))])
        power3_sp<-abs(spect_x$spec[which(abs(spect_x$freq-(1/3))==min(abs((spect_x$freq)-(1/3))))])
                
        return(c(freq3_fft,power3_fft,freq3_sp,power3_sp))
}


### This function plots the raw FFT periodogram, from http://stackoverflow.com/questions/3485456/useful-little-functions-in-r


plotFFT<-function(x, y, samplingFreq, shadeNyq=TRUE, showPeriod = TRUE)
{
        Nyq.Freq <- samplingFreq/2
        FFTFreqs <- getFFTFreqs(Nyq.Freq, y)
        
        FFT <- fft(y)
        modFFT <- Mod(FFT)
        FFTdata <- cbind(FFTFreqs, modFFT)
        plot(FFTdata[1:nrow(FFTdata)/2,], t="l", pch=20, lwd=2, cex=0.8, main="",
             xlab="Frequency (Hz)", ylab="Power")
        if (showPeriod == TRUE)
        {
                # Period axis on top        
                a <- axis(3, lty=0, labels=FALSE)
                axis(3, cex.axis=0.6, labels=format(1/a, digits=2), at=a)
        }
        if (shadeNyq == TRUE)
        {
                # Gray out lower frequencies
                rect(0, 0, 2/max(x), max(FFTdata[,2])*2, col="gray", density=30)
        }
        
        ret <- list("freq"=FFTFreqs, "FFT"=FFT, "modFFT"=modFFT)
}

### This function calculates the CSCPD as in Michel et al 2012 Gen Res


dual_take_CSCPDs<-function(tracks_to_analyze=all_tracks,index_tracks=all_tracks_index,exon_ids=all_tracks_index){
        unique_index<-unique(as.data.frame(exon_ids)[,"exon_id"])
        interpolation_mat1 = matrix(NA, nrow = length(unique_index), ncol=100)
        interpolation_mat2 = matrix(NA, nrow = length(unique_index), ncol=100)
        interpolation_mat3 = matrix(NA, nrow = length(unique_index), ncol=100)
        rownames(interpolation_mat1)<-unique_index
        rownames(interpolation_mat2)<-unique_index
        rownames(interpolation_mat3)<-unique_index
        for(i in 1:length(unique_index)){
                id<-unique_index[i]
                exon_track<-tracks_to_analyze[index_tracks[,1]==id]
                withsep<-strsplit(exon_track,split=" ")
                x<-t(data.frame(withsep))
                #rnames[i]<-x[1,1]
                strand<-x[1,2]
                if(length(grep("CCDS",id))>0){tracks_pre<-t(x[,-c(1:3)])} else {
                        tracks_pre<-t(x[,-c(1:2)])}
                if(strand=="-"){
                        tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
                } else if (strand=="+"){
                        tracks<-tracks_pre}
                colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
                mode(tracks)<-"numeric"
                length<-dim(tracks)[1]
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                FRAME_MAX_phase<-max.col(t(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))-1
                nts_toadd<-(3-FRAME_MAX_phase)%%3
                
                read_mat<-t(as.matrix(c(rep(x=0,nts_toadd),tracks[,1])))
                length<-(length(read_mat))
                
                upprop = matrix(NA,nrow = 1,ncol = length)
                downprop = matrix(NA,nrow = 1,ncol = length)
                seq_length1 = (length/3)-1
                seq_length2 = seq_length1 - 1
                
                #calculate the cumulative upstream proportions
                denom = cumsum(read_mat[1,])
                upprop[1,(1+3*(0:seq_length1))] = cumsum(read_mat[1,(1+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                upprop[1,(2+3*(0:seq_length1))] = cumsum(read_mat[1,(2+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                upprop[1,(3+3*(0:seq_length1))] = cumsum(read_mat[1,(3+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                #calculate the cumulative downstream proportions
                totalsDown = rep(NA, seq_length1)
                for (j in 1:seq_length1){
                        totalsDown[j] = sum(read_mat[1+3*(j:seq_length1)]+ read_mat[2+3*(j:seq_length1)]+
                                                    read_mat[3+3*(j:seq_length1)])}
                downprop[1,(1+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(1+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                downprop[1,(2+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(2+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                downprop[1,(3+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(3+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                #Calculate the CSCPD (absolute difference between cumulative upstream and downstream
                #proportions for sub-codon positions 1,2,3)
                third_full_seq_minus1 = (length/3)-1
                y1 = abs(upprop-downprop)[1+3*(0:third_full_seq_minus1)]
                y1_withoutNA = y1[which(y1 != "NA")] 
                y1_withoutNaN = y1_withoutNA[which(y1_withoutNA != "NaN")]
                
                y2 = abs(upprop-downprop)[2+3*(0:third_full_seq_minus1)]
                y2_withoutNA = y2[which(y2 != "NA")]
                y2_withoutNaN = y2_withoutNA[which(y2_withoutNA != "NaN")]
                
                y3 = abs(upprop-downprop)[3+3*(0:third_full_seq_minus1)]
                y3_withoutNA = y3[which(y3 != "NA")]
                y3_withoutNaN = y3_withoutNA[which(y3_withoutNA != "NaN")]
                
                
                
                if(length(y1_withoutNaN)==0 ){
                        y1_withoutNaN<-c(0.1,0.1,0.1,0.1,0.1)
                }
                if(length(y2_withoutNaN)==0 ){
                        y2_withoutNaN<-c(0.1,0.1,0.1,0.1,0.1)
                }
                if(length(y3_withoutNaN)==0 ){
                        y3_withoutNaN<-c(0.1,0.1,0.1,0.1,0.1)
                }
                
                
                
                
                if(length(y1_withoutNaN)<4 & length(y1_withoutNaN)>2 ){
                        y1_withoutNaN<-c(y1_withoutNaN[1],y1_withoutNaN,y1_withoutNaN[length(y1_withoutNaN)])
                }
                if(length(y2_withoutNaN)<4 & length(y2_withoutNaN)>2 ){
                        y2_withoutNaN<-c(y2_withoutNaN[1],y2_withoutNaN,y2_withoutNaN[length(y2_withoutNaN)])
                }
                if(length(y3_withoutNaN)<4 & length(y3_withoutNaN)>2 ){
                        y3_withoutNaN<-c(y3_withoutNaN[1],y3_withoutNaN,y3_withoutNaN[length(y3_withoutNaN)])
                }
                
                
                
                if(length(y1_withoutNaN)<3 & length(y1_withoutNaN)>0 ){
                        y1_withoutNaN<-c(y1_withoutNaN[1],y1_withoutNaN[1],y1_withoutNaN,y1_withoutNaN[length(y1_withoutNaN)],y1_withoutNaN[length(y1_withoutNaN)])
                }
                if(length(y2_withoutNaN)<3 & length(y2_withoutNaN)>0 ){
                        y2_withoutNaN<-c(y2_withoutNaN[1],y2_withoutNaN[1],y2_withoutNaN,y2_withoutNaN[length(y2_withoutNaN)],y2_withoutNaN[length(y2_withoutNaN)])
                }
                if(length(y3_withoutNaN)<3 & length(y3_withoutNaN)>0 ){
                        y3_withoutNaN<-c(y3_withoutNaN[1],y3_withoutNaN[1],y3_withoutNaN,y3_withoutNaN[length(y3_withoutNaN)],y3_withoutNaN[length(y3_withoutNaN)])
                }
                
                
                
                length_third_seq_withoutNaN = length(y1_withoutNaN)
                length_third_seq_withoutNaN_minus1= length(y1_withoutNaN)-1
                
                
                
                x1 = 1+3*(0:length_third_seq_withoutNaN_minus1)
                x2 = 2+3*(0:length_third_seq_withoutNaN_minus1)
                x3 = 3+3*(0:length_third_seq_withoutNaN_minus1)
                #Converting all coordinates in coding region to relative values between 0 and 1 and using a smoothing function
                ys1 = smooth.spline(x1/(length_third_seq_withoutNaN*3),y1_withoutNaN)
                ys2 = smooth.spline(x2/(length_third_seq_withoutNaN*3),y2_withoutNaN)
                ys3 = smooth.spline(x3/(length_third_seq_withoutNaN*3),y3_withoutNaN)
                #Sampling 100 equidistant CSCPD values between 0 and 1
                xout = 0.01*(1:100)
                yout1 = predict(ys1,xout)$y
                yout2 = predict(ys2,xout)$y
                yout3 = predict(ys3,xout)$y
                interpolation_mat1[i,] = yout1
                interpolation_mat2[i,] = yout2
                interpolation_mat3[i,] = yout3
        }
        interpolation_list<-list(interpolation_mat1,interpolation_mat2,interpolation_mat3)
        return(interpolation_list)
}


### This function calculates the PTS as in Michel et al 2012 Gen Res

dual_calculate_PTSs<-function(all_tracks,index,prev_percentiles){
        data_frame<-as.data.frame(prev_percentiles)
        attach(data_frame)
        unique_index<-unique(index)
        difference_mat1 = matrix(NA, nrow = dim(unique_index)[1], ncol=100)
        difference_mat2 = matrix(NA, nrow = dim(unique_index)[1], ncol=100)
        difference_mat3 = matrix(NA, nrow = dim(unique_index)[1], ncol=100)
        
        
        rnames = rep("",dim(unique_index)[1])
        cnames = c("PTS1", "PTS2", "PTS3", "PTS")
        PTS = matrix(0, nrow = dim(unique_index)[1], ncol=4, dimnames=list(rnames,cnames))
        
        
        for(i in 1:dim(unique_index)[1]){
                id<-unique_index[i,]
                exon_track<-all_tracks[index==id]
                withsep<-strsplit(exon_track,split=" ")
                x<-t(data.frame(withsep))
                rnames[i]<-x[1,1]
                strand<-x[1,2]
                if(length(grep("CCDS",id))>0){tracks_pre<-t(x[,-c(1:3)])} else {
                        tracks_pre<-t(x[,-c(1:2)])}
                if(strand=="-"){
                        tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
                } else if (strand=="+"){
                        tracks<-tracks_pre}
                colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
                mode(tracks)<-"numeric"
                length<-dim(tracks)[1]
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                FRAME_MAX_phase<-max.col(t(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))-1
                nts_toadd<-(3-FRAME_MAX_phase)%%3
                
                read_mat<-t(as.matrix(c(rep(x=0,nts_toadd),tracks[,1])))
                length<-(length(read_mat))
                
                upprop = matrix(NA,nrow = 1,ncol = length)
                downprop = matrix(NA,nrow = 1,ncol = length)
                seq_length1 = (length/3)-1
                seq_length2 = seq_length1 - 1
                
                #calculate the cumulative upstream proportions
                denom = cumsum(read_mat[1,])
                upprop[1,(1+3*(0:seq_length1))] = cumsum(read_mat[1,(1+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                upprop[1,(2+3*(0:seq_length1))] = cumsum(read_mat[1,(2+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                upprop[1,(3+3*(0:seq_length1))] = cumsum(read_mat[1,(3+3*(0:seq_length1))])/
                        denom[(3+3*(0:seq_length1))]
                #calculate the cumulative downstream proportions
                totalsDown = rep(NA, seq_length1)
                for (j in 1:seq_length1){
                        totalsDown[j] = sum(read_mat[1+3*(j:seq_length1)]+ read_mat[2+3*(j:seq_length1)]+
                                                    read_mat[3+3*(j:seq_length1)])}
                downprop[1,(1+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(1+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                downprop[1,(2+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(2+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                downprop[1,(3+3*(0:seq_length2))] = rev(cumsum(rev(read_mat[1,(3+3*(1:seq_length1))])))/
                        totalsDown[0:seq_length1]
                #Calculate the CSCPD (absolute difference between cumulative upstream and downstream
                #proportions for sub-codon positions 1,2,3)
                third_full_seq_minus1 = (length/3)-1
                y1 = abs(upprop-downprop)[1+3*(0:third_full_seq_minus1)]
                y1_withoutNA = y1[which(y1 != "NA")] 
                y1_withoutNaN = y1_withoutNA[which(y1_withoutNA != "NaN")]
                
                y2 = abs(upprop-downprop)[2+3*(0:third_full_seq_minus1)]
                y2_withoutNA = y2[which(y2 != "NA")]
                y2_withoutNaN = y2_withoutNA[which(y2_withoutNA != "NaN")]
                
                y3 = abs(upprop-downprop)[3+3*(0:third_full_seq_minus1)]
                y3_withoutNA = y3[which(y3 != "NA")]
                y3_withoutNaN = y3_withoutNA[which(y3_withoutNA != "NaN")]
                
                if(length(y1_withoutNaN)<4){
                        y1_withoutNaN<-c(y1_withoutNaN[1],y1_withoutNaN[1],y1_withoutNaN,y1_withoutNaN[length(y1_withoutNaN)],y1_withoutNaN[length(y1_withoutNaN)])
                }
                if(length(y2_withoutNaN)<4){
                        y2_withoutNaN<-c(y2_withoutNaN[1],y2_withoutNaN[1],y2_withoutNaN,y2_withoutNaN[length(y2_withoutNaN)],y2_withoutNaN[length(y2_withoutNaN)])
                }
                if(length(y3_withoutNaN)<4){
                        y3_withoutNaN<-c(y3_withoutNaN[1],y3_withoutNaN[1],y3_withoutNaN,y3_withoutNaN[length(y3_withoutNaN)],y3_withoutNaN[length(y3_withoutNaN)])
                }
                
                
                length_third_seq_withoutNaN = length(y1_withoutNaN)
                length_third_seq_withoutNaN_minus1= length(y1_withoutNaN)-1
                
                
                
                x1 = 1+3*(0:length_third_seq_withoutNaN_minus1)
                x2 = 2+3*(0:length_third_seq_withoutNaN_minus1)
                x3 = 3+3*(0:length_third_seq_withoutNaN_minus1)
                #Converting all coordinates in coding region to relative values between 0 and 1 and using a smoothing function
                ys1 = smooth.spline(x1/(length_third_seq_withoutNaN*3),y1_withoutNaN)
                ys2 = smooth.spline(x2/(length_third_seq_withoutNaN*3),y2_withoutNaN)
                ys3 = smooth.spline(x3/(length_third_seq_withoutNaN*3),y3_withoutNaN)
                #Sampling 100 equidistant CSCPD values between 0 and 1
                
                xout = 0.01*(1:100)
                yout1 = predict(ys1,xout)$y
                yout2 = predict(ys2,xout)$y
                yout3 = predict(ys3,xout)$y
                
                
                difference_mat1[i,] = yout1 - Percentile_P1
                difference_mat2[i,] = yout2 - Percentile_P2
                difference_mat3[i,] = yout3 - Percentile_P3
                gene_counter<-i
                
                for (k in 1:100)
                {
                        if (as.numeric(difference_mat1[gene_counter,k]) >= 0){
                                PTS[gene_counter,1] = PTS[gene_counter,1] + difference_mat1[gene_counter,k]
                        }
                        if (difference_mat2[gene_counter,k] >= 0){
                                PTS[gene_counter,2] = PTS[gene_counter,2] + difference_mat2[gene_counter,k]
                        }
                        if (difference_mat3[gene_counter,k] >= 0){
                                PTS[gene_counter,3] = PTS[gene_counter,3] + difference_mat3[gene_counter,k]                
                        }
                }
                
                PTS[gene_counter,4] = PTS[gene_counter,1] + PTS[gene_counter,2] + PTS[gene_counter,3]
                
        }       
        detach(data_frame)
        PTS<-as.data.frame(PTS,row.names=F,stringsAsFactors=F)
        PTS$exon_id<-unique_index[,1]
        PTS<-PTS[,c("exon_id","PTS1","PTS2","PTS3","PTS")]
        return(PTS)
}


### This function plots the PTS as in Michel et al 2012 Gen Res (outdated)



plot_CSCPDs<-function(x,y){
        one<-x[[1]][y,]
        two<-x[[2]][y,]
        three<-x[[3]][y,]
        exon_id<-rownames(x[[1]])[y]
        plot(one,type="l",col="red",ylim=c(0,1),ylab="CSCPDs")
        lines(two,type="l",col="green")
        lines(three,type="l",col="blue")
        legend("top",exon_id)
}

### This function plots data_tracks information (P-sites distribution, FFT etc...) (outdated)


plot_tracks_fig<-function(exon_id,complete_tracks=all_tracks,index=all_tracks_index){
        
        exon_track<-complete_tracks[index==exon_id]
        withsep<-strsplit(exon_track,split=" ")
        x<-t(data.frame(withsep))
        strand<-x[1,2]
        if(length(grep("CCDS",exon_id,))>0){tracks_pre<-t(x[,-c(1:3)])} else {
                tracks_pre<-t(x[,-c(1:2)])}
        if(strand=="-"){
                tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
        } else if (strand=="+"){
                tracks<-tracks_pre}
        colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
        mode(tracks)<-"numeric"
        length<-dim(tracks)[1]
        
        
        Phase0<-round(x=sum(tracks[seq(1,dim(tracks)[1],by=3),1])/sum(tracks[,1]),digits=4)
        Phase1<-round(x=sum(tracks[seq(2,dim(tracks)[1],by=3),1])/sum(tracks[,1]),digits=4)
        Phase2<-round(x=sum(tracks[seq(3,dim(tracks)[1],by=3),1])/sum(tracks[,1]),digits=4)
        
        Phase0_RNA<-round(x=sum(tracks[seq(1,dim(tracks)[1],by=3),4])/sum(tracks[,4]),digits=4)
        Phase1_RNA<-round(x=sum(tracks[seq(2,dim(tracks)[1],by=3),4])/sum(tracks[,4]),digits=4)
        Phase2_RNA<-round(x=sum(tracks[seq(3,dim(tracks)[1],by=3),4])/sum(tracks[,4]),digits=4)
        
        valuesribo<-c(min(tracks[,2]),tracks[,2],min(tracks[,2]))
        valuesrna<-c(min(tracks[,3]),tracks[,3],min(tracks[,3]))
        nucleot<-c(min(seq(1,dim(tracks)[1])),seq(1,dim(tracks)[1]),max(seq(1,dim(tracks)[1])))
        name_region<-exon_id
        x11(width=16,height=10)
        par(mar=c(4, 4, 1, 1))
        split.screen( figs = c( 2, 2 ) )
        split.screen( figs = c( 2, 1 ) ,screen=1)
        
        screen(5)   
        
        
        plot(tracks[,1],type="h",col=c("red","dark green","blue"),ylab="P_sites",xlab="nt")
        
        split.screen( figs = c( 1, 2 ) ,screen=6)
        screen(7)
        barplot(c(Phase0,Phase1,Phase2),xlab="Phases",ylim=c(0,1),ylab="%_Alignments",main="%Frames_RIBO",col=c("red","dark green","blue"))
        screen(8)
        plotFFT(x=seq(1,dim(tracks)[1]),y=tracks[,1],samplingFreq=1)
        
        screen(2)
        plot(tracks[,2],type="l",col="red",ylab="Ribo_cov",xlab="nt",main=x[1,1])
        polygon(x=nucleot,y=valuesribo,col="red")
        split.screen( figs = c( 2, 1 ) ,screen=3)
        screen(9)
        plot(tracks[,4],type="h",col=c("red","dark green","blue"),ylab="RNA_center",xlab="nt")
        split.screen( figs = c( 1, 2 ) ,screen=10)
        screen(11)
        barplot(c(Phase0_RNA,Phase1_RNA,Phase2_RNA),ylim=c(0,1),xlab="Phases",ylab="%_Alignments",main="%Frames_RNA",col=c("red","dark green","blue"))
        screen(12)
        plotFFT(x=seq(1,dim(tracks)[1]),y=tracks[,4],samplingFreq=1)
        screen(4)
        plot(tracks[,3],type="l",col="dark grey",ylab="RNA_cov",xlab="nt")
        polygon(x=nucleot,y=valuesrna,col="dark grey")  
        close.screen(all.screens=T)
}


### This function calculates the PTS from the CSPD, as in Michel et al 2012 Gen Res


calculate_PTS_from_CSPDs<-function(list_cscpds,quantile_value=0.95){
        
        cnames = c("Percentile_P1", "Percentile_P2", "Percentile_P3")
        rnames = rep("",100)
        quantiles = matrix(0, nrow = 100, ncol=3,dimnames=list(rnames,cnames))
        
        CSCPDs_1<-list_cscpds[[1]]
        CSCPDs_2<-list_cscpds[[2]]
        CSCPDs_3<-list_cscpds[[3]]
        
        difference_mat1 = matrix(NA, nrow = dim(CSCPDs_1)[1], ncol=100)
        difference_mat2 = matrix(NA, nrow = dim(CSCPDs_1)[1], ncol=100)
        difference_mat3 = matrix(NA, nrow = dim(CSCPDs_1)[1], ncol=100)
        
        for(j in 1:100){
                quantiles[j,1] = quantile(CSCPDs_1[,j],quantile_value,na.rm = T)
                quantiles[j,2] = quantile(CSCPDs_2[,j],quantile_value,na.rm = T)
                quantiles[j,3] = quantile(CSCPDs_3[,j],quantile_value,na.rm = T)
        }
        
        percentiles<-as.data.frame(quantiles,row.names=F)
        
        PTS = matrix(0, nrow = dim(CSCPDs_1)[1], ncol=4)
        
        

        difference_mat1 = t(apply(X=CSCPDs_1,MARGIN=1,FUN=function(x){x<-x-percentiles[,"Percentile_P1"]}))
        difference_mat2 = t(apply(X=CSCPDs_2,MARGIN=1,FUN=function(x){x<-x-percentiles[,"Percentile_P2"]}))
        difference_mat3 = t(apply(X=CSCPDs_3,MARGIN=1,FUN=function(x){x<-x-percentiles[,"Percentile_P3"]}))
        
        PTS[,1]<-t(apply(X=difference_mat1,MARGIN=1,FUN=function(x){sum(x[x>=0])}))
        PTS[,2]<-t(apply(X=difference_mat2,MARGIN=1,FUN=function(x){sum(x[x>=0])}))
        PTS[,3]<-t(apply(X=difference_mat3,MARGIN=1,FUN=function(x){sum(x[x>=0])}))
        PTS[,4]<-t(apply(X=PTS[,1:3],MARGIN=1,FUN=sum))
        PTS<-as.data.frame(PTS,row.names=NULL,stringsAsFactors=F)
        colnames(PTS)<-c("PTS1", "PTS2", "PTS3", "PTS")
        PTS$exon_id<-rownames(CSCPDs_1)
        PTS<-PTS[,c("exon_id","PTS1","PTS2","PTS3","PTS")]
        return(PTS)
}


### This function takes frequencies F-values and spectral coefficient for a data-track object.
### (you have to calculate slepian functions beforehand)

take_freqs_Fvalues_all_around_3nt_spec<-function(x,n_tapers,time_bw,slepians_values){
        if(length(x)<25){
                remain<-50-length(x)
                x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
        }
        if(length(x)<1024/2){padding<-1024}
        if(length(x)>=1024/2){padding<-"default"}
        resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
        resSpec2<-dropFreqs(resSpec1,0.1,0.45)
        freq_max<-resSpec2$freq[which(resSpec2$mtm$Ftest==max(resSpec2$mtm$Ftest))]
        Fmax<-resSpec2$mtm$Ftest[which(resSpec2$mtm$Ftest==max(resSpec2$mtm$Ftest))]
        P_all<-(pf(Fmax,df1=2,df2=(2*n_tapers)-2,lower.tail=F))
        
        resSpec2<-dropFreqs(resSpec2,0.29,0.39)
        
        freq_max_around_3nt<-resSpec2$freq[which(resSpec2$mtm$Ftest==max(resSpec2$mtm$Ftest))]
        
        Fmax_around_3nt<-resSpec2$mtm$Ftest[which(resSpec2$mtm$Ftest==max(resSpec2$mtm$Ftest))]
        P_around_3nt<-(pf(q=Fmax_around_3nt,df1=2,df2=(2*n_tapers)-2,lower.tail=F))
        
        
        freq_max_3nt<-resSpec1$freq[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
        
        Fmax_3nt<-resSpec1$mtm$Ftest[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
        P_3nt<-(pf(q=Fmax_3nt,df1=2,df2=(2*n_tapers)-2,lower.tail=F))
        Spec_3nt<-resSpec1$spec[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
        
        return(c(freq_max,P_all,freq_max_around_3nt,P_around_3nt,freq_max_3nt,P_3nt,Spec_3nt))
        
}


### This function calculates periodicity and other statistics on single exon tracks


make_analysis_exons<-function(x){
        strand<-x[1,2]
        tracks_pre<-t(x[,-c(1:2)])
        
        
        if(strand=="-"){
                tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]),rev(tracks_pre[,5]))
        } else if (strand=="+"){
                tracks<-tracks_pre}
        colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent","Seq")
        tracks<-tracks[,1:4]
        if(is.null(dim(tracks))){
                tracks<-t(as.matrix(tracks))
        }
        
        mode(tracks)<-"numeric"
        
        exon<-data.frame(exon_id=x[1,1],stringsAsFactors=F,row.names=NULL)
        exon$strand<-strand
        exon$frame_start_pred<-NA
        exon$frame_end_pred<-NA
        
        exon$length<-dim(tracks)[1]
        length<-dim(tracks)[1]
        
        
        P_sites_sum<-round(sum(tracks[,1]),digits=6)
        exon$P_sites_sum<-P_sites_sum
        
        
        Centered_sites_sum<-round(sum(tracks[,4]),digits=6)
        exon$RNA_sites_sum<-Centered_sites_sum
        exon$Ribocov_aver<-round(mean(tracks[,2]),digits=6)
        exon$RNAseqcov_aver<-round(mean(tracks[,3]),digits=6)
        exon$pctPhase_frame<-NA
        exon$pctPhase_frame_1<-NA
        exon$pctPhase_frame_2<-NA
        exon$pctPhaseCentered_frame<-NA
        exon$pctPhaseCentered_frame_1<-NA
        exon$pctPhaseCentered_frame_2<-NA
        if(length>2){
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                
                
                exon$pctPhase_frame<-Phase_P_sites_frame/P_sites_sum
                exon$pctPhase_frame_1<-Phase_P_sites_frame_1/P_sites_sum
                exon$pctPhase_frame_2<-Phase_P_sites_frame_2/P_sites_sum
                
                
                Phase_Centered_sites_frame<-sum(tracks[seq(1,length,by=3),4])
                Phase_Centered_sites_frame_1<-sum(tracks[seq(2,length,by=3),4])
                Phase_Centered_sites_frame_2<-sum(tracks[seq(3,length,by=3),4])
                
                
                exon$pctPhaseCentered_frame<-Phase_Centered_sites_frame/Centered_sites_sum
                exon$pctPhaseCentered_frame_1<-Phase_Centered_sites_frame_1/Centered_sites_sum
                exon$pctPhaseCentered_frame_2<-Phase_Centered_sites_frame_2/Centered_sites_sum
                
                
                MAXPhase_frame<-max(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2))
                FRAME_MAX_phase<-max.col(t(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2)))-1
                
                MAXPhaseCentered_frame<-max(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2))
                FRAME_MAX_phaseCentered<-max.col(t(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2)))-1
        }
        
        
        exon$multit_freq_best_ribo<-NA
        exon$pval_multit_3nt_ribo<-NA
        exon$spec_multit_3nt_ribo<-NA
        exon$fft_max_freq_ribo<-NA
        exon$fft_power_3_ribo<-NA
        exon$fft_aver_ribo<-NA
        exon$spec_max_freq_ribo<-NA
        exon$spec_power_3_ribo<-NA
        exon$spec_aver_power_ribo<-NA
        
        exon$multit_freq_best_rna<-NA
        exon$pval_multit_3nt_rna<-NA
        exon$spec_multit_3nt_rna<-NA
        exon$fft_max_freq_rna<-NA
        exon$fft_power_3_rna<-NA
        exon$fft_aver_rna<-NA
        exon$spec_max_freq_rna<-NA
        exon$spec_power_3_rna<-NA
        exon$spec_aver_power_rna<-NA
        
        exon$ORF_score_ribo<-NA
        exon$ORF_score_rna<-NA
        
        if(P_sites_sum>2 & length>5){
                if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                bestfreq_3ntpval_ribo<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,1],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                exon$multit_freq_best_ribo<-bestfreq_3ntpval_ribo[1]
                exon$pval_multit_3nt_ribo<-bestfreq_3ntpval_ribo[2]
                exon$spec_multit_3nt_ribo<-bestfreq_3ntpval_ribo[3]
                score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                exon$ORF_score_ribo<-log2(score1+score2+score3+1)
                
                if(max(tracks[,1])>(P_sites_sum*.7)){
                        new_track<-tracks
                        new_track[which(new_track[,1]==max(new_track[,1]))]<-0
                        exon$ORF_score_ribo<-NA
                        if(sum(new_track[,1])>2){
                                Phase_P_sites_frame_corr<-sum(new_track[seq(1,length,by=3),1])
                                Phase_P_sites_frame_1_corr<-sum(new_track[seq(2,length,by=3),1])
                                Phase_P_sites_frame_2_corr<-sum(new_track[seq(3,length,by=3),1])
                                score1<-((Phase_P_sites_frame_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                score2<-((Phase_P_sites_frame_1_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                score3<-((Phase_P_sites_frame_2_corr-sum(new_track[,1])/3)^2)/(sum(new_track[,1])/3)
                                exon$ORF_score_ribo<-log2(score1+score2+score3+1)
                        }
                }
                
                gino<-getFFTFreqs(Nyq.Freq=0.5,data=tracks[,1])
                modFFT <- Mod(fft(tracks[,1]))
                FFTdata <- cbind(gino, modFFT)
                exon$fft_aver_ribo<-mean(FFTdata[,2])
                exon$fft_power_3_ribo<-FFTdata[which(abs((gino-(1/3)))==min(abs((gino-(1/3))))),2]
                exon$fft_max_freq_ribo<-abs(gino[which(FFTdata==max((FFTdata[10:dim(FFTdata)[1]/2,2])),arr.ind=TRUE)[1]])[1]
                
                
                spect_P_sites<-spectrum(tracks[,1],plot=FALSE)
                exon$spec_max_freq_ribo<-spect_P_sites$freq[which(spect_P_sites$spec==max(spect_P_sites$spec),arr.ind=TRUE)][1]
                exon$spec_power_3_ribo<-spect_P_sites$spec[which(abs((spect_P_sites$freq-(1/3)))==min(abs((spect_P_sites$freq-(1/3)))))]
                exon$spec_aver_power_ribo<-mean(spect_P_sites$spec)        
                if(Centered_sites_sum>2){
                        
                        gino<-getFFTFreqs(Nyq.Freq=0.5,data=tracks[,4])
                        modFFT <- Mod(fft(tracks[,4]))
                        FFTdata <- cbind(gino, modFFT)
                        exon$fft_aver_rna<-mean(FFTdata[,2])
                        exon$fft_power_3_rna<-FFTdata[which(abs((gino-(1/3)))==min(abs((gino-(1/3))))),2]
                        exon$fft_max_freq_rna<-1/abs(gino[which(FFTdata==max((FFTdata[10:dim(FFTdata)[1]/2,2])),arr.ind=TRUE)[1]])[1]
                        
                        
                        
                        spect_rna<-spectrum(tracks[,4],plot=FALSE)
                        exon$spec_max_freq_rna<-spect_rna$freq[which(spect_rna$spec==max(spect_rna$spec),arr.ind=TRUE)][1]
                        exon$spec_power_3_rna<-spect_rna$spec[which(abs((spect_rna$freq-(1/3)))==min(abs((spect_rna$freq-(1/3)))))]
                        exon$spec_aver_power_rna<-mean(spect_rna$spec)
                        bestfreq_3ntpval_rna<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,4],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                        
                        exon$multit_freq_best_rna<-bestfreq_3ntpval_rna[1]
                        exon$pval_multit_3nt_rna<-bestfreq_3ntpval_rna[2]
                        exon$spec_multit_3nt_rna<-bestfreq_3ntpval_rna[3]
                        score_rna_1<-((Phase_Centered_sites_frame-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        score_rna_2<-((Phase_Centered_sites_frame_1-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        score_rna_3<-((Phase_Centered_sites_frame_2-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        exon$ORF_score_rna<-log2(score_rna_1+score_rna_2+score_rna_3+1)
                        if(max(tracks[,4])>(Centered_sites_sum*.7)){
                                new_track<-tracks
                                new_track[which(new_track[,4]==max(new_track[,4]))]<-0
                                exon$ORF_score_rna<-NA
                                if(sum(new_track[,4])>2){
                                        Phase_Centered_sites_frame_corr<-sum(new_track[seq(1,length,by=3),4])
                                        Phase_Centered_sites_frame_1_corr<-sum(new_track[seq(2,length,by=3),4])
                                        Phase_Centered_sites_frame_2_corr<-sum(new_track[seq(3,length,by=3),4])
                                        score1<-((Phase_Centered_sites_frame_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                        score2<-((Phase_Centered_sites_frame_1_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                        score3<-((Phase_Centered_sites_frame_2_corr-sum(new_track[,4])/3)^2)/(sum(new_track[,4])/3)
                                        exon$ORF_score_rna<-log2(score1+score2+score3+1)
                                }
                        }
                }
                
                
                
        }
        
        
        exon$chisq_ribo<-NA
        exon$chisq_rna<-NA
        
        
        if(P_sites_sum>15 & length>5){
                exon$chisq_ribo<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
        if(P_sites_sum<16 & P_sites_sum>0 & length>5){
                exon$chisq_ribo<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }        
        
        
        if(Centered_sites_sum>15 & length>5){
                exon$chisq_rna<-chisq.test(as.table(c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2)))$p.value}
        if(Centered_sites_sum<16 & Centered_sites_sum>0 & length>5){
                exon$chisq_rna<-xmulti(obs=c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }
                
        
        exon$max_notcov_ribo<-max((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len)))
        exon$coords_notcov_ribo<-max.col(t((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len))))-max((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len)))
        
        exon$max_notcov_rna<-max((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len)))
        exon$coords_notcov_rna<-max.col(t((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len))))-max((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len)))
        
        
        
        if(strand=="-"){
                exon$coords_notcov_ribo<-length-exon$coords_notcov_ribo
                exon$coords_notcov_rna<-length-exon$coords_notcov_rna
        }
        
        if(exon$max_notcov_ribo==0){
                exon$max_notcov_ribo<-"NA"
        }
        
        if(exon$max_notcov_rna==0){
                exon$max_notcov_rna<-"NA"
        }
        
        exon$notcovered_ribo<-sum(tracks[,2] == 0)
        exon$notcovered_rna<-sum(tracks[,3] == 0)
        
        
        
        
        if(length>2){
                exon$frame_start_pred<-FRAME_MAX_phase
                exon$frame_end_pred<-(length-(FRAME_MAX_phase+1))%%3
        }
        if(x[1,2]=="-" & length>2){
                
                exon$frame_end_pred<-FRAME_MAX_phase
                exon$frame_start_pred<-(length-(FRAME_MAX_phase+1))%%3
        }
        
        
        exon
        
        
}


### This function annotates exons based on their position relative to CCDS exons

annotate_exons<-function(x){
        annot_pos<-x[,c("type","start","end","length.x","P_sites_sum","RNA_sites_sum","notcovered_ribo","notcovered_rna","nt_more","nt_more_ribocovered","nt_more_P_sites","nt_more_rnacovered","nt_more_cent_sites","overlapping_ccds_start","overlapping_ccds_end")]
        
        ccdss<-which(annot_pos$type=="ccds")
        if(length(ccdss)>0){
                ccdss_coords<-annot_pos[ccdss,2:3]
                ccdss_all<-annot_pos[ccdss,]
                
                
                middle_ex<-which(annot_pos[,1]=="exon")
                middle_ex_coords<-annot_pos[which(annot_pos[,1]=="exon"),]
                listcoordsccds<-list()
                
                for(i in seq(1,dim(ccdss_coords)[1])){
                        listcoordsccds[[i]]<-seq(from=ccdss_coords[i,1],to=ccdss_coords[i,2])
                }
                
                if(length(middle_ex)>0){
                        for(y in seq(1,dim(middle_ex_coords)[1])){
                                a<-seq(from=middle_ex_coords[y,2],to=middle_ex_coords[y,3])
                                intersect<-c()
                                beginning<-c()
                                endpos<-c()
                                for(i in seq(1,length(listcoordsccds))){
                                        b<-listcoordsccds[[i]]
                                        if(sum(a%in%b)>0){
                                                intersect[i]<-TRUE
                                                beginning[i]<-(a%in%b)[1]
                                                endpos[i]<-(a%in%b)[length(a%in%b)]} else {
                                                        intersect[i]<-FALSE
                                                        beginning[i]<-FALSE
                                                        endpos[i]<-FALSE
                                                }
                                }
                                if(sum(intersect)>0 & sum(beginning)>0 & sum(endpos)>0){middle_ex_coords[y,1]<-"inside_ccds"}
                                if(sum(intersect)>0 & sum(beginning)==0 & sum(endpos)>0){middle_ex_coords[y,1]<-"overlapping_ccds"}
                                if(sum(intersect)>0 & sum(beginning)>0 & sum(endpos)==0){middle_ex_coords[y,1]<-"overlapping_ccds"}
                                if(sum(intersect)>0 & sum(beginning)==0 & sum(endpos)==0){middle_ex_coords[y,1]<-"containing_ccds"}
                                ccds_inters<-ccdss_all[intersect,]
                                if(dim(ccds_inters)[1]>1){middle_ex_coords[y,1]<-"overlapping_multiple_ccdss"}
                                if(dim(ccds_inters)[1]==1){middle_ex_coords[y,"nt_more"]<-middle_ex_coords[y,"length.x"]-ccds_inters[,"length.x"]
                                                           middle_ex_coords[y,"nt_more_ribocovered"]<-1-((middle_ex_coords[y,"notcovered_ribo"]-ccds_inters[,"notcovered_ribo"])/middle_ex_coords[y,"nt_more"])
                                                           middle_ex_coords[y,"nt_more_P_sites"]<-middle_ex_coords[y,"P_sites_sum"]-ccds_inters[,"P_sites_sum"]
                                                           middle_ex_coords[y,"nt_more_rnacovered"]<-1-((middle_ex_coords[y,"notcovered_rna"]-ccds_inters[,"notcovered_rna"])/middle_ex_coords[y,"nt_more"])
                                                           middle_ex_coords[y,"nt_more_cent_sites"]<-middle_ex_coords[y,"RNA_sites_sum"]-ccds_inters[,"RNA_sites_sum"]
                                                           middle_ex_coords[y,"overlapping_ccds_start"]<-ccds_inters[,"start"]
                                                           middle_ex_coords[y,"overlapping_ccds_end"]<-ccds_inters[,"end"]}
                        }
                }
                
                inside_ex<-middle_ex_coords[,1]=="inside_ccds"
                
                
                if(length(middle_ex)>0){
                        middle_ex_coords[middle_ex_coords[,2]%in%ccdss_coords[,1] & middle_ex_coords[,1]!="overlapping_multiple_ccdss",1]<-"exon_alt_donor"
                        middle_ex_coords[middle_ex_coords[,3]%in%ccdss_coords[,2] & middle_ex_coords[,1]!="overlapping_multiple_ccdss",1]<-"exon_alt_acceptor"
                        middle_ex_coords[middle_ex_coords[,2]%in%ccdss_coords[,1] & middle_ex_coords[,1]=="overlapping_multiple_ccdss",1]<-"overlapping_multiple_ccdss_alt_donor"
                        middle_ex_coords[middle_ex_coords[,3]%in%ccdss_coords[,2] & middle_ex_coords[,1]=="overlapping_multiple_ccdss",1]<-"overlapping_multiple_ccdss_alt_acceptor"
                        annot_pos[middle_ex,]<-middle_ex_coords
                }
                
                if(sum(inside_ex)>0){
                        middle_ex_coords[inside_ex & middle_ex_coords[,1]=="exon_alt_donor",1]<-"int_exon_alt_donor"
                        middle_ex_coords[inside_ex & middle_ex_coords[,1]=="exon_alt_acceptor",1]<-"int_exon_alt_acceptor"  
                }
                annot_pos[middle_ex,]<-middle_ex_coords
                
                
                
                annot_pos[1:(ccdss[1]-1),1]<-"5_utrs_ex"
                annot_pos[ccdss,1]<-"ccds"
                coords_start<-c(annot_pos[ccdss[1],2],annot_pos[ccdss[1],3])
                ccdss_start<-annot_pos[ccdss[1],]
                five_with_cds<-which(annot_pos[,2]<=coords_start[1] & annot_pos[,3]>=coords_start[2])
                five_with_cds<-five_with_cds[!five_with_cds%in%ccdss]
                annot_pos[five_with_cds,1]<-"5_utrs_st"
                annot_pos_fiveutr<-annot_pos[five_with_cds,]
                for(f in seq(1,dim(annot_pos_fiveutr)[1])){
                        annot_pos_fiveutr[f,"nt_more"]<-as.numeric(annot_pos_fiveutr[f,"length.x"]-ccdss_start[,"length.x"])
                        annot_pos_fiveutr[f,"nt_more_ribocovered"]<-1-((annot_pos_fiveutr[f,"notcovered_ribo"]-ccdss_start[,"notcovered_ribo"])/annot_pos_fiveutr[f,"nt_more"])
                        annot_pos_fiveutr[f,"nt_more_P_sites"]<-annot_pos_fiveutr[f,"P_sites_sum"]-ccdss_start[,"P_sites_sum"]
                        annot_pos_fiveutr[f,"nt_more_rnacovered"]<-1-((annot_pos_fiveutr[f,"notcovered_rna"]-ccdss_start[,"notcovered_rna"])/annot_pos_fiveutr[f,"nt_more"])
                        annot_pos_fiveutr[f,"nt_more_cent_sites"]<-annot_pos_fiveutr[f,"RNA_sites_sum"]-ccdss_start[,"RNA_sites_sum"]
                        annot_pos_fiveutr[f,"overlapping_ccds_start"]<-ccdss_start[1,2]
                        annot_pos_fiveutr[f,"overlapping_ccds_end"]<-ccdss_start[1,3]
                }
                annot_pos[five_with_cds,]<-annot_pos_fiveutr
                
                
                annot_pos[(1+(ccdss[length(ccdss)])):dim(annot_pos)[1],1]<-"3_utrs_ex"
                annot_pos[ccdss,1]<-"ccds"
                coords_stop<-c(annot_pos[tail(ccdss,1),2],annot_pos[tail(ccdss,1),3])
                ccdss_stop<-annot_pos[tail(ccdss,1),]
                three_with_cds<-which(annot_pos[,2]<=coords_stop[1] & annot_pos[,3]>=coords_stop[2])
                three_with_cds<-three_with_cds[!three_with_cds%in%ccdss]
                annot_pos[three_with_cds,1]<-"3_utrs_st"
                annot_pos_threeutr<-annot_pos[three_with_cds,]
                for(f in seq(1,dim(annot_pos_threeutr)[1])){
                        annot_pos_threeutr[f,"nt_more"]<-as.numeric(annot_pos_threeutr[f,"length.x"]-ccdss_stop[,"length.x"])
                        annot_pos_threeutr[f,"nt_more_ribocovered"]<-1-((annot_pos_threeutr[f,"notcovered_ribo"]-ccdss_stop[,"notcovered_ribo"])/annot_pos_threeutr[f,"nt_more"])
                        annot_pos_threeutr[f,"nt_more_P_sites"]<-annot_pos_threeutr[f,"P_sites_sum"]-ccdss_stop[,"P_sites_sum"]
                        annot_pos_threeutr[f,"nt_more_rnacovered"]<-1-((annot_pos_threeutr[f,"notcovered_rna"]-ccdss_stop[,"notcovered_rna"])/annot_pos_threeutr[f,"nt_more"])
                        annot_pos_threeutr[f,"nt_more_cent_sites"]<-annot_pos_threeutr[f,"RNA_sites_sum"]-ccdss_stop[,"RNA_sites_sum"]
                        annot_pos_threeutr[f,"overlapping_ccds_start"]<-ccdss_stop[,"start"]
                        annot_pos_threeutr[f,"overlapping_ccds_end"]<-ccdss_stop[,"end"]
                }
                annot_pos[ three_with_cds,]<-annot_pos_threeutr
                
                if(x$strand.x[1]=="-"){
                        int_don<-which(annot_pos[,1]=="int_exon_alt_donor")
                        int_acc<-which(annot_pos[,1]=="int_exon_alt_acceptor")
                        don<-which(annot_pos[,1]=="exon_alt_donor")
                        acc<-which(annot_pos[,1]=="exon_alt_acceptor")
                        multi_don<-which(annot_pos[,1]=="overlapping_multiple_ccdss_alt_donor")
                        multi_acc<-which(annot_pos[,1]=="overlapping_multiple_ccdss_alt_acceptor")
                        fiveex<-which(annot_pos[,1]=="5_utrs_ex")
                        fivest<-which(annot_pos[,1]=="5_utrs_st")
                        threeex<-which(annot_pos[,1]=="3_utrs_ex")
                        threest<-which(annot_pos[,1]=="3_utrs_st")
                        annot_pos[don,1]<-"exon_alt_acceptor"
                        annot_pos[acc,1]<-"exon_alt_donor"
                        annot_pos[int_don,1]<-"int_exon_alt_acceptor"
                        annot_pos[int_acc,1]<-"int_exon_alt_donor"
                        annot_pos[multi_don,1]<-"overlapping_multiple_ccdss_alt_acceptor"
                        annot_pos[multi_acc,1]<-"overlapping_multiple_ccdss_alt_donor"
                        annot_pos[fiveex,1]<-"3_utrs_ex"
                        annot_pos[fivest,1]<-"3_utrs_st"
                        annot_pos[threeex,1]<-"5_utrs_ex"
                        annot_pos[threest,1]<-"5_utrs_st"
                        
                }
                annot_pos<-annot_pos[!is.na(annot_pos[,"start"]),]
        }
        
        
        x[,c("type","start","end","length.x","P_sites_sum","RNA_sites_sum","notcovered_ribo","notcovered_rna","nt_more","nt_more_ribocovered","nt_more_P_sites","nt_more_rnacovered","nt_more_cent_sites","overlapping_ccds_start","overlapping_ccds_end")]<-annot_pos
        x
}


### This function calculates periodicity on NON-CCDS region of an exons


alt_exon_analysis<-function(x,sequences=seq_exons,tracks_exons=all_tracks,index_tracks=tracks_index){
        
        
        exon<-x
        names_exons<-names(sequences)
        
        seq_exon<-sequences[which(names_exons%in%exon["coords2"])][[1]]
        myexon_id<-exon[,"exon_id"]
        exon_track<-tracks_exons[index_tracks==myexon_id]
        withsep<-strsplit(exon_track,split=" ")
        x<-t(data.frame(withsep))
        
        strand<-x[1,2]
        tracks_pre<-t(x[,-c(1:2)])
        
        if(strand=="-"){
                tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]),rev(tracks_pre[,5]))
        } else if (strand=="+"){
                tracks<-tracks_pre}
        colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent","Seq")
        
        tracks<-tracks[,1:4]
        
        mode(tracks)<-"numeric"
        length<-dim(tracks)[1]
        exon$exon_id_noccds<-exon$exon_id
        if(exon$type=="exon_alt_acceptor"){
                
                tracks<-tracks[1:exon$nt_more,]
                seq_exon<-seq_exon[1:(exon$nt_more)]
                length<-exon$nt_more
                if(strand=="+"){exon$end=exon$overlapping_ccds_start-1}
                if(strand=="-"){exon$start=exon$overlapping_ccds_end+1}
                exon$exon_id_noccds<-paste(exon$chr,exon$start,exon$end,exon$type,exon$gene_id,sep="_")
        }
        
        
        if(exon$type=="exon_alt_donor"){
                tracks<-tracks[(length+1-exon$nt_more):length,]
                seq_exon<-seq_exon[(length+1-exon$nt_more):length]
                length<-exon$nt_more
                if(strand=="+"){exon$start=exon$overlapping_ccds_end+1}
                if(strand=="-"){exon$end=exon$overlapping_ccds_start-1}
                exon$exon_id_noccds<-paste(exon$chr,exon$start,exon$end,exon$type,exon$gene_id,sep="_")
        }
        
        exon<-data.frame(exon_id=exon$exon_id_noccds,exon_id_orig=exon$exon_id,type=exon$type,gene_id=exon$gene_id,annotation=exon$annotation)
        exon$strand<-strand
        exon$length<-dim(tracks)[1]
        exon$frame_start_pred<-NA
        exon$frame_end_pred<-NA
        
        
        length<-dim(tracks)[1]
        
        
        P_sites_sum<-round(sum(tracks[,1]),digits=6)
        exon$P_sites_sum<-P_sites_sum
        
        
        Centered_sites_sum<-round(sum(tracks[,4]),digits=6)
        exon$RNA_sites_sum<-Centered_sites_sum
        exon$Ribocov_aver<-round(mean(tracks[,2]),digits=6)
        exon$RNAseqcov_aver<-round(mean(tracks[,3]),digits=6)
        exon$pctPhase_frame<-NA
        exon$pctPhase_frame_1<-NA
        exon$pctPhase_frame_2<-NA
        exon$pctPhaseCentered_frame<-NA
        exon$pctPhaseCentered_frame_1<-NA
        exon$pctPhaseCentered_frame_2<-NA
        if(length>2){
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                
                
                exon$pctPhase_frame<-Phase_P_sites_frame/P_sites_sum
                exon$pctPhase_frame_1<-Phase_P_sites_frame_1/P_sites_sum
                exon$pctPhase_frame_2<-Phase_P_sites_frame_2/P_sites_sum
                
                
                Phase_Centered_sites_frame<-sum(tracks[seq(1,length,by=3),4])
                Phase_Centered_sites_frame_1<-sum(tracks[seq(2,length,by=3),4])
                Phase_Centered_sites_frame_2<-sum(tracks[seq(3,length,by=3),4])
                
                
                exon$pctPhaseCentered_frame<-Phase_Centered_sites_frame/Centered_sites_sum
                exon$pctPhaseCentered_frame_1<-Phase_Centered_sites_frame_1/Centered_sites_sum
                exon$pctPhaseCentered_frame_2<-Phase_Centered_sites_frame_2/Centered_sites_sum
                
                
                MAXPhase_frame<-max(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2))
                FRAME_MAX_phase<-max.col(t(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2)))-1
                
                MAXPhaseCentered_frame<-max(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2))
                FRAME_MAX_phaseCentered<-max.col(t(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2)))-1
        }
        
        
        exon$multit_freq_best_ribo<-NA
        exon$pval_multit_3nt_ribo<-NA
        exon$spec_multit_3nt_ribo<-NA
        exon$fft_max_freq_ribo<-NA
        exon$fft_power_3_ribo<-NA
        exon$fft_aver_ribo<-NA
        exon$spec_max_freq_ribo<-NA
        exon$spec_power_3_ribo<-NA
        exon$spec_aver_power_ribo<-NA
        
        exon$multit_freq_best_rna<-NA
        exon$pval_multit_3nt_rna<-NA
        exon$spec_multit_3nt_rna<-NA
        exon$fft_max_freq_rna<-NA
        exon$fft_power_3_rna<-NA
        exon$fft_aver_rna<-NA
        exon$spec_max_freq_rna<-NA
        exon$spec_power_3_rna<-NA
        exon$spec_aver_power_rna<-NA
        
        exon$ORF_score_ribo<-NA
        exon$ORF_score_rna<-NA
        
        if(P_sites_sum>2 & length>5){
                if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                bestfreq_3ntpval_ribo<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,1],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                exon$multit_freq_best_ribo<-bestfreq_3ntpval_ribo[1]
                exon$pval_multit_3nt_ribo<-bestfreq_3ntpval_ribo[2]
                exon$spec_multit_3nt_ribo<-bestfreq_3ntpval_ribo[3]
                score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                exon$ORF_score_ribo<-log2(score1+score2+score3+1)
                
                
                gino<-getFFTFreqs(Nyq.Freq=0.5,data=tracks[,1])
                modFFT <- Mod(fft(tracks[,1]))
                FFTdata <- cbind(gino, modFFT)
                exon$fft_aver_ribo<-mean(FFTdata[,2])
                exon$fft_power_3_ribo<-FFTdata[which(abs((gino-(1/3)))==min(abs((gino-(1/3))))),2]
                exon$fft_max_freq_ribo<-abs(gino[which(FFTdata==max((FFTdata[10:dim(FFTdata)[1]/2,2])),arr.ind=TRUE)[1]])[1]
                
                
                spect_P_sites<-spectrum(tracks[,1],plot=FALSE)
                exon$spec_max_freq_ribo<-spect_P_sites$freq[which(spect_P_sites$spec==max(spect_P_sites$spec),arr.ind=TRUE)][1]
                exon$spec_power_3_ribo<-spect_P_sites$spec[which(abs((spect_P_sites$freq-(1/3)))==min(abs((spect_P_sites$freq-(1/3)))))]
                exon$spec_aver_power_ribo<-mean(spect_P_sites$spec)        
                if(Centered_sites_sum>2){
                        
                        gino<-getFFTFreqs(Nyq.Freq=0.5,data=tracks[,4])
                        modFFT <- Mod(fft(tracks[,4]))
                        FFTdata <- cbind(gino, modFFT)
                        exon$fft_aver_rna<-mean(FFTdata[,2])
                        exon$fft_power_3_rna<-FFTdata[which(abs((gino-(1/3)))==min(abs((gino-(1/3))))),2]
                        exon$fft_max_freq_rna<-1/abs(gino[which(FFTdata==max((FFTdata[10:dim(FFTdata)[1]/2,2])),arr.ind=TRUE)[1]])[1]
                        
                        
                        
                        spect_rna<-spectrum(tracks[,4],plot=FALSE)
                        exon$spec_max_freq_rna<-spect_rna$freq[which(spect_rna$spec==max(spect_rna$spec),arr.ind=TRUE)][1]
                        exon$spec_power_3_rna<-spect_rna$spec[which(abs((spect_rna$freq-(1/3)))==min(abs((spect_rna$freq-(1/3)))))]
                        exon$spec_aver_power_rna<-mean(spect_rna$spec)
                        bestfreq_3ntpval_rna<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,4],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,7)]
                        
                        exon$multit_freq_best_rna<-bestfreq_3ntpval_rna[1]
                        exon$pval_multit_3nt_rna<-bestfreq_3ntpval_rna[2]
                        exon$spec_multit_3nt_rna<-bestfreq_3ntpval_rna[3]
                
                        score_rna_1<-((Phase_Centered_sites_frame-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        score_rna_2<-((Phase_Centered_sites_frame_1-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        score_rna_3<-((Phase_Centered_sites_frame_2-Centered_sites_sum/3)^2)/(Centered_sites_sum/3)
                        exon$ORF_score_rna<-log2(score_rna_1+score_rna_2+score_rna_3+1)
                }
                
                
                
        }
        
        
        exon$chisq_ribo<-NA
        exon$chisq_rna<-NA
        
        
        if(P_sites_sum>15 & length>5){
                exon$chisq_ribo<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
        if(P_sites_sum<16 & P_sites_sum>0 & length>5){
                exon$chisq_ribo<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }        
        
        
        if(Centered_sites_sum>15 & length>5){
                exon$chisq_rna<-chisq.test(as.table(c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2)))$p.value}
        if(Centered_sites_sum<16 & Centered_sites_sum>0 & length>5){
                exon$chisq_rna<-xmulti(obs=c(Phase_Centered_sites_frame,Phase_Centered_sites_frame_1,Phase_Centered_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
        }
        
        
        exon$max_notcov_ribo<-max((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len)))
        exon$coords_notcov_ribo<-max.col(t((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len))))-max((!tracks[,2]) * unlist(lapply(rle(tracks[,2])$lengths, seq_len)))
        
        exon$max_notcov_rna<-max((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len)))
        exon$coords_notcov_rna<-max.col(t((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len))))-max((!tracks[,3]) * unlist(lapply(rle(tracks[,3])$lengths, seq_len)))
        
        
        
        if(strand=="-"){
                exon$coords_notcov_ribo<-length-exon$coords_notcov_ribo
                exon$coords_notcov_rna<-length-exon$coords_notcov_rna
        }
        
        if(exon$max_notcov_ribo==0){
                exon$max_notcov_ribo<-"NA"
        }
        
        if(exon$max_notcov_rna==0){
                exon$max_notcov_rna<-"NA"
        }
        
        exon$notcovered_ribo<-sum(tracks[,2] == 0)
        exon$notcovered_rna<-sum(tracks[,3] == 0)
        
        
        
        
        if(length>2){
                exon$frame_start_pred<-FRAME_MAX_phase
                exon$frame_end_pred<-(length-(FRAME_MAX_phase+1))%%3
        }
        if(x[1,2]=="-" & length>2){
                
                exon$frame_end_pred<-FRAME_MAX_phase
                exon$frame_start_pred<-(length-(FRAME_MAX_phase+1))%%3
        }
        
        pept<-NA
        exon$transl_pept_notccds<-NA
        if(P_sites_sum>0){
                if(exon$strand=="-"){
                        pept<-unlist(getTrans(seq_exon,sens="F",frame=exon$frame_end_pred))
                } else {pept<-unlist(getTrans(seq_exon,sens="F",frame=exon$frame_start_pred))}
                exon$transl_pept_notccds<-paste(pept,sep="",collapse="")
        }
        
        return(exon)
}



### This function calculates coherence values for candidate regions with multi-frame translation


calculate_coherence<-function(x){
        strand<-x[1,2]
        tracks_pre<-t(x[,-c(1:2)])
        
        
        if(strand=="-"){
                tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
        } else if (strand=="+"){
                tracks<-tracks_pre}
        colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
        mode(tracks)<-"numeric"
        
        exon<-data.frame(exon_id=x[1,1],stringsAsFactors=F,row.names=NULL)
        exon$strand<-strand
        exon$frame_start_pred<-NA
        exon$frame_end_pred<-NA
        
        exon$length<-dim(tracks)[1]
        length<-dim(tracks)[1]
        
        
        P_sites_sum<-round(sum(tracks[,1]),digits=6)
        exon$P_sites_sum<-P_sites_sum
        
        
        Centered_sites_sum<-round(sum(tracks[,4]),digits=6)
        exon$RNA_sites_sum<-Centered_sites_sum
        exon$Ribocov_aver<-round(mean(tracks[,2]),digits=6)
        exon$RNAseqcov_aver<-round(mean(tracks[,3]),digits=6)
        exon$pctPhase_frame<-NA
        exon$pctPhase_frame_1<-NA
        exon$pctPhase_frame_2<-NA
        exon$pctPhaseCentered_frame<-NA
        exon$pctPhaseCentered_frame_1<-NA
        exon$pctPhaseCentered_frame_2<-NA
        if(length>2){
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                
                
                exon$pctPhase_frame<-Phase_P_sites_frame/P_sites_sum
                exon$pctPhase_frame_1<-Phase_P_sites_frame_1/P_sites_sum
                exon$pctPhase_frame_2<-Phase_P_sites_frame_2/P_sites_sum
                
                
                Phase_Centered_sites_frame<-sum(tracks[seq(1,length,by=3),4])
                Phase_Centered_sites_frame_1<-sum(tracks[seq(2,length,by=3),4])
                Phase_Centered_sites_frame_2<-sum(tracks[seq(3,length,by=3),4])
                
                
                exon$pctPhaseCentered_frame<-Phase_Centered_sites_frame/Centered_sites_sum
                exon$pctPhaseCentered_frame_1<-Phase_Centered_sites_frame_1/Centered_sites_sum
                exon$pctPhaseCentered_frame_2<-Phase_Centered_sites_frame_2/Centered_sites_sum
                
                
                MAXPhase_frame<-max(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2))
                FRAME_MAX_phase<-max.col(t(c(exon$pctPhase_frame,exon$pctPhase_frame_1,exon$pctPhase_frame_2)))-1
                
                MAXPhaseCentered_frame<-max(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2))
                FRAME_MAX_phaseCentered<-max.col(t(c(exon$pctPhaseCentered_frame,exon$pctPhaseCentered_frame_1,exon$pctPhaseCentered_frame_2)))-1
        }
        
        
        exon$multit_freq_best_ribo<-NA
        exon$pval_multit_3nt_ribo<-NA
        exon$spec_multit_3nt_ribo<-NA

        exon$coherence_1_2_ribo<-NA
        exon$coherence_1_3_ribo<-NA
        exon$coherence_2_3_ribo<-NA
        exon$min_coherence_ribo<-NA
        exon$multit_freq_best_rna<-NA
        exon$pval_multit_3nt_rna<-NA
        exon$spec_multit_3nt_rna<-NA
        exon$coherence_1_2_rna<-NA
        exon$coherence_1_3_rna<-NA
        exon$coherence_2_3_rna<-NA
        exon$min_coherence_rna<-NA
        if(P_sites_sum>10 & length>5){
                
                if(length<25){slepians<-dpss(n=length+(50-length),k=24,nw=12)}
                if(length>=25){slepians<-dpss(n=length,k=24,nw=12)}
                bestfreq_3ntpval_ribo<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,1],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,5,7)]
                exon$multit_freq_best_ribo<-bestfreq_3ntpval_ribo[1]
                exon$pval_multit_3nt_ribo<-bestfreq_3ntpval_ribo[2]
                exon$spec_multit_3nt_ribo<-bestfreq_3ntpval_ribo[4]
                
                y<-tracks[,1]
                
                if(length(y)<25){
                        remain<-50-length(y)
                        y<-c(rep(0,as.integer(remain/2)),y,rep(0,remain%%2+as.integer(remain/2)))
                }
                if(length(y)<1024/2){padding<-1024}
                if(length(y)>=1024/2){padding<-"default"}
                length<-length(y)
                
                y1<-rep(0,length)
                y2<-rep(0,length)
                y3<-rep(0,length)
                
                y1[seq(1,length,by=3)]<-y[seq(1,length,by=3)]
                y2[seq(2,length,by=3)]<-y[seq(2,length,by=3)]
                y3[seq(3,length,by=3)]<-y[seq(3,length,by=3)]
                
                
                
                spec_y1<-spec.mtm(timeSeries=as.ts(y1),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                spec_y2<-spec.mtm(timeSeries=as.ts(y2),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                spec_y3<-spec.mtm(timeSeries=as.ts(y3),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                
                coh1_2<-mtm.coh(spec_y1,spec_y2,plot=F)
                coh1_3<-mtm.coh(spec_y1,spec_y3,plot=F)
                coh2_3<-mtm.coh(spec_y2,spec_y3,plot=F)
                
                exon$coherence_1_2_ribo<-coh1_2$msc[which(coh1_2$freq==bestfreq_3ntpval_ribo[3])]
                exon$coherence_1_3_ribo<-coh1_3$msc[which(coh1_3$freq==bestfreq_3ntpval_ribo[3])]
                exon$coherence_2_3_ribo<-coh2_3$msc[which(coh2_3$freq==bestfreq_3ntpval_ribo[3])]
                exon$min_coherence_ribo<-min(c(exon$coherence_1_2_ribo,exon$coherence_1_3_ribo,exon$coherence_2_3_ribo),na.rm=T)
                if((Phase_Centered_sites_frame > 5 & Phase_Centered_sites_frame_1 > 5) | (Phase_Centered_sites_frame > 5 & Phase_Centered_sites_frame_2 > 5) | (Phase_Centered_sites_frame_1 > 5 & Phase_Centered_sites_frame_2 > 5)){
                        
                        bestfreq_3ntpval_rna<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,4],n_tapers=24,time_bw=12,slepians_values=slepians)[c(1,6,5,7)]
                        exon$multit_freq_best_rna<-bestfreq_3ntpval_rna[1]
                        exon$pval_multit_3nt_rna<-bestfreq_3ntpval_rna[2]
                        exon$spec_multit_3nt_rna<-bestfreq_3ntpval_rna[4]
                        
                        y<-tracks[,4]
                        
                        if(length(y)<25){
                                remain<-50-length(y)
                                y<-c(rep(0,as.integer(remain/2)),y,rep(0,remain%%2+as.integer(remain/2)))
                        }
                        if(length(y)<1024/2){padding<-1024}
                        if(length(y)>=1024/2){padding<-"default"}
                        length<-length(y)
                        
                        y1<-rep(0,length)
                        y2<-rep(0,length)
                        y3<-rep(0,length)
                        
                        y1[seq(1,length,by=3)]<-y[seq(1,length,by=3)]
                        y2[seq(2,length,by=3)]<-y[seq(2,length,by=3)]
                        y3[seq(3,length,by=3)]<-y[seq(3,length,by=3)]
                        
                        
                        
                        spec_y1<-spec.mtm(timeSeries=as.ts(y1),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                        spec_y2<-spec.mtm(timeSeries=as.ts(y2),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                        spec_y3<-spec.mtm(timeSeries=as.ts(y3),nw=12,k=24,dpssIN=slepians,returnInternals=T,plot=F,nFFT=padding)
                        
                        coh1_2<-mtm.coh(spec_y1,spec_y2,plot=F)
                        coh1_3<-mtm.coh(spec_y1,spec_y3,plot=F)
                        coh2_3<-mtm.coh(spec_y2,spec_y3,plot=F)
                        
                        exon$coherence_1_2_rna<-coh1_2$msc[which(coh1_2$freq==bestfreq_3ntpval_rna[3])]
                        exon$coherence_1_3_rna<-coh1_3$msc[which(coh1_3$freq==bestfreq_3ntpval_rna[3])]
                        exon$coherence_2_3_rna<-coh2_3$msc[which(coh2_3$freq==bestfreq_3ntpval_rna[3])]
                        exon$min_coherence_rna<-min(c(exon$coherence_1_2_rna,exon$coherence_1_3_rna,exon$coherence_2_3_rna),na.rm=T)
                        
                }
        }
        exon
}


### This function calculates exonic information on nonCCDS ORFs, to calculate multimapping information and CDS overlaps


pre_multi_nonCCDS_ORFs<-function(x,counter,all_exons_in_the_sign_transcr=exons_transcr_nonccds_sign,signif_exons=nonccds_res_sign){
        transcr<-x[,"transcript_id"]
        trascr_length<-x$length
        orf_strand<-x$strand
        transcr_data<-data.frame(transcript_id=transcr)
        
        exons_in_transcr<-all_exons_in_the_sign_transcr[all_exons_in_the_sign_transcr[,4]%in%transcr,"exon_id"]
        if(orf_strand=="-"){exons_in_transcr<-rev(exons_in_transcr)}
        
        exons_in_transcr_data<-nonccds_res[nonccds_res[,"exon_id"]%in%exons_in_transcr,]
        exons_in_transcr_data<-exons_in_transcr_data[match(exons_in_transcr,exons_in_transcr_data$exon_id),]
        
        orf_start<-x$start_pos
        orf_end<-x$st2vect
        cumsumexons<-cumsum(exons_in_transcr_data$length.x)
        
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
        if(x$strand=="+"){coord_start<-exons_in_transcr_data[st_ex,"start"] + (orf_start-rem_len)}
        if(x$strand=="-"){coord_start<-exons_in_transcr_data[st_ex,"end"] - (orf_start-rem_len)}
        
        if(length(in_betw_ex)==0){
                if(st_ex==end_ex){nt_to_rem<-0}
                if(st_ex!=end_ex){if(x$strand=="+"){
                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                }
                                  if(x$strand=="-"){
                                          nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                                  }
                }
        }
        
        if(length(in_betw_ex)>0){
                nt_in_betw<-sum(exons_in_transcr_data[in_betw_ex,"length.x"])
                if(x$strand=="+"){
                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                }
                if(x$strand=="-"){
                        nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                }
                nt_to_rem<-nt_to_rem+nt_in_betw
        }
        
        if(st_ex==end_ex & x$strand=="+"){coord_end<-coord_start+x$ORF_length+1}
        if(st_ex==end_ex & x$strand=="-"){coord_end<-coord_start-x$ORF_length+1}
        
        if(st_ex!=end_ex & x$strand=="+"){coord_end<-exons_in_transcr_data[end_ex,"start"] + (x$ORF_length-nt_to_rem)+1}
        if(st_ex!=end_ex & x$strand=="-"){coord_end<-exons_in_transcr_data[end_ex,"end"] - (x$ORF_length-nt_to_rem)+1}
        
        if(x$strand=="-"){
                coord_start2<-coord_start
                coord_start<-coord_end
                coord_end<-coord_start2
        }
        
        
        if(st_ex!=end_ex & x$strand=="+"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,exons_in_transcr_data[st_ex,"end"],"EXONnonCCDS",x$gene_id,x$strand,sep="_")
                                          to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],exons_in_transcr_data[end_ex,"start"],coord_end,"EXONnonCCDS",x$gene_id,x$strand,sep="_")
                                          to_check<-paste(to_check_st,to_check_end,sep=";")
                                          
        }
        if(st_ex!=end_ex & x$strand=="-"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],exons_in_transcr_data[st_ex,"start"],coord_end,"EXONnonCCDS",x$gene_id,x$strand,sep="_")
                                          to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],coord_start,exons_in_transcr_data[end_ex,"end"],"EXONnonCCDS",x$gene_id,x$strand,sep="_")
                                          to_check<-paste(to_check_st,to_check_end,sep=";")
        }
        
        if(st_ex==end_ex){to_check<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,"EXONnonCCDS",x$gene_id,x$strand,sep="_")}
        x$to_check<-to_check
        x$to_check_rem<-NA
        if(length(in_betw_ex)>0){
                x$to_check_rem<-paste(exon_inbetween_data$exon_id,collapse=";")
                
        }
        x$ORF_id_tr<-paste(transcr_data$transcript_id,orf_start,orf_end,sep="_")
        x$ORF_id_gen<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,sep="_")
        x
        
        
}

### This function calculates exonic information on CCDS ORFs, to calculate multimapping information and CDS overlaps



pre_multi_CCDS_ORFs<-function(x,counter,all_exons_in_the_sign_transcr=exons_transcr_nonccds_sign,signif_exons=nonccds_res){
        transcr<-x[,"transcript_id"]
        trascr_length<-x$length
        orf_strand<-x$strand
        transcr_data<-data.frame(transcript_id=transcr)
        
        exons_in_transcr<-all_exons_in_the_sign_transcr[all_exons_in_the_sign_transcr[,4]%in%transcr,"coords_id"]
        if(orf_strand=="-"){exons_in_transcr<-rev(exons_in_transcr)}
        
        exons_in_transcr_data<-signif_exons[signif_exons[,"coords"]%in%exons_in_transcr,]
        exons_in_transcr_data<-exons_in_transcr_data[match(exons_in_transcr,exons_in_transcr_data$coords),]
        
        orf_start<-x$start_pos
        orf_end<-x$st2vect
        cumsumexons<-cumsum(exons_in_transcr_data$length.x)
        
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
        if(x$strand=="+"){coord_start<-exons_in_transcr_data[st_ex,"start"] + (orf_start-rem_len)}
        if(x$strand=="-"){coord_start<-exons_in_transcr_data[st_ex,"end"] - (orf_start-rem_len)}
        
        if(length(in_betw_ex)==0){
                if(st_ex==end_ex){nt_to_rem<-0}
                if(st_ex!=end_ex){if(x$strand=="+"){
                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                }
                                  if(x$strand=="-"){
                                          nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                                  }
                }
        }
        
        if(length(in_betw_ex)>0){
                nt_in_betw<-sum(exons_in_transcr_data[in_betw_ex,"length.x"])
                if(x$strand=="+"){
                        nt_to_rem<-exons_in_transcr_data[st_ex,"end"]-coord_start
                }
                if(x$strand=="-"){
                        nt_to_rem<-coord_start-exons_in_transcr_data[st_ex,"start"]
                }
                nt_to_rem<-nt_to_rem+nt_in_betw
        }
        
        if(st_ex==end_ex & x$strand=="+"){coord_end<-coord_start+x$ORF_length+1}
        if(st_ex==end_ex & x$strand=="-"){coord_end<-coord_start-x$ORF_length+1}
        
        if(st_ex!=end_ex & x$strand=="+"){coord_end<-exons_in_transcr_data[end_ex,"start"] + (x$ORF_length-nt_to_rem)+1}
        if(st_ex!=end_ex & x$strand=="-"){coord_end<-exons_in_transcr_data[end_ex,"end"] - (x$ORF_length-nt_to_rem)+1}
        
        if(x$strand=="-"){
                coord_start2<-coord_start
                coord_start<-coord_end
                coord_end<-coord_start2
        }
        
        
        if(st_ex!=end_ex & x$strand=="+"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,exons_in_transcr_data[st_ex,"end"],"CCDS",x$gene_id,x$strand,sep="_")
                                          to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],exons_in_transcr_data[end_ex,"start"],coord_end,"CCDS",x$gene_id,x$strand,sep="_")
                                          to_check<-paste(to_check_st,to_check_end,sep=";")
                                          
        }
        if(st_ex!=end_ex & x$strand=="-"){to_check_st<-paste(exons_in_transcr_data[st_ex,"chr"],exons_in_transcr_data[st_ex,"start"],coord_end,"CCDS",x$gene_id,x$strand,sep="_")
                                          to_check_end<-paste(exons_in_transcr_data[end_ex,"chr"],coord_start,exons_in_transcr_data[end_ex,"end"],"CCDS",x$gene_id,x$strand,sep="_")
                                          to_check<-paste(to_check_st,to_check_end,sep=";")
        }
        
        if(st_ex==end_ex){to_check<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,"CCDS",x$gene_id,x$strand,sep="_")}
        x$to_check<-to_check
        x$to_check_rem<-NA
        if(length(in_betw_ex)>0){
                x$to_check_rem<-paste(exon_inbetween_data$exon_id,collapse=";")
                
        }
        x$ORF_id_tr<-paste(transcr_data$transcript_id,orf_start,orf_end,sep="_")
        x$ORF_id_gen<-paste(exons_in_transcr_data[st_ex,"chr"],coord_start,coord_end,sep="_")
        x
        
        
}




### This function calculates results for real and simulated exons for the multitaper analysis


take_simuls_multi<-function(x,tapers,bw,nsimul){
        unique_ex_id<-x[,"exon_id"]
        list_exons_tracks<-list()
        for(i in seq(1:length(unique_ex_id))){
                list_exons_tracks[[i]]<-all_tracks[index==unique_ex_id[i]]
        }
        simuls_eachexons<-list()
        for(s in 1:length(unique_ex_id)){
                withsep<-strsplit(list_exons_tracks[[s]],split=" ")
                x<-t(data.frame(withsep))
                id<-unique_ex_id[s]
                exon<-data.frame(exon_id=id,stringsAsFactors=F,row.names=NULL)
                strand<-x[1,2]
                tracks_pre<-t(x[,-c(1:2)])
                if(strand=="-"){
                        tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
                } else if (strand=="+"){
                        tracks<-tracks_pre}
                colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
                mode(tracks)<-"numeric"
                length<-dim(tracks)[1]
                exon$length<-length
                
                if(length<25){
                        slepians<-dpss(n=length+(50-length),k=tapers,nw=bw)
                }
                if(length>=25){
                        slepians<-dpss(n=length,k=tapers,nw=bw)
                }
                
                exon$pval_multi_ribo<-take_freqs_Fvalues_all_around_3nt_spec(x=tracks[,1],n_tapers=tapers,time_bw=bw,slepians_values=slepians)[6]
                ribo_covered_pos<-which(tracks[,2]>0)
                P_sites_sum<-sum(tracks[,1])
                exon$P_sites_sum<-P_sites_sum
                exon$RNA_sites_sum<-sum(tracks[,4])
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                
                simuls_results<-foreach(j=1:nsimul,.combine=c,.multicombine=T) %dopar%{
                        set.seed(j)
                        simtrack<-rep(0,length)
                        rand_pos<-sample(ribo_covered_pos,P_sites_sum,replace=T)
                        for(i in rand_pos){
                                simtrack[i]<-simtrack[i]+1
                        }
                        
                        simul_Pval_multi_3nt<-take_freqs_Fvalues_all_around_3nt_spec(x=simtrack,n_tapers=tapers,time_bw=bw,slepians_values=slepians)[6]
                        
                        return(simul_Pval_multi_3nt)
                }
                exon$n_simul_sign_multi<-sum(simuls_results<0.05)
                exon$pct_simul_sign_multi<-sum(simuls_results<0.05)/length(simuls_results)
                simuls_eachexons[[s]]<-exon
                
        }
        results_simuls<-do.call(args=simuls_eachexons,what=rbind.data.frame)
        results_simuls
}


### This function calculates results for real and simulated exons, for Chi-square and ORFscore


take_simuls_chisq_ORFscore<-function(x,nsimul,cutoff_ORFscore=quantile85_ORFscore){
        unique_ex_id<-x[,"exon_id"]
        list_exons_tracks<-list()
        for(i in seq(1:length(unique_ex_id))){
                list_exons_tracks[[i]]<-all_tracks[index==unique_ex_id[i]]
        }
        simuls_eachexons<-list()
        for(s in 1:length(unique_ex_id)){
                withsep<-strsplit(list_exons_tracks[[s]],split=" ")
                x<-t(data.frame(withsep))
                id<-unique_ex_id[s]
                exon<-data.frame(exon_id=id,stringsAsFactors=F,row.names=NULL)
                strand<-x[1,2]
                tracks_pre<-t(x[,-c(1:2)])
                if(strand=="-"){
                        tracks<-cbind(rev(tracks_pre[,1]),rev(tracks_pre[,2]),rev(tracks_pre[,3]),rev(tracks_pre[,4]))
                } else if (strand=="+"){
                        tracks<-tracks_pre}
                colnames(tracks)<-c("Psites","RiboCov","RNACov","RNAcent")
                mode(tracks)<-"numeric"
                length<-dim(tracks)[1]
                
                
                ribo_covered_pos<-which(tracks[,2]>0)
                P_sites_sum<-sum(tracks[,1])
                Phase_P_sites_frame<-sum(tracks[seq(1,length,by=3),1])
                Phase_P_sites_frame_1<-sum(tracks[seq(2,length,by=3),1])
                Phase_P_sites_frame_2<-sum(tracks[seq(3,length,by=3),1])
                score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                exon$ORF_score<-log2(score1+score2+score3+1)
                if(P_sites_sum>15){
                        exon$chisq<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
                if(P_sites_sum<16 & P_sites_sum>0){
                        exon$chisq<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
                }        
                simuls_results<-foreach(j=1:nsimul,.combine=rbind,.multicombine=T) %dopar%{       
                        set.seed(j)
                        simtrack<-rep(0,length)
                        rand_pos<-sample(ribo_covered_pos,P_sites_sum,replace=T)
                        for(i in rand_pos){
                                simtrack[i]<-simtrack[i]+1
                        }
                        
                        Phase_P_sites_frame<-sum(simtrack[seq(1,length,by=3)])
                        Phase_P_sites_frame_1<-sum(simtrack[seq(2,length,by=3)])
                        Phase_P_sites_frame_2<-sum(simtrack[seq(3,length,by=3)])
                        
                        score1<-((Phase_P_sites_frame-P_sites_sum/3)^2)/(P_sites_sum/3)
                        score2<-((Phase_P_sites_frame_1-P_sites_sum/3)^2)/(P_sites_sum/3)
                        score3<-((Phase_P_sites_frame_2-P_sites_sum/3)^2)/(P_sites_sum/3)
                        simul_ORF_score<-log2(score1+score2+score3+1)
                        if(P_sites_sum>15){
                                simul_Chisq<-chisq.test(as.table(c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2)))$p.value}
                        if(P_sites_sum<16 & P_sites_sum>0){
                                simul_Chisq<-xmulti(obs=c(Phase_P_sites_frame,Phase_P_sites_frame_1,Phase_P_sites_frame_2),expr=c(1,1,1),statName="Prob",detail=0)$pProb
                        }        
                        return(c(simul_Chisq,simul_ORF_score))
                }
                colnames(simuls_results)<-c("simul_Chisq","simul_ORF_score")
                exon$n_simul_sign_Chiq<-sum(simuls_results[,"simul_Chisq"]<0.05)
                exon$n_simul_sign_ORFscore<-sum(simuls_results[,"simul_ORF_score"]>cutoff_ORFscore)
                exon$pct_simul_sign_Chiq<-sum(simuls_results[,"simul_Chisq"]<0.05)/dim(simuls_results)[1]
                exon$pct_simul_sign_ORFscore<-sum(simuls_results[,"simul_ORF_score"]>6)/dim(simuls_results)[1]
                simuls_eachexons[[s]]<-exon
                
        }
        results_simuls<-do.call(args=simuls_eachexons,what=rbind.data.frame)
        results_simuls
}


# Multiple plot function, from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        require(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}
