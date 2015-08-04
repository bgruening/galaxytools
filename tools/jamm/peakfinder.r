########################################################################
# JAMMv1.0.7rev1 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2015  Mahmoud Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: mahmoud.ibrahim@mdc-berlin.de
########################################################################




# ======================= 
# User-defined variables
# ======================= 
samplingSeed = 1011414 #the seed makes the random sampling that JAMM does "deterministic". You can change this to your integer seed of choice to reproduce exact same results or randomize it this way: ceiling(rnorm(1,50,300)). 
reportNoClust = "n" #report windows for which clustering failed? Reported windows will be marked by "NoClust" in the peak name (4th column). 
cutoff = NA #To enforce an SNR cutoff ratio for bin enrichment calling, delete NA and enter the number you want.
strict = 1 #To make bin enrichment calling more / less strict, increase / decrease this number.
meanAdjust = "n" #Adjust the initialization mean vector for each window before clustering? If you want to do so, change this to "y".
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000).
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('mclust', installed.packages()[,1]) == FALSE) || (is.element('signal', installed.packages()[,1]) == FALSE) || (is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package 'mclust', 'signal' and 'parallel' are required. Please install them!")
}
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("signal"))
suppressPackageStartupMessages(library("parallel"))
#=======================> DONE! 





# ================= 
# Custom Functions
# =================
#Get per-row Geometric mean (takes list, returns vectors, not lists!)
geomeanL <- function(mat){
	n = length(mat)
	if (n > 1) {
		mult = (mat[[1]])*(mat[[2]])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[[i]])
			}
		}
		mult = mult^(1/n)
		mat = mult
		return(mat)
	} else {
		return(mat[[1]])
	}
}



#Get per-row Geometric mean (takes matrix, returns vectors, not lists!)
geomean <- function(mat){
	n = NCOL(mat)
	if (n > 1) {
		mult = (mat[,1])*(mat[,2])
		if (n > 2) {
			for (i in 3:n) {
				mult = mult*(mat[,i])
			}
		}
	mult = mult^(1/n)
	mat = mult
	}
return(mat)
}



#Read in bed(st2) file
parsein = function(bedfile) {
	l = read.table(bedfile, header = FALSE)[[1]]
	l = l + 1
	return(l)
}


#Read in bedpe(st2) file
parseinpe = function(bedfile) {
	l = read.table(bedfile, header = FALSE)
	l = cbind((l[[1]] + 1), l[[2]])
	return(l)
}


#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreads = function(bedfile, reads, frag, chromsize, filelist, chrcount) {
	
	o = which(filelist == bedfile)
	
	counts = vector(mode = "numeric", length = chromsize)

	for (j in 1:length(reads[[o]])) {
		if ((reads[[o]][j]+frag[o]-1) <= chromsize) {
			counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] =  counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] + 1
		}
	}
	
	mCount = mean(counts)
	
	if (chrcount == 1) {
		counts = counts/mCount
		write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
	} else {
		meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
		if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
			mCount = meanCounts
		} else {
			write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
		}
		counts = counts/mCount
	}

	return(counts)
}



#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreadspe = function(bedfile, reads, chromsize, filelist, chrcount) {
	
	o = which(filelist == bedfile)
	
	counts = vector(mode = "numeric", length = chromsize)

	for (j in 1:length(reads[[o]][,1])) {
		counts[(reads[[o]][j,1]):(reads[[o]][j,2])] = counts[(reads[[o]][j,1]):(reads[[o]][j,2])] + 1
	}
	
	mCount = mean(counts)
	
	if (chrcount == 1) {
		counts = counts/mCount
		write(paste(mCount), file = paste0(out, "/norma.", o, ".info"))
	} else {
		meanCounts = mean(as.numeric(read.table(paste0(out, "/norma.", o, ".info"))[[1]]))
		if ((mCount >  (5*meanCounts)) || (mCount <  (0.2*meanCounts))) {
			mCount = meanCounts
		} else {
			write(paste(mCount), file = paste0(out, "/norma.", o, ".info"), append = TRUE)
		}
		counts = counts/mCount
	}

	return(counts)
}



#find enriched bins
pickbins = function(winStart, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd) {

	if ((winStart + binSize) <= chromSize) {
		winEnd = winStart + binSize
	} else {
		winEnd = chromSize
	}
	binSizeTemp = winEnd - winStart
	tempend = winEnd - 1

	#extract subset of the background
	if (bkgd != "None") {
		Cs = counts[[numdup+1]][winStart:tempend]
		mCs = mean(Cs)
		dCs = sd(Cs)
	}

	go = rep(0, numdup)
	for (g in 1:numdup) {
		mS = (mean(counts[[g]][winStart:tempend])) 
		ratio = mS/dCs
		if ((mS > (mCs * strict)) && (ratio > cutoff)) {
			go[g] = 1
		}
	}
	veep = sum(go)
	return(veep)
}




#find enriched wins
pickwins = function(winStart, coffeeshopSud, counts, numdup, startlist, winSize) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	rWinSize = winEnd - winStart + 1


	if(rWinSize >= winSize) {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (mean(counts[[g]][winStart:winEnd])) 
		}
		veep = mean(mS)
	} else {
		veep = FALSE
	}
	
	return(veep)
}



#score windows for fast analysis
scorewindow = function(winStart, coffeeshopSud, numdup, C, bkgd, counts, startlist) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#will store peak information
	writethis = list()

	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	for (j in 1:numdup) {
		Rsr[,j] = counts[[j]][winStart:winEnd]
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
	}
	#extract subset of the background
	if (bkgd != "None") {
		Cs = counts[[numdup+1]]
		Cmin = min(Cs[Cs > 0])
		Cs = Cs[winStart:winEnd]
		Cs = filtfilt(rep(1,80)/80,1,Cs) + Cmin #gets rid of Inf in the fold change
	} else {
		set.seed(samplingSeed)
		Cs = sample(C, rWinSizeTemp, replace = TRUE)
		Cs = filtfilt(rep(1,80)/80,1,Cs)
	}
			
	#start scoring
	signal = (geomean(Rs))
	cairo = (mean(signal)) / (mean(Cs))
	return(cairo)
}


#Initialize MClust clustering parameters
smoothcounts = function(winStart, coffeeshopSud, numdup, counts, startlist) { #helper function1

	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	#extract subset of the IP
	Rs = matrix(0, nrow = (winEnd - winStart + 1), ncol = numdup)
	for (j in 1:numdup) {
		Rs[,j] = counts[[j]][winStart:winEnd]
	}
	#smooth extended read counts
	for (j in 1:numdup) {
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rs[,j])
	}	
	return(Rs)
}
cluster = function(model, sig, init, clustnummer, noise) { #helper function2
	set.seed(samplingSeed)
	noisy = sample(noise, length(sig[,1]), replace = TRUE)
	clust = me(model, sig+noisy, init)
	bicc =  bic(model, clust$loglik, length(sig[,1]), length(sig[1,]), clustnummer)
	out = list(bicc = bicc, param = clust$parameters)
	return(out)
}
initparam = function(coffeeshopNord, coffeeshopSud, numdup, counts, cornum, clustnummer, modelnames, noise) { #main function
	
	n = length(coffeeshopNord)
	#smooth extended read counts
	if (cornum > 1) {
		sig = mclapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, counts, startlist = coffeeshopNord, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		sig = lapply(coffeeshopNord, smoothcounts, coffeeshopSud, numdup, counts, startlist = coffeeshopNord)
	}
	sig = do.call(rbind, sig) 

	#kmeans initialization
	set.seed(samplingSeed)
	init = kmeans(sig, clustnummer, nstart = 20)
	init = unmap(init$cluster)
	

	if (cornum > 1) {
		param = mclapply(modelnames, cluster, sig, init, clustnummer, noise, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		param = lapply(modelnames, cluster, sig, init, clustnummer, noise)
	}
	
	bicc = vector(mode = "numeric", length = length(modelnames))
	for (i in 1:length(modelnames)) {
		bicc[i] = as.numeric(param[[i]]$bicc)
	}
	bicc = which.max(bicc)
	
	out = list(initparam = param[[bicc]]$param, modelname = modelnames[bicc])
	return(out) 
}




#find peaks
findpeak = function(winStart, coffeeshopSud, numdup, C, param, bkgd, resol, counts, noise, startlist, meanAdjust, clustnummer) {
	
		
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	
	
	#will store peak information
	writethis = list()
	ccx = 1 #default is clustering didNOT work

	
	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	Rsr = Rs
	if (meanAdjust == "y") {
		for (j in 1:numdup) {
			Rsr[,j] = counts[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
			kabel = which.max(param$init$mean[j,])
			param$initparam$mean[j,kabel] = mean(Rs[,j])
		}
	} else {
		for (j in 1:numdup) {
			Rsr[,j] = counts[[j]][winStart:winEnd]
			Rs[,j] = filtfilt(rep(1,80)/80,1,Rsr[,j])
		}
	}
	
	if (resol != "window") {
		
		#clustering (take 1)
		take = 1
		set.seed(samplingSeed)
		noisy = sample(noise, rWinSizeTemp, replace = TRUE)
		clust = em(param$modelname, Rs+noisy, param$initparam)
		clust$classification = map(clust$z)
		if (!((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification)))))) { #clustering didn't work, take1
			
			#repeat clustering from scratch, take 2!
			set.seed(samplingSeed)
			init = kmeans(Rs, clustnummer, nstart = 20)
			init = unmap(init$cluster)
			set.seed(samplingSeed)
			noisy = sample(noise, rWinSizeTemp, replace = TRUE)
			clust = me(param$modelname, Rs+noisy, init)
			clust$classification = map(clust$z)
			if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) {
				ccx = 0 #clustering worked, take2
				take = 2
			}
		} else {ccx = 0}  #clustering worked, take1
		 
		if (ccx != 1) { #clustering worked either in take1 or take2
					
			if (numdup > 1) { #check whether all components replicates agreed on clustering assignments
				cc = vector(mode = "numeric", length = numdup)
				for (g in 1:numdup) {
					cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
				}
				ccx = sum(diff(cc))
				cluster = cc[1]
				rm(cc)
			
				if ((ccx != 0) && (take == 1)) { #not all replicates agreed? Repeat the clustering with from scratch if not already done!
					set.seed(samplingSeed)
					init = kmeans(Rs, clustnummer, nstart = 20)
					init = unmap(init$cluster)
					set.seed(samplingSeed)
					noisy = sample(noise, rWinSizeTemp, replace = TRUE)
					clust = me(param$modelname, Rs+noisy, init)
					clust$classification = map(clust$z)
					if ((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) { #clustering worked? check whether replicates agreed take 3
						cc = vector(mode = "numeric", length = numdup)
						for (g in 1:numdup) {
							cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
						}
						ccx = sum(diff(cc))
						cluster = cc[1]
						rm(cc)
						take = 3
					}
				}
			} else { #no replicates!
				cluster = which.max(clust$parameters$mean) #which cluster has the largest mean (this is the peak cluster, hopefully!)
			}
		}
		
		if ((ccx != 0) && (reportNoClust=="y")) { resol = "window" } #clustering did not work and windows should be reported
			
		
		if (ccx == 0) { #clustering worked and all replicates agree on the cluster assignments
			
			#extract subset of the background
			if (bkgd != "None") {
				Cs = counts[[numdup+1]][winStart:winEnd]
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			} else {
				set.seed(samplingSeed)
				Cs = sample(C, rWinSizeTemp, replace = TRUE)
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			}
			
			#find region boundaries
			loc = 1:length(clust$classification)
			gmclass = cbind(loc, clust$classification)
			locPeak = gmclass[gmclass[,2] == cluster,,drop=FALSE]
			rStart = locPeak[1] #start position of the region
			rEnd = locPeak[length(locPeak[,1]),1] #end position of the region
		
			#peak resolution check
			if (resol == "region") {
				pSize = rEnd - rStart
				signal = (geomean(Rs[rStart:rEnd,])) 
				signal2 = (signal) - (Cs[rStart:rEnd])
				gm = mean(signal2)
				summit = which.max(geomean(Rsr[rStart:rEnd,])) - 1
				will2k = wilcox.test(signal, Cs[rStart:rEnd])
				
				#Is there signal in the region above background
				if (gm > 0) {
					writethis[[1]] = rStart + winStart - 1
					writethis[[2]] = rEnd + winStart
					writethis[[3]] = paste0(chromName, ".", rStart+winStart -1)
					writethis[[4]] = "1000"
					writethis[[5]] = "."
					writethis[[6]] = gm
					writethis[[7]] = will2k$p.value
					writethis[[8]] = "-1"
					writethis[[9]] = summit
				}
			} else if (resol == "peak") {
				#find out where separate peaks are
				d = diff(locPeak[,1]) 
				d[length(d)+1] = 0		
				locPeak = cbind(locPeak, d)
				bound1 = which(locPeak[,3] > 1, arr.in=TRUE)
				bound2 = bound1 + 1
				bound = locPeak[sort(c(bound1,bound2))]
				bound = c(rStart, bound, rEnd)
				w = 1
				warum = 0
				while (w < length(bound)) {
					pStart = bound[w] + winStart - 1 
					pEnd = bound[w+1] + winStart
					pSize = pEnd - pStart
					signal = (geomean(Rs[(bound[w]):(bound[w+1]),])) 
					signal2 = (signal) - (Cs[bound[w]:bound[w+1]])
					gm = mean(signal2)
					summit = which.max(geomean(Rsr[(bound[w]):(bound[w+1]),])) - 1
					will2k = wilcox.test(signal, Cs[(bound[w]):(bound[w+1])])
										
					weil = warum * 9
					#Is there signal in the region above background
					if (gm > 0) {
						writethis[[1+weil]] = pStart
						writethis[[2+weil]] = pEnd
						writethis[[3+weil]] = paste0(chromName, ".", pStart)
						writethis[[4+weil]] = "1000"
						writethis[[5+weil]] = "."
						writethis[[6+weil]] = gm
						writethis[[7+weil]] = will2k$p.value
						writethis[[8+weil]] = "-1"
						writethis[[9+weil]] = summit
					}
					w = w + 2
					warum = warum + 1			
				}
			} #peak resolution check
		} #clustering worked and all replicates agree on clustering assignments?
	} #window resolution check
	
	if (resol == "window") {
	
		#extract subset of the background
		if (bkgd != "None") {
			Cs = counts[[numdup+1]][winStart:winEnd]
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		} else {
			set.seed(samplingSeed)
			Cs = sample(C, rWinSizeTemp, replace = TRUE)
			Cs = filtfilt(rep(1,80)/80,1,Cs)
		}
		
		#calculate scores
		pSize = rWinSizeTemp
		signal = geomean(Rs) 
		signal2 = (signal) - (Cs)
		gm = mean(signal2)
		summit = which.max(geomean(Rsr)) - 1
		will2k = wilcox.test(signal, Cs)
		
		#Is there signal in the region above background
		if (gm > 0) {
			writethis[[1]] = winStart - 1
			writethis[[2]] = winEnd
			writethis[[3]] = paste0(chromName, ".", winStart -1, ".NoClust")
			writethis[[4]] = "1000"
			writethis[[5]] = "."
			writethis[[6]] = gm
			writethis[[7]] = will2k$p.value
			writethis[[8]] = "-1"
			writethis[[9]] = summit
		}
	} #window reporting

return(writethis)
}


#filter return value of findpeak()
processPeaks = function(peaks) {
	peaks = matrix(unlist(peaks), ncol=9, byrow=TRUE)
	peaks = peaks[peaks[,1] != FALSE,,drop=FALSE]
	peaks = data.frame(peaks)
	return(peaks)
}
#=======================> DONE!









# ========================== 
# Parse-in System Variables
# ==========================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Parsing arguments and storing values
for (each.arg in args) {
	#chormosome size file
	if (grepl('-sfile=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				size.file <- arg.split[2]
		} else {
			stop('No genome size file')
		} 
	}
	#bed file names
	if (grepl('-bednames=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bednames <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
	#bed files directory
	if (grepl('-frag=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				frag <- arg.split[2]
		} else {
			stop('No fragment length given')
		} 
	}	
	#bakcground files directory
	if (grepl('-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			stop('No background file')
		} 
	}	
	#bakcground files directory
	if (grepl('-out=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				out <- arg.split[2]
		} else {
			stop('No output directory given')
		} 
	}
	#Cluster number
	if (grepl('-clustnummer=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				clustnummer <- as.numeric(arg.split[2])
		} 
	}
	#resolution
	if (grepl('-resolution=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				resol <- arg.split[2]
		} 
	}
	#processor cores
	if (grepl('-p=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				cornum <- as.numeric(arg.split[2])
		} 
	}
	#minimum window size
	if (grepl('-window=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			winSize <- arg.split[2]
		} 
	}
	#window size
	if (grepl('-bin=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			binsize <- arg.split[2]
		} 
	}
	#type (paired / single)
	if (grepl('-type=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			type <- arg.split[2]
		} 
	}
	#chromosome number
	if (grepl('-chrcount=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			chrcount <- as.numeric(arg.split[2])
		} 
	}
	#window enrichment cutoff
	if (grepl('-windowe=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			windowe <- arg.split[2]
		} 
	}
	#initialize
	if (grepl('-initModel=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			initialize <- arg.split[2]
		} 
	}
	
}

##Parse in variables
chromosomes = read.table(size.file, header=FALSE)
chromName = chromosomes$V1; #which chromosome
chromSize = chromosomes$V2; #chromosomes size
rm(chromosomes)

if (chrcount == 1) {
	write(paste(samplingSeed), file = paste0(out, "/seed.info"))
} else {
	samplingSeed = as.numeric(read.table(paste0(out, "/seed.info"), header = FALSE)[[1]])
}
readsFiles = as.list(strsplit(bednames, ",", fixed = TRUE)[[1]])
numdup = length(readsFiles) #number of replicates
if (bkgd != "None") {
	readsFiles[[numdup+1]] = bkgd
}
winSize = as.numeric(winSize)
binSize = as.numeric(binsize)
winSize = binSize * winSize

if (type == "single") {
	frags = as.numeric(strsplit(frag, ",", fixed = TRUE)[[1]])
}
rm(bednames)

if (numdup > 1) {
	modelnames = c("VVV","VEV")
} else {
	modelnames = "V"
}

if (windowe != "auto") {
	windowe = as.numeric(windowe)
}

nothing = FALSE #default is we found some peaks
options(stringsAsFactors = FALSE)
#=======================> DONE!





# ======================= 
# Some preliminary stuff
# =======================
if (type == "single") {
	if (cornum > 1) {
		datain = mclapply(readsFiles, parsein, mc.cores = cornum, mc.preschedule = TRUE) #read in all bed files (samples and control)
	} else {
		datain = lapply(readsFiles, parsein) #read in all bed files (samples and control)
	}
}

if (type == "paired") {
	if (cornum > 1) {
		datain = mclapply(readsFiles, parseinpe, mc.cores = cornum, mc.preschedule = TRUE) #read in all bed files (samples and control)
	} else {
		datain = lapply(readsFiles, parseinpe) #read in all bed files (samples and control)
	}
}

#minimum peak size (only a recommendation)
minpeak = floor(binSize / 4)

#make bins vector
bins = seq(from = 1, to = (chromSize - 1), by = binSize)
#=======================> DONE!




# =============== 
# Counting Reads
# ===============
if (type == "single") {
	if (cornum > 1) {
		counts = mclapply(readsFiles, countreads, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		counts = lapply(readsFiles, countreads, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount)
	}
}

if (type == "paired") {
	if (cornum > 1) {
		counts = mclapply(readsFiles, countreadspe, reads = datain, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		counts = lapply(readsFiles, countreadspe, reads = datain, chromsize = chromSize, filelist = readsFiles, chrcount = chrcount)
	}
}

rm(datain)
#=======================> DONE!




# ============================ 
# Estimating Background Model
# ============================ 
if (chrcount == 1){ #first chromosome, estimate bkgd (includes SNR cutoff)
	if (is.na(cutoff)) {
		if (bkgd != "None") {
			cutoff = vector(length = numdup)
			sdC = sd(counts[[numdup+1]])
			for (x in 1:numdup) {
				cutoff[x] = (mean(counts[[x]]))/(sdC)
			}
			cutoff = max(cutoff)
			C = NULL
			mCs = NULL
			write(paste(c(cutoff,NA,NA)), file = paste0(out, "/bkgd.info"), append = TRUE)
		} else {
			cutoff = vector(length = numdup)
			mmV = var(geomeanL(counts))
			mmM = mean(geomeanL(counts))
			sigma = log(1+((mmV) / ((mmM)^2)))
			mu = (log(mmM)) - (0.5 * (sigma))
			set.seed(samplingSeed)
			C = rlnorm(100000, mu, sqrt(sigma))
			for (x in 1:numdup) {
				cutoff[x] = (mean(counts[[x]]))/(sd(C))
			}
			cutoff = max(cutoff)
			set.seed(samplingSeed)
			snow = sample(C, binSize*5, replace = TRUE)
			mCs = mean(snow)
			dCs = sd(snow)
			write(paste(c(cutoff,sigma,mu)), file = paste0(out, "/bkgd.info"), append = TRUE)
		}
	}
} else { #bkgd estiamted from before
	bkgdInfo = read.table(paste0(out, "/bkgd.info"), header = FALSE)
	if (is.na(cutoff)) {
		cutoff = as.numeric(bkgdInfo[[1]][1])
	}

	if (bkgd != "None") {
		C = NULL
		mCs = NULL
	} else {
		sigma = as.numeric(bkgdInfo[[1]][2])
		mu = as.numeric(bkgdInfo[[1]][3])
		set.seed(samplingSeed)
		C = rlnorm(100000, mu, sqrt(sigma))
		set.seed(samplingSeed)
		snow = sample(C, binSize*5, replace = TRUE)
		mCs = mean(snow)
		dCs = sd(snow)
	}
}
#=======================> DONE!




# ======================== 
# Picking Enriched Windows
# ========================
if (cornum > 1) {
	coffeeshop = mclapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd, mc.cores = cornum, mc.preschedule = TRUE)
} else {
	coffeeshop = lapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, strict, mCs, dCs, bkgd)
}
coffeeshop = as.numeric(unlist(coffeeshop))
coffeeshop[coffeeshop != numdup] = 0

if (sum(coffeeshop) != 0) { #Any enriched bins?

coffeeshop = c(0, diff(coffeeshop))
coffeeshop = cbind(coffeeshop, bins)
coffeeshopNord = coffeeshop[coffeeshop[,1] == numdup,,drop=FALSE]
coffeeshopSud = coffeeshop[coffeeshop[,1] == -numdup,,drop=FALSE]
coffeeshopNord = coffeeshopNord[,2]
coffeeshopSud = coffeeshopSud[,2] - 1
if (length(coffeeshopSud) < length(coffeeshopNord)) {
	coffeeshopSud = c(coffeeshopSud, chromSize) 
} else if (length(coffeeshopSud) > length(coffeeshopNord)) {
	coffeeshopNord = c(1, coffeeshopNord)
}
if (coffeeshopSud[length(coffeeshopSud)] > chromSize) {
	coffeeshopSud[length(coffeeshopSud)] = chromSize
}


if (cornum > 1) {
	coffeeshop = mclapply(coffeeshopNord, pickwins, coffeeshopSud, counts, numdup, startlist = coffeeshopNord, winSize, mc.cores = cornum, mc.preschedule = TRUE)
} else {
	coffeeshop = lapply(coffeeshopNord, pickwins, coffeeshopSud, counts, numdup, startlist = coffeeshopNord, winSize)
}
coffeeshop = as.numeric(unlist(coffeeshop))
coffeeshop = cbind(coffeeshopNord, coffeeshopSud, coffeeshop)
coffeeshop = coffeeshop[coffeeshop[,3] != FALSE,,drop=FALSE]
coffeeshop = coffeeshop[order(coffeeshop[,3], decreasing = TRUE),]
rm(bins)
#=======================> DONE!




# =================================== 
# Initializing Clustering Parameters
# ===================================
if (length(coffeeshop[,1]) > 0) { #any enriched windows detected?

if (initialize == "deterministic") {
	yummy = ceiling(length(coffeeshop[,1]) / 1000)
	if (yummy == 0) {
		yummy = 1
	}
}
if (initialize == "stochastic") {
	yummy = ceiling(length(coffeeshop[,1]) / 4)
	if (yummy > 20) {
		set.seed(samplingSeed)
		yummy = sample(1:yummy, 20)
	} else if (yummy > 0) {
		yummy = 1:yummy
	} else {
		yummy = 1
	}
}
coffeeshopNord = coffeeshop[yummy,1]
coffeeshopSud = coffeeshop[yummy,2]
set.seed(samplingSeed)
noise = rnorm(100000, mean=0, sd=0.1)
param = initparam(coffeeshopNord, coffeeshopSud, numdup, counts, cornum, clustnummer, modelnames, noise)
#=======================> DONE!





# ========================== 
# Enriched Window Filtering
# ==========================
if (windowe != 1) { #do it only if window fold enrichment filtering is required
	if (cornum > 1) {
		scores = mclapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgd, counts, startlist = coffeeshop[,1], mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		scores = lapply(coffeeshop[,1], scorewindow, coffeeshop[,2], numdup, C, bkgd, counts, startlist = coffeeshop[,1])
	}
	scores = unlist(scores)

	if (windowe == "auto") {
		lscores = log(scores)
		if (length(scores) > 0) {
			if (chrcount == 1) {
				cutthisTEMP = ((mean(lscores)) + (sd(lscores)*1))
				write(paste(cutthisTEMP), file = paste0(out, "/bkgd.info"), append = TRUE)
			} else {
				cutthisTEMP = as.numeric(bkgdInfo[[1]][4])
			}
			finalwins = which(lscores > cutthisTEMP)
			cutthisW = min(scores[finalwins])
			coffeeshop = coffeeshop[finalwins,]
		} else {
			if (chrcount == 1) {
				cutthisTEMP = 0
				cutthisW = "Not Applicable, All Windows Analyzed!"
				write(paste(cutthisTEMP), file = paste0(out, "/bkgd.info"), append = TRUE)
			}
		}
	} else {
		cutthisW = windowe
		if (length(scores) > 0) {
			finalwins = which(scores >= windowe)
			coffeeshop = coffeeshop[finalwins,]
		}
	}
}
else { cutthisW = 1 }
#=======================> DONE!





# ============== 
# Finding Peaks
# ==============
if (length(coffeeshop[,1]) > 0) { #any enriched windows left after filtering?
coffeeshop = coffeeshop[,-3]
if (cornum > 1) {
	peaks = mclapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgd, resol, counts, noise, startlist = coffeeshop[,1], meanAdjust, clustnummer, mc.cores = cornum, mc.preschedule = TRUE)
} else {
	peaks = lapply(coffeeshop[,1], findpeak, coffeeshop[,2], numdup, C, param, bkgd, resol, counts, noise, startlist = coffeeshop[,1], meanAdjust, clustnummer)
}
if (!(is.null(peaks))) { #any peaks discovered?
writethis = processPeaks(peaks)
#=======================> DONE!





# =========================
# Writing Peak Information
# =========================
} else { nothing = TRUE } #no peaks
} else { nothing = TRUE } #no enriched windows left after filtering
} else { nothing = TRUE } #no enriched widnows discovered
} else { nothing = TRUE; cutthisW = windowe } #no enriched bins discovered

if (isTRUE(nothing)) {
	file.create(paste0(out, "/", chromName, ".peaks.bed"))
	write(paste(chromName, minpeak, sep = "	"), file = paste0(out, "/min.peaksize"), append=TRUE)

	if (chrcount == 1) {
		message(paste0("No peaks found! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
	} else {
		message("No peaks found!")
	}
} else {
	write(paste(chromName, writethis$X1, writethis$X2, writethis$X3, writethis$X4, writethis$X5, writethis$X6, writethis$X7, writethis$X8, writethis$X9, minpeak, sep = "	"), file = paste0(out, "/", chromName, ".peaks.bed"), ncolumns = 1)
	write(paste(chromName, minpeak, sep = "	"), file = paste0(out, "/min.peaksize"), append=TRUE)
	
	
	if (chrcount == 1) {
		message(paste0("Done! - Window Fold Enrichment: ", cutthisW, " - Seed: ", samplingSeed))
	} else {
		message("Done!")
	}
}
#=======================> DONE!


