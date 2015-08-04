###############################
#### Peak finder for NGS Data
#### R script		
###############################




# ======================= 
# User-defined variables
# ======================= 
cutoff = NA #To enforce an SNR cutoff ratio, delete NA and enter the number you want
strict = 1 #to make peak calling more / less strict, increase / decrease this number - Only use this if you are getting an unusually large or small number of peaks.
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000) 
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('mclust', installed.packages()[,1]) == FALSE) || (is.element('signal', installed.packages()[,1]) == FALSE) || (is.element('sqldf', installed.packages()[,1]) == FALSE) || (is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package 'Mclust', 'signal', 'sqldf' and 'parallel' are required. Please install them!")
}
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("signal"))
suppressPackageStartupMessages(library("sqldf"))
suppressPackageStartupMessages(library("tcltk"))
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
	l = sqldf()
	readsFile = file(bedfile)
	return(sqldf("SELECT V1 + 1 As 'start' FROM readsFile", file.format=list(header = FALSE, sep = "\t"))$start)
	l = sqldf()
	close(readsFile)
}


#Read in bedpe(st2) file
parseinpe = function(bedfile) {
	l = sqldf()
	readsFile = file(bedfile)
	sebastian = sqldf("SELECT * FROM readsFile", file.format=list(header = FALSE, sep = "\t"))
	return(cbind(sebastian[[1]]+1, sebastian[[2]]))
	l = sqldf()
	close(readsFile)
}


#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreads = function(bedfile, reads, frag, chromsize, filelist) {
	
	o = which(filelist == bedfile)
	
	counts = vector(mode = "numeric", length = chromsize)

	for (j in 1:length(reads[[o]])) {
		if ((reads[[o]][j]+frag[o]-1) <= chromsize) {
			counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] =  counts[(reads[[o]][j]):(reads[[o]][j]+frag[o]-1)] + 1
		}
	}
	counts = counts/(mean(counts))
	return(counts)
}



#Produces normalized extended read counts (takes output of parsein(), return a vector of floats)
countreadspe = function(bedfile, reads, chromsize, filelist) {
	
	o = which(filelist == bedfile)
	
	counts = vector(mode = "numeric", length = chromsize)

	for (j in 1:length(reads[[o]][,1])) {
		counts[(reads[[o]][j,1]):(reads[[o]][j,2])] = counts[(reads[[o]][j,1]):(reads[[o]][j,2])] + 1
	}
	counts = counts/(mean(counts))
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
	} else {
		Cs = sample(C, binSizeTemp, replace = TRUE)
	}

	
	#find out whether it's enriched
	go = rep(0, numdup)
	for (g in 1:numdup) {
		mS = (mean(counts[[g]][winStart:tempend])) 
		ratio = mS/dCs
		if ((mS > (mCs * strict)) && (ratio > cutoff)) {
			go[g] = 1
		}
	}
	veep = sum(go)

	#get counts when enriched	
	if (veep == numdup) {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (mean(counts[[g]][winStart:winEnd])) 
		}
		cairo = mean(mS)
	} else {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (mean(counts[[g]][winStart:winEnd])) 
		}
		cairo = -(mean(mS))
	}
	return(cairo)
}



#find enriched wins
pickwins = function(winStart, coffeeshopSud, counts, numdup, startlist, winSize) {
	
	plz = which(startlist == winStart)
	winEnd = coffeeshopSud[plz]
	rWinSize = winEnd - winStart + 1


	if(rWinSize >= winSize) {
		mS = rep(0, numdup)
		for (g in 1:numdup) {
			mS[g] = (sum(counts[[g]][winStart:winEnd])) 
		}
		veep = sum(mS)
	} else {
		veep = FALSE
	}
	
	return(veep)
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
	init = kmeans(sig, clustnummer)
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
findpeak = function(winStart, coffeeshopSud, numdup, C, param, bkgd, resol, counts, noise, startlist) {
	
	
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

	
	####CLUSTERING START HERE#####
	noisy = sample(noise, rWinSizeTemp, replace = TRUE)
	clust = em(param$modelname, Rs+noisy, param$initparam)
	clust$classification = map(clust$z)
	####CLUSTERING DONE####
	
	#check whether clutering succeeded
	if((any(diff(clust$classification)) != 0) && (!(any(is.na(clust$classification))))) {
	
		
		#check whether all components replicates agreed on clustering assignments
		ccx = 1
		if (numdup > 1) {
			cc = vector(mode = "numeric", length = numdup)
			for (g in 1:numdup) {
				cc[g] = which.max(clust$parameters$mean[g,]) #which cluster has the largest mean (this is the peak cluster, hopefully!)
			}
			ccx = sum(diff(cc))
			cluster = cc[1]
			rm(cc)
		} else {
			cluster = which.max(clust$parameters$mean) #which cluster has the largest mean (this is the peak cluster, hopefully!)
			ccx = 0
		}
		if (ccx == 0) {
			
			#extract subset of the background
			if (bkgd != "None") {
				Cs = counts[[numdup+1]][winStart:winEnd]
				Cs = filtfilt(rep(1,80)/80,1,Cs)
			} else {
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
		} #all replicates agree on clustering assignments?
	} #clustering worked?
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
	if (grepl('^-sfile=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				size.file <- arg.split[2]
		} else {
			stop('No genome size file')
		} 
	}
	#bed file names
	if (grepl('^-bednames=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				bednames <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
	#bed files directory
	if (grepl('^-frag=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				frag <- arg.split[2]
		} else {
			stop('No fragment length given')
		} 
	}	
	#bakcground files directory
	if (grepl('^-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			stop('No background file')
		} 
	}	
	#bakcground files directory
	if (grepl('^-out=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				out <- arg.split[2]
		} else {
			stop('No output directory given')
		} 
	}
	#Cluster number
	if (grepl('^-clustnummer=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				clustnummer <- as.numeric(arg.split[2])
		} 
	}
	#resolution
	if (grepl('^-resolution=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				resol <- arg.split[2]
		} 
	}
	#processor cores
	if (grepl('^-p=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				cornum <- as.numeric(arg.split[2])
		} 
	}
	#minimum window size
	if (grepl('^-window=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			winSize <- arg.split[2]
		} 
	}
	#window size
	if (grepl('^-bin=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			binsize <- arg.split[2]
		} 
	}
	#type (paired / single)
	if (grepl('^-type=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			type <- arg.split[2]
		} 
	}
}

##Parse in variables
chromosomes = read.table(size.file, header=FALSE)
chromName = chromosomes$V1; #which chromosome
chromSize = chromosomes$V2; #chromosomes size
rm(chromosomes)

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
		counts = mclapply(readsFiles, countreads, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		counts = lapply(readsFiles, countreads, reads = datain, frag = frags, chromsize = chromSize, filelist = readsFiles)
	}
}

if (type == "paired") {
	if (cornum > 1) {
		counts = mclapply(readsFiles, countreadspe, reads = datain, chromsize = chromSize, filelist = readsFiles, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		counts = lapply(readsFiles, countreadspe, reads = datain, chromsize = chromSize, filelist = readsFiles)
	}
}

rm(datain)

#get total counts
mS = vector(mode = "numeric", length = numdup)
for (g in 1:numdup) {
	mS[g] = sum(counts[[g]]) 
}
totalCounts = sum(mS)
rm(mS)
#=======================> DONE!



# =================== 
# Estimating Cutoff
# ===================
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
	} else {
		cutoff = vector(length = numdup)
		mmV = var(geomeanL(counts))
		mmM = mean(geomeanL(counts))
		sigma = log(1+((mmV) / ((mmM)^2)))
		mu = (log(mmM)) - (0.5 * (sigma))
		C = rlnorm(100000, mu, sqrt(sigma))
		for (x in 1:numdup) {
			cutoff[x] = (mean(counts[[x]]))/(sd(C))
		}
		cutoff = max(cutoff)
		mCs = (mean(sample(C, binSize*5, replace = TRUE)))
		dCs = (sd(sample(C, binSize*5, replace = TRUE)))
	}	
}
#=======================> DONE!




# ======================== 
# Picking Enriched Windows
# ========================
binNum = vector(mode = "numeric", length = 41)
binTotalSize = vector(mode = "numeric", length = 41)
binTotalSizeRatio = vector(mode = "numeric", length = 41)
binTotalCount = vector(mode = "numeric", length = 41)
binTotalCountNO = vector(mode = "numeric", length = 41)
binTotalCountRatio = vector(mode = "numeric", length = 41)
maxiFun = vector(mode = "numeric", length = 41)
ggg = seq(1,5, by = 0.1)
for (i in 1:41) {
	if (cornum > 1) {
		coffeeshop = mclapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, ggg[i], mCs, dCs, bkgd, mc.cores = cornum, mc.preschedule = TRUE)
	} else {
		coffeeshop = lapply(bins, pickbins, counts, binSize, chromSize, numdup, C, cutoff, ggg[i], mCs, dCs, bkgd)
	}
	coffeeshop = as.numeric(unlist(coffeeshop))
	coffeeshopYES = coffeeshop[coffeeshop > 0]
	coffeeshopNO = -(coffeeshop[coffeeshop <= 0])
	
	binNum[i] = length(coffeeshopYES)	
	binTotalSize[i] = binNum[i] * binSize
	binTotalSizeRatio[i] = binTotalSize[i] / chromSize
	
	binTotalCount[i] = sum(coffeeshopYES)
	binTotalCountNO[i] = sum(coffeeshopNO)
	binTotalCountRatio[i] = binTotalCount[i] / binTotalCountNO[i]
	
	maxiFun[i] = binTotalCountRatio[i] / binTotalSizeRatio[i]		
}
rm(bins)
#=======================> DONE!


write(paste(ggg, binNum, binTotalSize, binTotalSizeRatio, binTotalCount, binTotalCountNO, binTotalCountRatio, maxiFun, sep = "	"), file = paste0(file = paste0("/home/mibrahim/Documents/", chromName, ".GOYEAH"), ncolumns = 1))
