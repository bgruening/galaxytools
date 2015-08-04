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
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000) 
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('signal', installed.packages()[,1]) == FALSE) || (is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R packages 'signal'and 'parallel' are required. Please install them!")
}
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
	reads = read.table(bedfile, header = FALSE)[[1]]
	return(reads)
}

#Read in bedpe(st2)
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




#Produces normalized extended read counts (takes output of parsein(), return a vector of floats / PAIRED END)
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


#generate signal
getsig = function(winStart, winEnds, startList, stepSize, chromName) {
	
	
	plz = which(startList == winStart)
	winEnd = winEnds[plz]
	rWinSizeTemp = winEnd - winStart + 1
	
	#extract subset of the IP
	Rs = matrix(nrow = rWinSizeTemp, ncol = numdup)
	for (j in 1:numdup) {
		Rs[,j] = counts[[j]][winStart:winEnd]
		Rs[,j] = filtfilt(rep(1,80)/80,1,Rs[,j])
	}

	#extract subset of the background and initial signal
	if (bkgd != "None") {
		Cs = counts[[numdup+1]][winStart:winEnd]
		Cs = filtfilt(rep(1,80)/80,1,Cs)
		signal = (geomean(Rs)) - (Cs)
	} else {
		signal = geomean(Rs)
	}
			

	#get binned signal
	signal = round(colMeans(matrix(c(signal, rep(0, stepSize - length(signal) %% stepSize)), stepSize)), digits = 2)
	steps = seq(winStart, winEnd, by = stepSize)
	starts = steps - 1
	ends = steps + stepSize - 1
	if (ends[length(ends)] > winEnd) {
		ends[length(ends)] = winEnd
	}
	
	#write to file
	write(paste(chromName, starts, ends, signal, sep = "	"), file = paste0(out, "/", chromName, plz, ".bedSignal"), ncolumns = 1)

	
return(NULL)
}
#=======================> DONE!





# ========================== 
# Parse-in System Variables
# ==========================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Parsing arguments and storing values
for (each.arg in args) {
	#bed file names
	if (grepl('^-bednames=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bednames <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
	#bed file names
	if (grepl('^-chromo=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				chromName <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
	#bed files directory
	if (grepl('^-frag=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			frag <- arg.split[2]
		}
	}	
	#bakcground files directory
	if (grepl('^-bkgd=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				bkgd <- arg.split[2]
		} else {
			stop('No background file')
		} 
	}	
	#bakcground files directory
	if (grepl('^-out=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				out <- arg.split[2]
		} else {
			stop('No output directory given')
		} 
	}
	#Cluster number
	if (grepl('^-clustnummer=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				clustnummer <- as.numeric(arg.split[2])
		} 
	}
	#regions
	if (grepl('^-regions=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				regions <- arg.split[2]
		} 
	}
	#processor cores
	if (grepl('^-p=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				cornum <- as.numeric(arg.split[2])
		} 
	}
	#chromSize
	if (grepl('^-chromoS=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				chromSize <- as.numeric(arg.split[2])
		} 
	}
	#stepsize
	if (grepl('^-sig=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				sigstep <- as.numeric(arg.split[2])
		} 
	}
	#chrcount
	if (grepl('^-chrcount=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				chrcount <- as.numeric(arg.split[2])
		} 
	}
	#alignment type
	if (grepl('^-t=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				type <- arg.split[2]
		} 
	}
}

##Parse in variables
readsFiles = as.list(strsplit(bednames, ",", fixed = TRUE)[[1]])
numdup = length(readsFiles) #number of replicates
if (bkgd != "None") {
	readsFiles[[numdup+1]] = bkgd
}
if (type == "single") {
	frags = as.numeric(strsplit(frag, ",", fixed = TRUE)[[1]])
}
rm(bednames)

options(stringsAsFactors = FALSE)
#=======================> DONE!




# ======================= 
# Some preliminary stuff
# =======================
#import regions
regions = read.table(regions, header = FALSE)
regions = regions[regions[[1]] == chromName, , drop = FALSE]
regions = cbind(regions[[2]] + 1, regions[[3]])

#import read data
if (type == "single") {
	if (cornum > 1) {
		datain = mclapply(readsFiles, parsein, mc.cores = cornum, mc.preschedule = TRUE) #read in all bed files (samples and control)
	} else {
		datain = lapply(readsFiles, parsein) #read in all bed files (samples and control)
	}
}
#import read data
if (type == "paired") {
	if (cornum > 1) {
		datain = mclapply(readsFiles, parseinpe, mc.cores = cornum, mc.preschedule = TRUE) #read in all bed files (samples and control)
	} else {
		datain = lapply(readsFiles, parseinpe) #read in all bed files (samples and control)
	}
}
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



# ============== 
# Getting signal
# ==============
if (cornum > 1) {
	sig = mclapply(regions[,1], getsig, winEnds = regions[,2], startList = regions[,1], stepSize = sigstep, chromName, mc.cores = cornum, mc.preschedule = TRUE)
} else {
	sig = lapply(regions[,1], getsig, winEnds = regions[,2], startList = regions[,1], stepSize = sigstep, chromName)
}
message("Done!")
#=======================> DONE!
