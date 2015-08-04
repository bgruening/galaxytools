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
defaultBins = seq(50, 50*15, by = 50) # default binsize search space, used when fragment length is equal or less than read length 
#=======================> DONE! 




# ================================ 
# Required Libraries check & load
# ================================ 
if ((is.element('parallel', installed.packages()[,1]) == FALSE)) {
	stop("R package'parallel' is required. Please install it!")
}
suppressPackageStartupMessages(library("parallel"))
#=======================> DONE! 




# ================= 
# Custom Functions
# =================
#Implements the Shimazaki procedure
shimazaki = function(bedfile, rl, bins, maxIter, filelist, chromSize, type) {
		
	#read in the file data
	if (type == "single") {
		reads = read.table(bedfile, header = FALSE)
		reads = cbind(as.character(reads[[1]]), as.character(reads[[2]]))
	}
	if (type == "paired") {
		#reads = sqldf("SELECT V1 As 'start', V2 As 'end' FROM readsFile", file.format=list(header = FALSE, sep = "\t"))
		reads = read.table(bedfile, header = FALSE)
		reads = cbind(as.numeric(reads[[1]]), as.numeric(reads[[2]]))
	}
	readnum = length(reads[[1]])
	o = which(filelist == bedfile)
	readlen = rl[o]
	jack = o - 1
	bins = bins[(1+(jack*15)):(15+(jack*15))]

	costs = vector(mode = "numeric", length = length(bins))
	#Shimazaki procedure
	for (i in 1:length(bins)) {
		#construct the counting breaks vector
		genomevec = seq(0, chromSize, by = bins[i]);
		if (max(genomevec) < chromSize) {
			genomevec = append(genomevec, chromSize);
		}

		#create a vector of read counts
		if (type == "single") {
			ameirah = sort(c((as.numeric(reads[reads[,2] == "+",,drop = FALSE][,1])), ((as.numeric(reads[reads[,2] == "-",,drop = FALSE][,1])) + readlen - 1)))
		}
		if (type == "paired") {
			ameirah = sort(c((reads[[1]]), (reads[[2]])))
		}
		ameirah = hist(ameirah, breaks = genomevec, plot = FALSE)
		ameirah = ameirah$counts
	
		#get cost function
		m = mean(ameirah)
		v = (sum((ameirah - m)^2)) / (length(ameirah))
		num = ((2*m) - v)
		den = ((bins[i]) * readnum)^2
		cost = -(log(abs(num)) - log(den))
		costs[i] = cost
	}

index = which.min(costs)
finbin = bins[index]
	
return(finbin)
}
#=======================> DONE! 




# ========================== 
# Parse-in System Variables
# ==========================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Set arguments to default values
ibed = NA # input bed file
sFile = NA # chromosome size
storeFile = NA # file to store result
cornum = 1 # number of processors to use
rl = NA # read length
frags = NA # fragment lengths
bins = NA

#Parsing arguments and storing values
for (each.arg in args) {
	if (grepl('^-ibed=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
				ibed <- arg.split[2]
		} else {
			stop('No input bed file')
		} 
	}
	if (grepl('^-s=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			sFile <- arg.split[2] 
		} else {
			stop('No chromosome size file')
		}
	}
	if (grepl('^-rl=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			rl <- arg.split[2] 
		}
	}
	if (grepl('^-d=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			storeFile <- arg.split[2] 
		} else {
			stop('No file to store result')
		}
	}
	if (grepl('^-p=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			cornum <- as.numeric(arg.split[2]) 
		} else {
			stop('No number of cores given')
		}
	}
	if (grepl('^-f=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			frags <- arg.split[2] 
		} else {
			stop('No Fragment lengths given')
		}
	}
	if (grepl('^-type=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] 
		if (! is.na(arg.split[2]) ) {
			type <- arg.split[2] 
		} else {
			stop('No type given')
		}
	}
}

#Read in variables
chromosomes = read.table(sFile, header=FALSE)
chromSize = as.numeric(chromosomes$V2) #chromosome size
chromSize = max(chromSize) #get maximum chrom size
rm(chromosomes)

ibed = strsplit(ibed, ",", fixed = TRUE)[[1]]
rl = as.numeric(strsplit(rl, ",", fixed = TRUE)[[1]])
frags = as.numeric(strsplit(frags, ",", fixed = TRUE)[[1]])
#=======================> DONE! 



# ===================================================
# Shimazaki Procedure (Shimazaki and Shinomoto 2007)
# ===================================================
for (i in 1:length(ibed)) {
	if (frags[i] > rl[i]) {
		minbin = floor(frags[i] / 2)
		bins = c(bins, seq(minbin, minbin*15, by = minbin)) 
	} else {
		bins = c(bins, defaultBins)
	}
}
bins = bins[!is.na(bins)]
bins = mclapply(ibed, shimazaki, rl, bins, maxIter, ibed, chromSize, type = type, mc.cores = cornum)
bins = min(unlist(bins))
#=======================> DONE! 


# ================== 
# Write Information
# ==================
write(paste0(bins), file = paste0(storeFile, "/binsize.txt"))
message(bins)
#=======================> DONE!
