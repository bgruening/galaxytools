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




# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line


#Parsing arguments and storing values
for (each.arg in args) {
	#bed file names
	if (grepl('^-filelist=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				bednames <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
}
options(stringsAsFactors = FALSE)
#=======================> DONE! 



# ================= 
# Custom Functions
# =================
score <- function(x){
	x  = ((x-min(x))/(max(x)-min(x)))*10000
	return(x)
}

fixnpf = function(bedfile) {

	writethis = read.table(bedfile, header=FALSE)

	
	message(paste0("Minimum Peak Width Used to Produce Filtered List: ", writethis$V11[1]))
	writethis$V9 = p.adjust(writethis$V8, method = "BH")
	writethis$V8[writethis$V8 == 0] = 10
	writethis$V8[writethis$V8 == 10] = min(writethis$V8)
	writethis$V8 = -(log10(writethis$V8))
	writethis$V9[writethis$V9 == 0] = 10
	writethis$V9[writethis$V9 == 10] = min(writethis$V9)
	writethis$V9 = -(log10(writethis$V9))	
	writethis$V7 = as.numeric(writethis$V7) * writethis$V9
	writethis$V5 = score(as.numeric(writethis$V7))
	
	
	#find peak score cutoff (pareto distribution)
	sizes = writethis$V3 - writethis$V2
	scores = writethis$V7
	ss = cbind(sizes, scores)
	ss = ss[ss[,1] > writethis$V11[1],,drop=FALSE]
	scores = ss[,2]
	scores = scores[scores > 0]
	num = length(scores)
	if (length(scores) > 300000) {
		scores = scores[1:300000]
		num = 300000
	}
	den = log(min(scores))
	den = sum(log(scores) - den)
	alpha = num / den	
	geom = (exp(1/alpha)) * (min(scores)) 
	message(paste0("Minimum Peak Score Used to Produce Filtered List: ", geom))
	
	#done now write something
	write(paste(writethis$V1, writethis$V2, writethis$V3, writethis$V4, writethis$V5, ".", writethis$V7, writethis$V8, writethis$V9, writethis$V10, writethis$V11, geom, sep = "\t"), file = bedfile, ncolumns = 1)
	return(bedfile)
}
#=======================> DONE! 



# ==================
# Fixing narrowpeak
# ==================
bedfiles = fixnpf(bednames) #read in all bed files (samples and control)
#=======================> DONE!
