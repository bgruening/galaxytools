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


# ============================================== 
# Parsing Arguments (source: phantom SPP script)
# ============================================== 
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line
#nargs = length(args) # number of arguments

#Set arguments to default values
infile = NA # input file
out = NA # output directory

#Parsing arguments and storing values
for (each.arg in args) {
	#input bed file
	if (grepl('^-infile=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
				infile <- arg.split[2]
		} else {
			stop('No input file')
		} 
	}
	if (grepl('^-out=',each.arg)) {			
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			out <- arg.split[2] 
		} else {
			stop('No output directory')
		}
	}
}

#=======================> DONE! 


# =================== 
# Do the processing
# ===================

#Read in variables
if (file.exists(infile)) {

message("@", infile)
xcorrs = read.table(infile, header=FALSE)
shifts = as.numeric(xcorrs$V2) 


#how many chromosomes
num = length(xcorrs$V1)

if (length(xcorrs$V2) > 0) {
if (length(xcorrs$V2) > 1)  {
#shift stats
avg = mean(shifts)
std = sd(shifts)
mod = as.numeric(names(sort(-table(shifts)))[1])

#find out fragment length
fraglength = 100 #default
default = 1 #was default shift assumed
if (length(shifts) != 0) {
	fraglength = ceiling(mean(shifts))
	default = 0
} 
} else {
avg = xcorrs$V2
std = xcorrs$V2
mod = xcorrs$V2
fraglength = xcorrs$V2
default = 0
}
} else {
fraglength = 100 #default
}
} else {
fraglength = 100
default = 1
}
#=======================> DONE! 





# ===================== 
# Write result to File
# =====================
if (default == 0) {
writethis = paste0("Number of chromosomes used:", num, "\nAverage strand shift:", avg, "\nStandard Deviation:", std, "\nMode:", mod, "\nFragment Length:", fraglength)
write(writethis, file = paste0(out, "/xcorrsummary.txt"))

if (length(xcorrs$V2) > 1)  {
	i = pdf(file=paste0(out, "/xcorrsummary.pdf"))
	i = plot(density(shifts), main="Strand Shifts Distribution over Chromosomes")
	i = dev.off()
}
	message("Fragment Length:", fraglength, "\n")
} else if (default == 1) {

writethis = paste0("Number of chromosomes used:", "0", "\nAverage strand shift:", "NA", "\nStandard Deviation:", "NA", "\nMode:", "NA", "\nFragment Length:", fraglength)
write(writethis, file = paste0(out, "/xcorrsummary.txt"))

	message("Fragment length estimation failed, default fragment length assumed:", fraglength, "\n")
}
#=======================> DONE!
