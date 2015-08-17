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



use feature qw(say);


# ================== 
# Parsing Arguments
# ================== 

#initialize to NULL
my $bed_file = NULL; #bed file
my $shift_size = NULL; #shift size
my $read_length = NULL; #read length

#Parse the arguments
$bed_file = $ARGV[0]; #bed file
$shift_size = $ARGV[1];	#shift size
$read_length = $ARGV[2]; #read length

#=======================> DONE! 



# ========================================
# Parse the bed file and extend the reads
# ========================================

#open the file
open(DATA, $bed_file) || die("Can't open the bed file, probably you gave me the wrong path!");

#loop through the rest of the file line by line (To do: look for a faster way)
while (<DATA>) {
    my ($start, $strand) = split(/\t/,$_,2);
    $strand =~ s/\015?\012?$//;


    
	#plus strand
	if ($strand eq '+') {
		#do nothing
	}
	#minus strand
	elsif ($strand eq '-') {
		$start = $start - $shift_size + $read_length; #difference between $start and $end should be equal to fragment length ($shift_size)
	}
	#bad format 
	else {
		die("It appears the bed file is not formatted properly. Specifically, I expect to find either + or - in the second column, but I found something else!");	
	}
	if ($start >= 0) {
		#Now write the new line
		say join "\t", $start, $strand;
	}
}
#=======================> DONE! 
