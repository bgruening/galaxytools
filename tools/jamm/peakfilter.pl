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

#Parse the arguments
$bed_file = $ARGV[0]; #bed file
#=======================> DONE! 




# ========================================
# Parse the bed file and extend the reads
# ========================================
#open the file
open(DATA, $bed_file) || die("Can't open the bed file, probably you gave me the wrong path!");

while (<DATA>) {
    my ($chr, $start, $end, $name, $score, $strand, $signal, $pvalue, $qvalue, $summit, $minpeak, $geom) = split(/\t/,$_,12);
    
	my $size = $end - $start;
	
	if($size >= $minpeak) {
		if($signal > $geom) {
			#Now write the new line
			$summit =~ s/\015?\012?$//;
			say join "\t", $chr, $start, $end, $name, $score, $strand, $signal, $pvalue, $qvalue, $summit;
		}
	}
}
#=======================> DONE! 
