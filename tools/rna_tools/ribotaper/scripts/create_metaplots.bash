#!/bin/bash


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


##This script creates aggregate plots around start-stop codons, it uses as arguments a bam file, a bed file for start-stop positions (e.g. the one produced by the create_annotation.bash script), a name as an appendix for further analysis, the RiboTaper scripts directory


if [ $# -ne 5 ]; then
	echo "Usage: create_metaplots.bash <ribo.bam> <bedfile> <name> <bedtools_dir> <scripts_dir>"
	exit 1
fi
if ! [[ -f "$1" ]]; then
     echo "!!!!!   ribo.bam file not found!."
     exit 1
   fi
if ! [[ -f "$2" ]]; then
     echo "!!!!!   start_stop bed file not found!."
     exit 1
   fi

bedtools_dir=$4

echo "Downsampling to 10%..."

samtools view -s 1.1 $1 > sample_to_metapl.sam

cat <( samtools view -H $1 ) <(cat sample_to_metapl.sam) | samtools view - -bS > sample_to_metapl.bam
rm sample_to_metapl.sam

echo "Intersecting alignments with start/stop sites ..."

$bedtools_dir"/bamToBed" -i sample_to_metapl.bam -bed12 -split | $bedtools_dir"/windowBed" -w 100 -sm -b stdin -a $2 | awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18}' | sort -k1,1 -k2,2g | $bedtools_dir"/closestBed" -s -t "last" -a stdin -b $2 > $3

if !  [[ -s $3 ]]; then
     echo "!!!!!   no intersections found, check input files"
     exit 1
   fi



echo "Creating metaplots..."

scripts_dir=$5

Rscript $scripts_dir"/metag.R" $3

mkdir metaplots

mv $3*.pdf metaplots/



echo "Done !!! "
