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


##This script counts the multi-mapping/unique reads ratio per region, including coverage information, it uses as arguments a bed file, a name as an appendix for further analysis, the bedtools exec directory

if [ $# -ne 3 ]; then  
	echo "Usage: analyze_multi.bash <bed_file> <name> <bedtools_dir>"
	exit 1
fi
if ! [[ -f "$1" ]]; then
     echo "!!!!!   bed file not found!."
     exit 1
   fi

bedtools_dir=$3

echo "-----Intersecting with unique/best alignments-----"

$bedtools_dir"/coverageBed" -s -split -abam RIBO_unique.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RIBO_unique_counts"_$2"

$bedtools_dir"/coverageBed" -s -split -abam RNA_unique.bam -b $1 | sort -k1,1 -k2,2g |  sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RNA_unique_counts"_$2"

$bedtools_dir"/coverageBed" -s -split -abam RIBO_best.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' |  awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RIBO_best_counts"_$2"

$bedtools_dir"/coverageBed" -s -split -abam RNA_best.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RNA_best_counts"_$2"

echo "-----Done !!!-----"



