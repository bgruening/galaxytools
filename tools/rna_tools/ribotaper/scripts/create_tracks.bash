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


##This script creates the data_tracks files, it uses as arguments a bed file, a name as an appendix for further analysis, the bedtools exec directory

if [ $# -ne 4 ]; then  
	echo "Usage: create_tracks.bash <bed_file> <fasta_file> <name> <bedtools_dir>"
	exit 1
fi
if ! [[ -f "$1" ]]; then
     echo "!!!!!   bed file not found!."
     exit 1
   fi

if ! [[ -f "$2" ]]; then
     echo "!!!!!   fasta file not found!."
     exit 1
   fi
bedtools_dir=$4

mkdir -p data_tracks

echo "-----Calculating coverage tracks for each exon-----"

$bedtools_dir"/coverageBed" -s -split -abam RIBO_unique.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RIBO_unique_counts"_$3"

$bedtools_dir"/coverageBed" -s -split -abam RNA_unique.bam -b $1 | sort -k1,1 -k2,2g |  sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RNA_unique_counts"_$3"

$bedtools_dir"/coverageBed" -s -split -abam RIBO_best.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' |  awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RIBO_best_counts"_$3"

$bedtools_dir"/coverageBed" -s -split -abam RNA_best.bam -b $1 | sort -k1,1 -k2,2g | sed 's/_//g' | awk '{ print $1 "_" $2 "_" $3 "_" $4 "_" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > RNA_best_counts"_$3"


$bedtools_dir"/coverageBed" -s -d -a P_sites_all -b $1 |  awk '{ print $1 ";" $2 ";" $3 ";" $4 ";" $5";" $6 "\t" "_" $8}' | awk -F"\t" '{if(a[$1])a[$1]=a[$1]" "$NF; else a[$1]=$NF}END{for (i in a)print i "\t" a[i]}'| sed 's/_//g' | sed 's/;/\t/g' | sort -k1,1 -k4,4 -k2,2g > data_tracks/P_sites_all_tracks"_$3"


$bedtools_dir"/coverageBed" -s -d -split -abam RIBO_best.bam -b $1 | awk '{ print $1 ";" $2 ";" $3 ";" $4 ";" $5";" $6 "\t" "_" $8}' | awk -F"\t" '{if(a[$1])a[$1]=a[$1]" "$NF; else a[$1]=$NF}END{for (i in a)print i "\t" a[i]}'| sed 's/_//g' | sed 's/;/\t/g' | sort -k1,1 -k4,4 -k2,2g > data_tracks/RIBO_tracks"_$3"

$bedtools_dir"/coverageBed" -s -d -split -abam RNA_best.bam -b $1 | awk '{ print $1 ";" $2 ";" $3 ";" $4 ";" $5";" $6 "\t" "_" $8}' | awk -F"\t" '{if(a[$1])a[$1]=a[$1]" "$NF; else a[$1]=$NF}END{for (i in a)print i "\t" a[i]}'| sed 's/_//g' | sed 's/;/\t/g' | sort -k1,1 -k4,4 -k2,2g > data_tracks/RNA_tracks"_$3"

$bedtools_dir"/coverageBed" -s -d -a Centered_RNA -b $1 |  awk '{ print $1 ";" $2 ";" $3 ";" $4 ";" $5";" $6 "\t" "_" $8}' | awk -F"\t" '{if(a[$1])a[$1]=a[$1]" "$NF; else a[$1]=$NF}END{for (i in a)print i "\t" a[i]}'| sed 's/_//g' | sed 's/;/\t/g' | sort -k1,1 -k4,4 -k2,2g > data_tracks/Centered_RNA_tracks"_$3"

echo "-----Merging tracks together-----"


cat data_tracks/P_sites_all_tracks"_$3" data_tracks/RIBO_tracks"_$3" data_tracks/RNA_tracks"_$3" data_tracks/Centered_RNA_tracks"_$3" $2 | tr '\t' '_' | sed 's/_/\t/6' | awk -F"\t" '{a[$1]=a[$1]"\n" $1 "\t" $2}END{for (i in a)print i "\t" a[i];}' | awk -F"\t" '{if($2>=0)print $0}' | sed 's/_/ /5' | sed 's/\t/ /1' > data_tracks/Psit_Ribo_Rna_Cent_tracks"_$3"

cut -f 1 data_tracks/Psit_Ribo_Rna_Cent_tracks"_$3" -d" " > data_tracks/index_tracks"_$3"



