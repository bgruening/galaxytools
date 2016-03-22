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

### RiboTaper master File


set -e

if [ $# -ne 8 ]; then
	echo "Usage: ./Ribotaper.sh <Ribo_bamfile> <RNA_bamfile> <annotation_dir> <comma-sep_read_lenghts_ribo> <comma-sep_cutoffs> <scripts_dir> <bedtools_dir> <n_cores> "
	exit 1
fi
if ! [[ -f "$1" ]]; then
     echo "!!!!!   ribo_bam file not found!."
     exit 1
   fi

if ! [[ -f "$2" ]]; then
     echo "!!!!!   ribo_bam not found!."
     exit 1
   fi

if [ ! -d "$3" ]; then
     echo "!!!!!   annotation_directory not found!."
    exit 1
fi

if [ ! -d "$6" ]; then
     echo "!!!!!   scripts_directory not found!."
    exit 1
fi

if [ ! -d "$7" ]; then
     echo "!!!!!   bedtools_directory not found!."
    exit 1
fi


re='^[0-9]+$'
if ! [[ "$8" =~ $re ]] ; then
   echo "!!!!!   n of cores not valid"
   exit 1
fi


if [ "$8" == 1 ]; then
     echo "!!!!!   n of cores required >1."
    exit 1
fi



ribo_bam="`readlink -f $1`"
rna_bam="`readlink -f $2`"
annot_dir="`readlink -f $3`"
read_len=$4
cutoffs=$5
scripts_dir="`readlink -f $6`"
bedtools_dir="`readlink -f $7`"
n_of_cores=$8


echo "Parameters used:"
echo ""

echo "<Ribo_bamfile> $ribo_bam"
echo "<RNA_bamfile> $rna_bam"
echo "<annotation_dir> $annot_dir"
echo "<comma-sep_read_lenghts_ribo> $read_len"
echo "<comma-sep_cutoffs> $cutoffs"
echo "<scripts_dir> $scripts_dir"
echo "<bedtools_dir> $bedtools_dir"
echo "<n_cores> $8"
echo ""
echo "---------------"
echo ""



echo "CCDS ORF finding..."

Rscript $scripts_dir"/CCDS_orf_finder.R" $annot_dir $scripts_dir $bedtools_dir $n_of_cores

echo "NONCCDS ORF finding..."

Rscript $scripts_dir"/NONCCDS_orf_finder.R" $annot_dir $scripts_dir $bedtools_dir $n_of_cores

# Groups ORFs and creates BED files + protein fasta database

echo "Grouping ORFs and creating protein fasta database..."

Rscript $scripts_dir"/create_protein_db.R"

# makes summary plot for the found ORFs

echo "Summarizing ORF finding results"

Rscript $scripts_dir"/ORF_final_results.R"

echo "RiboTaper analysis finished !!!"
