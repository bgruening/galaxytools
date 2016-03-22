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


##This script creates annotation files to be used in the RiboTaper pipeline, it uses as arguments a gtf file, a genome fasta file, a logical value for using the CCDS annotation (true or false) , a logical value for using the APPRIS annotation (true or false), a destination folder, the bedtools executables directory,  the RiboTaper scripts directory




if [ $# -ne 7 ]; then
	echo "Usage: ./create_annotation_files.bash <gtf_file> <genome_fasta_file(indexed)> <use_ccdsid?> <use_appris?> <dest_folder> <bedtools_path> <scripts_dir>"
	exit 1
fi
if ! [[ -f "$1" ]]; then
     echo "!!!!!   gtf_file not found!."
     exit 1
   fi

if ! [[ -f "$2" ]]; then
     echo "!!!!!   genome fasta not found!."
     exit 1
   fi

if ! [ "$3" = true ]  ; then
     if ! [ "$3" = false ]; then
          echo "use_ccdsid = "true" or "false""
          exit 1
          fi
fi

if ! [ "$4" = true ]  ; then
     if ! [ "$4" = false ]; then
          echo "use_appris = "true" or "false""
          exit 1
          fi
fi




gencode_ann=$1
genc_full="`readlink -e $gencode_ann`"


genome=$2
genome_full=`readlink -e $genome`

scripts_dir=$7
scripts_dir_full=`readlink -e $scripts_dir`

dest_folder=$5
dest_folder_full=`readlink -f $dest_folder`


bedtools_path=$6
bedtools_path_full=`readlink -e $bedtools_path`

echo "Parameters used:"
echo ""


echo "<gtf_file> $genc_full"
echo "<genome_fasta_file(indexed)> $genome_full"
echo "<use_ccdsid?> $3"
echo "<use_appris?> $4"
echo "<dest_folder> $dest_folder_full"
echo "<bedtools_path> $bedtools_path_full"
echo "<scripts_dir> $scripts_dir_full"
echo ""
echo "---------------"
echo ""




echo "creating directory..."
mkdir -p $dest_folder_full

cd $dest_folder_full

echo "Extracting gene names + biotypes from gtf..."
awk '{ for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++)  sub(/_/,"-", $(x+1))}1' $genc_full | awk '{ for (x=1;x<=NF;x++) if ($x~"^transcript_id") for (y=1;y<=NF;y++)  sub(/_/,"-", $(x+1))}1'   > gtf_corr

#less $genome_full | sed 's/_/-/g' > genome_corr.fasta

#genome_full=genome_corr.fasta

genc_full=gtf_corr

#grep out at each step!

awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~/gene_type|gene_biotype/) for (z=1;z<=NF;z++) if ($z~"gene_name") print $(x+1) "\t" $(y+1) "\t" $(z+1)}' $genc_full  | sort | uniq | sed 's/;//g' | sed 's/"//g' > gene_name_type
less gene_name_type | cut -f 1 | grep -Fvf - $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (z=1;z<=NF;z++) if ($z~"gene_name") print $(x+1) "\t" "no_biotype" "\t" $(z+1)}' | sort | uniq | sed 's/;//g' | sed 's/"//g' > gene_name_notype
less gene_name_type | cut -f 1 | grep -Fvf - $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~/gene_type|gene_biotype/)  print $(x+1) "\t" $(y+1) "\t" "no_name"}'  | sort | uniq | sed 's/;//g' | sed 's/"//g' > gene_noname_type


cat gene_name_type gene_name_notype gene_noname_type > gene_annot_name_pre
less gene_annot_name_pre | cut -f 1 | grep -Fvf - $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") print $(x+1) "\t" "no_biotype" "\t" "no_name"}'  | sort | uniq | sed 's/;//g' | sed 's/"//g' > gene_noname_notype

cat gene_annot_name_pre gene_noname_notype | sort | uniq > gene_annot_names


rm gene_name_type gene_name_notype gene_noname_type gene_annot_name_pre gene_noname_notype



echo "creating bed_files..."

#TAKE CDS OF CCDS REGIONS

if [ "$3" = true ] ; then
    awk '{if($3=="CDS") print $0}' $genc_full | grep ccdsid | awk '{ for (x=1;x<=NF;x++) if ($x~"^gene_id") print $1 "\t" $4-1 "\t" $5 "\t" "CCDS" "\t" $(x+1) "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > unique_ccds.bed
    less $genc_full | grep ccdsid | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1)"\t" $(x+1) "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_ccds_ccdsid.bed

fi

if [ "$3" = false ] ; then
    awk '{if($3=="CDS") print $0}' $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") print $1 "\t" $4-1 "\t" $5 "\t" "CCDS" "\t" $(x+1)  "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > unique_ccds.bed
fi

#STORE CCDS GENES
less unique_ccds.bed | cut -f 5 | sort | uniq > genes_ccds

#STORE COORDINATES CDS CCDS
less unique_ccds.bed | cut -f 1-3 | tr '\t' '_'  > coords_ccds

#TAKE ALL EXONS OF CCDS GENES
grep -Ff genes_ccds $genc_full | awk '{if($3=="exon") print $0}' | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") print $1 "\t" $4-1 "\t" $5 "\t" "EXONCCDS" "\t" $(x+1) "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > unique_exons_ccds.bed

#STORE COORDINATES EXONS CCDS
less unique_exons_ccds.bed | awk '{print $1"_"$2"_"$3 "\t" $0}' > coords_unique_exons_ccds.bed
#TAKE OUT CDS CCDS FROM EXONS CCDS
grep -Fvf coords_ccds coords_unique_exons_ccds.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" "EXONCCDS" "\t" $6 "\t" $7}' > unique_exons_ccds.bed
#REMOVE COORDS
rm coords_ccds coords_unique_exons_ccds.bed




#TAKE EXONS OF NONCCDS GENES
grep -Fvf genes_ccds $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" "EXONnonCCDS" "\t" $(x+1) "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > unique_nonccds.bed

#TAKE SEQUENCES, STRANDED INFO

echo "creating fasta sequences..."

fastaFb="$bedtools_path_full/fastaFromBed"

$fastaFb -s -fi $genome_full -bed unique_ccds.bed  -fo unique_ccds_seq.fa

awk '{$7=$1; gsub( "_","",$7 ) ; $4=$7"_"$2"_"$3"_"$4"_"$5"_"$6; print $0}' OFS="\t" unique_ccds.bed  | cut -f  1-6| $fastaFb -fi $genome_full -name -bed - -tab -fo unique_ccds_seq_name_tab
paste <(cut -f 1 unique_ccds_seq_name_tab |  tr '_' '\t') <(cut -f 2 unique_ccds_seq_name_tab | sed 's/[A-Z]/& /g') > sequences_ccds

awk '{$7=$1; gsub( "_","",$7 ) ; $4=$7"_"$2"_"$3"_"$4"_"$5"_"$6; print $0}' OFS="\t" unique_exons_ccds.bed  | cut -f  1-6| $fastaFb -fi $genome_full -name -bed - -tab -fo unique_exons_ccds_seq_name_tab
paste <(cut -f 1 unique_exons_ccds_seq_name_tab   |  tr '_' '\t') <(cut -f 2 unique_exons_ccds_seq_name_tab | sed 's/[A-Z]/& /g') > sequences_exonsccds

awk '{$7=$1; gsub( "_","",$7 ) ; $4=$7"_"$2"_"$3"_"$4"_"$5"_"$6; print $0}' OFS="\t" unique_nonccds.bed  | cut -f  1-6|  $fastaFb -fi $genome_full -name -bed - -tab -fo unique_nonccds_seq_name_tab
paste <(cut -f 1 unique_nonccds_seq_name_tab   | tr '_' '\t') <(cut -f 2 unique_nonccds_seq_name_tab | sed 's/[A-Z]/& /g') > sequences_nonccds


$fastaFb -s -fi $genome_full -bed unique_exons_ccds.bed  -fo unique_exons_ccds_seq.fa
$fastaFb -s -fi $genome_full -bed unique_nonccds.bed -fo unique_exons_nonccds_seq.fa

#CAT SEQUENCES TOGETHER FOR ORF FINDING
cat unique_ccds_seq.fa unique_exons_ccds_seq.fa > unique_ccds_exonccds_seq.fa

#make all CDS regions
less $genc_full | awk '{if($3=="CDS") print $0}' | awk '{for (x=1;x<=NF;x++) if ($x~"^transcript_id") print $1 "\t" $4-1 "\t" $5 "\t" "cds" "\t" $(x+1) "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > all_cds.bed

echo "assembling transcript information..."

#TAKE TRANSCR CCDS
grep -Ff genes_ccds $genc_full | awk '{ for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1) "\t" $(x+1)  "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_ccds.bed


#TAKE TRANSCR APPRIS CCDS

if [ "$4" = true ] ; then

grep -Ff genes_ccds $genc_full | grep appris_prin | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1)"\t" $(x+1) "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_ccds_appris_prin.bed
cut -f 5 transcr_exons_ccds_appris_prin.bed | sort | uniq > genes_appris_prin
grep -Ff genes_ccds $genc_full | grep -Fvf genes_appris_prin - | grep appris | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1)"\t" $(x+1) "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_ccds_appris_noprin.bed
cut -f 5 transcr_exons_ccds_appris_noprin.bed | sort | uniq > genes_appris_noprin
grep -Ff genes_ccds $genc_full | grep -Fvf genes_appris_prin - | grep -Fvf genes_appris_noprin | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1)"\t" $(x+1) "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_ccds_noappris_noprin.bed
cut -f 5 transcr_exons_ccds_noappris_noprin.bed | sort | uniq > genes_noappris_noprin
cat transcr_exons_ccds_appris_prin.bed transcr_exons_ccds_appris_noprin.bed transcr_exons_ccds_noappris_noprin.bed > transcr_exons_ccds_appris.bed
cat genes_appris_prin genes_appris_noprin genes_noappris_noprin > genes_ccds_appris

fi

#TAKE TRANSCR NONCCDS
grep -Fvf genes_ccds $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") for (y=1;y<=NF;y++) if ($y~"^transcript_id") if($3=="exon") print $1 "\t" $4-1 "\t" $5 "\t" $(y+1)"\t" $(x+1) "\t" $7}'| sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k2,2n | uniq > transcr_exons_nonccds.bed

#start_stop_cds

less $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") if($3=="start_codon" || $3=="stop_codon") print $1 "\t"$4-1 "\t"$5 "\t" $3 "\t" $(x+1) "\t"$7}' | sed 's/;//g' | sed 's/"//g' | sort -k1,1 -k2,2g | uniq | awk 'p{print $0 "\t" $2-p}{p=$2}' | tac | awk 'p{print $0 "\t" $2-p}{p=$2}' | tac | awk '{if($NF<-100 || $(NF-1)>100) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > start_stops_FAR.bed
Rscript $scripts_dir_full"/write_startstops.R"
#make cds transcript coords
echo "Creating transcript cds coordinates from gtf..."

awk '{ for (x=1;x<=NF;x++) if ($x~"^transcript_id") if ( $3=="exon" || $3=="CDS" ) print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $(x+1)}' $genc_full | sed 's/"//g' | sed 's/;//g' | sort -k1,1 -k3,3g > exons_cds_all
Rscript $scripts_dir_full"/gtf_to_start_stop_tr.R"

#make cds frames

less $genc_full | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") if($3=="CDS") print $1 "_" $4-1 "_" $5 "_" "CCDS" "_" $(x+1) "\t" $8 "\t"$7 "\t" $5-($4-1)}' | sed 's/;//g' | sed 's/"//g' | sort -k1,1 -k2,2 | uniq > frames_ccds


#take all exonic regions
less $genc_full | awk '{if($3=="exon") print $0}' | awk '{for (x=1;x<=NF;x++) if ($x~"^gene_id") print $1 "\t" $4-1 "\t" $5 "\t" "exon" "\t" $(x+1) "\t" $7 }' | sort -k1,1 -k2,2n | uniq | sed 's/;//g' | sed 's/"//g' > all_exons.bed

Rscript $scripts_dir_full"/genes_coor.R"
echo "Done!"
