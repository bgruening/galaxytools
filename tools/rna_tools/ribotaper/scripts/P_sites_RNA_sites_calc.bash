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


##This script calculates P-sites and RNA sites positions, it uses as arguments a comma-separated list of read lengths to be used, a comma-separated list of offsets, the bedtools exec directory

if [ $# -ne 3 ]; then  
        echo  "------------Usage: P_sites_RNA_sites_calc.bash <read_lengths> <offsets> <bedtools_dir> "
        exit 1
fi

lengths=$1

lengths=(${lengths//,/ })

bedtools_dir=$3

offsets=$2
offsets=(${offsets//,/ })

if [ ${#lengths[@]} -ne ${#offsets[@]} ]; then
	echo ${#lengths[@]}
	echo ${#offsets[@]}
        echo  "------------The number of read lengths and offsets differ! Insert comma-separated value for read lengths and offsets e.g. 28,29 11,12 Usage: P_sites_RNA_sites_calc.bash <read_lengths> <offsets> <bedtools_dir> "
        exit 1

fi
n_frag=${#lengths[@]}

for (( i=0; i<${n_frag}; i++ ));
do
  len=${lengths[$i]}
  offs=${offsets[$i]}
  echo "------------processing" $len "nt reads with offset of +" $offs
  $bedtools_dir"/bamToBed" -cigar -bed12 -i RIBO_best.bam | awk -v env_var=$len -F"\t" '{split($11,c,","); if((c[1]+c[2]+c[3])==env_var) print $0 "\t" c[1]+c[2]+c[3] > "tmp_align_len"}'
  less tmp_align_len | awk -F"\t" '{split($11,c,","); print $0 "\t" c[1]+c[2]+c[3]+c[4]+c[5] }' | awk -v env_var=$offs '{split($11,c,","); split($12,d,","); if($6=="+" && c[1]<=env_var && (c[1]+c[2])>env_var) $2=($2)+d[2]+(env_var-c[1]); if($6=="+" && c[1]>env_var) $2=($2)+env_var; if($6=="+" && c[1]<=env_var && (c[1]+c[2])<12 && (c[1]+c[2]+c[3])>env_var) $2=($2)+d[3]+(env_var-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<env_var) $2=($2)+d[4]+(env_var-c[1]-c[2]-c[3]); if($6=="-" && c[1]>=$NF-env_var) $2=$2+$NF-env_var-1; if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])>($NF-env_var)) $2=($2)+d[2]+(($NF-env_var-1)-c[1]); if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])<($NF-env_var) && (c[1]+c[2]+c[3])>=($NF-env_var)) $2=($2)+d[3]+(($NF-env_var-1)-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<($NF-env_var)) $2=($2)+d[4]+(($NF-env_var)-c[1]-c[2]-c[3]); print $0}' OFS="\t" | awk '{$3=$2+1 ; print $0 }' OFS="\t" | awk '{if($2>0) print $0}' OFS="\t" > P_sites_len
  rm tmp_align_len
  mv P_sites_len tmp_P_sites_"$len"

  
done

cat tmp_P_sites_* > P_sites_all
rm tmp_P_sites_*
echo "------------Done!"

echo "------------processing RNA-seq with offset of + 25"

$bedtools_dir"/bamToBed"  -cigar -bed12 -i RNA_best.bam |  awk -F"\t" '{split($11,c,","); print $0 "\t" c[1]+c[2]+c[3]+c[4]+c[5] }' | awk -v env_var=25 '{split($11,c,","); split($12,d,","); if($6=="+" && c[1]<env_var && (c[1]+c[2])>env_var) $2=($2)+d[2]+(env_var-c[1]); if($6=="+" && c[1]>=env_var) $2=($2)+env_var; if($6=="+" && c[1]<env_var && (c[1]+c[2])<12 && (c[1]+c[2]+c[3])>env_var) $2=($2)+d[3]+(env_var-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<env_var) $2=($2)+d[4]+(env_var-c[1]-c[2]-c[3]); if($6=="-" && c[1]>=$NF-env_var) $2=$2+$NF-env_var-1; if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])>=($NF-env_var)) $2=($2)+d[2]+(($NF-env_var-1)-c[1]); if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])<($NF-env_var) && (c[1]+c[2]+c[3])>=($NF-env_var)) $2=($2)+d[3]+(($NF-env_var-1)-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<($NF-env_var)) $2=($2)+d[4]+(($NF-env_var)-c[1]-c[2]-c[3]); print $0}' OFS="\t" | awk '{$3=$2+1 ; print $0 }' OFS="\t" | awk '{if($2>0) print $0}' OFS="\t" > Centered_RNA

echo "------------P_sites and RNA_sites calculated !!!"

