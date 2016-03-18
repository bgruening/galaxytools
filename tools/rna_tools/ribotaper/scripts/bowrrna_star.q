#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=15G
#$ -e "error_mapp_bowrrnastar"
#$ -o "out_mapp_bowrrnastar"
#$ -cwd


##fastq, rRNA ref, star_in, start_stop bed file

fastq=$1

full_fastq="`readlink -f $fastq`"

name_exp="`echo $fastq | sed 's/\.fastq//g'`"

full_name_exp="`echo $full_fastq | sed 's/\.fastq//g'`"

/data/ohler/Lorenzo/bins/bowtie1/bowtie --best -S -p 4 --al $name_exp"_rRNA.fastq" --un $name_exp"_notrRNA.fastq" $2 $1 > /dev/null


mkdir "starmapp_star_"$name_exp/

cd "starmapp_star_"$name_exp/

/data/ohler/Lorenzo/STAR_2.3.1z1/STAR --genomeDir $3 --alignEndsType EndToEnd --readFilesIn $full_name_exp"_notrRNA.fastq" --runThreadN 4 --outFilterMismatchNmax 4 --outFilterMultimapNmax 8 --chimScoreSeparation 10 --chimScoreMin 20 --chimSegmentMin 15 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJoverhangMin 500 --outFileNamePrefix "star_"$name_exp"_" --outReadsUnmapped Fastx 
samtools view -bS "star_"$name_exp"_"Aligned.out.sam | samtools sort - "star_"$name_exp"_"Aligned.out.sorted
samtools index "star_"$name_exp"_"Aligned.out.sorted.bam


/data/ohler/website/files/RiboTaper/Version_1.2/create_metaplots.bash "star_"$name_exp"_"Aligned.out.sorted.bam $4 $name_exp"_metaplots"

echo "done"$name_exp"!!!"

