#flye 2.9.3
threads=8
#Test 1
flye --pacbio-raw ecoli_{01..07}.fastq.gz --threads ${threads}  --out-dir test1

#Test 2
flye --nano-raw nanopore.fasta.gz --out-dir test2 --iterations 0 --threads ${threads}

#Test 3
flye --nano-raw ecoli_hifi_* --genome-size 3980000 --threads ${threads} --asm-coverage 30 --out-dir test3

#Test 4
flye --pacbio-raw ecoli_{01..07}.fastq.gz --threads ${threads} --out-dir test4 --meta

#Test 5
flye --nano-hq ecoli_hifi_* --threads ${threads} --out-dir test5 --min-overlap 1000

#Test 6
flye --pacbio-hifi ecoli_hifi_* --threads ${threads}  --hifi-error 0.21 --out-dir test6

#Test 7
flye --pacbio-corr ecoli_hifi_* --threads ${threads}  --hifi-error 0.21 --out-dir test7 --keep-haplotypes

#Test 8
flye --nano-hq ecoli_hifi_* --threads ${threads}  --out-dir test8 --min-overlap 1000 --scaffold

#Test 9
flye --nano-raw nanopore.fasta.gz --threads ${threads}  --out-dir test9 --iterations 0 --no-alt-contigs