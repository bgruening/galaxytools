This file will contain example commands for newly introduced options. (If you want to copy and paste the commands, do not copy the $-sign.)
For more examples and information, take a look at the Manual (http://prinseq.sourceforge.net/manual.html).


(1) Graph data
==============
To generate the graphs from the web version or the HTML report, you can use the -graph_data option:
$ perl prinseq-lite.pl -verbose -fastq example1.fastq -graph_data example1.gd -out_good null -out_bad null

The verbose mode shows the progress and the "-out_good null -out_bad null" prevents PRINSEQ from generating any other output files than the specified test.gd file containing the graphs data.

To generate the graph as PNG files, you can use the prinseq-graphs -png_all option:
$ perl prinseq-graphs.pl -i example1.gd -png_all -o example1

To generate the HTML report containing all the tables and figures from the web version, you can use the prinseq-graphs -html_all option:
$ perl prinseq-graphs.pl -i example1.gd -html_all -o example1


(2) Consider exact duplicates only
==================================
When you process large files, the duplicate removal will require a lot of memory. To reduce the amount of memory required and speed up the process (at the cost of only removing forward and reverse exact duplicates), you can use the option -exact_only when generating the graphs data:
$ perl prinseq-lite.pl -verbose -fastq example1.fastq -graph_data example1.gd -out_good null -out_bad null -exact_only

Note that for processing the data, if you specify -derep 1, -derep 4, or -derep 14 then the exact_only option will be used automatically.


(3) Duplicate threshold and no quality header information
=========================================================
Process the data (-fastq example1.fastq) with status report (-verbose), remove exact sequence duplicates (-derep 1) that occur more than 2 times (-derep_min 3) and save the sequences passing the filter in example1_good.fastq (-out_good example1_good) without the quality header (-no_qual_header) and the filtered sequences in example1_bad.fastq (-out_bad example1_bad):
$ perl prinseq-lite.pl -verbose -fastq example1.fastq -derep 1 -derep_min 3 -out_good example1_good -out_bad example1_bad -no_qual_header


(4) Paired-end data
===================
Paired-end data is processed similar to single read data. The only difference is that two input files are required (either two FASTA or two FASTQ files). The second file is specified either with "-fasta2 file.fa" or "-fastq2 file.fq".

