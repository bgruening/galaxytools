Important notes
----------------

- At least 2 cores should be used for running ``ribotaper part 3: ribosome profiling``,  therefore please configure **job_conf.xml** accordingly.

- We ran the RiboTaper analysis on an SGE cluster, using 7 cores and h_vmem 8G. For each dataset, the complete RiboTaper workflow (from the bam files to final results) took ~ 1 day.

- The current RiboTaper framework is not designed to identify and quantify ORFs on different transcripts. This means the transcript annotation is crucial.

- Be careful about using scaffolds, both in the genome and GTF files, which may slow the whole pipeline.
