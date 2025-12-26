cwlVersion: v1.2
class: Workflow

inputs:
        rna_reads_human: File
	ref_genome: Directory
        annotations: File

steps:
 quality_control:
 run: bio-cwl-tools/fastqc/fastqc_2.cwl
 in:
 reads_file: rna_reads_human
 out: [html_file]

mapping_reads:
 requirements:
 ResourceRequirement:
 ramMin: 9000
 run: bio-cwl-tools/STAR/STAR-Align.cwl
 in:
 RunThreadN: {default: 4}
 GenomeDir: ref_genome
 ForwardReads: rna_reads_human
 OutSAMtype: {default: BAM}
 SortedByCoordinate: {default: true}
 OutSAMunmapped: {default: Within}
 out: [alignment]

index_alignment:
 run: bio-cwl-tools/samtools/samtools_index.cwl
 in:
 bam_sorted: mapping_reads/alignment
 out: [bam_sorted_indexed]

count_reads:
 requirements:
 ResourceRequirement:
 ramMin: 500
 run: bio-cwl-tools/subread/featureCounts.cwl
 in:
 mapped_reads: index_alignment/bam_sorted_indexed
 annotations: annotations
 out: [featurecounts]

outputs:
qc_html:
 type: File
 outputSource: quality_control/html_file
 bam_sorted_indexed:
 type: File
 outputSource: index_alignment/bam_sorted_indexed
 featurecounts:
 type: File
 outputSource: count_reads/featurecounts

