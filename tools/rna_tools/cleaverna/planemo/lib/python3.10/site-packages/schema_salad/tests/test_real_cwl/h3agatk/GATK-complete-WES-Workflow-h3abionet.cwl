#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement

doc: |
  # H3ABioNet GATK Germline Workflow

  # Overview
  A [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS) germline workflow designed to work with GATK 3.5  (Van der Auwera et al., 2013).

  For more information see our [GitHub](https://github.com/h3abionet/h3agatk) site.

  # Workflow Summary

  ![pipeline](https://raw.githubusercontent.com/h3abionet/h3agatk/master/workflows/GATK/gatk_germline_small.png)

  # Workflow Tool Details

  ## FastQC
  FastQC is used as an initial QC step where the input files are checked for usual metrics such as:
  	- Read length
  	- Reads distribution
  	- GC content
  	- ...

  ## Trimmomatic
  Trimmomatic is the entry point of the pipeline, it is used to cleanup the reads in the input fastq files from any sequencing adaptors.

  ## BWA
  [BWA](http://bio-bwa.sourceforge.net) is used to align the reads from the the input fastq files -paired-ends- (Li, 2013). We use specifically `bwa mem` as recommended by the [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS). BWA produces a SAM file containing the aligned reads against the human reference genome (hg19, GATK bundle build 2.8).

  As GATK tools downstream requires properly formatted Read Group information. We add by default 'toy' Read Group information while processing the alignment to the output SAM file. we specifically use the flag `-R '@RG\tID:foo\tSM:bar\tLB:library1'`.

  ## SAMtools
  [SAMtools](http://www.htslib.org) (Li et al., 2009) are used few times in the pipeline:
  	1. Convert BWA's output from a SAM format to a BAM format
  	2. Sort the reads in the generated BAM file in step 1 (above)
  	3. Indexing the BAM file for the following tools to use

  ## Picard
  [Picard tools](https://broadinstitute.github.io/picard/) are used to mark duplicated reads in the aligned and sorted BAM file, making thus the files lighter and less prone to errors in the downstream steps of the pipeline.

  ## GATK
  [Genome Analysis Tool Kit](https://software.broadinstitute.org/gatk) refered to as GATK (DePristo et al., 2011) is used to process the data throught multiple steps as described by the [GATK best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS) (i.e. figure bellow).
  ![GATK best-practices pipeline](https://raw.githubusercontent.com/h3abionet/h3agatk/master/workflows/GATK/gatk.png)
  The GATK steps are the following:
  	1. Indel Realignment:
  		1. [Realign Target Creator](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php)
  		2. [Indel Realigner](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)
  	2. Mark Duplicates (a picard step)
  	3. Base Quality Score Recalibration (BQSR):
  		1. [Base Recalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
  		2. [Print Reads](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php)
  	4. [Haplotype Caller](https://software.broadinstitute.org/gatk/documentation/tooldocs/)
  	5. Variant Quality Score Recalibration (VQSR):
  		1. [Variant Recalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php)
  		2. [Apply Recalibration](https://software.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.php)

  ## SnpEff
  SNPEff is used in this pipeline to annotate the variant calls (Cingolani et al., 2012). The annotation is extensive and uses multi-database approach to provide the user with as much information about the called variants as possible.

  ## BAMStat
  [BAMStats](http://bamstats.sourceforge.net), is a simple software tool built on the Picard Java API (2), which can calculate and graphically display various metrics derived from SAM/BAM files of value in QC assessments.

inputs:
  reference:
    type: File
    doc: reference human genome file

  uncompressed_reference:
    type: File
    doc: reference human genome file

  reads:
    type: File[]?
    doc: files containing the paired end reads in fastq format required for bwa-mem

  bwa_output_name:
    type: string
    doc: name of bwa-mem output file

  bwa_read_group:
    type: string
    doc: read group

  bwa_threads:
    type: int
    doc: number of threads

  gatk_threads:
    type: int
    doc: number of threads

  samtools_threads:
    type: int
    doc: number of threads

  output_RefDictionaryFile:
    type: string
    doc: output file name for picard create dictionary command from picard toolkit

  samtools-view-isbam:
    type: boolean
    doc: boolean set to output bam file from samtools view

  samtools-index-bai:
    type: boolean
    doc: boolean set to output bam file from samtools view

  output_samtools-view:
    type: string
    doc: output file name for bam file generated by samtools view

  output_samtools-sort:
    type: string
    doc: output file name for bam file generated by samtools sort

  outputFileName_MarkDuplicates:
    type: string
    doc: output file name generated as a result of Markduplicates command from picard toolkit

  metricsFile_MarkDuplicates:
    type: string
    doc: metric file generated by MarkDupicates command listing duplicates

  readSorted_MarkDuplicates:
    type: string
    doc: set to be true showing that reads are sorted already

  removeDuplicates_MarkDuplicates:
    type: string
    doc: set to be true

  createIndex_MarkDuplicates:
    type: string
    doc: set to be true to create .bai file from Picard Mark Duplicates

  outputFileName_RealignTargetCreator:
    type: string
    doc: name of realignTargetCreator output file

  known_variant_db:
    type: File[]?
    doc: array of known variant files for realign target creator

  outputFileName_IndelRealigner:
    type: string
    doc: name of indelRealigner output file

  outputFileName_BaseRecalibrator:
    type: string
    doc: name of BaseRecalibrator output file

  outputFileName_PrintReads:
    type: string
    doc: name of PrintReads command output file

  outputFileName_HaplotypeCaller:
    type: string
    doc: name of Haplotype caller command output file

  dbsnp:
    type: File
    doc: vcf file containing SNP variations used for Haplotype caller

  tmpdir:
    type: string?
    doc: temporary dir for picard

  samtools-view-sambam:
    type: string?
    doc: temporary dir for picard

  covariate:
    type: string[]?
    doc: required for base recalibrator

  depth_omitIntervalStatistics:
    type: boolean?
    doc: Do not calculate per-interval statistics

  depth_omitDepthOutputAtEachBase:
    type: boolean?
    doc: Do not output depth of coverage at each base

  depth_outputfile_DepthOfCoverage:
    type: string?
    doc: name of the output report basename

  filter_expression:
    type: string
    default: "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

  snpf_genome:
    type: string

  snpf_nodownload:
    type: boolean

  snpf_data_dir:
    type: Directory

  resource_mills:
    type: File

  resource_hapmap:
    type: File

  resource_omni:
    type: File

  resource_dbsnp:
    type: File

outputs:

  output_bamstat:
    type: File
    outputSource: HaplotypeCaller/output_bamstat

  output_printReads:
    type: File
    outputSource: HaplotypeCaller/output_printReads

  output_HaplotypeCaller:
    type: File
    outputSource: HaplotypeCaller/output_HaplotypeCaller

  output_SnpVQSR_recal_File:
    type: File
    outputSource: SnpVQSR/recal_File

  output_SnpVQSR_annotated_snps:
    type: File
    outputSource: SnpVQSR/annotated_snps

  output_IndelFilter_annotated_indels:
    type: File
    outputSource: IndelFilter/annotated_indels

steps:

  HaplotypeCaller:
    run: GATK-Sub-Workflow-Workflow-h3abionet-haplotype.cwl
    in:
      reference: reference
      uncompressed_reference: uncompressed_reference
      reads: reads
      bwa_output_name: bwa_output_name
      bwa_read_group: bwa_read_group
      bwa_threads: bwa_threads
      gatk_threads: gatk_threads
      samtools_threads: samtools_threads
      output_RefDictionaryFile: output_RefDictionaryFile
      samtools-view-isbam: samtools-view-isbam
      samtools-index-bai: samtools-index-bai
      output_samtools-view: output_samtools-view
      output_samtools-sort: output_samtools-sort
      outputFileName_MarkDuplicates: outputFileName_MarkDuplicates
      metricsFile_MarkDuplicates: metricsFile_MarkDuplicates
      readSorted_MarkDuplicates: readSorted_MarkDuplicates
      removeDuplicates_MarkDuplicates: removeDuplicates_MarkDuplicates
      createIndex_MarkDuplicates: createIndex_MarkDuplicates
      outputFileName_RealignTargetCreator: outputFileName_RealignTargetCreator
      known_variant_db: known_variant_db
      outputFileName_IndelRealigner: outputFileName_IndelRealigner
      outputFileName_BaseRecalibrator: outputFileName_BaseRecalibrator
      outputFileName_PrintReads: outputFileName_PrintReads
      outputFileName_HaplotypeCaller: outputFileName_HaplotypeCaller
      dbsnp: dbsnp
      tmpdir: tmpdir
      samtools-view-sambam: samtools-view-sambam
      covariate: covariate
      depth_omitIntervalStatistics: depth_omitIntervalStatistics
      depth_omitDepthOutputAtEachBase: depth_omitDepthOutputAtEachBase
      depth_outputfile_DepthOfCoverage: depth_outputfile_DepthOfCoverage
    out: [ output_bamstat, output_printReads, output_HaplotypeCaller ]

  SnpVQSR:
    run: GATK-Sub-Workflow-h3abionet-snp.cwl
    in:
      reference: reference
      snpf_genome: snpf_genome
      snpf_nodownload: snpf_nodownload
      snpf_data_dir: snpf_data_dir
      resource_mills: resource_mills
      haplotest_vcf: HaplotypeCaller/output_HaplotypeCaller
      resource_hapmap: resource_hapmap
      resource_omni: resource_omni
      resource_dbsnp: resource_dbsnp
    out: [ recal_File, annotated_snps ]

  IndelFilter:
    run: GATK-Sub-Workflow-h3abionet-indel-no-vqsr.cwl
    in:
      reference: reference
      snpf_genome: snpf_genome
      snpf_nodownload: snpf_nodownload
      snpf_data_dir: snpf_data_dir
      resource_mills: resource_mills
      haplotest_vcf: HaplotypeCaller/output_HaplotypeCaller
      resource_hapmap: resource_hapmap
      resource_omni: resource_omni
      resource_dbsnp: resource_dbsnp
      filter_expression: filter_expression
    out: [ annotated_indels ]
