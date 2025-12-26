#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  gatk_jar:
    type: File
    doc: Jar executable of the GATK tool

  reference:
    type: File
    doc: reference human genome file

  indel_mode:
    type: string
    default: 'INDEL'

  snp_mode:
    type: string
    default: 'SNP'

  snpf_genome:
    type: string

  snpf_nodownload:
    type: boolean

  snpf_data_dir:
    type: Directory

  resource_mills:
    type: File

  resource_1kg:
    type: File

  haplotest_vcf:
    type: File

  resource_hapmap:
    type: File

  resource_omni:
    type: File

  resource_dbsnp:
    type: File


  snpeff_java_mem_opts:
      type: string[]?

outputs:
  recal_File:
    type: File
    outputSource: vqsr_snps/recal_File

  annotated_snps:
    type: File
    outputSource: snpeff_snps/annotated_vcf

steps:

  vqsr_snps:
    run: ../../tools/GATK-VariantRecalibrator-SNPs.cwl
    in:
      gatk_jar: gatk_jar
      haplotypecaller_snps_vcf: haplotest_vcf
      reference: reference
      resource_dbsnp: resource_dbsnp
      resource_omni: resource_omni
      resource_1kg: resource_1kg
      resource_hapmap: resource_hapmap
    out: [tranches_File, recal_File]

  apply_recalibration_snps:
    run: ../../tools/GATK-ApplyRecalibration.cwl
    in:
      gatk_jar: gatk_jar
      #raw_vcf: HaplotypeCaller/output_HaplotypeCaller
      raw_vcf: haplotest_vcf
      reference: reference
      recal_file: vqsr_snps/recal_File
      tranches_file: vqsr_snps/tranches_File
      mode: snp_mode
    out: [ vqsr_vcf ]

  snpeff_snps:
    run: ../../tools/snpEff.cwl
    in:
      genome: snpf_genome
      variant_calling_file: apply_recalibration_snps/vqsr_vcf
      nodownload: snpf_nodownload
      data_dir: snpf_data_dir
      java_mem_opts: snpeff_java_mem_opts
    out: [ annotated_vcf ]
