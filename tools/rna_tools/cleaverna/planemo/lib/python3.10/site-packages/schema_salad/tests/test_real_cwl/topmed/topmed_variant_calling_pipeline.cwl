#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/9
label: TOPMed Variant Calling Pipeline CWL1
dct:creator:
 foaf:name: Seven Bridges
 foaf:mbox: "mailto:support@sbgenomics.com"
inputs:
  - id: reference
    'sbg:fileTypes': FA
    type: File
    'sbg:x': -460
    'sbg:y': -266.8571472167969
  - id: bam_cram_file
    'sbg:fileTypes': 'BAM, CRAM'
    type: 'File[]'
    'sbg:x': -457.8571472167969
    'sbg:y': -143.42857360839844
  - id: bai_crai_file
    'sbg:fileTypes': 'BAI, CRAI'
    type: 'File[]'
    'sbg:x': -456.4285888671875
    'sbg:y': -17.285715103149414
  - id: reference_file
    'sbg:fileTypes': TGZ
    type: File
    'sbg:x': -119
    'sbg:y': -629.1428833007812
  - id: pedigree_file
    'sbg:fileTypes': PED
    type: File
    'sbg:x': -119.71428680419922
    'sbg:y': -515.8571166992188
  - id: num_of_jobs
    type: int?
    'sbg:x': -121.57142639160156
    'sbg:y': -401.71429443359375
  - id: genotype_unit
    type: int
    'sbg:x': -127
    'sbg:y': -50.85714340209961
  - id: discover_unit
    type: int
    'sbg:x': -133
    'sbg:y': 73
  - id: chromosomes
    type: 'string[]'
    'sbg:x': -128.14285278320312
    'sbg:y': 198.2857208251953
  - id: reference_genome_1
    type:
      type: enum
      symbols:
        - hg19
        - hg38
      name: reference_genome
    'sbg:x': -459.2857360839844
    'sbg:y': -385.71429443359375
outputs:
  - id: called_variant_sites
    outputSource:
      - topmed_freeze3_calling/called_variant_sites
    type: File
    'sbg:x': 400.1481018066406
    'sbg:y': -54.76030349731445
  - id: genotypes
    outputSource:
      - topmed_freeze3_calling/genotypes
    type: File
    'sbg:x': 425.2412109375
    'sbg:y': -194.8717498779297
  - id: makefile_log
    outputSource:
      - topmed_freeze3_calling/makefile_log
    type: File?
    'sbg:x': 402.1481018066406
    'sbg:y': -341.23968505859375
steps:
  - id: verifybamid_cwl1
    in:
      - id: bai_crai_file
        source:
          - bai_crai_file
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: reference
        source:
          - reference
      - id: reference_genome
        source:
          - reference_genome_1
    out:
      - id: output_index_file
    run: verifybamid/verifybamid.cwl
    label: VerifyBamID_CWL1
    scatter:
      - bam_cram_file
    'sbg:x': -233.57144165039062
    'sbg:y': -197.14285278320312
  - id: topmed_freeze3_calling
    in:
      - id: bai_crai_files
        source:
          - bai_crai_file
      - id: bam_cram_files
        source:
          - bam_cram_file
      - id: chromosomes
        default: []
        source:
          - chromosomes
      - id: discover_unit
        source:
          - discover_unit
      - id: genotype_unit
        source:
          - genotype_unit
      - id: index_files
        source:
          - verifybamid_cwl1/output_index_file
      - id: num_of_jobs
        source:
          - num_of_jobs
      - id: pedigree_file
        source:
          - pedigree_file
      - id: reference_file
        source:
          - reference_file
      - id: reference_genome
        source:
          - reference_genome_1
    out:
      - id: called_variant_sites
      - id: genotypes
      - id: makefile_log
    run: topmed_freeze3_calling/topmed_freeze3_calling.cwl
    label: Topmed_freeze3_CWL1
    'sbg:x': 157.14285278320312
    'sbg:y': -198
requirements:
  - class: ScatterFeatureRequirement
$namespaces:
  sbg: 'https://sevenbridges.com'
'sbg:appVersion':
  - v1.0
'sbg:contributors':
  - vladimir_obucina
'sbg:createdBy': vladimir_obucina
'sbg:createdOn': 1526996458
'sbg:id': >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/9
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/9.png
'sbg:latestRevision': 9
'sbg:modifiedBy': vladimir_obucina
'sbg:modifiedOn': 1527071783
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:publisher': sbg
'sbg:revision': 9
'sbg:revisionNotes': >-
  UPDATE: changed VerifyBamID, error was lack od \n sign at the end of each
  output index file.
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1526996458
    'sbg:revision': 0
    'sbg:revisionNotes': null
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1526996882
    'sbg:revision': 1
    'sbg:revisionNotes': 'Firste Version with CWL1 tools, scatter on VerifyBAMId is of type none'
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1526997137
    'sbg:revision': 2
    'sbg:revisionNotes': Exposed ports
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1526997206
    'sbg:revision': 3
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1526997265
    'sbg:revision': 4
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527001305
    'sbg:revision': 5
    'sbg:revisionNotes': Separated bam and bai inputs for VerifyBAMId and Topmed_freeze3
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527007399
    'sbg:revision': 6
    'sbg:revisionNotes': Added output to VerifyBamID
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527060490
    'sbg:revision': 7
    'sbg:revisionNotes': 'UPDATE: Removed symlinks'
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527060551
    'sbg:revision': 8
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527071783
    'sbg:revision': 9
    'sbg:revisionNotes': >-
      UPDATE: changed VerifyBamID, error was lack od \n sign at the end of each
      output index file.
'sbg:sbgMaintained': false
'sbg:validationErrors': []
