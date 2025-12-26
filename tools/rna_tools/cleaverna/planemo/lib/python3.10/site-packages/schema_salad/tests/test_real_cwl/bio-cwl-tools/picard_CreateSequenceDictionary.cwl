#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

baseCommand:
- picard
- CreateSequenceDictionary

doc: |-
  Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
   SAMRecords, and the header contains only sequence records.

requirements:
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.REFERENCE)
  InlineJavascriptRequirement:
    expressionLib:
    - |
      function generateGATK4BooleanValue(){
          /**
           * Boolean types in GATK 4 are expressed on the command line as --<PREFIX> "true"/"false",
           * so patch here
           */
          if(self === true || self === false){
              return self.toString()
          }

          return self;
      }
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/picard:2.22.2--0
inputs:
- doc: Input reference fasta or fasta.gz [synonymous with -R]
  id: REFERENCE
  type: File
  format: edam:format_1929  # FASTA
  inputBinding:
    valueFrom: REFERENCE=$(self.basename)
- doc: Put into AS field of sequence dictionary entry if supplied [synonymous with
    -AS]
  id: GENOME_ASSEMBLY
  type: string?
  inputBinding:
    prefix: GENOME_ASSEMBLY=
    separate: false
- doc: Put into UR field of sequence dictionary entry.  If not supplied, input reference
    file is used [synonymous with -UR]
  id: URI
  type: string?
  inputBinding:
    prefix: URI=
    separate: false
- doc: Put into SP field of sequence dictionary entry [synonymous with -SP]
  id: SPECIES
  type: string?
  inputBinding:
    prefix: SPECIES=
    separate: false
- doc: Make sequence name the first word from the > line in the fasta file.  By default
    the entire contents of the > line is used, excluding leading and trailing whitespace.
  id: TRUNCATE_NAMES_AT_WHITESPACE
  type: boolean?
  inputBinding:
    prefix: TRUNCATE_NAMES_AT_WHITESPACE=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: Stop after writing this many sequences.  For testing.
  id: NUM_SEQUENCES
  type: int?
  inputBinding:
    prefix: NUM_SEQUENCES=
    separate: false
- doc: "Optional file containing the alternative names for the contigs. Tools may\
    \ use this information to consider different contig notations as identical (e.g:\
    \ 'chr1' and '1'). The alternative names will be put into the appropriate @AN\
    \ annotation for each contig. No header. First column is the original name, the\
    \ second column is an alternative name. One contig may have more than one alternative\
    \ name. [synonymous with -AN]"
  id: ALT_NAMES
  type: File?
  inputBinding:
    prefix: ALT_NAMES=
    separate: false
- doc: Control verbosity of logging.
  id: VERBOSITY
  type:
  - 'null'
  - type: enum
    symbols:
    - ERROR
    - WARNING
    - INFO
    - DEBUG
  inputBinding:
    prefix: VERBOSITY=
    separate: false
- doc: Whether to suppress job-summary info on System.err.
  id: QUIET
  type: boolean?
  inputBinding:
    prefix: QUIET=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: Validation stringency for all SAM files read by this program.  Setting stringency
    to SILENT can improve performance when processing a BAM file in which variable-length
    data (read, qualities, tags) do not otherwise need to be decoded.
  id: VALIDATION_STRINGENCY
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - LENIENT
    - SILENT
  inputBinding:
    prefix: VALIDATION_STRINGENCY=
    separate: false
- doc: Compression level for all compressed files created (e.g. BAM and VCF).
  id: COMPRESSION_LEVEL
  type: int?
  inputBinding:
    prefix: COMPRESSION_LEVEL=
    separate: false
- doc: When writing files that need to be sorted, this will specify the number of
    records stored in RAM before spilling to disk. Increasing this number reduces
    the number of file handles needed to sort the file, and increases the amount of
    RAM needed.
  id: MAX_RECORDS_IN_RAM
  type: int?
  inputBinding:
    prefix: MAX_RECORDS_IN_RAM=
    separate: false
- doc: Use the JDK Deflater instead of the Intel Deflater for writing compressed output
    [synonymous with -use_jdk_deflater]
  id: USE_JDK_DEFLATER
  type: boolean?
  inputBinding:
    prefix: USE_JDK_DEFLATER=
    separate: false
    valueFrom: $(generateGATK4BooleanValue())
- doc: Use the JDK Inflater instead of the Intel Inflater for reading compressed input
    [synonymous with -use_jdk_inflater]
  id: USE_JDK_INFLATER
  type: boolean?
  inputBinding:
    prefix: USE_JDK_INFLATER=
    separate: false
    valueFrom: $(generateGATK4BooleanValue())
- doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
  id: CREATE_INDEX
  type: boolean?
  inputBinding:
    prefix: CREATE_INDEX=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: 'Whether to create an MD5 digest for any BAM or FASTQ files created.  '
  id: CREATE_MD5_FILE
  type: boolean?
  inputBinding:
    prefix: CREATE_MD5_FILE=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: Google Genomics API client_secrets.json file path.
  id: GA4GH_CLIENT_SECRETS
  type: File?
  inputBinding:
    prefix: GA4GH_CLIENT_SECRETS=
    separate: false

arguments:
 - TMP_DIR=$(runtime.tmpdir)
 - OUTPUT=$(inputs.REFERENCE.nameroot).dict

outputs:
  sequences_with_dictionary:
    type: File
    format: edam:format_2573  # SAM
    secondaryFiles:
      - ^.dict
      - .fai?
    outputBinding:
      glob: $(inputs.REFERENCE.basename)

  sequence_dictionary:
    type: File
    outputBinding:
      glob: $(inputs.REFERENCE.nameroot).dict

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
