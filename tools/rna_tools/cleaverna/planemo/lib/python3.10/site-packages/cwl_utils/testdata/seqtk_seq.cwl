#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: "Convert to FASTA (seqtk)"

baseCommand: ["seqtk", "seq"]

inputs:
  sequences:
    type: File
    inputBinding:
      prefix: "-a"

outputs:
  result: stdout

hints:
  SoftwareRequirement:
    packages:
    - package: seqtk
      version:
      - r93
