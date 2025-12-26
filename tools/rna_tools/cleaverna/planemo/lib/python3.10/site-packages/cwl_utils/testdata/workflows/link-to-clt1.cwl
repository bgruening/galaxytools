#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: inputs.txt
      entry: |
        $(inputs.in1.file.path)
        $(inputs.in1.meta.species)
  SchemaDefRequirement:
    types:
    - $import: ../types/recursive.yml
    - $import: ../types/array.yml
    - $import: ../types/singletype.yml
    # - $import: ../types/singletype2.yml

inputs:
  in1: ../types/recursive.yml#file_with_sample_meta
  in2:
    type: ../types/array.yml#study_meta_too
  in3:
    type: ../types/singletype.yml#simple_record
#   in4:
#     type: ../types/singletype2.yml#simple_record2
  in4:
    type: [string, ../types/recursive.yml#sample_meta]
  in5:
    type: Any?

outputs:
  out1:
    type: File
    outputBinding:
      glob: '*.txt'
  out2:
    type: ../types/array.yml#study_meta_too
    outputBinding:
      outputEval: $(inputs.in2)
  out3: stdout

baseCommand: [echo]
arguments: [hello world]
