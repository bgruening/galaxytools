#!/usr/bin/env cwl-runner
# Checks symbolic links on github

class: Workflow
cwlVersion: v1.2
inputs:
  in1:
    type: ../types/recursive.yml#file_with_sample_meta
  in2:
    type: ../types/array.yml#study_meta_too
  in3:
    type: ../types/singletype.yml#simple_record
  in4:
    type: [string, ../types/recursive.yml#sample_meta]

steps:
  s1:
    run: link-to-clt1.cwl
    in:
      in1: "#in1"   # This should be normalized out
      in2: in2
      in3: in3
      in4: in4
    out: [out2]

outputs:
  - id: out1
    type: ../types/array.yml#study_meta_too
    outputSource: "#s1/out2"

requirements:
  SchemaDefRequirement:
    types:
    - $import: ../types/recursive.yml
    - $import: ../types/array.yml
    - $import: ../types/singletype.yml
    - name: user_type1  # For tools/clt2.cwl
      type: record
      fields:
          - name: prop
            type: string
