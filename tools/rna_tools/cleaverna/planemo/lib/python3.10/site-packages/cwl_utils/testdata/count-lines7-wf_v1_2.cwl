#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  - class: MultipleInputFeatureRequirement

inputs:
    file1: File[]
    file2: File[]

outputs:
    count_output:
      type: int
      outputSource: step1/output

steps:
  step1:
    run: wc3-tool_v1_2.cwl
    in:
      file1:
        source: [file1, file2]
        linkMerge: merge_flattened
    out: [output]