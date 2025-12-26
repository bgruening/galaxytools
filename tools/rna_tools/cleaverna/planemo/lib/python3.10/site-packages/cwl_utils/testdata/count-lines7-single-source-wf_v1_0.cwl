#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

inputs:
    file1: File[]

outputs:
    count_output:
      type: int
      outputSource: step1/output

steps:
  step1:
    run: wc3-tool_v1_0.cwl
    in:
      file1:
        source: file1
        linkMerge: merge_flattened
    out: [output]