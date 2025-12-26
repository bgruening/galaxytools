#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.1

requirements:
  - class: ScatterFeatureRequirement

inputs:
    file1: File[]

outputs:
    count_output:
      type: int
      outputSource: step1/output

steps:
  step1:
    run: wc3-tool_v1_1.cwl
    scatter: file1
    in:
      file1:
        source: file1
        linkMerge: merge_nested
    out: [output]