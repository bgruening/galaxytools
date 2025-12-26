#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
inputs:
  val: int

steps:

  step1:
    in:
      in1: val
      a_new_var: val
    run: foo-array.cwl
    when: $(inputs.in1 < 1)
    out: [out1]

outputs:
  out1:
    type: string
    outputSource: step1/out1
    pickValue: first_non_null

requirements:
  InlineJavascriptRequirement: {}
