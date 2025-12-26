#!/usr/bin/env cwl-runner
# We have this tool to test both local and remote packing

class: CommandLineTool
cwlVersion: v1.2
inputs:
  in1: 
    type: string
    inputBinding: 
      position: 1
      valueFrom: A_$(inputs.in1)_B_${return inputs.in1}_C_$(inputs.in1)
baseCommand: echo
arguments:
  - valueFrom: $(runtime)
outputs:
  out1:
    type: string
    outputBinding:
      glob: out.txt
      loadContents: true
      outputEval: $(self)_D_$(runtime)

stdout: out.txt
requirements:
  InlineJavascriptRequirement:
    expressionLib:
      - $include: ../lib.js
