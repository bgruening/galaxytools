#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
inputs:
  in1: int
baseCommand: [echo]
outputs:
  out1:
    type: string[]
    outputBinding:
      outputEval: ${ return [foo, inputs.in1]; }
requirements:
  InlineJavascriptRequirement: {}
