#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
requirements:
  - class: InlineJavascriptRequirement
inputs:
  d: Directory
outputs:
  out:
    type: boolean
    outputBinding:
      outputEval: '$(inputs.d.listing.length === 1 && inputs.d.listing[0].listing.length === 1)'
baseCommand: "true"
