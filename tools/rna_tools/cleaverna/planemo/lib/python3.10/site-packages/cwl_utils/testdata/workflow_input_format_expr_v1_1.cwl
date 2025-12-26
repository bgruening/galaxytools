#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: Workflow

requirements:
  InlineJavascriptRequirement: {}

inputs:
  first:
    type: File
    format: |
      ${ return "http://edamontology.org/format_3016"; }
    default:
      class: File
      basename: default
      format: http://edamontology.org/format_3016
      contents: "42"

steps:
  format_extract:
    in:
      target: first
    run:
      class: ExpressionTool
      inputs:
        target: File
      expression: |
        ${ return {"result": inputs.target.format}; }
      outputs:
        result: string
    out: [ result ]

outputs:
  format_uri:
    type: string
    outputSource: format_extract/result
