#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}

inputs:
  first:
    type: File
    secondaryFiles: |
      ${ return self.basename+".2"; }
    default:
      class: File
      basename: default.txt
      contents: "42"
      secondaryFiles:
       - class: File
         basename: default.txt.2
         contents: "23"

steps:
  sf_extract:
    in:
      target: first
    run:
      class: ExpressionTool
      inputs:
        target: File
      expression: |
        ${ return {"result": inputs.target.secondaryFiles[0].basename}; }
      outputs:
        result: string
    out: [ result ]

outputs:
  sf_name:
    type: string
    outputSource: sf_extract/result
