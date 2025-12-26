#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
inputs:
 09first_input: string
 05second_input: int
 01third_input: File
steps:
 zz_step_one:
   run:
     class: ExpressionTool
     inputs: []
     outputs: []
     expression: ${return {}; }
     requirements:
       InlineJavascriptRequirement: {}
   in: []
   out: []
 00_step_two:
   out: []
   run:
     inputs: []
     requirements:
       InlineJavascriptRequirement: {}
     outputs: []
     expression: ${return {}; }
     class: ExpressionTool
   in: []
outputs:
  zz_first_output:
    type: File
    outputSource: 01third_input
  ll_second_output:
    type: string
    outputSource: 09first_input
  aa_third_output:
    type: int
    outputSource: 05second_input
