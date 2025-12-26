#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
inputs:
  input_file: File?
steps:
  md5:
    run: md5.cwl
    out:
    - report
    in:
      input_file: 'input_file'
  validatefiles:
    run: validate.cwl
    out:
    - report
    in:
      input_file: 'input_file'
      type: {}
outputs:
  validatefiles_report:
    type: File?
    outputSource: validatefiles/report
  md5_report:
    type: File?
    outputSource: md5/report
