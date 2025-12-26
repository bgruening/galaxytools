#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

inputs:
  input_file: File

outputs:
  output_file:
    type: File
    outputSource: md5sum/output_file

steps:
  md5sum:
    run: dockstore-tool-md5sum_v12.cwl
    in:
      input_file: input_file
    out: [output_file]

