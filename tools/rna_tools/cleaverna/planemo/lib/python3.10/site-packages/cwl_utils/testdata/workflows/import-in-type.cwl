#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.2

inputs:
 file1:
  type:
    $import: type-import.yaml

outputs:
  output:
    type: File
    outputBinding: { glob: output }

baseCommand: [wc, -l]

stdin: $(inputs.file1.path)
stdout: output
