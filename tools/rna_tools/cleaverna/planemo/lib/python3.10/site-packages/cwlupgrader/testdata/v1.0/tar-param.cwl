#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, --extract]
inputs:
  tarfile:
    type: File
    inputBinding:
      prefix: --file
  extractfile:
    type: string
    inputBinding:
      position: 1
outputs:
  extracted_file:
    type: File
    outputBinding:
      glob: $(inputs.extractfile)
