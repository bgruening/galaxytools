#!/usr/bin/env cwl-runner
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - EDAM_subset.owl
class: CommandLineTool
cwlVersion: v1.2
doc: "Reverse each line using the `rev` command"
hints:
  DockerRequirement:
    dockerPull: docker.io/debian:stable-slim

inputs:
  input:
    type: File
    inputBinding: {}
    format: edam:format_2330

outputs:
  output:
    type: File
    outputBinding:
      glob: output.txt
    format: $(inputs.input.format)

baseCommand: rev
stdout: output.txt
